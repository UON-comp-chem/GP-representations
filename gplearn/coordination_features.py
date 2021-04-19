#!/usr/bin/env python
# coding: utf-8
# GNU Lesser General Public License v3.0

"""
This source code is licensed under the GNU Lesser General Public License v3.0
found in the LICENSE file in the root directory of this source tree.

---
This code borrows heavily from
https://github.com/ulissigroup/GASpy_regressions/blob/master/gaspy_regress/fingerprinters.py
---

"""


__author__ = ''
__email__ = 'ktran@andrew.cmu.edu'

from collections import defaultdict
import warnings
from abc import ABC, abstractmethod
import pickle
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', message='numpy.dtype size changed')
    import mendeleev
from pymatgen.ext.matproj import MPRester
from .defaults import get_median_energies, group, period
from ase.data import atomic_numbers as anumbers


CACHE_LOCATION = "../data/"

class Fingerprinter(ABC, BaseEstimator, TransformerMixin):
    '''
    This is a template fingerprinter that is meant to be extended before using.
    It is meant to be extended, and needs a `transform` method.
    The especially useful attributes that this class has are `dummy_fp_`,
    `max_num_species_`, `median_adsorption_energies_`, and `mendeleev_data_`.
    For more details on what they are, refer to the respective methods.
    Refer to Tran & Ulissi (Nature Catalysis, 2018) for even more details.
    '''
    def __init__(self, features = ["atomic_number", "median_energy", "electronegativity", 'count'],  elements = 'auto', adsorbate = 'H'):
        self.features = features
        self.elements = elements
        self.adsorbate = adsorbate


    def fit(self, X, y=None):
        '''
        This method behaves as a `fit` method from SKLearn would:  given
        training data, it calculates and assigns the attributes needed
        to run the `transform` method.
        Args:
            X   A list of dictionaries that should be fetched using the
                `gaspy.gasdb.get_adsorption_docs` function. You can change the
                filter used to get the documents, but you probably shouldn't
                remove any fingerprints unless you know what you're doing.
            y   Argument required to be able to fit within the SKLearn
                framework. You should ignore this argument and not pass
                anything to it.
        Returns:
            self
        '''
        # Get the data that we need to calculate the prerequisite information
        self.adsorption_docs = X

        # Calculate the information we need to make a fingerprint
        self._calculate_dummy_fp()
        self._get_max_num_species()

        # Delete some data to save memory
        del self.adsorption_docs

        return self


    def fit_transform(self, X, y=None):
        '''
        This method behaves as a `fit_transform` method from SKLearn would:
        given training data, it calculates and assigns the attributes needed
        to run the `transform` method. It then returns the transformation
        of the training data.
        Args:
            X   A list of dictionaries that should be fetched using the
                `gaspy.gasdb.get_adsorption_docs` function. You can
                change the filter used to get the documents, but you
                probably shouldn't remove any fingerprints unless you
                know what you're doing.
            y   Argument required to be able to fit within
                the SKLearn framework. You should ignore this argument
                and not pass anything to it.
        Returns:
            fingerprints    A numpy.array object that is a numerical representation the
                            documents that you gave this method, as per the docstring of
                            this class.
        '''
        self.fit(X, y=y)
        fingerprints = self.transform(X)
        return fingerprints


    def transform(self, docs):
        fingerprints = np.array([self.fingerprint_doc(doc) for doc in docs])
        return fingerprints


    def fingerprint_doc(self, doc):
        '''
        Convert a document into a numerical fingerprint.
        Inputs:
            doc     A dictionary that should have the keys 'mpid' and 'coordination'.
                    The value for 'mpid' should be in the form 'mpid-23' and the value
                    for 'coordination' should be in the form 'Cu-Cu-Cu'.
                    Should probably come from the `gaspy.gasdb.get_catalog_docs` function.
        Returns:
            fingerprint A numpy.array object that is a numerical representation the
                        document that you gave this method, as per the docstring of
                        this class. Note that the array is actually a flattened,
                        1-dimensional object.
        '''
        fingerprint = []
        shell_atoms = self._concatenate_shell(doc)
        adsorbate = self.adsorbate
        # Add and sort the elemental information for each element present
        
        for element in set(shell_atoms):
            rep = []
            element_data = self.mendeleev_data_[element]
            for feature in self.features:
                if feature == 'electronegativity':
                    v = element_data.electronegativity(scale='pauling')
                    rep.append(v)
                elif feature == 'median_energy':
                    try:
                        median_energies_by_element = self.median_adsorption_energies_[adsorbate]
                        try:
                            median_energy = median_energies_by_element[element]
                        except KeyError:
                            raise RuntimeError('You did not initialize the '
                                               'fingerprinter with the %s element, and '
                                               'so we cannot make predictions on it.'
                                               % element)
                    except KeyError:
                        raise RuntimeError('You did not initialize the fingerprinter '
                                           'with the %s adsorbate, and so we cannot make'
                                           'predictions on it.' % adsorbate)

                    rep.append(median_energy)
                elif feature == 'count':
                    count = shell_atoms.count(element)
                    rep.append(count)
                else:
                    v = getattr(element_data, feature)
                    rep.append(v)
            rep = tuple(rep)
            fingerprint.append(rep)
        
                
        fingerprint = sorted(fingerprint)

        # Fill in the dummy fingerprints
        for _ in range(len(fingerprint), self.max_num_species_):
            fingerprint.append(self.dummy_fp_[adsorbate])

        return np.array(fingerprint).flatten()


    @abstractmethod
    def _concatenate_shell(self, doc):
        pass


    def _calculate_dummy_fp(self):
        '''
        This method establishes a "dummy" value for a 1x4 vector of
        information. The numbers in the 1x4 vector is average of all the median
        adsorption energies we have between a given adsorbate and various
        monometallics, the average atomic number of all elements we're
        considering, their average Pauling electronegativity, and an atomic
        count of zero.
        This dummy value is useful when using variable number of features;
        reference Davie et al (Kriging atomic properties with a variable number
        of inputs, J Chem Phys 2016). The out-of-bounds feature we choose is
        the atomic count.
        Resulting attributes:
            dummy_fp_   A dictionary whose keys are adsorbates and whose values
                        are 4-tuples that represents a single fingerprint, but
                        has the "dummy" values
        '''
        # Prerequisite calculations
        self._get_compositions_by_mpid()
        self._get_elements_in_scope()
        self._get_mendeleev_data()
        self._calculate_median_adsorption_energies()
        elements = self.elements_

        # Calculate `dummy_fp_` elements that are element-based
        dummy_count = 0

        # Calculate the adsorbate-specific average-median binding energies
        self.dummy_fp_ = {}
        temp_dict = {}
        for feature in self.features:
            if feature == 'electronegativity':
                    v = np.average([self.mendeleev_data_[element].electronegativity(scale='pauling')
                                     for element in elements])
                    temp_dict[feature] = v
            elif feature == 'median_energy':
                pass
            elif feature == 'count':
                count = 0
                temp_dict[feature]  = count
            else:
                avg = np.average([[getattr(self.mendeleev_data_[element], feature)
                                     for element in elements]])
                temp_dict[feature]  = avg
                
        for adsorbate, median_energies_by_element in self.median_adsorption_energies_.items():
            self.dummy_fp_[adsorbate] = []
            for feature in self.features:
                if feature == 'median_energy':
                    energies = list(median_energies_by_element.values())
                    self.dummy_fp_[adsorbate].append(np.mean(energies))
                else:
                    self.dummy_fp_[adsorbate].append(temp_dict[feature])
            self.dummy_fp_[adsorbate] = tuple(self.dummy_fp_[adsorbate])



    def _get_compositions_by_mpid(self):
        '''
        We use the Materials Project's python API to find the composition of
        various materials given their MPIDs. This can take awhile though, so we
        also cache the results and modify the cache as necessary.
        Resulting attribute:
            compositions_by_mpid_   A dictionary whose keys are MPIDs and whose
                                    values are lists of strings for each
                                    element that is present in the
                                    corresponding material. This object is
                                    cached and therefore may have extra
                                    key:value pairings that you may not need.
        '''
        if 'slab_name' in self.adsorption_docs[0].keys():
            compositions_by_mpid = {}
            # method come from gaspy, in here slab_name is unique which could work as mpid
            mpids = {doc['slab_name']:[doc['metalA'],doc['metalB']] for doc in self.adsorption_docs}

            # If metalB is None, then delete it
            for mpid in mpids.keys():
                if None in mpids[mpid]:
                    mpids[mpid] = [mpids[mpid][0]]
            compositions_by_mpid_ = mpids
            
        else:
            with open(CACHE_LOCATION + 'mp_comp_data.pkl', 'rb') as file_handle:
                compositions_by_mpid = pickle.load(file_handle)
            mpids = np.unique([n['mpid'] for n in self.adsorption_docs])
            compositions_by_mpid_ = {}
            for mpid in mpids:
                compositions_by_mpid_[mpid] = compositions_by_mpid[mpid]
        self.compositions_by_mpid_ = compositions_by_mpid_


    def _get_elements_in_scope(self):
        '''
        This class has quite a few attributes that use elements as keys.
        We set the scope of the elements (i.e., the keys for these attributes)
        here by figuring out what elements are included in the MPIDs that
        we are considering.
        Resulting attribute:
            elements    A set of strings, where each string is the 2-letter
                        symbol for each element that shows up in the MPIDs
                        that we are considering. The MPIDs we are considering
                        are all of the MPIDs we can find in the catalog
                        and adsorption collections.
        '''
        if isinstance(self.elements, str):
            if self.elements =='auto':
                elements = set()

                # First grab all of the things in our cache of bulks
                for _, composition in self.compositions_by_mpid_.items():
                    elements = elements | set(composition)

                # Then add in anything in our training data
                for doc in self.adsorption_docs:
                    for neighbor in doc['neighborcoord']:
                        next_neighbors = neighbor.split(':')[-1].split('-')
                        elements = elements | set(next_neighbors)

                self.elements_ = elements
        else:
            self.elements_ = self.elements


    def _get_mendeleev_data(self):
        '''
        This method will get all of the Mendeleev data for the substrates
        included in a set of aggregated Mongo documents.
        Resulting attribute:
            mendeleev_data  A dictionary whose keys are the elements present in
                            `docs` and whose values are the Mendeleev data
        '''
        # Get the Mendeleev data for each element
        mendeleev_data = dict.fromkeys(self.elements_)
        for element in mendeleev_data:
            mendeleev_data[element] = getattr(mendeleev, element)
        self.mendeleev_data_ = mendeleev_data


    def _calculate_median_adsorption_energies(self):
        '''
        This method calculates the median adsorption energies on each
        monometallic bulk.
        Resulting attribute:
            median_adsorption_energies_ A nested dictionary first layer of keys
                                        are the adsorbate and whose second
                                        layer are the element. The values are
                                        the median adsorption energy of their
                                        respective adsorbate-element pairing.
        '''
        # Different with GASpy, here the median energy is defined in defaults.py
        adsorption_energies = defaultdict(dict)
        median_adsorption_energies = defaultdict(dict)
        adsorbate = self.adsorbate
        median_adsorption_energies[adsorbate] = get_median_energies(adsorbate)
        
        self.median_adsorption_energies_ = median_adsorption_energies


    def _get_max_num_species(self):
        '''
        When populating "dummy fingerprints", we need to know how many of them to make.
        We set this number equal to the maximum number of elements present in any one
        alloy in the catalog, and we find this number here.
        Resulting attributes:
            max_num_species_    An integer for the maximum number of elements/species
                                present in any single mpid we are looking at. This is useful
                                for figuring out how many dummy features you need to add.
        '''
        mpids = set(self.compositions_by_mpid_.keys())
        num_species_per_mpid = [len(self.compositions_by_mpid_[mpid]) for mpid in mpids]
        self.max_num_species_ = max(num_species_per_mpid)


class InnerShellFingerprinter(Fingerprinter):
    '''
    This fingerprinter converts the "inner shell" atoms---i.e., the coordinated
    atoms---into a Nx4 array of numbers, where N is the maximum number of
    elements seen in any of the inner shells all sites in the catalog. Each 1x4
    vector corresponds to one of the elements present in this inner shell. The
    numbers in the 1x4 vectors are the element's median adsorption energy, its
    atomic number, its Pauling electronegativity, and the number of those
    elements that show up in the coordination. We also sort the 1x4 vectors
    such that the first 1x4 vector that shows up is the one with the lowest
    median adsorption energy.
    We also use "dummy" vectors to address the fact that we will have a
    variable number of features/elements present in the inner shell; reference
    Davie et al (Kriging atomic properties with a variable number of inputs, J
    Chem Phys 2016). The out-of-bounds feature we choose is the atomic count.
    '''

    @staticmethod
    def _concatenate_shell(doc):
        shell_atoms = doc['coordination'].split('-')

        # Sometimes there is no coordination. If this happens, then hackily reformat it
        if shell_atoms == ['']:
            shell_atoms = []

        return shell_atoms


class OuterShellFingerprinter(Fingerprinter):
    '''
    This fingerprinter converts the "outer shell" atoms---i.e., the next neighbor
    atoms---into a Nx4 array of numbers, where N is the maximum number of
    elements seen in any of the outer shells all sites in the catalog. Each 1x4
    vector corresponds to one of the elements present in this outer shell. The
    numbers in the 1x4 vectors are the element's median adsorption energy, its
    atomic number, its Pauling electronegativity, and the sum of the number
    of times that the element shows up as being coordinated with a binding atom.
    We also sort the 1x4 vectors such that the first 1x4 vector that shows up is
    the one with the lowest median adsorption energy.
    We also use "dummy" vectors to address the fact that we will have a
    variable number of features/elements present in the outer shell; reference
    Davie et al (Kriging atomic properties with a variable number of inputs, J
    Chem Phys 2016). The out-of-bounds feature we choose is the atomic count.
    '''

    @staticmethod
    def _concatenate_shell(doc):
        '''
        This is a helper method to parse a neighborcoord string and
        concatenate all of the neighbors of the binding atoms together. Note
        that the counting that we do here allows for redundant counting of
        atoms. In other words:  If an atom is bound to three different binding
        atoms, then it will show up in this method's output three times.
        Arg:
            doc     A dictionary with the 'neighborcoord' string, whose contents
                    should look like:
                        ['Cu:Cu-Cu-Cu-Cu-Cu-Al',
                         'Al:Cu-Cu-Cu-Cu-Cu-Cu']
        Returns:
            second_shell_atoms  An extended list of the coordinations of all
                                binding atoms. Continiuing from the example
                                shown in the description for the `doc` argument,
                                we would get:
                                ['Cu', 'Cu', 'Cu', 'Cu', 'Cu', Al, 'Cu', 'Cu', 'Cu', 'Cu', 'Cu', 'Cu']
        '''
        second_shell_atoms = []
        for neighbor_coordination in doc['neighborcoord']:
            _, coordination = neighbor_coordination.split(':')
            coordination = coordination.split('-')
            second_shell_atoms.extend(coordination)
        return second_shell_atoms

class NextNearestFingerprinter(Fingerprinter):
    '''
    This fingerprinter converts Next Nearest atoms
    '''

    @staticmethod
    def _concatenate_shell(doc):
        shell_atoms = doc['nextnearest'].split('-')

        # Sometimes there is no coordination. If this happens, then hackily reformat it
        if shell_atoms == ['']:
            shell_atoms = []

        return shell_atoms

class StackedFingerprinter(object):
    '''
    If you have multiple fingerprinters that you want to feed into
    the same pipeline, then you can use this class to stack/concatenate
    the fingerprinters into one object. This new, stacked fingerprinter
    will stack the results from every fingerprinter that you provide it.
    '''
    def __init__(self, *fingerprinters):
        '''
        Args:
            *fingerprinters A sequence of instances of the fingerprinter
                            class objects that you want stacked together.
                            We do not recommend fitting them before
                            passing them to this stacker.
        '''
        self.fingerprinters = fingerprinters


    def fit(self, X, y=None):
        '''
        Executes the `fit` methods of all the fingerprinters you passed
        during instantiation.
        Args:
            X   A list of dictionaries that should be fetched using the
                `gaspy.gasdb.get_adsorption_docs` function. You can
                change the filter used to get the documents, but you
                probably shouldn't remove any fingerprints unless you
                know what you're doing.
            y   Argument required to be able to fit within
                the SKLearn framework. You should ignore this argument
                and not pass anything to it.
        Returns:
            self
        '''
        for fingerprinter in self.fingerprinters:
            _ = fingerprinter.fit(X)    # noqa: F841
        return self


    def fit_transform(self, X, y=None):
        '''
        This method behaves as a `fit_transform` method from SKLearn would:
        given training data, it calculates and assigns the attributes needed
        to run the `transform` method. It then returns the transformation
        of the training data.
        Args:
            X   A list of dictionaries that should be fetched using the
                `gaspy.gasdb.get_adsorption_docs` function. You can
                change the filter used to get the documents, but you
                probably shouldn't remove any fingerprints unless you
                know what you're doing.
            y   Argument required to be able to fit within
                the SKLearn framework. You should ignore this argument
                and not pass anything to it.
        Returns:
            features    A numpy.array object that is a numerical representation the
                        documents that you gave this method, as per the docstring of
                        this class.
        '''
        self.fit(X, y=y)
        features = self.transform(X)
        return features


    def transform(self, docs):
        '''
        Convert a list of documents into a list of numerical fingerprints.
        Inputs:
            docs    A list of dictionaries that contain information you need for
                    fingerprinting. The required contents of these dictionaries
                    inherit the requirements of the fingerprinters that you are
                    stacking. Should probably come from the
                    `gaspy.gasdb.get_catalog_docs` function.
        Returns:
            fingerprints    A list of numpy.array objects. Each numpy array is a
                            numerical representation of each document that you gave this
                            method, as per the docstrings of the fingerprinters you used
                            to initialize this class. Note that the array is actually
                            a flattened, 1-dimensional object.
        '''
        tupled_fingerprint = tuple(fingerprinter.transform(docs)
                                   for fingerprinter in self.fingerprinters)
        stacked_fingerprint = np.concatenate(tupled_fingerprint, axis=1)
        return stacked_fingerprint
