# GP-representations
Scripts used to generate group and period-based Coordination Atom Fingerprints and FCHL19 representations (GP-CAF and GP-FCHL19).
As this is an improvement based on 
[Coordination Atom Fingerprints](https://doi.org/10.1038/s41929-018-0142-1) and 
[FCHL19](https://aip.scitation.org/doi/abs/10.1063/1.5126701). 

It contains lots code from the repos of above publications, which are

[GASpy](https://github.com/ulissigroup/GASpy),
[GASpy_regressions](https://github.com/ulissigroup/GASpy_regressions),
[qml](https://github.com/qmlcode/qml)

As well as the repo used to get the dataset

[CatHub](https://github.com/SUNCAT-Center/CatHub)

## Installation

**Requirements**

**1. Python environment**

We recommend using Conda package manager

```bash
conda create -n gprep python=3.7.3
source activate gprep
```

**2. Essential python packages**
  - ase=3.18.1
  - pymatgen=2019.11.11
  - scikit-learn=0.22.1
  - tpot=0.11.0
  - xgboost=0.90
  - mendeleev==0.5.1
Thes can be installed through
```bash
conda env update --file gprep.yml
```

**3. Installing QML, CatHub**

As the GP-FCHL19 based on the Fortran/Python interface used in qml, qml is essential to the installation of this code. Please go to http://www.qmlcode.org/installation.html to get more information about QML installation. If QML cannot be installed, this code cannot as well.

Install QML
```bash
git clone -b develop https://github.com/qmlcode/qml.git
cd qml
python setup.py install
cd ../
```

Install CatHub
```bash
https://github.com/SUNCAT-Center/CatHub.git
cd CatHub
pip install .
cd ../
```

**4. Installing GPREP**

```bash
python setup.py install
```

## Usage 
```bash
python get_gpcaf.py
python get_gpfchl19.py
```

## Acknowledgements
- The GP-CAF is based on the CAF implemented in [GASpy_regressions](https://github.com/ulissigroup/GASpy_regressions).
- The GP-FCHL19 is based on the FCHL19 implemented in [qml](https://github.com/qmlcode/qml).
