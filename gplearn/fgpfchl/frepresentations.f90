module representations

    implicit none

contains

function decay(r, invrc, natoms) result(f)

    implicit none

    double precision, intent(in), dimension(:,:) :: r
    double precision, intent(in) :: invrc
    integer, intent(in) :: natoms
    double precision, dimension(natoms, natoms) :: f

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    ! Decaying function reaching 0 at rc
    f = 0.5d0 * (cos(pi * r * invrc) + 1.0d0)


end function decay

function calc_angle(a, b, c) result(angle)

    implicit none

    double precision, intent(in), dimension(3) :: a
    double precision, intent(in), dimension(3) :: b
    double precision, intent(in), dimension(3) :: c

    double precision, dimension(3) :: v1
    double precision, dimension(3) :: v2

    double precision :: cos_angle
    double precision :: angle

    v1 = a - b
    v2 = c - b

    v1 = v1 / norm2(v1)
    v2 = v2 / norm2(v2)

    cos_angle = dot_product(v1,v2)

    ! Clipping
    if (cos_angle > 1.0d0) cos_angle = 1.0d0
    if (cos_angle < -1.0d0) cos_angle = -1.0d0

    angle = acos(cos_angle)

end function calc_angle

function calc_cos_angle(a, b, c) result(cos_angle)

    implicit none

    double precision, intent(in), dimension(3) :: a
    double precision, intent(in), dimension(3) :: b
    double precision, intent(in), dimension(3) :: c

    double precision, dimension(3) :: v1
    double precision, dimension(3) :: v2

    double precision :: cos_angle

    v1 = a - b
    v2 = c - b

    v1 = v1 / norm2(v1)
    v2 = v2 / norm2(v2)

    cos_angle = dot_product(v1,v2)

end function calc_cos_angle

end module representations


subroutine fgenerate_local_bob(cent_pos, cent_charges, atomic_charges, coordinates, cutoff, nuclear_charges, lmax, nmax, cm)
    implicit none
    double precision, dimension(:), intent(in) :: cent_pos
    integer, intent(in) :: cent_charges
    double precision, dimension(:), intent(in) :: atomic_charges
    double precision, dimension(:,:), intent(in) :: coordinates
    double precision, intent(in) :: cutoff
    integer, dimension(:), intent(in) :: nuclear_charges
    integer, intent(in) :: lmax
    integer, intent(in) :: nmax
    double precision, dimension(lmax, nmax), intent(out) :: cm
    
    integer :: m, n, i, j, natoms, nelements
    integer, allocatable, dimension(:) :: element_types
    integer, allocatable, dimension(:) :: located
    
    double precision :: pair_norm
    double precision, allocatable, dimension(:) :: distance_matrix
    double precision, allocatable, dimension(:) :: bag
    

    if (size(coordinates, dim=1) /= size(atomic_charges, dim=1)) then
        write(*,*) "ERROR: Bag of Bonds generation"
        write(*,*) size(coordinates, dim=1), "coordinates, but", &
            & size(atomic_charges, dim=1), "atom_types!"
        stop
    else
        natoms = size(atomic_charges, dim=1)
    endif
    
    ! Allocate temporary
    allocate(distance_matrix(natoms))
    !$OMP PARALLEL DO PRIVATE(pair_norm)
    do i = 1, natoms
        pair_norm = sqrt(sum((coordinates(i,:) - cent_pos)**2))
        distance_matrix(i) = pair_norm
    enddo  
    !$OMP END PARALLEL DO
    
    nelements = size(nuclear_charges, dim=1)
    allocate(element_types(natoms))
    
    !$OMP PARALLEL DO
    do i = 1, natoms
        do j = 1, nelements
            if (atomic_charges(i) .eq. nuclear_charges(j)) then
                element_types(i) = j
                continue
            endif
        enddo
    enddo
    !$OMP END PARALLEL DO

    allocate(located(nelements))
    located = 1
    !$OMP PARALLEL DO PRIVATE (m, n)
    do i =1, natoms
        if (distance_matrix(i) .lt. cutoff ) then
            m = element_types(i)
            n=located(m)
            cm(m,n) = cent_charges * atomic_charges(i)/ distance_matrix(i) ** 2
            located(m) = located(m) + 1
        endif
    enddo
    !$OMP END PARALLEL DO

    allocate(bag(nmax))

    do i =1, lmax
        bag(:) = cm(i, :)
        do j = 1, nmax
            m = maxloc(bag(:), dim=1)
            cm(i, j) = bag(m)
            bag(m) = -1.
        end do
    end do

    
    ! Clean up
    deallocate(distance_matrix)
    deallocate(element_types)
    deallocate(bag)
    deallocate(located)
end subroutine fgenerate_local_bob

subroutine fgenerate_local_fchl_acsf(idx, coordinates, nuclear_charges, elements, &
                          & Rs2, Rs3, Ts, eta2, eta3, zeta, rcut, acut, natoms, rep_size, &
                          & two_body_decay, three_body_decay, three_body_weight, rep)
    use representations, only: decay, calc_angle, calc_cos_angle
    implicit none
    integer, intent(in) :: idx
    double precision, intent(in), dimension(:, :) :: coordinates
    integer, intent(in), dimension(:) :: nuclear_charges
    integer, intent(in), dimension(:) :: elements
    double precision, intent(in), dimension(:) :: Rs2
    double precision, intent(in), dimension(:) :: Rs3
    double precision, intent(in), dimension(:) :: Ts
    double precision, intent(in) :: eta2
    double precision, intent(in) :: eta3
    double precision, intent(in) :: zeta
    double precision, intent(in) :: rcut
    double precision, intent(in) :: acut
    integer, intent(in) :: natoms
    integer, intent(in) :: rep_size
    double precision, intent(in) :: two_body_decay
    double precision, intent(in) :: three_body_decay
    double precision, intent(in) :: three_body_weight
    
    double precision, intent(out), dimension(rep_size) :: rep
    
    integer :: i, j, k, l, n, m, o, p, q, s, z, nelements, nbasis2, nbasis3, nabasis
    integer, allocatable, dimension(:) :: element_types
    double precision :: rij, rik, angle, cos_1, cos_2, cos_3, invcut
    ! double precision :: angle_1, angle_2, angle_3
    double precision, allocatable, dimension(:) :: radial, angular, a, b, c
    double precision, allocatable, dimension(:, :) :: distance_matrix, rdecay 
    double precision, allocatable, dimension(:) :: rep3
    double precision :: mu, sigma, ksi3
    
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
    if (natoms /= size(nuclear_charges, dim=1)) then
        write(*,*) "ERROR: Atom Centered Symmetry Functions creation"
        write(*,*) natoms, "coordinates, but", &
            & size(nuclear_charges, dim=1), "atom_types!"
        stop
    endif
    ! number of element types
    nelements = size(elements)
    ! Allocate temporary
    allocate(element_types(natoms))
    ! Store element index of every atom
    ! !$OMP PARALLEL DO
    do i = 1, natoms
        do j = 1, nelements
            if (nuclear_charges(i) .eq. elements(j)) then
                element_types(i) = j
                continue
            endif
        enddo
    enddo
    ! !$OMP END PARALLEL DO
    ! Get distance matrix
    ! Allocate temporary
    allocate(distance_matrix(natoms, natoms))
    distance_matrix = 0.0d0
    !  !$OMP PARALLEL DO PRIVATE(rij)
    do i = 1, natoms
        do j = i+1, natoms
            rij = norm2(coordinates(j,:) - coordinates(i,:))
            distance_matrix(i, j) = rij
            distance_matrix(j, i) = rij
        enddo
    enddo
    ! !$OMP END PARALLEL DO
    ! number of basis functions in the two body term
    nbasis2 = size(Rs2)
    ! Inverse of the two body cutoff
    invcut = 1.0d0 / rcut
    ! pre-calculate the radial decay in the two body terms
    rdecay = decay(distance_matrix, invcut, natoms)
    ! Allocate temporary
    allocate(radial(nbasis2))
    
    rep = 0.0d0
    radial = 0.0d0
    ! 
    i = idx +1
    ! index of the element of atom i
    m = element_types(i)
    ! !$OMP PARALLEL DO PRIVATE(n,rij,radial) REDUCTION(+:rep)
    do j = 1, natoms
        ! index of the element of atom j
        n = element_types(j)
        ! distance between atoms i and j
        rij = distance_matrix(i,j)
        if ((rij <= rcut) .AND. (rij >= 0.1)) then
            ! two body term of the representation
            mu    = log(rij / sqrt(1.0d0 + eta2  / rij**2))
            sigma = sqrt(log(1.0d0 + eta2  / rij**2))
            radial(:) = 0.0d0
            do k = 1, nbasis2 
               radial(k) = 1.0d0/(sigma* sqrt(2.0d0*pi) * Rs2(k))&
                          & * rdecay(i,j)  * exp( - (log(Rs2(k)) - mu)**2 / (2.0d0 * sigma**2) ) &
                          & / rij**two_body_decay
            enddo

            rep((n-1)*nbasis2 + 1:n*nbasis2) = rep((n-1)*nbasis2 + 1:n*nbasis2) + radial
        endif
    enddo
    ! !$OMP END PARALLEL DO
    
    deallocate(radial)
    ! number of radial basis functions in the three body term
    nbasis3 = size(Rs3)
    ! number of radial basis functions in the three body term
    nabasis = size(Ts)
    ! Inverse of the three body cutoff
    invcut = 1.0d0 / acut
    ! pre-calculate the radial decay in the three body terms
    rdecay = decay(distance_matrix, invcut, natoms)
    ! Allocate temporary
    allocate(rep3(rep_size))
    allocate(a(3))
    allocate(b(3))
    allocate(c(3))
    allocate(radial(nbasis3))
    allocate(angular(nabasis))
    rep3 = 0.0d0
    ! This could probably be done more efficiently if it's a bottleneck
    ! Also the order is a bit wobbly compared to the tensorflow implementation
    ! !$OMP PARALLEL DO PRIVATE(rij, n, rik, m, a, b, c, angle, radial, angular, &
    ! !$OMP cos_1, cos_2, cos_3, mu, sigma, o, ksi3, &
    ! !$OMP p, q, s, z) REDUCTION(+:rep3) COLLAPSE(2) SCHEDULE(dynamic)
    i = idx + 1
    do j = 1, natoms - 1
        if (i .eq. j) cycle
        ! distance between atoms i and j
        rij = distance_matrix(i,j)
        if (rij > acut)  cycle
        ! index of the element of atom j
        n = element_types(j)
        do k = j + 1, natoms
            if (i .eq. k) cycle
            if (j .eq. k) cycle
            ! distance between atoms i and k
            rik = distance_matrix(i,k)
            if (rik > acut) cycle
            ! index of the element of atom k
            m = element_types(k)
            ! coordinates of atoms j, i, k
            a = coordinates(j,:)
            b = coordinates(i,:)
            c = coordinates(k,:)
            ! angle between atoms i, j and k centered on i
            angle   = calc_angle(a,b,c)
            cos_1 = calc_cos_angle(a,b,c)
            cos_2 = calc_cos_angle(a,c,b)
            cos_3 = calc_cos_angle(b,a,c)

            ! The radial part of the three body terms including decay
            radial = exp(-eta3*(0.5d0 * (rij+rik) - Rs3)**2) * rdecay(i,j) * rdecay(i,k)
           
            ksi3 = (1.0d0 + 3.0d0 * cos_1 * cos_2 * cos_3) &
                 & / (distance_matrix(i,k) * distance_matrix(i,j) * distance_matrix(j,k) &
             & )**three_body_decay * three_body_weight

            angular = 0.0d0 
            do l = 1, nabasis/2

                o = l*2-1
                angular(2*l-1) = angular(2*l-1) + 2*cos(o * angle) &
                    & * exp(-(zeta * o)**2 /2)
                
                angular(2*l) = angular(2*l) + 2*sin(o * angle) &
                    & * exp(-(zeta * o)**2 /2)

            enddo
            
            ! The lowest of the element indices for atoms j and k
            p = min(n,m) - 1
            ! The highest of the element indices for atoms j and k
            q = max(n,m) - 1
            ! calculate the indices that the three body terms should be added to
            s = nelements * nbasis2 + nbasis3 * nabasis * (-(p * (p + 1))/2 + q + nelements * p) + 1

            do l = 1, nbasis3
                ! calculate the indices that the three body terms should be added to
                z = s + (l-1) * nabasis
                ! Add the contributions from atoms i,j and k
                rep3(z:z + nabasis - 1) = rep3(z:z + nabasis - 1) + angular * radial(l) * ksi3
            enddo
        enddo
    enddo
    ! !$OMP END PARALLEL DO

    rep = rep + rep3
    deallocate(element_types)
    deallocate(rdecay)
    deallocate(distance_matrix)
    deallocate(rep3)
    deallocate(a)
    deallocate(b)
    deallocate(c)
    deallocate(radial)
    deallocate(angular)

end subroutine fgenerate_local_fchl_acsf

subroutine fgenerate_local_two_body_fchl_acsf(idx, coordinates, nuclear_charges, elements, &
                          & Rs2,  eta2, rcut, natoms, rep_size, two_body_decay, rep)
    use representations, only: decay, calc_angle, calc_cos_angle
    implicit none
    integer, intent(in) :: idx
    double precision, intent(in), dimension(:, :) :: coordinates
    integer, intent(in), dimension(:) :: nuclear_charges
    integer, intent(in), dimension(:) :: elements
    double precision, intent(in), dimension(:) :: Rs2
    double precision, intent(in) :: eta2
    double precision, intent(in) :: rcut
    integer, intent(in) :: natoms
    integer, intent(in) :: rep_size
    double precision, intent(in) :: two_body_decay  
    double precision, intent(out), dimension(rep_size) :: rep
    
    integer :: i, j, k, n, m, nelements, nbasis2
    integer, allocatable, dimension(:) :: element_types
    double precision :: rij,  invcut
    ! double precision :: angle_1, angle_2, angle_3
    double precision, allocatable, dimension(:) :: radial
    double precision, allocatable, dimension(:, :) :: distance_matrix, rdecay
    double precision :: mu, sigma
    
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
    if (natoms /= size(nuclear_charges, dim=1)) then
        write(*,*) "ERROR: Atom Centered Symmetry Functions creation"
        write(*,*) natoms, "coordinates, but", &
            & size(nuclear_charges, dim=1), "atom_types!"
        stop
    endif
    ! number of element types
    nelements = size(elements)
    ! Allocate temporary
    allocate(element_types(natoms))
    ! Store element index of every atom
    ! !$OMP PARALLEL DO
    do i = 1, natoms
        do j = 1, nelements
            if (nuclear_charges(i) .eq. elements(j)) then
                element_types(i) = j
                continue
            endif
        enddo
    enddo
    ! !$OMP END PARALLEL DO
    ! Get distance matrix
    ! Allocate temporary
    allocate(distance_matrix(natoms, natoms))
    distance_matrix = 0.0d0
    !  !$OMP PARALLEL DO PRIVATE(rij)
    do i = 1, natoms
        do j = i+1, natoms
            rij = norm2(coordinates(j,:) - coordinates(i,:))
            distance_matrix(i, j) = rij
            distance_matrix(j, i) = rij
        enddo
    enddo
    ! !$OMP END PARALLEL DO
    ! number of basis functions in the two body term
    nbasis2 = size(Rs2)
    ! Inverse of the two body cutoff
    invcut = 1.0d0 / rcut
    ! pre-calculate the radial decay in the two body terms
    rdecay = decay(distance_matrix, invcut, natoms)
    ! Allocate temporary
    allocate(radial(nbasis2))
    
    rep = 0.0d0
    radial = 0.0d0
    
    i = idx + 1
    m = element_types(i)
    ! !$OMP PARALLEL DO PRIVATE(n,rij,radial) REDUCTION(+:rep)
    do j = 1, natoms         
        ! index of the element of atom j
        n = element_types(j)
        ! distance between atoms i and j
        rij = distance_matrix(i,j)
        if ((rij <= rcut) .AND. (rij >= 0.1)) then
                ! two body term of the representation
                mu    = log(rij / sqrt(1.0d0 + eta2  / rij**2))
                sigma = sqrt(log(1.0d0 + eta2  / rij**2))
                radial(:) = 0.0d0
                do k = 1, nbasis2 
                   radial(k) = 1.0d0/(sigma* sqrt(2.0d0*pi) * Rs2(k)) * rdecay(i,j) &
                              & * exp( - (log(Rs2(k)) - mu)**2 / (2.0d0 * sigma**2) ) / rij**two_body_decay
                enddo
                rep((n-1)*nbasis2 + 1:n*nbasis2) = rep((n-1)*nbasis2 + 1:n*nbasis2) + radial
            endif
    enddo
    ! !$OMP END PARALLEL DO

    deallocate(element_types)
    deallocate(rdecay)
    deallocate(distance_matrix)
    deallocate(radial)
end subroutine fgenerate_local_two_body_fchl_acsf

subroutine fgenerate_local_gp_fchl(idx, coordinates, nuclear_charges, periods, groups, ugroups, &
                          & Rs2, Rs3, Ts, eta2, eta3, zeta, rcut, acut, natoms, rep_size, &
                          & two_body_decay, three_body_decay, three_body_weight, rep)

    use representations, only: decay, calc_angle, calc_cos_angle

    implicit none
    integer, intent(in) :: idx
    double precision, intent(in), dimension(:, :) :: coordinates
    integer, intent(in), dimension(:) :: nuclear_charges
    integer, intent(in), dimension(:) :: periods
    integer, intent(in), dimension(:) :: groups
    integer, intent(in), dimension(:) :: ugroups
    
    double precision, intent(in), dimension(:) :: Rs2
    double precision, intent(in), dimension(:) :: Rs3
    double precision, intent(in), dimension(:) :: Ts
    double precision, intent(in) :: eta2
    double precision, intent(in) :: eta3
    double precision, intent(in) :: zeta
    double precision, intent(in) :: rcut
    double precision, intent(in) :: acut
    integer, intent(in) :: natoms
    integer, intent(in) :: rep_size
    double precision, intent(in) :: two_body_decay
    double precision, intent(in) :: three_body_decay
    double precision, intent(in) :: three_body_weight    
    double precision, intent(out), dimension(rep_size) :: rep
    
    integer :: i, j, k, l, n, m, o, p, q, s, z, nelements, nbasis2, nbasis3, nabasis
    integer, allocatable, dimension(:) :: element_types
    double precision :: rij, rik, angle, cos_1, cos_2, cos_3, invcut
    ! double precision :: angle_1, angle_2, angle_3
    double precision, allocatable, dimension(:) :: radial, angular, a, b, c
    double precision, allocatable, dimension(:, :) :: distance_matrix, rdecay
    double precision, allocatable, dimension(:) :: rep3
    double precision :: mu, sigma, ksi3
    
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)


    if (natoms /= size(nuclear_charges, dim=1)) then
        write(*,*) "ERROR: Atom Centered Symmetry Functions creation"
        write(*,*) natoms, "coordinates, but", &
            & size(nuclear_charges, dim=1), "atom_types!"
        stop
    endif


    ! number of element types
    nelements = size(ugroups)
    ! Allocate temporary
    allocate(element_types(natoms))

    ! Store element index of every atom
    ! !$OMP PARALLEL DO
    do i = 1, natoms
        do j = 1, nelements
            if (groups(i) .eq. ugroups(j)) then
                element_types(i) = j
                continue
            endif
        enddo
    enddo
    ! !$OMP END PARALLEL DO


    ! Get distance matrix
    ! Allocate temporary
    allocate(distance_matrix(natoms, natoms))
    distance_matrix = 0.0d0


    !  !$OMP PARALLEL DO PRIVATE(rij)
    do i = 1, natoms
        do j = i+1, natoms
            rij = norm2(coordinates(j,:) - coordinates(i,:))
            distance_matrix(i, j) = rij
            distance_matrix(j, i) = rij
        enddo
    enddo
    ! !$OMP END PARALLEL DO

    ! number of basis functions in the two body term
    nbasis2 = size(Rs2)

    ! Inverse of the two body cutoff
    invcut = 1.0d0 / rcut

    ! pre-calculate the radial decay in the two body terms
    rdecay = decay(distance_matrix, invcut, natoms)

    ! Allocate temporary
    allocate(radial(nbasis2))
    
    rep = 0.0d0
    radial = 0.0d0
    ! 
    i = idx +1
    ! index of the element of atom i
    m = element_types(i)
    ! !$OMP PARALLEL DO PRIVATE(n,rij,radial) REDUCTION(+:rep)
    do j = 1, natoms
        ! index of the element of atom j
        n = element_types(j)
        ! distance between atoms i and j
        rij = distance_matrix(i,j)
        if ((rij <= rcut) .AND. (rij >= 0.1)) then
            ! two body term of the representation
            mu    = log(rij / sqrt(1.0d0 + eta2  / rij**2))
            sigma = sqrt(log(1.0d0 + eta2  / rij**2))
            radial(:) = 0.0d0
            do k = 1, nbasis2 
               radial(k) = periods(j)*1.0d0/(sigma* sqrt(2.0d0*pi) * Rs2(k))&
                          & * rdecay(i,j)  * exp( - (log(Rs2(k)) - mu)**2 / (2.0d0 * sigma**2) ) &
                          & / rij**two_body_decay
            enddo

            rep((n-1)*nbasis2 + 1:n*nbasis2) = rep((n-1)*nbasis2 + 1:n*nbasis2) + radial
        endif
    enddo
    ! !$OMP END PARALLEL DO

    deallocate(radial)

    ! number of radial basis functions in the three body term
    nbasis3 = size(Rs3)
    ! number of radial basis functions in the three body term
    nabasis = size(Ts)

    ! Inverse of the three body cutoff
    invcut = 1.0d0 / acut
    ! pre-calculate the radial decay in the three body terms
    rdecay = decay(distance_matrix, invcut, natoms)

    ! Allocate temporary
    allocate(rep3(rep_size))
    allocate(a(3))
    allocate(b(3))
    allocate(c(3))
    allocate(radial(nbasis3))
    allocate(angular(nabasis))

    rep3 = 0.0d0

    ! This could probably be done more efficiently if it's a bottleneck
    ! Also the order is a bit wobbly compared to the tensorflow implementation
    ! !$OMP PARALLEL DO PRIVATE(rij, n, rik, m, a, b, c, angle, radial, angular, &
    ! !$OMP cos_1, cos_2, cos_3, mu, sigma, o, ksi3, &
    ! !$OMP p, q, s, z) REDUCTION(+:rep3) COLLAPSE(2) SCHEDULE(dynamic)
    i = idx + 1
    do j = 1, natoms - 1
        if (i .eq. j) cycle
        ! distance between atoms i and j
        rij = distance_matrix(i,j)
        if (rij > acut)  cycle
        ! index of the element of atom j
        n = element_types(j)
        do k = j + 1, natoms
            if (i .eq. k) cycle
            if (j .eq. k) cycle
            ! distance between atoms i and k
            rik = distance_matrix(i,k)
            if (rik > acut) cycle
            ! index of the element of atom k
            m = element_types(k)
            ! coordinates of atoms j, i, k
            a = coordinates(j,:)
            b = coordinates(i,:)
            c = coordinates(k,:)
            ! angle between atoms i, j and k centered on i
            angle   = calc_angle(a,b,c)
            cos_1 = calc_cos_angle(a,b,c)
            cos_2 = calc_cos_angle(a,c,b)
            cos_3 = calc_cos_angle(b,a,c)
            ! The radial part of the three body terms including decay
            radial = exp(-eta3*(0.5d0 * (rij+rik) - Rs3)**2) * rdecay(i,j) * rdecay(i,k)
           
            ksi3 = (1.0d0 + 3.0d0 * cos_1 * cos_2 * cos_3) &
                 & / (distance_matrix(i,k) * distance_matrix(i,j) * distance_matrix(j,k) &
             & )**three_body_decay * three_body_weight
            angular = 0.0d0 
            do l = 1, nabasis/2
                o = l*2-1
                angular(2*l-1) = angular(2*l-1) + 2*cos(o * angle) &
                    & * exp(-(zeta * o)**2 /2)
                
                angular(2*l) = angular(2*l) + 2*sin(o * angle) &
                    & * exp(-(zeta * o)**2 /2)
            enddo
            
            ! The lowest of the element indices for atoms j and k
            p = min(n,m) - 1
            ! The highest of the element indices for atoms j and k
            q = max(n,m) - 1
            ! calculate the indices that the three body terms should be added to
            s = nelements * nbasis2 + nbasis3 * nabasis * (-(p * (p + 1))/2 + q + nelements * p) + 1
            do l = 1, nbasis3
                ! calculate the indices that the three body terms should be added to
                z = s + (l-1) * nabasis
                ! Add the contributions from atoms i,j and k
                rep3(z:z + nabasis - 1) = rep3(z:z + nabasis - 1) + periods(j) &
                                             & * periods(k) * angular * radial(l) * ksi3
            enddo
        enddo
    enddo
    ! !$OMP END PARALLEL DO
    rep = rep + rep3
    deallocate(element_types)
    deallocate(rdecay)
    deallocate(distance_matrix)
    deallocate(rep3)
    deallocate(a)
    deallocate(b)
    deallocate(c)
    deallocate(radial)
    deallocate(angular)
end subroutine fgenerate_local_gp_fchl

subroutine fgenerate_local_two_body_gp_fchl(idx, coordinates, nuclear_charges, periods, groups, ugroups, &
                          & Rs2, eta2, rcut, natoms, rep_size, two_body_decay, rep)
    use representations, only: decay, calc_angle, calc_cos_angle
    implicit none
    integer, intent(in) :: idx
    double precision, intent(in), dimension(:, :) :: coordinates
    integer, intent(in), dimension(:) :: nuclear_charges
    integer, intent(in), dimension(:) :: periods
    integer, intent(in), dimension(:) :: groups
    integer, intent(in), dimension(:) :: ugroups
    double precision, intent(in), dimension(:) :: Rs2
    double precision, intent(in) :: eta2
    double precision, intent(in) :: rcut
    integer, intent(in) :: natoms
    integer, intent(in) :: rep_size
    double precision, intent(in) :: two_body_decay  
    double precision, intent(out), dimension(rep_size) :: rep
    
    integer :: i, j, k, n, m, nelements, nbasis2
    integer, allocatable, dimension(:) :: element_types
    double precision :: rij,  invcut
    ! double precision :: angle_1, angle_2, angle_3
    double precision, allocatable, dimension(:) :: radial
    double precision, allocatable, dimension(:, :) :: distance_matrix, rdecay
    double precision :: mu, sigma
    
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    if (natoms /= size(nuclear_charges, dim=1)) then
        write(*,*) "ERROR: Atom Centered Symmetry Functions creation"
        write(*,*) natoms, "coordinates, but", &
            & size(nuclear_charges, dim=1), "atom_types!"
        stop
    endif

    ! number of element types
    nelements = size(ugroups)
    ! Allocate temporary
    allocate(element_types(natoms))

    ! Store element index of every atom
    ! !$OMP PARALLEL DO
    do i = 1, natoms
        do j = 1, nelements
            if (groups(i) .eq. ugroups(j)) then
                element_types(i) = j
                continue
            endif
        enddo
    enddo
    ! !$OMP END PARALLEL DO


    ! Get distance matrix
    ! Allocate temporary
    allocate(distance_matrix(natoms, natoms))
    distance_matrix = 0.0d0


    !  !$OMP PARALLEL DO PRIVATE(rij)
    do i = 1, natoms
        do j = i+1, natoms
            rij = norm2(coordinates(j,:) - coordinates(i,:))
            distance_matrix(i, j) = rij
            distance_matrix(j, i) = rij
        enddo
    enddo
    ! !$OMP END PARALLEL DO

    ! number of basis functions in the two body term
    nbasis2 = size(Rs2)

    ! Inverse of the two body cutoff
    invcut = 1.0d0 / rcut

    ! pre-calculate the radial decay in the two body terms
    rdecay = decay(distance_matrix, invcut, natoms)

    ! Allocate temporary
    allocate(radial(nbasis2))
    
    rep = 0.0d0
    radial = 0.0d0
    ! 
    i = idx + 1
    ! index of the element of atom i
    m = element_types(i)
    ! !$OMP PARALLEL DO PRIVATE(n,rij,radial) REDUCTION(+:rep)
    do j = 1, natoms
        ! index of the element of atom j
        n = element_types(j)
        ! distance between atoms i and j
        rij = distance_matrix(i,j)
        if ((rij <= rcut) .AND. (rij >= 0.1)) then
            ! two body term of the representation
            mu    = log(rij / sqrt(1.0d0 + eta2  / rij**2))
            sigma = sqrt(log(1.0d0 + eta2  / rij**2))
            radial(:) = 0.0d0
            do k = 1, nbasis2 
               radial(k) = periods(j)*1.0d0/(sigma* sqrt(2.0d0*pi) * Rs2(k))&
                          & * rdecay(i,j)  * exp( - (log(Rs2(k)) - mu)**2 / (2.0d0 * sigma**2) ) &
                          & / rij**two_body_decay
            enddo

            rep((n-1)*nbasis2 + 1:n*nbasis2) = rep((n-1)*nbasis2 + 1:n*nbasis2) + radial
        endif
    enddo
    ! !$OMP END PARALLEL DO

    deallocate(element_types)
    deallocate(rdecay)
    deallocate(distance_matrix)
    deallocate(radial)
    
end subroutine fgenerate_local_two_body_gp_fchl

