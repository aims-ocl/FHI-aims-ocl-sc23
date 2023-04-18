!****h* FHI-aims/hartree_potential_recip
!  NAME
!   hartree_potential_recip
!  SYNOPSIS

module hartree_potential_recip

  !  PURPOSE
  !   This module contains the Fourier part of the far-distance Hartree potential routines.
  !   The method and the equations are taken from:
  !   Bernard Delley, J. Phys. Chem 1996, 100, 6107-6110.
  !
  !  USES
  implicit none
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE

  private

  integer :: n_k_points_hartree                           ! Total number of G vectors for reciprocal space Ewald sum
  integer, dimension(:,:), allocatable :: k_points_hartree  ! List of G-vectors
  logical :: not_allocated                                ! Is the G-vector list already initialized ?
  integer, dimension(3) :: k_points_max                   ! Value of largest G-vector index in recipr. direction 1, 2, 3
  complex*16, dimension(:), allocatable :: hartree_coef     ! Parts of the potential calculations which do not
  ! depend on place
  complex*16, dimension(:,:,:), allocatable :: hartree_c_grad   ! Parts of the gradient  calculations which do not
  ! depend on place

  integer, allocatable, dimension(:,:) :: number_z ! counts valid G_z points for a given G_x,G_y combination

  !VB: real*8,private:: pot_zero_level                                ! Shift of potential level.
  !    pot_zero_level is now commented out and no longer used, as it is fully superseded by the
  !    present treatment of charged periodic systems, using the average electrostatic potential
  !    as a reference level. If you want to restore pot_zero_level, let me know.

  public :: update_hartree_potential_recip_v2
  public :: update_hartree_potential_recip
  public :: update_hartree_gradient_recip
  public :: evaluate_hartree_recip_coef
  public :: evaluate_hartree_recip_gradient_coef_atom
  public :: evaluate_potential_at_vacuum
  public :: evaluate_dipole_correction
  public :: deallocate_hartree_p_recip
  public :: evaluate_dipole_correction_from_dipole
  public :: initialize_recip_hartree_potential

  private :: update_hartree_gradient_recip_test





contains






  !****s* hartree_potential_recip/update_hartree_potential_recip_v2
  !  NAME
  !    update_hartree_potential_recip_v2
  !  SYNOPSIS

  subroutine update_hartree_potential_recip_v2( coord_current, potential, sigma_lm )

    ! PURPOSE
    ! This routine calculates the Fourier part of the compensating Ewald potential
    ! of the periodic long-range Hartree potential at the point coord_current.
    ! The potential is added to the variable 'potential'.
    ! The routine uses coefficients calculated by  'evaluate_hartree_recip_coef', which
    ! should be called once after every electron density update.
    ! update_hartree_potential_recip only does the part of the work that depends
    ! on coord_current explicitly. - Paula Havu
    !
    ! _v2 is a new (faster) version that is kept in order to introduce a switch for testing,
    ! for now (Feb 2014).
    !
    ! _v2 will become the default if no problems arise.
    !
    use analytical_stress
    use geometry, only: recip_lattice_vector
    use runtime_choices, only: use_analytical_stress, potential_well_requested,&
        potential_well_start, potential_well_end, potential_well_depth
    implicit none
    !  ARGUMENTS

    real*8:: coord_current(3)
    real*8:: potential
    real*8,dimension(1:3,1:3),optional  :: sigma_lm

    !  INPUTS
    !    o coord_current -- coordinates where potential is calculated
    !
    !  OUTPUT
    !    o potential -- potential is added here.
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    real*8:: proj1, proj2, proj3

    complex*16:: dpot, h_pot
    complex*16:: W1, W2, W3

    integer:: i_atom
    integer:: n, nx, ny, nz

    complex*16,dimension(-k_points_max(1):k_points_max(1)):: Wp1
    complex*16,dimension(-k_points_max(2):k_points_max(2)):: Wp2
    complex*16 :: Wp3
    complex*16 :: Wp12

    real*8, dimension(-k_points_max(1):k_points_max(1),-k_points_max(2):k_points_max(2) ) :: ReWp12, ImWp12
    real*8, dimension(-k_points_max(3):k_points_max(3)) :: ReWp3, ImWp3

    real*8 :: ReFactor
    real*8 :: ImFactor

    real*8 :: ReFactors(1:3,1:3)
    real*8 :: ImFactors(1:3,1:3)

    real*8,external:: ddot

    real*8:: kx,ky,kz

    integer :: i_z

    ! Stress tensor
    integer :: i_l, i_m

    ! What we have to do, is the inverse Fourier transform.
    ! The equation is written in the exponential form. W**(n), where n is
    ! integer.


    proj1 = ddot(3,recip_lattice_vector(1:3,1),1, coord_current(1:3),1)
    proj2 = ddot(3,recip_lattice_vector(1:3,2),1, coord_current(1:3),1)
    proj3 = ddot(3,recip_lattice_vector(1:3,3),1, coord_current(1:3),1)

    W1= exp((0,1)*proj1)
    W2= exp((0,1)*proj2)
    W3= exp((0,1)*proj3)

    ! Different exponents W**n are calculated once and saved to the table.

    Wp1(-k_points_max(1)) = W1**(-k_points_max(1))

    do n=1-k_points_max(1),k_points_max(1)
       Wp1(n) = Wp1(n-1)*W1
    end do

    Wp2(-k_points_max(2)) = W2**(-k_points_max(2))

    do n=1-k_points_max(2),k_points_max(2)
       Wp2(n) = Wp2(n-1)*W2
    end do

    Wp3 = W3**(-k_points_max(3))
    ReWp3(-k_points_max(3)) = dble(Wp3)
    ImWp3(-k_points_max(3)) = dimag(Wp3)

    do n=1-k_points_max(3),k_points_max(3)
       Wp3 = Wp3*W3
       ReWp3(n)=dble(Wp3)
       ImWp3(n)=dimag(Wp3)
    end do

    ! The potential is calculated in this variable.
    dpot = 0

    ! CC: B/L matrices
    if ( (use_analytical_stress) .and. present(sigma_lm) ) then
      sigma_lm(:,:)       = 0.0d0
    end if


    if ( (use_analytical_stress) .and. present(sigma_lm) ) then
      ! Go through k-points list and construct the actual potential.
      n=1
      do while (n.le.n_k_points_hartree)

         nx = k_points_hartree(n, 1)
         ny = k_points_hartree(n, 2)

         Wp12 =  Wp1(nx)* Wp2(ny)

         Refactor = 0.d0
         Imfactor = 0.d0
         Refactors(:,:) = 0.d0
         Imfactors(:,:) = 0.d0
         do i_z = 1, number_z(nx,ny), 1
            nz = k_points_hartree(n, 3)
            ReFactor = ReFactor + dble (hartree_coef(n))*ReWp3(nz) &
                                - dimag(hartree_coef(n))*ImWp3(nz)
            ImFactor = ImFactor - dble (hartree_coef(n))*ImWp3(nz) &
                                - dimag(hartree_coef(n))*ReWp3(nz)

            do i_l=1,3
               do i_m = 1,3
                  ReFactors(i_l,i_m) = ReFactors(i_l,i_m) + dble (AS_hartree_coeff_stress(i_l,i_m,n))*ReWp3(nz) &
                                - dimag(AS_hartree_coeff_stress(i_l,i_m,n))*ImWp3(nz)
                  ImFactors(i_l,i_m) = ImFactors(i_l,i_m) - dble (AS_hartree_coeff_stress(i_l,i_m,n))*ImWp3(nz) &
                                - dimag(AS_hartree_coeff_stress(i_l,i_m,n))*ReWp3(nz)
               enddo
            enddo

            n = n+1
         enddo

         dpot = dpot + ReFactor*dble(Wp12) + ImFactor*dimag(Wp12)

         sigma_lm(:,:) = sigma_lm(:,:) + ReFactors(:,:)*dble(Wp12) + ImFactors(:,:)*dimag(Wp12)
      enddo

    else
      ! Go through k-points list and construct the actual potential.

      ! The following loop structure looks complicated but essentially uses
      ! the fact that the potential is real to take a factor out of a bracket
      ! and out of the innermost loop.

      ! The original, simpler look structure is further below.

      ! It is quite possible that this is still a factor 2 too much work
      ! because either the imaginary part is zero or the sum over the imaginary
      ! parts cancels. I cannot prove this momentarily, however. (may only hold
      ! for structures of certain symmetry)

      n=1
      do while (n.le.n_k_points_hartree)

         nx = k_points_hartree(n, 1)
         ny = k_points_hartree(n, 2)

         Wp12 =  Wp1(nx)* Wp2(ny)

         Refactor = 0.d0
         Imfactor = 0.d0
         do i_z = 1, number_z(nx,ny), 1
            nz = k_points_hartree(n, 3)
            ReFactor = ReFactor + dble (hartree_coef(n))*ReWp3(nz) &
                                - dimag(hartree_coef(n))*ImWp3(nz)
            ImFactor = ImFactor - dble (hartree_coef(n))*ImWp3(nz) &
                                - dimag(hartree_coef(n))*ReWp3(nz)
            n = n+1
         enddo

         dpot = dpot + ReFactor*dble(Wp12) + ImFactor*dimag(Wp12)

      enddo

       ! The following, commented-out block is what the loop above does.

!      do n = 1, n_k_points_hartree
!
!
!         nx = k_points_hartree(n, 1)
!         ny = k_points_hartree(n, 2)
!         nz = k_points_hartree(n, 3)
!
!!         dpot = dpot +  dble(hartree_coef(n) *Wp12(nx,ny) * Wp3(nz))
!
!      end do
    end if

    potential = potential +  dpot


    if (potential_well_requested) then
       if ((coord_current(3).gt.potential_well_start) .and. (coord_current(3).lt. potential_well_end)) then
          potential= potential+potential_well_depth
       endif
    endif


    ! If the total charge of the system is non-zero we have to shift the potential level
    ! by 'pot_zero_level' so that average potential is zero. If the system is charge neutral
    ! then 'pot_zero_level' should be zero.
    ! VB: Now commented out potential = potential + pot_zero_level

  end subroutine update_hartree_potential_recip_v2
  !******







  !****s* hartree_potential_recip/update_hartree_potential_recip
  !  NAME
  !    update_hartree_potential_recip
  !  SYNOPSIS

  subroutine update_hartree_potential_recip( coord_current, potential, sigma_lm )

    ! PURPOSE
    ! This routine calculates the Fourier part of the compensating Ewald potential
    ! of the periodic long-range Hartree potential at the point coord_current.
    ! The potential is added to the variable 'potential'.
    ! The routine uses coefficients calculated by  'evaluate_hartree_recip_coef', which
    ! should be called once after every electron density update.
    ! update_hartree_potential_recip only does the part of the work that depends
    ! on coord_current explicitly. - Paula Havu
    !
    ! Old (slower) version that is kept in order to introduce a switch for testing,
    ! for now (Feb 2014). The "_v2" version above is a faster version that _should_
    ! do exactly the same thing, only focusing just on the real part of the potential.
    !
    ! _v2 will become the default if no problems arise.
    !
    use analytical_stress
    use geometry, only: recip_lattice_vector
    use runtime_choices, only: use_analytical_stress, potential_well_requested,&
        potential_well_start, potential_well_end, potential_well_depth
    implicit none
    !  ARGUMENTS

    real*8:: coord_current(3)
    real*8:: potential
    real*8,dimension(1:3,1:3),optional  :: sigma_lm

    !  INPUTS
    !    o coord_current -- coordinates where potential is calculated
    !
    !  OUTPUT
    !    o potential -- potential is added here.
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    real*8:: proj1, proj2, proj3

    complex*16:: dpot, h_pot
    complex*16:: W1, W2, W3

    integer:: i_atom
    integer:: n, nx, ny, nz

    complex*16,dimension(-k_points_max(1):k_points_max(1)):: Wp1
    complex*16,dimension(-k_points_max(2):k_points_max(2)):: Wp2
    complex*16,dimension(-k_points_max(3):k_points_max(3)):: Wp3
    complex*16,dimension(-k_points_max(1):k_points_max(1),-k_points_max(2):k_points_max(2) ):: Wp12

    real*8,external:: ddot

    real*8:: kx,ky,kz

    ! Stress tensor
    integer :: i_l, i_m

    ! What we have to do, is the inverse Fourier transform.
    ! The equation is written in the exponential form. W**(n), where n is
    ! integer.


    proj1 = ddot(3,recip_lattice_vector(1:3,1),1, coord_current(1:3),1)
    proj2 = ddot(3,recip_lattice_vector(1:3,2),1, coord_current(1:3),1)
    proj3 = ddot(3,recip_lattice_vector(1:3,3),1, coord_current(1:3),1)

    W1= exp((0,1)*proj1)
    W2= exp((0,1)*proj2)
    W3= exp((0,1)*proj3)

    ! Different exponents W**n are calculated once and saved to the table.

    Wp1(-k_points_max(1)) = W1**(-k_points_max(1))

    do n=1-k_points_max(1),k_points_max(1)
       Wp1(n) = Wp1(n-1)*W1
    end do

    Wp2(-k_points_max(2)) = W2**(-k_points_max(2))

    do n=1-k_points_max(2),k_points_max(2)
       Wp2(n) = Wp2(n-1)*W2
    end do

    Wp3(-k_points_max(3)) = W3**(-k_points_max(3))

    do n=1-k_points_max(3),k_points_max(3)
       Wp3(n) = Wp3(n-1)*W3
    end do

    do nx=-k_points_max(1),k_points_max(1)
       do ny=-k_points_max(2),k_points_max(2)
          Wp12(nx,ny) =  Wp1(nx)* Wp2(ny)
       end do
    end do

    ! The potential is calculated in this variable.
    dpot = 0

    ! CC: B/L matrices
    if ( (use_analytical_stress) .and. present(sigma_lm) ) then
      sigma_lm(:,:)       = 0.0d0
    end if


    if ( (use_analytical_stress) .and. present(sigma_lm) ) then
      ! Go through k-points list and construct the actual potential.
      do n = 1, n_k_points_hartree


         nx = k_points_hartree(n, 1)
         ny = k_points_hartree(n, 2)
         nz = k_points_hartree(n, 3)


         dpot = dpot +  dble(hartree_coef(n) *Wp12(nx,ny) * Wp3(nz))

         !! CC: Sum over all indices
         do i_m=1,3
           do i_l=1,3
             sigma_lm(i_l,i_m) = sigma_lm(i_l,i_m) + dble( Wp12(nx,ny) * Wp3(nz) * AS_hartree_coeff_stress(i_l,i_m,n))
             !CC: Slower explicit Code : sigma_lm(i_l,i_m) = sigma_lm(i_l,i_m) + dble( Wp12(nx,ny) * Wp3(nz) * ( (hartree_coef(n) * AS_B_matrix(i_l,i_m,n)) - AS_Lambda_matrix(i_l,i_m,n) ))
           end do
         end do

      end do
    else
      ! Go through k-points list and construct the actual potential.
      do n = 1, n_k_points_hartree


         nx = k_points_hartree(n, 1)
         ny = k_points_hartree(n, 2)
         nz = k_points_hartree(n, 3)

         dpot = dpot +  dble(hartree_coef(n) *Wp12(nx,ny) * Wp3(nz))

!         dpot = dpot +  dble(hartree_coef(n) *Wp1(nx) * Wp2(ny) * Wp3(nz))

         !       write(use_unit,*) hartree_coef(n)
         !       write(use_unit,*)  Wp3(nz)
         !       write(use_unit,*)  ny,ny,nz

      end do
    end if

    potential = potential +  dpot

    if (potential_well_requested) then
       if ((coord_current(3).gt.potential_well_start) .and. (coord_current(3).lt. potential_well_end)) then
          potential= potential+potential_well_depth
       endif
    endif

    ! If the total charge of the system is non-zero we have to shift the potential level
    ! by 'pot_zero_level' so that average potential is zero. If the system is charge neutral
    ! then 'pot_zero_level' should be zero.
    ! VB: Now commented out potential = potential + pot_zero_level

  end subroutine update_hartree_potential_recip
  !******





  !****s* hartree_potential_recip/update_hartree_gradient_recip_test
  !  NAME
  !   update_hartree_gradient_recip_test
  !  SYNOPSIS

  subroutine update_hartree_gradient_recip_test( coord_current, gradient)

    !  PURPOSE
    !
    ! This routine calculates the Fourier part of the long-range Hartree potential gradient
    ! at point coord_current.
    ! The routine uses coefficients calculated by  'evaluate_hartree_recip_coef', which
    ! should be called once after every electron density update. The present routine only
    ! evaluates components that depend on the currect grid point.
    !
    ! Note: this is only a test routine for debugging reasons. The real one is below.
    !
    use geometry, only: recip_lattice_vector
    implicit none
    !  ARGUMENTS

    real*8:: coord_current(3)
    real*8:: gradient(3)

    !  INPUTS
    !   o coord_current -- coordinates
    !  OUTPUT
    !   o gradient -- gradient of the potential is added here
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    real*8:: kx,ky,kz
    real*8:: proj1, proj2, proj3

    complex*16:: dgrad(3), h_pot
    complex*16:: W1, W2, W3

    integer:: i_atom
    integer:: n, nx, ny, nz

    complex*16,dimension(-k_points_max(1):k_points_max(1)):: Wp1
    complex*16,dimension(-k_points_max(2):k_points_max(2)):: Wp2
    complex*16,dimension(-k_points_max(3):k_points_max(3)):: Wp3



    ! What we have to do, is the inverse Fourier transform.
    ! The equation is written to the exponential form. W**(n), where n is
    ! integer.

    proj1 = dot_product(recip_lattice_vector(1:3,1), coord_current(1:3))
    proj2 = dot_product(recip_lattice_vector(1:3,2), coord_current(1:3))
    proj3 = dot_product(recip_lattice_vector(1:3,3), coord_current(1:3))

    W1= exp((0,1)*proj1)
    W2= exp((0,1)*proj2)
    W3= exp((0,1)*proj3)

    ! Different exponents W**n are calculated once and saved to the table.

    Wp1(-k_points_max(1)) = W1**(-k_points_max(1))

    do n=1-k_points_max(1),k_points_max(1)
       Wp1(n) = Wp1(n-1)*W1
    end do

    Wp2(-k_points_max(2)) = W2**(-k_points_max(2))

    do n=1-k_points_max(2),k_points_max(2)
       Wp2(n) = Wp2(n-1)*W2
    end do

    Wp3(-k_points_max(3)) = W3**(-k_points_max(3))

    do n=1-k_points_max(3),k_points_max(3)
       Wp3(n) = Wp3(n-1)*W3
    end do



    ! The gradient is calculated to this variable.
    dgrad = 0


    ! Go through k-points list and construct the actual potential.

    do n = 1, n_k_points_hartree


       nx = k_points_hartree(n, 1)
       ny = k_points_hartree(n, 2)
       nz = k_points_hartree(n, 3)



       kx = nx * recip_lattice_vector(1,1) + ny * recip_lattice_vector(1,2) +  nz * recip_lattice_vector(1,3)
       ky = nx * recip_lattice_vector(2,1) + ny * recip_lattice_vector(2,2) +  nz * recip_lattice_vector(2,3)
       kz = nx * recip_lattice_vector(3,1) + ny * recip_lattice_vector(3,2) +  nz * recip_lattice_vector(3,3)


       dgrad(1) = dgrad(1) +  dble( (kx*(0,-1))* hartree_coef(n) *Wp1(nx) * Wp2(ny) * Wp3(nz))
       dgrad(2) = dgrad(2) +  dble( (ky*(0,-1))* hartree_coef(n) *Wp1(nx) * Wp2(ny) * Wp3(nz))
       dgrad(3) = dgrad(3) +  dble( (kz*(0,-1))* hartree_coef(n) *Wp1(nx) * Wp2(ny) * Wp3(nz))

       !       write(use_unit,*) hartree_coef(n)
       !       write(use_unit,*)  Wp3(nz)
       !       write(use_unit,*)  ny,ny,nz



    end do

    gradient = gradient + dgrad


  end subroutine update_hartree_gradient_recip_test
  !******







  !****s* hartree_potential_recip/update_hartree_gradient_recip
  !  NAME
  !    update_hartree_gradient_recip
  !  SYNOPSIS

  subroutine update_hartree_gradient_recip( coord_current, gradient, i_atom)

    !  PURPOSE
    ! This routine calculates the Fourier part of the long-range Hartree potential gradient
    ! at point coord_current.
    ! The routine uses coefficients calculated by  'evaluate_hartree_recip_coef', which
    ! should be called once after every electron density update. The present routine only
    ! evaluates components that depend on the currect grid point.
    !
    ! Note that before calling this one, the coefficients of Fourier transform have to be calculated:
    ! call  evaluate_hartree_recip_gradient_coef_atom
    !
    use geometry, only: recip_lattice_vector, length_of_lattice_vector
    implicit none
    !  ARGUMENTS

    real*8:: coord_current(3)
    real*8:: gradient(3)
    integer:: i_atom

    !  INPUTS
    !    o coord_current -- coordinates
    !    o i_atom -- force of this atom is calculated
    !
    !  OUTPUT
    !    o gradient -- the potential gradient term is added here
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    real*8:: kx,ky,kz
    real*8:: proj1, proj2, proj3

    complex*16:: dpot, h_pot
    complex*16:: W1, W2, W3

    integer:: i_center, i_coord
    integer:: n, nx, ny, nz

    complex*16,dimension(-k_points_max(1):k_points_max(1)):: Wp1
    complex*16,dimension(-k_points_max(2):k_points_max(2)):: Wp2
    complex*16,dimension(-k_points_max(3):k_points_max(3)):: Wp3



    ! What we have to do, is the inverse Fourier transform.
    ! The equation is writen to the expodential form. W**(n), where n is
    ! integer.

    proj1 = dot_product(recip_lattice_vector(1:3,1), coord_current(1:3))
    proj2 = dot_product(recip_lattice_vector(1:3,2), coord_current(1:3))
    proj3 = dot_product(recip_lattice_vector(1:3,3), coord_current(1:3))

    W1= exp((0,1)*proj1)
    W2= exp((0,1)*proj2)
    W3= exp((0,1)*proj3)

    ! Different exponents W**n are calculated once and saved to the table.

    Wp1(-k_points_max(1)) = W1**(-k_points_max(1))

    do n=1-k_points_max(1),k_points_max(1)
       Wp1(n) = Wp1(n-1)*W1
    end do

    Wp2(-k_points_max(2)) = W2**(-k_points_max(2))

    do n=1-k_points_max(2),k_points_max(2)
       Wp2(n) = Wp2(n-1)*W2
    end do

    Wp3(-k_points_max(3)) = W3**(-k_points_max(3))

    do n=1-k_points_max(3),k_points_max(3)
       Wp3(n) = Wp3(n-1)*W3
    end do



    ! The potential is calculated to this variable.
    dpot = 0


    ! Go through k-points list and construct the actual potential.

    !    do i_center = 1, n_centers_hartree_potential
    !       i_atom = center_to_atom(centers_hartree_potential(i_center))


    do i_coord = 1,3
       do n = 1, n_k_points_hartree


          nx = k_points_hartree(n, 1)
          ny = k_points_hartree(n, 2)
          nz = k_points_hartree(n, 3)

          dpot = dpot +  dble(hartree_c_grad(i_coord,n,1) *Wp1(nx) * Wp2(ny) * Wp3(nz))

       end do

       gradient(i_coord) = gradient(i_coord)  + dpot
       dpot = 0

    end do
    !    end do

  end subroutine update_hartree_gradient_recip
  !******







  !****s* hartree_potential_recip/initialize_recip_hartree_potential
  !  NAME
  !    initialize_recip_hartree_potential
  !  SYNOPSIS

  subroutine initialize_recip_hartree_potential

    !  PURPOSE
    ! Rewritten subroutine to determine the reciprocal-space vectors which enter the
    ! reciprocal space sum in the Ewald potential. The sum is chosen to reproduce the
    ! a set of smooth Gaussian charge densities, compensating the real-spaced charge multipoles
    ! (charges, dipole moments, quadrupole moments, etc.) situated on each dipole.
    !
    ! The sole important criterion is therefore that this artificial density be reproduced
    ! accurately enough. The criteria that enter are
    ! (1) the width parameter of the compensating charge density, r_0 [hartree_convergence_parameter in control.in]
    ! (2) the maximum angular momentum used in the Hartree potential [hartree_l_max in each species]
    ! (3) a threshold energy parameter below which potential components may be neglected
    !     [hartree_fourier_part_th in control.in]
    !
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    use constants, only: pi
    use dimensions, only: n_atoms
    use geometry, only: recip_lattice_vector, length_of_lattice_vector
    use localorb_io, only: localorb_info, use_unit, OL_norm
    use runtime_choices, only: hartree_potential_fourier_part_threshold, &
        ewald_radius
    use pbc_lists, only: species_center
    use species_data, only: l_hartree
    use synchronize_mpi, only: check_allocation
    implicit none

    ! General purpuse:
    character*100 :: info_str
    integer       :: info

    ! Section 1:
    real*8  :: delta_energy
    integer :: l_max, i_atom
    real*8  :: actual_threshold
    real*8  :: G_max, G_upper, G_lower
    real*8  :: dpot, dpot_upper, dpot_lower
    real*8  :: bifurc_threshold = 1e-10

    ! Section 2, 3 and 4:
    integer, dimension(3) :: nmax
    integer               :: n1, n2, n3
    real*8,  dimension(3) :: k


    ! 1. First, determine maximum momentum (cutoff momentum) in Ewald reciprocal space sum.
    !    The actual equation that we solve is
    !
    !             G_max^(l_max-2) * exp [-(r_0 * G_max/2)^2] > threshold / (20*4*pi)
    !
    !    (for the unusual case l_max < 2, we simply forget about the G_max^(l_max-2) term though!)
    !
    !    This equation has no analytical solution. To derive it, a simple look at Ewald's method
    !    and the parts that enter the reciprocal space part should be sufficient. Note that with our
    !    default "threshold", this is also a pretty conservative estimate.

    ! Rename the threshold parameter (for convenience only)
    delta_energy = hartree_potential_fourier_part_threshold

    ! Find maximal l_max in Hartree potential
    l_max = 0
    do i_atom = 1, n_atoms
       l_max = max(l_max,l_hartree(species_center(i_atom)))
    end do

    ! actual threshold
    actual_threshold = delta_energy / (4.d0*pi*20.d0)

    ! Solve equation by bifurcation
    ! Minimum threshold is 1.0d0 though.
    G_max = 1.d0
    dpot_upper = exp(-Ewald_radius**2 * G_max**2 / 4.)
    if (l_max.gt.2) then
      dpot_upper = G_max**(l_max-2) * dpot_upper
    end if

    ! find lower and upper bounds for G_max
    dpot_lower = 0.d0
    do while ( (dpot_upper.gt.actual_threshold) .or. (dpot_lower.le.dpot_upper) )
      G_lower = G_max
      dpot_lower = dpot_upper

      G_max = G_max + 1.d0
      G_upper = G_max
      dpot_upper = exp(-Ewald_radius**2 * G_max**2 / 4.)
      if (l_max.gt.2) then
        dpot_upper = G_max**(l_max-2) * dpot_upper
      end if
    end do
    ! once we have left the preceding loop, the correct value of G_max is between G_lower and G_upper.
    ! Can now start the bifurcation
    do while ( abs(G_upper-G_lower).gt.bifurc_threshold )
      G_max = G_lower + (G_upper-G_lower)/2.d0

      dpot = exp(-Ewald_radius**2 * G_max**2 / 4.)
      if (l_max.gt.2) then
        dpot = G_max**(l_max-2) * dpot
      end if

      if (dpot.gt.actual_threshold) then
        G_lower = G_max
      else if (dpot.lt.actual_threshold) then
        G_upper = G_max
      else if (dpot.eq.actual_threshold) then
        ! in this case, we're done.
        G_lower = G_max
        G_upper = G_max
      end if
    end do

    write(info_str,'(2X,A,F18.8,A)') '| Estimated reciprocal-space cutoff momentum G_max: ',G_max,' bohr^-1 .'
    call localorb_info( info_str, use_unit, '(A)', OL_norm )



    ! 2. All k-points inside a sphere with radius G_max have to be considered in
    !    the Ewald method (except for the Gamma point). In this section, a
    !    certain parallelepiped is defined that encloses the sphere. Later (in
    !    section 3 and 4), all discrete lattice points within the parallelepiped
    !    will be tested one after the other whether their magnitude is less than
    !    or equal to G_max.
    !      The following derivation for the parallelepiped is mentioned briefly
    !    in Richard Martin, "Electronic Structure: Basic Theory and Practical
    !    Methods", Cambridge Univ. Pr., 2004, Chapter 4.2. Note that the same
    !    method is used in "get_n_supercells" in bravais.f90, but for the case
    !    of points in _real_ space.
    !      First some notation: if the lattice vectors are denoted by a1, a2, a3
    !    and the reciprocal lattice vectors by b1, b2, b3, then the discrete
    !    points k in reciprocal space are given by:
    !
    !             k = n1*b1 + n2*b2 + n3*b3     (1)
    !
    !    Here, n1, n2, n3 are integers. In order to determine those k-points
    !    that lie in the G_max-sphere, a parallelepiped is defined that is
    !    aligned with the reciprocal lattice vectors and that contains the
    !    sphere. The following drawing illustrates the situation in two
    !    dimensions:
    !
    !                                  /\ b2
    !                                  /
    !                                 /M2
    !                     .--ooooooo-/----------.
    !                    /ooo       ooo        /
    !                   oo         /   oo     /
    !                  o          /      o   /
    !                 o          /        o /M1
    !            ----/o---------.__-------o/----> b1
    !               / o        /   \___   o
    !              /   o      /  G_max \_o_
    !             /     oo   /         oo d\___
    !            /        ooo       ooo/       \__\
    !           .----------/-ooooooo--.           / a1
    !                     /
    !                    /
    !
    !    In order to minimize the parallelepiped, its six faces should touch the
    !    sphere. Since the parallelepiped is aligned with the reciprocal lattice
    !    vectors, the two faces that belong to direction i (=1,2,3), are given
    !    by:
    !
    !             k = +- Mi*bi + mj*bj + mk*bk     (2)
    !
    !    Here, Mi is a positive fixed real number, while mj and mk vary across
    !    the face. The two integers j and k refer to the remaining two
    !    directions (for a given direction i). If a face touches the sphere,
    !    there is vector d from the origin to the face that satisfies two
    !    conditions: d is orthogonal to the face and has magnitude G_max. Since
    !    a surface normal of (2) is given by the normalized real space vector
    !    ai/|ai|, we can write d = +- G_max ai/|ai|. The sign is the same as in
    !    (2), because ai and bi enclose an angle of less than pi/2. Inserting d
    !    in (2) yields:
    !
    !             G_max ai/|ai| = Mi*bi + mj*bj + mk*bk
    !
    !    Multiplying with ai/|ai|:
    !
    !             G_max         = Mi*bi*ai/|ai|
    !                           = Mi*2*pi /|ai|
    !    Finally:
    !             Mi = G_max * |ai| / (2*pi)
    !
    !    Concerning the discreteness of the lattice points in (1), it is
    !    sufficient to consider only ni = -nmaxi, ..., nmaxi with nmaxi = |_Mi_|
    !    (floor!), since rounding _up_ would yield a face of the parallelepiped
    !    that does not intersect with the sphere. Finally, we add 10^-10 to
    !    account for limited numerical representation.

    nmax(:)    =    floor(   G_max * length_of_lattice_vector(:) / (2*pi)  +  1E-10   )



    ! 3. Now go through the parallelepiped to identify the k-points within the
    !    G_max-sphere. In the first pass, we merely count, to determine the
    !    appropriate array dimension.

    n_k_points_hartree = 0

    do n1 = -nmax(1), nmax(1)
       do n2 = -nmax(2), nmax(2)
          do n3 = -nmax(3), nmax(3)

             ! Implementing (1) with k(:) = (kx, ky, kz):
             k(:)  =   n1 * recip_lattice_vector(:,1)  &
                     + n2 * recip_lattice_vector(:,2)  &
                     + n3 * recip_lattice_vector(:,3)

             ! Find k-points with magnitude <= G_max, but exclude the Gamma point:
             if (  .not. ( n1 == 0 .and. n2 == 0 .and. n3 == 0 )  &
                   .and. sqrt( sum(k(:)**2) ) <= G_max             ) then
                n_k_points_hartree = n_k_points_hartree + 1
             end if

          end do
       end do
    end do

    write(info_str,'(2X,A,I8)') '| Reciprocal lattice points for long-range Hartree potential:',n_k_points_hartree
    call localorb_info( info_str, use_unit, '(A)', OL_norm )

    allocate( k_points_hartree(n_k_points_hartree,3), stat=info )
    call check_allocation(info, 'k_points_hartree              ')

    allocate( hartree_coef(n_k_points_hartree),       stat=info )
    call check_allocation(info, 'hartree_coef                  ')

    allocate(number_z(-nmax(1):nmax(1),-nmax(2):nmax(2)))
    call check_allocation(info, 'number_z                      ')

    not_allocated = .false.



    ! 4. In the second pass, we store the actually needed k-vectors in a list
    !    (now allocated).

    n_k_points_hartree = 0
    k_points_max(:) = 0

    do n1 = -nmax(1), nmax(1)
       do n2 = -nmax(2), nmax(2)
          number_z(n1,n2)=0
          do n3 = -nmax(3), nmax(3)

             ! Implementing (1) with k(:) = (kx, ky, kz):
             k(:)  =   n1 * recip_lattice_vector(:,1)  &
                     + n2 * recip_lattice_vector(:,2)  &
                     + n3 * recip_lattice_vector(:,3)

             ! Find k-points with magnitude <= G_max, but exclude the Gamma point:
             if (  .not. ( n1 == 0 .and. n2 == 0 .and. n3 == 0 )  &
                   .and.  sqrt( sum(k(:)**2) ) <= G_max             ) then

                n_k_points_hartree = n_k_points_hartree + 1

                k_points_hartree(n_k_points_hartree,:)  =  (/ n1, n2, n3 /)

                ! Fast update routine needs maximum k-points value to the different directions.
                k_points_max(1) = max(k_points_max(1),n1)
                k_points_max(2) = max(k_points_max(2),n2)
                k_points_max(3) = max(k_points_max(3),n3)

                number_z(n1,n2)=number_z(n1,n2)+1

             end if

          end do
       end do
    end do

  end subroutine initialize_recip_hartree_potential
  !******







  !****s* hartree_potential_recip/evaluate_hartree_recip_coef
  !  NAME
  !   evaluate_hartree_recip_coef
  !  SYNOPSIS

  subroutine evaluate_hartree_recip_coef( index_lm, multipole_moments,l_hartree_max_far_distance, flag_stress )

    use analytic_multipole_coefficients, only: n_dd_lm_ijk, index_dd, dd
    use analytical_stress
    use constants
    use dimensions, only: n_atoms, l_pot_max
    use geometry, only: unit_cell_shape_changed, recip_lattice_vector, &
        cell_volume
    use mpi_tasks, only: myid, n_tasks
    use pbc_lists, only: coords_center, center_to_atom
    use runtime_choices, only: ewald_radius
    use synchronize_mpi, only: sync_vector_complex

    !  PURPOSE
    ! This routine calculates all the math of the Fourier part of the Hartree potential
    ! which do not depend on the coordinate point of the potential. Routine is called
    ! once in the iteration at begin on 'Sum of whole potential', before going on every
    ! grid point.
    !
    implicit none
    !  ARGUMENTS

    integer index_lm(-l_pot_max:l_pot_max, 0:l_pot_max)
    real*8, dimension( ( l_pot_max + 1)**2, n_atoms) ::  multipole_moments
    integer:: l_hartree_max_far_distance(n_atoms)
    logical,optional :: flag_stress

    !  INPUTS
    !    o index_lm -- order of l and m components in multipole_moments table
    !    o multipole_moments -- multipole moments
    !    o l_hartree_max_far_distance -- maximum l component of Hartree potential in this part
    !
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    integer :: l_max
    integer :: i_atom
    integer :: i_l
    integer :: i_m
    integer :: ii, jj, kk, nx,ny,nz,i_ml_ijk
    integer :: n

    real*8:: kx,ky,kz,k2

    integer, parameter :: ik8 = selected_int_kind(18) ! kind for integer*8 variables
    integer(kind=ik8) :: n_work, n_work_s, n_work_e, n_atom_s, n_atom_e, n_k_point_s, n_k_point_e

    real*8 :: prefactor
    complex*16 :: cmplx_prefactor

    !! CC: Added for stress tensor
    complex*16                 :: hartree_coef_wo_k
    real*8,dimension(1:3)      :: vec_k
    integer                    :: i_coord,i_coord_2
    complex*16,dimension(1:3)  :: der_k_hartree_coef
    logical                    :: analytical_stress_on



    ! If this routine is called for the first time, or after a geometry step
    ! that changed the unit cell shape, perform initialization of
    ! the G vectors:

    if (unit_cell_shape_changed) then ! This variable is memorized in predict_new_geometry()
       call deallocate_hartree_p_recip()
       unit_cell_shape_changed = .false.
    end if

    if (not_allocated) then
       call initialize_recip_hartree_potential
    end if

    ! CC: Decide if the B/L matrix for the stress need to be computed
    if ( present(flag_stress) ) then
      analytical_stress_on = flag_stress
    else
      analytical_stress_on = .false.
    end if

    ! CC: Allocate stress B/L matrix
    if ( analytical_stress_on ) then
      if (.not.allocated(AS_B_matrix))             allocate(AS_B_matrix(1:3,1:3,1:n_k_points_hartree))
      if (.not.allocated(AS_Lambda_matrix))        allocate(AS_Lambda_matrix(1:3,1:3,1:n_k_points_hartree))
      if (.not.allocated(AS_hartree_coeff_stress)) allocate(AS_hartree_coeff_stress(1:3,1:3,1:n_k_points_hartree))
      AS_B_matrix(:,:,:)             =  0.0d0
      AS_Lambda_matrix(:,:,:)        = (0.0d0, 0.0d0)
      AS_hartree_coeff_stress(:,:,:) = (0.0d0, 0.0d0)
    end if




    ! If no k-points in reciprocal space have to be considered, hartree_coef(:)
    ! is of size zero and does not have to be calculated or even initialized.
    ! We can thus exit this subroutine to avoid division by n_k_points_hartree.

    if ( n_k_points_hartree == 0 ) return


    hartree_coef = (0.d0,0.d0)

    ! For efficient parallelization on massively parallel systems, we must
    ! distribute work among n_atoms AND n_k_points_hartree
    ! For safety we have to use INTEGER*8 arithmetic here!!!

    n_work = int(n_atoms,ik8) * int(n_k_points_hartree,ik8) ! Total number of k-points to work on
    n_work_s = (n_work*myid)/n_tasks+1   ! first k-point (within total number) for my task
    n_work_e = (n_work*(myid+1))/n_tasks ! last k-point  (within total number) for my task

    n_atom_s = (n_work_s-1)/n_k_points_hartree + 1 ! first atom with work for current task
    n_atom_e = (n_work_e-1)/n_k_points_hartree + 1 ! last  atom with work for current task


    ! JW: Should efficiency be of importance for this subroutine, this could
    !     be done along the lines of the charge density update in
    !     [Wieferink, Krueger, Pollmann, PRB 74, 205311 (2006)].
    !     The first step would be to perform the inner loop calculations only
    !     once for each ii,jj,kk triple.  But for a regular grid, things can
    !     be done even more efficient.

    do i_atom = n_atom_s, n_atom_e


       ! get start/end K-point for current atom

       n_k_point_s = max(n_work_s-int(i_atom-1,ik8)*int(n_k_points_hartree,ik8),int(1,ik8))
       n_k_point_e = min(n_work_e-int(i_atom-1,ik8)*int(n_k_points_hartree,ik8),int(n_k_points_hartree,ik8))

       l_max =    l_hartree_max_far_distance(i_atom)

       do n = n_k_point_s, n_k_point_e

          ! k-points were calculated beforehand.

          nx = k_points_hartree(n, 1)
          ny = k_points_hartree(n, 2)
          nz = k_points_hartree(n, 3)

          kx = nx * recip_lattice_vector(1,1) + ny * recip_lattice_vector(1,2) +  nz * recip_lattice_vector(1,3)
          ky = nx * recip_lattice_vector(2,1) + ny * recip_lattice_vector(2,2) +  nz * recip_lattice_vector(2,3)
          kz = nx * recip_lattice_vector(3,1) + ny * recip_lattice_vector(3,2) +  nz * recip_lattice_vector(3,3)

          k2 = kx**2 + ky**2 + kz**2

          prefactor = 4*pi / cell_volume /k2 * exp(-Ewald_radius**2 * k2 / 4.)
          cmplx_prefactor = prefactor * &
             exp ( (0,1)*(kx* ( - coords_center(1, i_atom)) + ky* (  - coords_center(2, i_atom)) +  &
                          kz* ( - coords_center(3, i_atom)) ) )

          ! We go through the list of so called D coefficients.

          do i_ml_ijk = 1, n_dd_lm_ijk

             ! l and m points of the multipole belong for particular d component
             i_l = index_dd(i_ml_ijk, 1)
             i_m = index_dd(i_ml_ijk, 2)

             ! Exponential terms for x, y and z direction  belong for particular d component
             ii  = index_dd(i_ml_ijk, 3)
             jj  = index_dd(i_ml_ijk, 4)
             kk  = index_dd(i_ml_ijk, 5)

             if(i_l <= l_max)then

                ! This is equation (6) in Delleys paper.

                !             hartree_coef(n)  =     hartree_coef(n) &
                !                  + 4*pi / cell_volume /k2 &
                !                  *  multipole_moments(  index_lm(i_m, i_l), center_to_atom(i_atom))  &
                !                  * exp(-Ewald_radius**2 * k2 / 4.) &
                !                  * dd(i_ml_ijk) *(0,1)**i_l  * kx**ii * ky**jj * kz**kk &
                !                  !!
                !                  * exp( (0,1)*( &
                !                  kx* ( - coords_center(1, i_atom)) + ky* (  - coords_center(2, i_atom)) +  &
                !                  kz* ( - coords_center(3, i_atom))))

                !! CC: Original version
                hartree_coef(n)  =     hartree_coef(n) &
                     +  cmplx_prefactor &
                     *  multipole_moments(  index_lm(i_m, i_l), center_to_atom(i_atom))  &
                     * dd(i_ml_ijk) *(0,1)**i_l  * kx**ii * ky**jj * kz**kk

                !! CC: We construct Lambda matrix
                if ( analytical_stress_on ) then

                  !! CC: I split up the orig formulation into hartree_coef_wo_k and k-factor for non-redundant implementation of Lambda_matrix
                  !!     This however leads to small numerical errors in the forces that break the symmetry. FIXME: Check if similar errors arise
                  !!     for the stress
                  hartree_coef_wo_k = 4*pi / cell_volume /k2 &
                     *  multipole_moments(  index_lm(i_m, i_l), center_to_atom(i_atom))  &
                     * dd(i_ml_ijk) *(0,1)**i_l &
                     * exp( (-Ewald_radius**2 * k2 / 4.) + &
                     (0,1)*( &
                     kx* ( - coords_center(1, i_atom)) + ky* (  - coords_center(2, i_atom)) +  &
                     kz* ( - coords_center(3, i_atom))))
                  !! hartree_coef(n) = hartree_coef(n) + &
                  !!      ( hartree_coef_wo_k * kx**ii * ky**jj * kz**kk )

                  !! CC: we basically computed [ d/dk_m kx**ii ky**jj kz**kk ] * hartree_coef_wo_k * k_l
                  der_k_hartree_coef(:) = (0.0d0, 0.0d0)
                  if ( ii .gt. 0 ) then
                    der_k_hartree_coef(1) = ( ii*kx**(ii-1) *    ky**jj     *    kz**kk     ) * hartree_coef_wo_k
                  end if
                  if ( jj .gt. 0 ) then
                    der_k_hartree_coef(2) = (    kx**ii     * jj*ky**(jj-1) *    kz**kk     ) * hartree_coef_wo_k
                  end if
                  if ( kk .gt. 0 ) then
                    der_k_hartree_coef(3) = (    kx**ii     *    ky**jj     * kk*kz**(kk-1) ) * hartree_coef_wo_k
                  end if
                  do i_coord=1,3
                    AS_Lambda_matrix(1,i_coord,n) = AS_Lambda_matrix(1,i_coord,n) + ( kx * der_k_hartree_coef(i_coord) )
                    AS_Lambda_matrix(2,i_coord,n) = AS_Lambda_matrix(2,i_coord,n) + ( ky * der_k_hartree_coef(i_coord) )
                    AS_Lambda_matrix(3,i_coord,n) = AS_Lambda_matrix(3,i_coord,n) + ( kz * der_k_hartree_coef(i_coord) )
                  end do

                end if

             end if
          end do  ! i_ml_ijk
          !             write(use_unit,*)  hartree_coef(n),kx,ky,kz,'------------'
       end do  ! n
    end do  ! i_atom

    !! CC: calculate so called B_matrix
    if ( analytical_stress_on ) then

      do n = 1,n_k_points_hartree

         ! k-points were calculated beforehand.
         nx = k_points_hartree(n, 1)
         ny = k_points_hartree(n, 2)
         nz = k_points_hartree(n, 3)
         kx = nx * recip_lattice_vector(1,1) + ny * recip_lattice_vector(1,2) +  nz * recip_lattice_vector(1,3)
         ky = nx * recip_lattice_vector(2,1) + ny * recip_lattice_vector(2,2) +  nz * recip_lattice_vector(2,3)
         kz = nx * recip_lattice_vector(3,1) + ny * recip_lattice_vector(3,2) +  nz * recip_lattice_vector(3,3)
         k2 = kx**2 + ky**2 + kz**2
         vec_k(1) = kx
         vec_k(2) = ky
         vec_k(3) = kz

         ! CC: B-Matrix as defined in Heyes
         do i_coord=1,3
           AS_B_matrix(i_coord,i_coord,n) = -1.0d0
           do i_coord_2=1,3
             AS_B_matrix(i_coord,i_coord_2,n) = AS_B_matrix(i_coord,i_coord_2,n) + &
               ( ( (2.0d0 / k2) + (0.5d0 * Ewald_radius**2.0d0) ) * vec_k(i_coord) * vec_k(i_coord_2) )
           end do
         end do

      end do
    end if

    call sync_vector_complex(hartree_coef, n_k_points_hartree)
    if (analytical_stress_on) then
      do i_coord=1,3
        do i_coord_2=1,3
          call sync_vector_complex( AS_Lambda_matrix(i_coord,i_coord_2,:), n_k_points_hartree )
          ! CC: Sum up everything here to avoid recomputation during d/de v evaluation
          do n = 1,n_k_points_hartree
            AS_hartree_coeff_stress(i_coord,i_coord_2,n) =  &
               (hartree_coef(n) * AS_B_matrix(i_coord,i_coord_2,n)) &
               - AS_Lambda_matrix(i_coord,i_coord_2,n)
          end do
        end do
      end do
    end if


    ! VB: Legacy only: This would be the integration constant (self-term) for the charged Ewald summation
    !
    ! If we have non-neutral system we have to shift the potential level, so that the averace of
    ! the final is zero. This potential is calculated here and saved for use later.
    !
    !    pot_zero_level = 0.0
    !    do i_atom = 1, n_atoms
    !
    !   ! VB: orig incorrect version  pot_zero_level = pot_zero_level - multipole_moments(1, i_atom)*pi*Ewald_radius**2/cell_volume
    !   ! VB: This would have been correct:
    !   !     pot_zero_level = pot_zero_level - sqrt(pi4) * multipole_moments(1, i_atom)*pi*Ewald_radius**2/cell_volume
    !
    !    end do

  end subroutine evaluate_hartree_recip_coef
  !******







  !****s* hartree_potential_recip/evaluate_hartree_recip_gradient_coef_atom
  !  NAME
  !    evaluate_hartree_recip_gradient_coef_atom
  !  SYNOPSIS

  subroutine evaluate_hartree_recip_gradient_coef_atom( index_lm, multipole_moments,l_hartree_max_far_distance, i_atom)

    !  PURPOSE
    ! This routine calculates all the math of the Fourier part of the Hartree potential gradient
    ! which do not depend on the coordinate point of the potential for single atom i_atom. Routine is called
    ! once for every atom in the iteration at begin on 'Sum of whole potential', before going on every
    ! grid point.
    !
    ! This is needed in the periodic forces.
    !
    use analytic_multipole_coefficients, only: n_dd_lm_ijk, index_dd, dd
    use constants, only: pi
    use dimensions, only: n_atoms, l_pot_max
    use geometry, only: recip_lattice_vector, cell_volume
    use pbc_lists, only: center_to_atom, coords_center
    use runtime_choices, only: ewald_radius
    use synchronize_mpi, only: check_allocation
    implicit none
    !  ARGUMENTS

    integer index_lm(-l_pot_max:l_pot_max, 0:l_pot_max)
    real*8, dimension( ( l_pot_max + 1)**2, n_atoms) ::  multipole_moments
    integer:: l_hartree_max_far_distance(n_atoms)
    integer:: i_atom

    !  INPUTS
    !    o index_lm -- order of l and m components in multipole_moments table
    !    o multipole_moments -- multipole moments
    !    o l_hartree_max_far_distance -- maximum l component of Hartree potential in this part
    !    o i_atom -- atom which components are calculated
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE



    integer :: l_max
    integer :: i_l
    integer :: i_m, info
    integer :: ii, jj, kk, nx,ny,nz,i_ml_ijk
    integer :: n, i_coord
    real*8  :: kx,ky,kz,k2


    ! If this is the very first time routine is called perform initialization
    ! of the k-points
    if(.not. allocated(hartree_c_grad))then
       allocate(hartree_c_grad(3,n_k_points_hartree,1),stat=info)
       call check_allocation(info, 'hartree_c_grad                ')
    end if



    hartree_c_grad = (0.d0,0.d0)


    l_max =    l_hartree_max_far_distance(i_atom)

    do n = 1, n_k_points_hartree


       ! k-points were calculated beforehand.

       nx = k_points_hartree(n, 1)
       ny = k_points_hartree(n, 2)
       nz = k_points_hartree(n, 3)

       kx = nx * recip_lattice_vector(1,1) + ny * recip_lattice_vector(1,2) +  nz * recip_lattice_vector(1,3)
       ky = nx * recip_lattice_vector(2,1) + ny * recip_lattice_vector(2,2) +  nz * recip_lattice_vector(2,3)
       kz = nx * recip_lattice_vector(3,1) + ny * recip_lattice_vector(3,2) +  nz * recip_lattice_vector(3,3)

       k2 = kx**2 + ky**2 + kz**2


       ! We go trought list of so called D coefficients.

       do i_ml_ijk = 1, n_dd_lm_ijk

          ! l and m points of the multipole belong for particular d component
          i_l = index_dd(i_ml_ijk, 1)
          i_m = index_dd(i_ml_ijk, 2)

          ! Expodential terms for x, y and z direction  belong for particular d component
          ii  = index_dd(i_ml_ijk, 3)
          jj  = index_dd(i_ml_ijk, 4)
          kk  = index_dd(i_ml_ijk, 5)

          if(i_l <= l_max)then

             hartree_c_grad(1,n,1)  =     hartree_c_grad(1,n,1) &
                  + 4*pi / cell_volume /k2 &
                  *  multipole_moments(  index_lm(i_m, i_l), center_to_atom(i_atom))  &
                  * dd(i_ml_ijk) *(0,1)**i_l  * kx**ii * ky**jj * kz**kk &
                  * exp( (-Ewald_radius**2 * k2 / 4.) + &
                  (0,1)*( &
                  kx* ( - coords_center(1, i_atom)) + ky* (  - coords_center(2, i_atom)) +  &
                  kz* ( - coords_center(3, i_atom)))) &
                  *  ((0,1) * kx)


             hartree_c_grad(2,n,1)  =     hartree_c_grad(2,n,1) &
                  + 4*pi / cell_volume /k2 &
                  *  multipole_moments(  index_lm(i_m, i_l), center_to_atom(i_atom))  &
                  * dd(i_ml_ijk) *(0,1)**i_l  * kx**ii * ky**jj * kz**kk &
                  * exp( (-Ewald_radius**2 * k2 / 4.) + &
                  (0,1)*( &
                  kx* ( - coords_center(1, i_atom)) + ky* (  - coords_center(2, i_atom)) +  &
                  kz* ( - coords_center(3, i_atom)))) &
                  *  ((0,1) * ky)

             hartree_c_grad(3,n,1)  =     hartree_c_grad(3,n,1) &
                  + 4*pi / cell_volume /k2 &
                  *  multipole_moments(  index_lm(i_m, i_l), center_to_atom(i_atom))  &
                  * dd(i_ml_ijk) *(0,1)**i_l  * kx**ii * ky**jj * kz**kk &
                  * exp( (-Ewald_radius**2 * k2 / 4.) + &
                  (0,1)*( &
                  kx* ( - coords_center(1, i_atom)) + ky* (  - coords_center(2, i_atom)) +  &
                  kz* ( - coords_center(3, i_atom)))) &
                  *  ((0,1) * kz)
          end if

       end do
    end do


    !    write(use_unit,*) ' done'

  end subroutine evaluate_hartree_recip_gradient_coef_atom
  !******







  !****s* hartree_potential_recip/evaluate_potential_at_vacuum
  !  NAME
  !    evaluate_potential_at_vacuum
  !  SYNOPSIS

  subroutine evaluate_potential_at_vacuum &
    (z_coord, numx, numy, gradient, dip_origin, dip_length, potential_shift)

    !  PURPOSE
    !  The subroutine calculates and prints out Hartree potential at vacuum (= on the region
    !  where is only Fourier part of the potential) at given z-axis. Potential is averaged over
    !  given grid points.
    !
    !  USES

    use constants, only: bohr, hartree
    use geometry
    use localorb_io, only: localorb_info, use_unit, OL_norm
    use mpi_tasks, only: myid, n_tasks
    use runtime_choices, only: use_dipole_correction, vacuum_z_level
    use synchronize_mpi, only: sync_real_number
    implicit none

    !  ARGUMENTS

    real*8:: z_coord
    integer:: numx, numy
    real*8:: gradient
    real*8:: dip_origin
    real*8:: dip_length
    real*8 :: potential_shift

    !  INPUTS
    !   o z_coord -- z-coordinate of the place where potential is calculated
    !   o numx -- number of grid points averaged in x-direction
    !   o numy -- number of grid points averaged in y-direction
    !   o gradient -- gradient of dipole correction
    !   o dip_origin -- the point where dipole correction potential is zero
    !   o dip_length -- lenght of the supercell in z-direction
    !   o potential_shift --- The potential zero for periodic systems is shifted to the
    !       average real-space electrostatic potential of the unit cell, and this
    !       must be reflected in the work function
    !
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8:: coord_current(3), potential_average
    integer:: i_x, i_y, mul_int
    character*100 :: info_str

    coord_current(3) = z_coord
    potential_average = 0.d0

    do i_x = 1, numx

       if(mod(i_x-1,n_tasks) /= myid) cycle

       do i_y = 1, numy


          coord_current(1) = lattice_vector(1,1) * (i_x-1) / (numx-1) +  lattice_vector(1,2) * (i_y-1) / (numy-1)

          coord_current(2) = lattice_vector(2,1) * (i_x-1) / (numx-1) +  lattice_vector(2,2) * (i_y-1) / (numy-1)

          call  update_hartree_potential_recip( coord_current, potential_average)

       end do
    end do

    call sync_real_number(potential_average)

    potential_average = potential_average/ (numx * numy)

    potential_average = potential_average - potential_shift

    if( use_dipole_correction)then

       mul_int = int(floor((coord_current(3) - vacuum_z_level)/dip_length))
       coord_current(3) = coord_current(3) - mul_int * dip_length

       if( coord_current(3) <  vacuum_z_level)then
          coord_current(3) = coord_current(3) + dip_length
       end if

       potential_average = potential_average -  (coord_current(3)-dip_origin) * gradient

    end if

    write(info_str,*) z_coord*bohr,  potential_average*hartree,  -(coord_current(3)-dip_origin) * gradient*hartree
    call localorb_info(info_str,use_unit,'(A)',OL_norm)


  end subroutine evaluate_potential_at_vacuum
  !******







  !****s*  hartree_potential_recip/evaluate_dipole_correction
  !  NAME
  !   evaluate_dipole_correction
  !  SYNOPSIS

  subroutine evaluate_dipole_correction(potential_shift, gradient, dip_origin, dip_length, verbose)

    !  PURPOSE
    !  The subroutine calculates the size of dipole correction gradient and
    !  the work functions of surfaces.
    !
    !  USES

    use constants, only: hartree_over_bohr, bohr, hartree
    use geometry
    use physics, only: chemical_potential, pot_jump, vacuum_level_upper, &
        vacuum_level_lower
    use localorb_io, only: localorb_info, use_unit, OL_norm
    use runtime_choices, only: vacuum_z_level, use_dipole_correction
    implicit none

    !  ARGUMENTS

    real*8 :: potential_shift
    real*8:: gradient
    real*8:: dip_origin
    real*8:: dip_length
    logical :: verbose

    !  INPUTS
    !   o potential_shift --- The potential zero for periodic systems is shifted to the
    !       average real-space electrostatic potential of the unit cell, and this
    !       must be reflected in the work function
    !   o verbose -- determines whether here should be text output or not
    !  OUTPUT
    !   o gradient -- gradient of dipole correction
    !   o dip_origin -- the point where dipole correction potential is zero
    !   o dip_length -- lenght of the supercell in z-direction
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !    Computer Physics Communications (2008), submitted.
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE



    real*8:: z_coord
    integer, parameter:: numx = 10
    integer, parameter:: numy = 10


    real*8:: coord_current(3), potential_average_1, potential_average_2,  potential_average_3
    integer:: i_x, i_y, mul_int
    character*100 :: info_str
    !real*8:: pot_jump !OTH, 30Jun2011 --> moved to physics
    real*8 ::  gradient_1, gradient_2
    real*8 :: out_jump


    dip_length =  abs((maxval(lattice_vector(3,:)) - minval(lattice_vector(3,:))))

    mul_int = int(floor((vacuum_z_level)/dip_length))
    coord_current(3) = vacuum_z_level - mul_int * dip_length

    dip_origin = vacuum_z_level + dip_length * 0.5d0

    gradient = 0.d0
    pot_jump = 0.d0

    if(use_dipole_correction)then

       coord_current(3) = vacuum_z_level - 1

       potential_average_1 = 0.d0

       do i_x = 1, numx

          do i_y = 1, numy


             coord_current(1) = lattice_vector(1,1) * (i_x-1) / (numx-1) +  lattice_vector(1,2) * (i_y-1) / (numy-1)

             coord_current(2) = lattice_vector(2,1) * (i_x-1) / (numx-1) +  lattice_vector(2,2) * (i_y-1) / (numy-1)

             call  update_hartree_potential_recip( coord_current, potential_average_1)

          end do
       end do

       potential_average_1 = potential_average_1/ (numx * numy)

    end if

    coord_current(3) = vacuum_z_level

    potential_average_2 = 0.d0

    do i_x = 1, numx

       do i_y = 1, numy


          coord_current(1) = lattice_vector(1,1) * (i_x-1) / (numx-1) +  lattice_vector(1,2) * (i_y-1) / (numy-1)

          coord_current(2) = lattice_vector(2,1) * (i_x-1) / (numx-1) +  lattice_vector(2,2) * (i_y-1) / (numy-1)

          call  update_hartree_potential_recip( coord_current, potential_average_2)

       end do
    end do

    potential_average_2 = potential_average_2/ (numx * numy)

    if(use_dipole_correction)then
       coord_current(3) = vacuum_z_level +1

       potential_average_3 = 0.d0

       do i_x = 1, numx

          do i_y = 1, numy


             coord_current(1) = lattice_vector(1,1) * (i_x-1) / (numx-1) +  lattice_vector(1,2) * (i_y-1) / (numy-1)

             coord_current(2) = lattice_vector(2,1) * (i_x-1) / (numx-1) +  lattice_vector(2,2) * (i_y-1) / (numy-1)

             call  update_hartree_potential_recip( coord_current, potential_average_3)

          end do
       end do

       potential_average_3 = potential_average_3/ (numx * numy)

       gradient_1 = (potential_average_2 - potential_average_1)

       gradient_2 = (potential_average_3 - potential_average_2)



       if(abs(gradient_1) > 1e-9)then
          if(abs( (gradient_1 - gradient_2) / gradient_1) <  0.1)then

            gradient =  (gradient_1 + gradient_2) *0.5d0
            pot_jump =  - dip_length * gradient * 0.5d0

           if (verbose) then
             write(info_str,'(2X,A,F15.8,A)') &
                               '| Dipole correction gradient                  :', &
                               gradient * hartree_over_bohr, ' eV/Angstrom'
             call localorb_info(info_str,use_unit,'(A)',OL_norm)

             write(info_str,'(2X,A,F15.8,A)') &
                               '| Dipole correction potential jump            :', &
                               pot_jump * hartree * 2.d0, ' eV'
             call localorb_info(info_str,use_unit,'(A)',OL_norm)
           endif

          else

             out_jump = - dip_length * (gradient_1 + gradient_2) * 0.5d0 * 0.5d0

           if (verbose) then
             write(info_str,'(2X,A,F15.8,A)') &
                               '| Dipole correction jump would be         :', &
                               out_jump * hartree_over_bohr * 2.d0, ' eV/Angstrom'
             call localorb_info(info_str,use_unit,'(A)',OL_norm)

             write(info_str,'(2X,A,F10.2,A)') &
                  '| Gradient for dipole correction is not accurate enough. Error:', &
                  abs( (gradient_1 - gradient_2) / gradient_1)*100,' %'
             call localorb_info(info_str,use_unit,'(A)',OL_norm)

             write(info_str,'(2X,A)') &
                  '| Dipole correction not performed.'
             call localorb_info(info_str,use_unit,'(A)',OL_norm)
           endif

          end if
       end if

    end if




   !This work function output is done  before SCF
   !save chemical potential and vacuum levels for outputing work function after SCF:
   vacuum_level_lower = potential_average_2 - potential_shift - pot_jump
   vacuum_level_upper = potential_average_2 - potential_shift + pot_jump

!   if(verbose) then
   if(.false.) then
    write(info_str,'(2X,A,F15.8,A)') &
                               '| Potential vacuum level, "upper" slab surface:', &
                               (potential_average_2 - potential_shift + pot_jump)*hartree, ' eV'
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A,F15.8,A)') &
                               '| Potential vacuum level, "lower" slab surface:', &
                               (potential_average_2 - potential_shift - pot_jump)*hartree, ' eV'
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A,F15.8,A)') &
                               '| Work function ("upper" slab surface)        :', &
                      - (chemical_potential- (potential_average_2 - potential_shift + pot_jump))*hartree, ' eV'
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A,F15.8,A)') &
                               '| Work function ("lower" slab surface)        :', &
                      - (chemical_potential- (potential_average_2 - potential_shift - pot_jump))*hartree, ' eV'
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A)') &
      '*** Attention. The work function is given relative to the Fermi level ("chemical potential").'
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A)') &
      '*** The Fermi level is numerically determined and written to the output of FHI-aims, but'
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A)') &
      '*** for insulating systems, it may not be very meaningful. In that case, you may wish '
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A)') &
      '*** to correct the work function to reflect the highest occupied state '
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A)') &
      '*** (valence band maximum) of your system instead.'
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
   endif


    end subroutine evaluate_dipole_correction
    !******







    !****s* hartree_potential_recip/deallocate_hartree_p_recip
    !  NAME
    !   deallocate_hartree_p_recip
    !  SYNOPSIS

    subroutine deallocate_hartree_p_recip

    !  PURPOSE
    !  The subroutine deallocates modulevariables in hartree_potential_recip.

    use analytical_stress
    implicit none
    !
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

        if (allocated(k_points_hartree))   deallocate(k_points_hartree)
        if (allocated(hartree_coef))       deallocate(hartree_coef)
        if (allocated(hartree_c_grad))     deallocate(hartree_c_grad)

        if (allocated(number_z))           deallocate(number_z)

        if (allocated(AS_B_matrix))             deallocate(AS_B_matrix)
        if (allocated(AS_Lambda_matrix))        deallocate(AS_Lambda_matrix)
        if (allocated(AS_hartree_coeff_stress)) deallocate(AS_hartree_coeff_stress)

        n_k_points_hartree = 0
        k_points_max = 0
        not_allocated = .true.

    end subroutine deallocate_hartree_p_recip
    !******

    subroutine evaluate_dipole_correction_from_dipole(potential_shift, dip_grad, dip_origin, dip_length, verbose)
    !  PURPOSE
    !  The subroutine calculates the size of dipole correction gradient and
    !  the work functions of surfaces based on the dipole of the system
    !  (rather than the electrostatic potential at the VL)
    !
    !  USES

    use constants, only: hartree_over_bohr, pi, hartree
    use geometry
    use localorb_io, only: localorb_info, use_unit, OL_norm
    use physics, only: chemical_potential, pot_jump, vacuum_level_upper, &
        vacuum_level_lower
    use runtime_choices, only: vacuum_z_level
    implicit none

    !  ARGUMENTS


    real*8 :: potential_shift
    real*8:: dip_grad
    real*8:: dip_origin
    real*8:: dip_length
    logical :: verbose


    !  INPUTS
    !   o verbose -- determines whether here should be text output or not
    !  OUTPUT
    !   o gradient -- gradient of dipole correction
    !   o dip_origin -- the point where dipole correction potential is zero
    !   o dip_length -- lenght of the supercell in z-direction
    !   o verbose -- give output
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !    Computer Physics Communications (2008), submitted.
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    !auxilary variables
    real*8 :: dip(3) ! dipole moment
    real*8 :: Area ! surface area of the unit cell
    character(LEN=100) :: info_str
    integer :: i_x, i_y
    real*8 :: potential_average_2
    real*8 :: coord_current(3)

    !Here starts the work
    dip(:)=0.0

    Area=lattice_vector(1,1)*lattice_vector(2,2)-lattice_vector(2,1)*lattice_vector(1,2)


    call output_periodic_dipole_moment(dip)
    dip_length =  abs((maxval(lattice_vector(3,:)) - minval(lattice_vector(3,:))))
    dip_origin = vacuum_z_level + dip_length * 0.5d0
    pot_jump=-dip(3)/Area*2*pi
    dip_grad=-pot_jump/dip_length

!OTH Testing purpose
    potential_average_2 = 0.d0
    coord_current(3) = vacuum_z_level
    do i_x = 1, 10
       do i_y = 1, 10
          coord_current(1) = lattice_vector(1,1) * (i_x-1) / (10-1) +  lattice_vector(1,2) * (i_y-1) / (10-1)
          coord_current(2) = lattice_vector(2,1) * (i_x-1) / (10-1) +  lattice_vector(2,2) * (i_y-1) / (10-1)
          call  update_hartree_potential_recip( coord_current, potential_average_2)
       end do
    end do
    potential_average_2 = potential_average_2/(10 * 10)

   !save chemical potential and vacuum levels for outputing work function after SCF:
   vacuum_level_lower = potential_average_2 - potential_shift - pot_jump
   vacuum_level_upper = potential_average_2 - potential_shift + pot_jump
!OTH End testing




    if(verbose) then
       write(info_str,*) 'Area of unit cell: ', Area, ' Bohr**2'
       call localorb_info(info_str,use_unit,'(A)',OL_norm)

       write(info_str,'(2X,A,F15.8,A)') &
                          '| Dipole correction gradient                  :', &
                           2*dip_grad * hartree_over_bohr, ' eV/Angstrom'
       call localorb_info(info_str,use_unit,'(A)',OL_norm)

       write(info_str,'(2X,A,F15.8,A)') &
                         '| Dipole correction potential jump            :', &
                         2*pot_jump * hartree, ' eV'
       call localorb_info(info_str,use_unit,'(A)',OL_norm)


    endif


    end subroutine evaluate_dipole_correction_from_dipole




end module hartree_potential_recip
!******
