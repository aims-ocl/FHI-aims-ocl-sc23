!****h* FHI-aims/cartesian_ylm
!  NAME
!    cartesian_ylm 
!  SYNOPSIS

module cartesian_ylm

!  PURPOSE 
!    This module encapsulates everything needed to expand the spherical harmonics
!    into cartesian gaussian according to 
!
!    H. Bernhard Schlegel and Michael J. Frisch, 
!    "Transformation Between Cartesian and Pure Spherical Harmonic Gaussians",
!    International Journal of Quantum Chemistry, 54, 83-87 (1995)
!
!    Needed for the evaluation of the hessian of the basis functions which would 
!    be tedious if done in a straight forward way with spherical harmonics expressed 
!    in spherical coordinates.
!
!      Y_lm(theta, phi) = r^{-l} sum_{l_x l_y l_z with l_x + l_y + l_z = l} *
!                        coeff (l, m, l_x, l_y, l_z) * x^l_x * y^l_y * z^l_z
!
!    derivative terms needed for the gradient and hessian of a basis function
!      phi(r,theta,phi) = psi(r) / r^{l} sum_{l_x l_y l_z with l_x + l_y + l_z = l} *
!                       coeff (l, m, l_x, l_y, l_z) * x^l_x * y^l_y * z^l_z
!    are provided
!
!    individual subroutines:
!    * initialize_cartesian_ylm -- all necessary initializations
!    * cleanup_cartesian_ylm    -- deallocations
!    * evaluate_cartesians     -- needs to be called first to tabulize cartesian products x^l_x * y^l_y * z^l_z
!    * tab_ylm_cartesian     -- evaluates ylm-functions based upon the tabulized cartesians
!    * evaluate_cartesian_gradient_terms   -- required terms for the first derivative of a basis functions 
!    * evaluate_cartesian_hessian_and_gradient_terms    -- required terms for the hessian of a basis functions
!    * evaluate_cartesian_hessian_and_gradient_terms_p2 
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
!



!  Declaration of variables

  private

  integer :: max_l_prepared = -1
  ! l(1:3, i_cartesian, i_l) -> exponent l_x, l_y, l_z in cartesians
  integer, dimension(:,:,:), allocatable :: l
  ! coefficients of expansion of ylms into cartesians; use like in:
  ! ylm_from_cart(-L:L, 1:n_tree(L)) = 0.d0
  ! do M = -L, L
  !    i_lm = index_lm(M, L)
  !    do i_coeff = 1, n_cartesian(i_lm)
  !       i_cart = which_cartesian(i_coeff, i_lm)
  !       ylm_from_cart(M, i_cart) = ylm_from_cart(M, i_cart) &
  !       &                        + coeff(i_coeff, i_lm)
  !    end do
  ! end do
  real*8, dimension(:,:), allocatable :: coeff
  ! n_tree(i_l) number of all possible cartesians for a given l
  integer, dimension(:), allocatable :: n_tree
  integer, dimension(:,:), allocatable :: tree ! describes how cartesians for l are calculated from the cartesians for l-1
  integer, dimension(:,:), allocatable :: which_cartesian ! maps coefficients to the corresponding cartesians
  integer, dimension(:,:,:), allocatable :: which_cartesian_gradient ! maps coefficients to the corresponding cartesians for gradient terms
  integer, dimension(:,:,:), allocatable :: which_cartesian_hessian ! maps coefficients to the corresponding cartesians for hessian terms
  integer, dimension(:), allocatable :: n_cartesian ! i_lm -> number of cartesians with non-zero coefficients for given Y_lm.
  integer, dimension(:,:), allocatable, private :: index_lm
  integer, public :: n_max_cartesian


  public :: initialize_cartesian_ylm, cleanup_cartesian_ylm, &
       evaluate_cartesians, & 
       tab_ylm_cartesian, evaluate_cartesian_gradient_terms, &
       evaluate_cartesian_hessian_and_gradient_terms, &
       evaluate_cartesian_hessian_and_gradient_terms_p2, &
       evaluate_onecenter_cartesians, & 
       tab_ylm_onecenter_cartesian, evaluate_onecenter_cartesian_gradient_terms, &
       generate_cartesian_rotation_element, generate_cartesian_rotation, &
       test_cartesian_rotation, test_ylm_rotation

  contains


!******
!----------------------------------------------------------------------------------------------
!****s* cartesian_ylm/initialize_cartesian_ylm
!  NAME
!    initialize_cartesian_ylm
!  SYNOPSIS
    subroutine initialize_cartesian_ylm(max_l)
!  PURPOSE 
!    All required working arrays are allocated.
!    Lookuptables are initialized:
!    * tree, ntree 
!        used to calculate cartesians x^l_x * y^l_y * z^l_z for l_x + l_y + l_z = l 
!        efficiently from x^l_x * y^l_y * z^l_z for l_x + l_y + l_z = l-1
!    * which_cartesian 
!        gives for each Y_lm functions the indizes for the required cartesians
!    * which_cartesian_gradient 
!        gives for each gradient term of a basis function the indizes of the required cartesians
!    * which_cartesian_hessian 
!        gives for each hessian term of a basis function the indizes for the required cartesians
!    * l()  exponent l_x, l_y, l_z in cartesians
!    
!    Expansion coefficients coeff() are evaluated.
!
!  USES
      use runtime_choices, only: coeff_thr
      use mpi_tasks, only: check_allocation
      implicit none
!  ARGUMENTS
!  INPUTS
!    max_l -- maximum angular momentum for which to prepare arrays
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      integer, intent(IN) :: max_l


      ! local variables

      real*8 :: temp_coeff
      integer, dimension(3) :: target_l
      integer :: info

      ! functions
      real*8, external :: binomial_coefficient
      
      ! counter
      integer :: i_coords
      integer :: i_coords_2
      integer :: i_coords_3
      integer :: i_tree
      integer :: i_tree_2
      integer :: i_coeff
      integer :: i_l
      integer :: i_m
      integer :: i_index
      integer :: i_counter

      if (max_l <= max_l_prepared) then
         return
      else if (max_l_prepared >= 0) then
         call cleanup_cartesian_ylm()
      end if
      max_l_prepared = max_l

!      pi = (4.0D+0)*ATAN(1.0D+0)
!      pi4_inv = 1. / (4 * pi)
      ! n_max_cartesians is the number of combinations with repetition binomial_coefficient(n + k - 1    k)
      n_max_cartesian = int(binomial_coefficient(3 + max_l - 1, max_l))
      if (.not.allocated(l)) then
         allocate(l(3, n_max_cartesian, 0:max_l),stat=info)
         call check_allocation(info, 'l                             ')
      end if
      if (.not.allocated(coeff)) then
         allocate(coeff(n_max_cartesian, (max_l+1) ** 2),stat=info)
         call check_allocation(info, 'coeff                         ')
      end if
      if (.not.allocated(index_lm)) then
         allocate(index_lm(-max_l:max_l, 0:max_l),stat=info)
         call check_allocation(info, 'index_lm                      ')
      end if
      if (.not.allocated(n_tree)) then
         allocate(n_tree(0:max_l),stat=info)
         call check_allocation(info, 'n_tree                        ')
      end if
      if (.not.allocated(tree)) then
         allocate(tree(n_max_cartesian, 0:max_l),stat=info)
         call check_allocation(info, 'tree                          ')
      end if
      if (.not.allocated(which_cartesian)) then
         allocate(which_cartesian(n_max_cartesian, (max_l+1) ** 2),stat=info)
         call check_allocation(info, 'which_cartesian               ')
      end if
      if (.not.allocated(which_cartesian_gradient)) then
         allocate(which_cartesian_gradient(3, n_max_cartesian, (max_l+1) ** 2),stat=info)
         call check_allocation(info, 'which_cartesian_gradient      ')
      end if
      if (.not.allocated(which_cartesian_hessian)) then
         allocate(which_cartesian_hessian(6, n_max_cartesian, (max_l+1) ** 2),stat=info)
         call check_allocation(info, 'which_cartesian_hessian       ')
      end if
      if (.not.allocated(n_cartesian)) then
         allocate(n_cartesian((max_l+1) ** 2),stat=info)
         call check_allocation(info, 'n_cartesian                   ')
      end if
      i_index = 0
      do i_l = 0, max_l, 1
         do i_m = -i_l, i_l
            i_index = i_index + 1
            index_lm(i_m, i_l) = i_index
         enddo
      enddo
      ! initialize lookup-tables tree, n_tree, l(:)
      n_tree(0) = 1
      l(:, 1, 0) = 0
      if (max_l .gt. 0) then
         do i_coords = 1, 3, 1
            tree(i_coords, 1) = i_coords 
         end do
         do i_counter = 1, 3, 1
            l(:, i_counter, 1) = 0
         end do
         do i_counter = 1, 3, 1
            l(i_counter, i_counter, 1) = 1
         end do
         n_tree(1) = 3
         do i_l = 2, max_l, 1
            i_counter = 0
            do i_tree = 1, n_tree(i_l - 1), 1
               do i_coords = tree(i_tree, i_l - 1), 3, 1
                  i_counter = i_counter + 1
                  tree(i_counter, i_l) = i_coords
                  l(:, i_counter, i_l) = l(:, i_tree, i_l - 1)
                  l(i_coords, i_counter, i_l) = l(i_coords, i_tree, i_l - 1) + 1
               end do
            end do
            n_tree(i_l) = i_counter
         end do
      end if
      ! evaluate coefficients and initialize which_cartesian
      do i_l = 0, max_l, 1
         do i_m = -i_l, i_l, 1
            i_coeff = 0
            do i_tree = 1, n_tree(i_l), 1
               temp_coeff = aims_real_cartesian_coefficient(m=i_m, &
                                 lx=l(1,i_tree,i_l), ly=l(2,i_tree,i_l), &
                                 lz=l(3,i_tree,i_l) )

               if (abs(temp_coeff) .gt. coeff_thr) then
                  i_coeff = i_coeff + 1
                  coeff(i_coeff, index_lm(i_m, i_l)) = temp_coeff

                  ! make assignment in which_cartesian
                  which_cartesian(i_coeff, index_lm(i_m, i_l)) = i_tree

                  ! make assignment in which_cartesian_gradient
                  if (i_l .gt. 0) then
                     which_cartesian_gradient(:, i_coeff, index_lm(i_m, i_l)) = -1
                     do i_tree_2 = 1, n_tree(i_l - 1), 1
                        do i_coords = 1, 3, 1
                           do i_coords_2 = 1, 3, 1
                              target_l(i_coords_2) = l(i_coords_2, i_tree, i_l)
                           end do
                           target_l(i_coords) = target_l(i_coords) - 1
                           if ((l(1, i_tree_2, i_l - 1) .eq. target_l(1)) .and. & 
                                (l(2, i_tree_2, i_l - 1) .eq. target_l(2)) .and. &
                                (l(3, i_tree_2, i_l - 1) .eq. target_l(3))) then
                              which_cartesian_gradient(i_coords, i_coeff, index_lm(i_m, i_l)) = i_tree_2
                           end if
                        end do
                     end do
                  end if

                  ! make assignment in which_cartesian_hessian
                  if (i_l .gt. 1) then
                     which_cartesian_hessian(:, i_coeff, index_lm(i_m, i_l)) = -1
                     do i_tree_2 = 1, n_tree(i_l - 2), 1
                        i_counter = 0
                        do i_coords = 1, 3, 1
                           do i_coords_2 = i_coords, 3, 1
                              i_counter = i_counter + 1
                              do i_coords_3 = 1, 3, 1
                                 target_l(i_coords_3) = l(i_coords_3, i_tree, i_l)
                              end do
                              target_l(i_coords) = target_l(i_coords) - 1
                              target_l(i_coords_2) = target_l(i_coords_2) - 1
                              if ( (l(1, i_tree_2, i_l - 2) .eq. target_l(1)) .and. & 
                                   (l(2, i_tree_2, i_l - 2) .eq. target_l(2)) .and. &
                                   (l(3, i_tree_2, i_l - 2) .eq. target_l(3))) then
                                 which_cartesian_hessian(i_counter, i_coeff, index_lm(i_m, i_l)) = i_tree_2
                              end if
                           end do
                        end do
                     end do
                  end if

               end if ! abs(coeff) > threshold
            end do ! i_tree
            n_cartesian(index_lm(i_m, i_l)) = i_coeff
         end do ! i_m
      end do ! i_l

    end subroutine initialize_cartesian_ylm

!******
!------------------------------------------------------------------------------------------------------------
!****s* cartesian_ylm/cleanup_cartesian_ylm
!  NAME
!   cleanup_cartesian_ylm 
!  SYNOPSIS   
    subroutine cleanup_cartesian_ylm()
!  PURPOSE
!   deallocations of module cartesian ylm  
!
!  USES
!  ARGUMENTS


!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      max_l_prepared = -1
      if (allocated(n_cartesian)) then
         deallocate(n_cartesian)
      end if
      if (allocated(which_cartesian_hessian)) then
         deallocate(which_cartesian_hessian)
      end if
      if (allocated(which_cartesian_gradient)) then
         deallocate(which_cartesian_gradient)
      end if
      if (allocated(which_cartesian)) then
         deallocate(which_cartesian)
      end if
      if (allocated(tree)) then
         deallocate(tree)
      end if
      if (allocated(n_tree)) then
         deallocate(n_tree)
      end if
      if (allocated(index_lm)) then
         deallocate(index_lm)
      end if
      if (allocated(coeff)) then
         deallocate(coeff)
      end if
      if (allocated(l)) then
         deallocate(l)
      end if


    end subroutine cleanup_cartesian_ylm

!******
!------------------------------------------------------------------------------------------------------
!****s* cartesian_ylm/evaluate_onecenter_cartesians
!  NAME
!   evaluate_onecenter_cartesians
!  SYNOPSIS

    subroutine evaluate_onecenter_cartesians(relvec, l_max, cartesians)

!  PURPOSE
!    Evaluate and tabulates all cartesian terms 
!
!      x^l_x * y^l_y * z^l_z for l_x + l_y + l_z = l 
!
!    for the relvec == (/x, y, z/).
!
!  USES

      use dimensions, only: n_species
      use mpi_tasks, only: aims_stop
      use pbc_lists, only: species_center
      implicit none

!  ARGUMENTS

      real*8, dimension(3), intent(in) :: relvec
      integer, intent(in) :: l_max
      real*8, dimension(n_max_cartesian, 0:l_max), intent(out) :: cartesians

!  INPUTS
!   o relvec -- direction to atoms
!   o l_max -- maximum l index
!
!  OUTPUT
!   o cartesians -- cartesian terms x^l_x * y^l_y * z^l_z
      !             for l_x + l_y + l_z = l
! 
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      ! counter
      integer :: i_coords
      integer :: i_l
      integer :: i_tree
      integer :: i_counter
      character(*), parameter :: func = 'evaluate_onecenter_cartesians'

      if (l_max > max_l_prepared) then
         call aims_stop('Need to initialize_cartesian_ylm with higher l', func)
      end if
      cartesians(1, 0) = 1.d0
      if (l_max > 0) then
         do i_coords = 1, 3, 1
            cartesians(i_coords, 1) = relvec(i_coords)
         end do
         do i_l = 2, l_max
            i_counter = 1
            do i_tree = 1, n_tree(i_l)
               cartesians(i_tree, i_l) = cartesians(i_counter, i_l - 1) &
               &                         * relvec(tree(i_tree, i_l))
               if (tree(i_tree, i_l) .eq. 3) then
                  i_counter = i_counter + 1
               end if
            end do
         end do
      end if
      
    end subroutine evaluate_onecenter_cartesians

!******
!------------------------------------------------------------------------------------------------------
!****s* cartesian_ylm/evaluate_cartesians
!  NAME
!   evaluate_cartesians
!  SYNOPSIS

    subroutine evaluate_cartesians(dir_tab, l_max, l_ylm_max, atom_index, n_compute_atoms, cartesians)

!  PURPOSE
!    Evaluate and tabulates all cartesian terms 
!
!      x^l_x * y^l_y * z^l_z for l_x + l_y + l_z = l 
!
!    for all atoms for a given configurations described in dir_tab.
!
!  USES

      use dimensions, only: n_species
      use pbc_lists, only: species_center
      implicit none

!  ARGUMENTS

      real*8, dimension(3, n_compute_atoms), intent(in) :: dir_tab
      integer, dimension(n_species), intent(in) :: l_max
      integer, intent(in) :: l_ylm_max
      integer, dimension(n_compute_atoms), intent(in) :: atom_index
      integer, intent(in) :: n_compute_atoms
      real*8, dimension(n_max_cartesian, 0:l_ylm_max, n_compute_atoms), intent(out) :: cartesians

!  INPUTS
!   o dir_tab -- direction to atoms
!   o l_max -- maximum l index
!   o l_ylm_max -- maximum l index in Y_lm functions
!   o atom_index -- list of relevant atoms
!   o n_compute_atoms -- number of relevant atoms
!
!  OUTPUT
!   o cartesians --  cartesian terms x^l_x * y^l_y * z^l_z for l_x + l_y + l_z = l
! 
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      ! counter
      integer :: i_compute
      integer :: i_atom
      integer :: i_species

      ! evaluate all cartesian gaussians with lookup-tables
      do i_compute = 1, n_compute_atoms, 1
         i_atom = atom_index(i_compute)
         i_species = species_center(i_atom)
         call evaluate_onecenter_cartesians( &
         & dir_tab(:, i_compute), l_max(i_species), cartesians(:,:, i_compute))
      end do
      
    end subroutine evaluate_cartesians

!******
!---------------------------------------------------------------------------------------------------------------
!****s* cartesian_ylm/tab_ylm_onecenter_cartesian
!  NAME
!    tab_ylm_onecenter_cartesian
!  SYNOPSIS   

    subroutine tab_ylm_onecenter_cartesian(l_max, cartesians, rlylm)

!  PURPOSE
!    Evaluates ylm_functions assuming that cartesians are already evaluated !!!
!    -> subroutine evaluate_onecenter_cartesians()
!
!  USES

    use constants, only: pi4_inv
    use mpi_tasks, only: aims_stop
    implicit none

!  ARGUMENTS

    integer, intent(in) :: l_max
    real*8, dimension(n_max_cartesian, 0:l_max), intent(in) :: cartesians
    real*8, dimension((l_max+1)*(l_max+1)), intent(out) :: rlylm

!  INPUTS
!   o l_max -- maximum l index
!   o cartesians -- cartesian terms x^l_x * y^l_y * z^l_z for l_x + l_y + l_z = !
!  OUTPUT
!   o rlylm -- r^l * Y_lm
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      ! local variables
      integer :: lm_index

      ! counter
      integer :: i_l
      integer :: i_m
      integer :: i_coeff
      integer :: i_atom
      character(*), parameter :: func = 'tab_ylm_onecenter_cartesian'

      if (l_max > max_l_prepared) then
         call aims_stop('Need to initialize_cartesian_ylm with higher l', func)
      end if
      rlylm = 0.d0
      rlylm(1) = sqrt(pi4_inv)
      do i_l = 1, l_max
         do i_m = -i_l, i_l, 1
            lm_index = index_lm(i_m, i_l)
            do i_coeff = 1, n_cartesian(lm_index), 1
               rlylm(lm_index) = rlylm(lm_index) + &
               & coeff(i_coeff, lm_index) &
               & * cartesians(which_cartesian(i_coeff, lm_index), i_l)
            end do
         end do
      end do

    end subroutine tab_ylm_onecenter_cartesian

!******
!---------------------------------------------------------------------------------------------------------------
!****s* cartesian_ylm/tab_ylm_cartesian
!  NAME
!    tab_ylm_cartesian
!  SYNOPSIS   

    subroutine tab_ylm_cartesian(l_max, l_ylm_max, cartesians, atom_index, n_compute_atoms, ylm_tab)

!  PURPOSE
!    Evaluates ylm_functions assuming that cartesians are already evaluated !!!
!    -> subroutine evaluate_cartesians()
!
!  USES

    use dimensions, only: n_species, n_atoms
    use pbc_lists, only: species_center
    implicit none

!  ARGUMENTS

    integer, dimension(n_species), intent(in) :: l_max
    integer, intent(in) :: l_ylm_max
    real*8, dimension(n_max_cartesian, 0:l_ylm_max, n_atoms), intent(in) :: cartesians
    integer, dimension(n_atoms), intent(in) :: atom_index
    integer, intent(in) :: n_compute_atoms
    real*8, dimension((l_ylm_max+1)*(l_ylm_max+1), n_atoms), intent(out) :: ylm_tab

!  INPUTS
!   o l_max -- maximum l index
!   o l_ylm_max -- maximum l index in Y_lm functions
!   o cartesians -- cartesian terms x^l_x * y^l_y * z^l_z for l_x + l_y + l_z = l
!   o atom_index -- list of relevant atoms
!   o n_compute_atoms -- number of relevant atoms
!
!  OUTPUT
!   o ylm_tab -- r^l * Y_lm
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      ! local variables
      integer :: lm_index
      integer :: i_atom
      integer :: i_species
      integer :: i_compute

      ! evaluate spherical harmonic with lookup-tables
      do i_compute = 1, n_compute_atoms, 1
         i_atom = atom_index(i_compute)
         i_species = species_center(i_atom)
         call tab_ylm_onecenter_cartesian( &
         & l_max(i_species), cartesians(:,:, i_compute), ylm_tab(:, i_compute))
      end do

    end subroutine tab_ylm_cartesian

!******
!------------------------------------------------------------------------------------------------------------
!****s* cartesian_ylm/evaluate_onecenter_cartesian_gradient_terms
!  NAME
!   evaluate_onecenter_cartesian_gradient_terms 
!  SYNOPSIS

    subroutine evaluate_onecenter_cartesian_gradient_terms(l_max, cartesians, &
    &                                                      rlylm, drlylm)

!  PURPOSE
!    Evaluates terms that appear in connection with the first derivative of a basis function,
!    like e.g. 
!
!      sum_{l_x l_y l_z with l_x + l_y + l_z = l} * 
!              coeff (l, m, l_x, l_y, l_z) * l_x * x^{l_x - 1} * y^l_y * z^l_z
!
!  USES

      use mpi_tasks, only: aims_stop
      implicit none

!  ARGUMENTS

      integer, intent(in) :: l_max
      real*8, dimension(n_max_cartesian, 0:l_max), intent(in) :: cartesians
      real*8, dimension((l_max+1) ** 2) :: rlylm
      real*8, dimension(3, (l_max+1) ** 2) :: drlylm

!  INPUTS
!   o l_max -- maximum l index
!   o cartesians -- cartesian terms x^l_x * y^l_y * z^l_z for l_x + l_y + l_z = l
!  OUTPUTS
!   o rlylm -- r^l * Y_lm
!   o drlylm -- drlylm(i, lm) := d (r^l * Y_lm) / dx_i
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      ! counter
      integer :: i_atom
      integer :: i_l
      integer :: i_m
      integer :: i_lm
      integer :: i_coeff
      integer :: i_coords, i_power
      integer :: i_cart, i_cart_grad
      character(*), parameter :: func = 'evaluate_onecenter_cartesian_gradient_terms'

      if (l_max > max_l_prepared) then
         call aims_stop('Need to initialize_cartesian_ylm with higher l', func)
      end if
      do i_l = 0, l_max
         do i_m = -i_l, i_l, 1
            i_lm = index_lm(i_m, i_l)
            rlylm(i_lm) = 0.d0
            drlylm(:, i_lm) = 0.d0
            do i_coeff = 1, n_cartesian(i_lm)
               i_cart = which_cartesian(i_coeff, i_lm)
               rlylm(i_lm) = rlylm(i_lm) + &
               & coeff(i_coeff, i_lm) * cartesians(i_cart, i_l)
               do i_coords = 1, 3
                  i_power = l(i_coords, i_cart, i_l)
                  if (i_power .gt. 0) then
                     i_cart_grad = which_cartesian_gradient(i_coords, i_coeff, i_lm)
                     drlylm(i_coords, i_lm) &
                     & = drlylm(i_coords, i_lm) + &
                     & i_power * coeff(i_coeff, i_lm) * cartesians(i_cart_grad, i_l-1)
                  end if
               end do
            end do
         end do
      end do
    end subroutine evaluate_onecenter_cartesian_gradient_terms

!******
!--------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!****s* cartesian_ylm/evaluate_cartesian_gradient_terms
!  NAME
!   evaluate_cartesian_gradient_terms 
!  SYNOPSIS

    subroutine evaluate_cartesian_gradient_terms(l_max, l_ylm_max, cartesians, atom_index, n_compute_atoms, &
         rlylm, drlylm)

!  PURPOSE
!    Evaluates terms that appear in connection with the first derivative of a basis function,
!    like e.g. 
!
!      sum_{l_x l_y l_z with l_x + l_y + l_z = l} * 
!              coeff (l, m, l_x, l_y, l_z) * l_x * x^{l_x - 1} * y^l_y * z^l_z
!
!  USES

      use dimensions, only: n_species, n_atoms
      use pbc_lists, only: species_center
      implicit none

!  ARGUMENTS

      integer, dimension(n_species), intent(in) :: l_max
      integer :: l_ylm_max
      real*8, dimension(n_max_cartesian, 0:l_ylm_max, n_atoms), intent(in) :: cartesians
      integer, dimension(n_atoms), intent(in) :: atom_index
      integer, intent(in) :: n_compute_atoms
      real*8, dimension((l_ylm_max+1) ** 2, n_atoms) :: rlylm
      real*8, dimension(3, (l_ylm_max+1) ** 2, n_atoms) :: drlylm

!  INPUTS
!   o l_max -- maximum l index
!   o l_ylm_max -- maximum l index in Y_lm functions
!   o cartesians -- cartesian terms x^l_x * y^l_y * z^l_z for l_x + l_y + l_z = l
!   o atom_index -- list of relevant atoms
!   o n_compute_atoms -- number of relevant atoms
!
!  OUTPUTS
!   o rlylm -- r^l * Y_lm
!   o drlylm -- drylm(i, lm, i_compute) d (r^l Y_lm) / dx_i
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      ! counter
      integer :: i_compute
      integer :: i_atom
      integer :: i_species

      do i_compute = 1, n_compute_atoms, 1
         i_atom = atom_index(i_compute)
         i_species = species_center(i_atom)
         call evaluate_onecenter_cartesian_gradient_terms( &
         & l_max(i_species), cartesians(:,:, i_compute), &
         & rlylm(:, i_compute), drlylm(:,:, i_compute))
      end do

    end subroutine evaluate_cartesian_gradient_terms

!******
!--------------------------------------------------------------------------------------------------------------
!****s* cartesian_ylm/evaluate_cartesian_hessian_and_gradient_terms
!  NAME
!   evaluate_cartesian_hessian_and_gradient_terms
!  SYNOPSIS    

    subroutine evaluate_cartesian_hessian_and_gradient_terms(l_max, l_ylm_max, cartesians, atom_index, n_compute_atoms, &
         ylm_tab, sum_gradient, sum_hessian)

!  PURPOSE
!    Evaluates terms that appear in connection with the first derivative and hessian of a basis function,
!    like e.g. 
!        sum_{l_x l_y l_z with l_x + l_y + l_z = l} 
!                * coeff (l, m, l_x, l_y, l_z) * l_x * l_y * x^{l_x - 1} * y^{l_y - 1} * z^l_z
!
!  USES

    use dimensions, only: n_species
    use pbc_lists, only: species_center
    implicit none

!  ARGUMENTS

      integer, dimension(n_species), intent(in) :: l_max
      integer, intent(in) :: l_ylm_max
      real*8, dimension(n_max_cartesian, 0:l_ylm_max, n_compute_atoms), intent(in) :: cartesians
      integer, dimension(n_compute_atoms), intent(in) :: atom_index
      integer, intent(in) :: n_compute_atoms
      real*8, dimension((l_ylm_max+1) ** 2, n_compute_atoms) :: ylm_tab
      real*8, dimension(3, (l_ylm_max+1) ** 2, n_compute_atoms) :: sum_gradient
      real*8, dimension(6, (l_ylm_max+1) ** 2, n_compute_atoms) :: sum_hessian

!  INPUTS
!   o l_max -- maximum l index
!   o l_ylm_max -- maximum l index in Y_lm functions
!   o cartesians -- cartesian terms x^l_x * y^l_y * z^l_z for l_x + l_y + l_z = l
!   o atom_index -- list of relevant atoms
!   o n_compute_atoms -- number of relevant atoms
!
!  OUTPUTS
!   o ylm_tab -- Y_lm functions
!   o sum_gradient -- ???????
!   o sum_hessian -- ??????????
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



      ! local variables
      integer :: lm_index
      real*8 :: coeff_cartesian
      integer :: current_atom

      ! counter
      integer :: i_atom
      integer :: i_l
      integer :: i_m
      integer :: i_coeff
      integer :: i_coords
      integer :: i_coords_2
      integer :: i_counter

      do i_atom = 1, n_compute_atoms, 1
         current_atom = atom_index(i_atom)
         do i_l = 0, l_max(species_center(current_atom)), 1
            do i_m = -i_l, i_l, 1
               lm_index = index_lm(i_m, i_l)
               ylm_tab(lm_index, i_atom)    = 0.d0
               sum_gradient(:, lm_index, i_atom) = 0.d0
               sum_hessian(:, lm_index, i_atom)      = 0.d0
               do i_coeff = 1, n_cartesian(lm_index), 1
                  ylm_tab(lm_index, i_atom) = ylm_tab(lm_index, i_atom) + &
                       coeff(i_coeff, lm_index) * cartesians(which_cartesian(i_coeff, lm_index), i_l, i_atom)
                  i_counter = 0
                  do i_coords = 1, 3, 1
                     ! gradient terms
                     if (l(i_coords, which_cartesian(i_coeff, lm_index), i_l) .gt. 0) then
                        sum_gradient(i_coords, lm_index, i_atom) = sum_gradient(i_coords, lm_index, i_atom) + &
                             l(i_coords, which_cartesian(i_coeff, lm_index), i_l) * &
                             coeff(i_coeff, lm_index) * & 
                             cartesians(which_cartesian_gradient(i_coords, i_coeff, lm_index), i_l-1, i_atom)
                     end if

                     ! diagonal terms for hessian
                     i_counter = i_counter + 1
                     if (l(i_coords, which_cartesian(i_coeff, lm_index), i_l) .gt. 1) then
                        sum_hessian(i_counter, lm_index, i_atom) = sum_hessian(i_counter, lm_index, i_atom) + &
                             l(i_coords, which_cartesian(i_coeff, lm_index), i_l) * & 
                             (l(i_coords, which_cartesian(i_coeff, lm_index), i_l) - 1) * &
                             coeff(i_coeff, lm_index) * & 
                             cartesians(which_cartesian_hessian(i_counter, i_coeff, lm_index), i_l-2, i_atom)
                     end if

                     ! non-diagonal terms for hessian
                     do i_coords_2 = i_coords + 1, 3, 1
                        i_counter = i_counter + 1
                        if ((l(i_coords, which_cartesian(i_coeff, lm_index), i_l) .gt. 0) .and. &
                             (l(i_coords_2, which_cartesian(i_coeff, lm_index), i_l) .gt. 0)) then
                           sum_hessian(i_counter, lm_index, i_atom) = sum_hessian(i_counter, lm_index, i_atom) + &
                                l(i_coords, which_cartesian(i_coeff, lm_index), i_l) * &
                                l(i_coords_2, which_cartesian(i_coeff, lm_index), i_l) * &
                                coeff(i_coeff, lm_index) * & 
                                cartesians(which_cartesian_hessian(i_counter, i_coeff, lm_index), i_l-2, i_atom)
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do
    end subroutine evaluate_cartesian_hessian_and_gradient_terms
!******
!------------------------------------------------------------------------------------------------------------
!****s* cartesian_ylm/evaluate_cartesian_hessian_and_gradient_terms_p2
!  NAME
!   evaluate_cartesian_hessian_and_gradient_terms_p2
!  SYNOPSIS    

    subroutine evaluate_cartesian_hessian_and_gradient_terms_p2(l_max, l_ylm_max, cartesians, atom_index, n_compute_atoms, &
         ylm_tab, sum_gradient, sum_hessian)

!  PURPOSE
!    Evaluates terms that appear in connection with the first derivative and hessian of a basis function,
!    like e.g. 
!
!       sum_{l_x l_y l_z with l_x + l_y + l_z = l} 
!           * coeff (l, m, l_x, l_y, l_z) * l_x * l_y * x^{l_x - 1} * y^{l_y - 1} * z^l_z
!
!    VB: Modified version that changes the order of indices in sum_gradient!!
!

!  USES
    use dimensions, only: n_species
    use pbc_lists, only: species_center
    implicit none
!  ARGUMENTS

    integer, dimension(n_species), intent(in) :: l_max
    integer, intent(in) :: l_ylm_max
    real*8, dimension(n_max_cartesian, 0:l_ylm_max, n_compute_atoms), intent(in) :: cartesians
    integer, dimension(n_compute_atoms), intent(in) :: atom_index
    integer, intent(in) :: n_compute_atoms
    real*8, dimension((l_ylm_max+1) ** 2, n_compute_atoms) :: ylm_tab
    real*8, dimension((l_ylm_max+1) ** 2, 3, n_compute_atoms) :: sum_gradient
    real*8, dimension(6, (l_ylm_max+1) ** 2, n_compute_atoms) :: sum_hessian

!  INPUTS
!   o l_max -- maximum l index
!   o l_ylm_max -- maximum l index in Y_lm functions
!   o cartesians -- cartesian terms x^l_x * y^l_y * z^l_z for l_x + l_y + l_z = l
!   o atom_index -- list of relevant atoms
!   o n_compute_atoms -- number of relevant atoms
!
!  OUTPUTS
!   o ylm_tab -- Y_lm functions
!   o sum_gradient -- ???????
!   o sum_hessian -- ??????????
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE





      ! local variables
      integer :: lm_index
      real*8 :: coeff_cartesian
      integer :: current_atom

      ! counter
      integer :: i_atom
      integer :: i_l
      integer :: i_m
      integer :: i_coeff
      integer :: i_coords
      integer :: i_coords_2
      integer :: i_counter

      do i_atom = 1, n_compute_atoms, 1
         current_atom = atom_index(i_atom)
         do i_l = 0, l_max(species_center(current_atom)), 1
            do i_m = -i_l, i_l, 1
               lm_index = index_lm(i_m, i_l)
               ylm_tab(lm_index, i_atom)    = 0.d0
               sum_gradient(lm_index, :, i_atom) = 0.d0
               sum_hessian(:, lm_index, i_atom)      = 0.d0
               do i_coeff = 1, n_cartesian(lm_index), 1
                  ylm_tab(lm_index, i_atom) = ylm_tab(lm_index, i_atom) + &
                       coeff(i_coeff, lm_index) * cartesians(which_cartesian(i_coeff, lm_index), i_l, i_atom)
                  i_counter = 0
                  do i_coords = 1, 3, 1
                     ! gradient terms
                     if (l(i_coords, which_cartesian(i_coeff, lm_index), i_l) .gt. 0) then
                        sum_gradient(lm_index, i_coords, i_atom) = sum_gradient(lm_index, i_coords, i_atom) + &
                             l(i_coords, which_cartesian(i_coeff, lm_index), i_l) * &
                             coeff(i_coeff, lm_index) * & 
                             cartesians(which_cartesian_gradient(i_coords, i_coeff, lm_index), i_l-1, i_atom)
                     end if

                     ! diagonal terms for hessian
                     i_counter = i_counter + 1
                     if (l(i_coords, which_cartesian(i_coeff, lm_index), i_l) .gt. 1) then
                        sum_hessian(i_counter, lm_index, i_atom) = sum_hessian(i_counter, lm_index, i_atom) + &
                             l(i_coords, which_cartesian(i_coeff, lm_index), i_l) * & 
                             (l(i_coords, which_cartesian(i_coeff, lm_index), i_l) - 1) * &
                             coeff(i_coeff, lm_index) * & 
                             cartesians(which_cartesian_hessian(i_counter, i_coeff, lm_index), i_l-2, i_atom)
                     end if

                     ! non-diagonal terms for hessian
                     do i_coords_2 = i_coords + 1, 3, 1
                        i_counter = i_counter + 1
                        if ((l(i_coords, which_cartesian(i_coeff, lm_index), i_l) .gt. 0) .and. &
                             (l(i_coords_2, which_cartesian(i_coeff, lm_index), i_l) .gt. 0)) then
                           sum_hessian(i_counter, lm_index, i_atom) = sum_hessian(i_counter, lm_index, i_atom) + &
                                l(i_coords, which_cartesian(i_coeff, lm_index), i_l) * &
                                l(i_coords_2, which_cartesian(i_coeff, lm_index), i_l) * &
                                coeff(i_coeff, lm_index) * & 
                                cartesians(which_cartesian_hessian(i_counter, i_coeff, lm_index), i_l-2, i_atom)
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do
    end subroutine evaluate_cartesian_hessian_and_gradient_terms_p2
!******
  !----------------------------------------------------------------------------
  !****s* cartesian_ylm/generate_cartesian_rotation
  !  NAME
  !    generate_cartesian_rotation
  !  SYNOPSIS

  subroutine generate_cartesian_rotation(l_min, l_max, rotmat, cartrotmats)

    !  PURPOSE
    !
    !    Generate Cartesian rotation matrix
    !       [rotmat rvec]^l == cartrot * [rvec]^powers
    !    where [rvec]^powers is the vector of all powers of a random realspace
    !    vector.
    !
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: l_min, l_max
    real*8, intent(IN) :: rotmat(3, 3)
    real*8, intent(OUT) :: cartrotmats(n_max_cartesian, n_max_cartesian, &
    &                                  l_min:l_max)

    !  INPUTS
    !    o l_max -- Generate rotation matrices for powers up to ...
    !    o rotmat -- Rotation for which to generate ...
    !  OUTPUTS
    !    o cartrotmats -- Rotation matrix for cartesian terms
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    real*8 :: testmat(3, 3), fac
    integer :: i, i_l
    integer :: i_cart, j_cart
    integer :: l_lhs(3), l_rhs(3)
    integer :: l11, l12, l13, l21, l22, l23, l31, l32, l33
    real*8, external :: factorial
    character(*), parameter :: func = 'generate_cartesian_rotation'

    ! Check orthogonality
    testmat = matmul(rotmat, transpose(rotmat))
    do i = 1, 3
       testmat(i, i) = testmat(i, i) - 1.d0
    end do
    if (any(abs(testmat) > 1d-10)) call aims_stop('Not unitary', func)

    ! Generate cartesian rotation matrix
    do i_l = l_min, l_max
       do i_cart = 1, n_tree(i_l)
          l_lhs = l(:, i_cart, i_l)
          do j_cart = 1, n_tree(i_l)
             l_rhs = l(:, j_cart, i_l)
             call generate_cartesian_rotation_element(l_lhs, l_rhs, rotmat, &
             &                                cartrotmats(i_cart, j_cart, i_l))
          end do
       end do
    end do

  end subroutine generate_cartesian_rotation
  !******
  !----------------------------------------------------------------------------
  !****s* cartesian_ylm/generate_cartesian_rotation_element
  !  NAME
  !    generate_cartesian_rotation_element
  !  SYNOPSIS

  subroutine generate_cartesian_rotation_element(l_lhs, l_rhs, rotmat, element)

    !  PURPOSE
    !
    !    Generate one single Cartesian rotation matrix element.
    !
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: l_lhs(3), l_rhs(3)
    real*8, intent(IN) :: rotmat(3, 3)
    real*8, intent(OUT) :: element

    !  INPUTS
    !    o l_lhs, l_rhs -- Powers (l_x, l_y, l_z) of left/right hand side
    !    o rotmat -- Rotation for which to generate ...
    !  OUTPUTS
    !    o element -- Rotation matrix element for cartesian terms
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: l11, l12, l13, l21, l22, l23, l31, l32, l33
    real*8 :: fac
    real*8, external :: factorial
    character(*), parameter :: func = 'generate_cartesian_rotation_element'

    if (sum(l_lhs) /= sum(l_rhs)) call aims_stop('sum rule violation', func)

    ! JW: The actual formula for the matrix elements is ... involved.
    ! However, the central point, here, is to write out the expression for
    ! [rotmat rvec]^lhs == product(matmul(rotmat, rvec)**lhs) and use the
    ! multinomial expansion
    !     (a1 + a2 + a3)**n
    !       = \sum_{l1+l2+l3=n} n!/(l1! l2! l3!) a1^l1 a2^l2 a3^l3.
    ! The rest is straight-forward if tedious administration of terms.

    ! Loop over all combinations of lij with
    !   sum_j lij == l_lhs(i)  .and. sum_i lij == l_rhs(j)

    element = 0.d0
    do l11 = 0, min(l_lhs(1), l_rhs(1))
       do l12 = 0, min(l_lhs(1)-l11, l_rhs(2))
          l13 = l_lhs(1)-l11-l12
          do l21 = 0, min(l_lhs(2), l_rhs(1)-l11)
             l31 = l_rhs(1)-l11-l21
             do l22 = 0, min(l_lhs(2)-l21, l_rhs(2)-l12)
                l23 = l_lhs(2)-l21-l22
                l32 = l_rhs(2)-l12-l22
                l33 = l_lhs(3)-l31-l32
                if (l33 < 0) cycle
                if (l33 /= l_rhs(3)-l13-l23) call aims_stop('l mismatch', func)
                if (l11+l12+l13 /= l_lhs(1)) call aims_stop('1x mis', func)
                if (l21+l22+l23 /= l_lhs(2)) call aims_stop('2x mis', func)
                if (l31+l32+l33 /= l_lhs(3)) call aims_stop('3x mis', func)
                if (l11+l21+l31 /= l_rhs(1)) call aims_stop('x1 mis', func)
                if (l12+l22+l32 /= l_rhs(2)) call aims_stop('x2 mis', func)
                if (l13+l23+l33 /= l_rhs(3)) call aims_stop('x3 mis', func)

                fac = 1.d0
                fac = fac * factorial(l_lhs(1)) &
                &         * rotmat(1,1)**l11 / factorial(l11) &
                &         * rotmat(1,2)**l12 / factorial(l12) &
                &         * rotmat(1,3)**l13 / factorial(l13)
                fac = fac * factorial(l_lhs(2)) &
                &         * rotmat(2,1)**l21 / factorial(l21) &
                &         * rotmat(2,2)**l22 / factorial(l22) &
                &         * rotmat(2,3)**l23 / factorial(l23)
                fac = fac * factorial(l_lhs(3)) &
                &         * rotmat(3,1)**l31 / factorial(l31) &
                &         * rotmat(3,2)**l32 / factorial(l32) &
                &         * rotmat(3,3)**l33 / factorial(l33)
                element = element + fac
                !write(0,"('add',6I3,':',F10.6,'  ->',I6,'|',3(3I2,', '))") &
                !& l_lhs, l_rhs, fac, i_entry, &
                !& l11, l21, l31, l12, l22, l32, l13, l23, l33
             end do
          end do
       end do
    end do

  end subroutine generate_cartesian_rotation_element
  !******
  !----------------------------------------------------------------------------
  !****s* cartesian_ylm/generate_ylm_rotation
  !  NAME
  !    generate_ylm_rotation
  !  SYNOPSIS

  subroutine generate_ylm_rotation(L, rotmat, l_ylm_max, ylm_rotmat)

    !  PURPOSE
    !
    !    Generate rotation with Ylm(rotmat rvec) == ylm_rotmat Ylm(rvec).
    !
    !  USES

    use numerical_utilities
    use mpi_tasks, only: check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: L
    real*8, intent(IN) :: rotmat(3, 3)
    integer, intent(IN) :: l_ylm_max
    real*8, intent(OUT) :: ylm_rotmat(-l_ylm_max:l_ylm_max, &
    &                                 -l_ylm_max:l_ylm_max)

    !  INPUTS
    !    o L -- angular momentum
    !    o rotmat -- rotation matrix in R^3
    !    o l_ylm_max -- dimension of ylm_rotmat (needs to be >= L)
    !  OUTPUTS
    !    o ylm_rotmat -- rotation matrix for 'Y_lm's.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: M, i_lm, i_coeff, i_cart
    real*8, allocatable :: cart_from_ylm(:,:)
    real*8, allocatable :: ylm_from_cart(:,:)
    real*8, allocatable :: cartrotmat(:,:)
    integer :: info
    character(*), parameter :: func = 'generate_ylm_rotation'

    ! JW: The idea of this subroutine is to generate a matrix which first
    ! transfers a covariant ylm vector (as returned e.g. by ylm_real) into
    ! cartesian terms \sum_ijk x^k y^j z^k, to rotate these and transfer it
    ! back.

    allocate(ylm_from_cart(-L:L, n_tree(L)), stat=info)
    call check_allocation(info, 'ylm_from_cart', func)

    ! Set up explicit cart -> ylm transfer matrix
    ylm_from_cart = 0.d0
    do M = -L, L
       i_lm = index_lm(M, L)
       do i_coeff = 1, n_cartesian(i_lm)
          i_cart = which_cartesian(i_coeff, i_lm)
          ylm_from_cart(M, i_cart) = ylm_from_cart(M, i_cart) &
          &                        + coeff(i_coeff, i_lm)
       end do
    end do

    ! Obtain explicit ylm -> cart transfer matrix
    allocate(cart_from_ylm(n_tree(L), -L:L), stat=info)
    call check_allocation(info, 'ylm_from_cart', func)
    ! JW: This step is (strictly speaking) not uniquely defined; however, this
    ! does not really matter as long as the back-transform gives back the
    ! original vector, which the pseudo-inverse ensures.
    call pseudo_inverse(func, 2*L+1, n_tree(L), &
    &                   ylm_from_cart, cart_from_ylm, 1d-10)

    ! Obtain cartesian rotation
    allocate(cartrotmat(n_max_cartesian, n_max_cartesian), stat=info)
    call check_allocation(info, 'cartrotmat', func)
    call generate_cartesian_rotation(L, L, rotmat, cartrotmat)

    ! Matrix multiply: (cart->ylm) cartrotmat (cart->ylm)^(-1)
    ylm_rotmat = 0.d0
    ylm_rotmat(-L:L, -L:L) &
    & = matmul(ylm_from_cart, &
    &          matmul(cartrotmat(1:n_tree(L), 1:n_tree(L)), cart_from_ylm))
    deallocate(ylm_from_cart, cart_from_ylm, cartrotmat)

  end subroutine generate_ylm_rotation
  !******
  !----------------------------------------------------------------------------
  !****s* cartesian_ylm/test_cartesian_rotation
  !  NAME
  !    test_cartesian_rotation
  !  SYNOPSIS

  subroutine test_cartesian_rotation(rotmat, rvec)

    !  PURPOSE
    !
    !    Test generate_cartesian_rotation() for the given rotmat and a test
    !    vector generated by evaluating the cartesian powers at rvec(1:3).
    !    The module is assumed to be initialized.  All powers up to
    !    max_l_prepared are tested.
    !
    !  USES

    use debug_output
    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: rotmat(3, 3)
    real*8, intent(IN) :: rvec(3)

    !  INPUTS
    !    o rotmat -- Rotation matrix to test
    !    o rvec -- Realspace position to generate test vector
    !  OUTPUTS
    !    none [bails out on error]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: l_max, i_l, n_cart
    real*8 :: rotrvec(3)
    real*8, allocatable :: cartrotmats(:,:,:)
    real*8, allocatable :: cartvec(:,:), refcartvec(:,:), rescartvec(:,:)
    real*8, allocatable :: tmpvecs(:,:)
    integer :: info
    character(*), parameter :: func = 'test_cartesian_rotation'

    l_max = max_l_prepared

    ! Generate test vector
    allocate(cartvec(n_max_cartesian, 0:l_max), stat=info)
    call check_allocation(info, 'cartvec', func)
    cartvec = 0.d0
    call evaluate_onecenter_cartesians(rvec, l_max, cartvec)
    ! Generate reference result
    allocate(refcartvec(n_max_cartesian, 0:l_max), stat=info)
    call check_allocation(info, 'refcartvec', func)
    refcartvec = 0.d0
    rotrvec = matmul(rotmat, rvec)
    call evaluate_onecenter_cartesians(rotrvec, l_max, refcartvec)

    ! Generate rotation
    allocate(cartrotmats(n_max_cartesian, n_max_cartesian, 0:l_max), stat=info)
    call check_allocation(info, 'cartrotmats', func)
    call generate_cartesian_rotation(0, l_max, rotmat, cartrotmats)

    ! Apply rotation
    allocate(rescartvec(n_max_cartesian, 0:l_max), stat=info)
    call check_allocation(info, 'rescartvec', func)
    rescartvec = 0.d0
    do i_l = 0, l_max
       n_cart = n_tree(i_l)
       rescartvec(1:n_cart, i_l) &
       & = matmul(cartrotmats(1:n_cart, 1:n_cart, i_l), &
       &          cartvec(1:n_cart, i_l))
    end do

    ! --- Compare result

    if (any(abs(rescartvec - refcartvec) > 1d-10)) then
       call debug_array(0, rotmat, 'rot: ', 'F10.6')
       call debug_array(0, reshape(rvec, (/3, 1/)), 'vec: ', 'F10.6')
       do i_l = 0, l_max
          n_cart = n_tree(i_l)
          allocate(tmpvecs(n_cart, 3))
          tmpvecs(:, 1) = cartvec(1:n_cart, i_l)
          tmpvecs(:, 2) = refcartvec(1:n_cart, i_l)
          tmpvecs(:, 3) = rescartvec(1:n_cart, i_l)
          write(0, "(5X,3A10,'; l =',I5)") 'orig', 'ref', 'res', i_l
          call debug_array(0, tmpvecs, 'vec: ', 'F10.6')
          deallocate(tmpvecs)
       end do
       call aims_stop('Test failed', func)
    end if

    deallocate(cartvec, refcartvec, rescartvec)
    deallocate(cartrotmats)

  end subroutine test_cartesian_rotation
  !******
  !----------------------------------------------------------------------------
  !****s* cartesian_ylm/test_ylm_rotation
  !  NAME
  !    test_ylm_rotation
  !  SYNOPSIS

  subroutine test_ylm_rotation(rotmat, rvec)

    !  PURPOSE
    !
    !    Test generate_ylm_rotation() for the given rotmat and a test
    !    vector generated by evaluating the cartesian powers at rvec(1:3).
    !    The module is assumed to be initialized.  All powers up to
    !    max_l_prepared are tested.
    !
    !  USES

    use debug_output
    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: rotmat(3, 3)
    real*8, intent(IN) :: rvec(3)

    !  INPUTS
    !    o rotmat -- Rotation matrix to test
    !    o rvec -- Realspace position to generate test vector
    !  OUTPUTS
    !    none [bails out on error]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: l_max, L, M, i_lm, n_cart, n_ylm
    real*8 :: rotrvec(3)
    real*8, allocatable :: ylmrotmats(:,:,:)
    real*8, allocatable :: ylm(:), refylm(:), resylm(:)
    real*8, allocatable :: tmpvecs(:,:)
    integer :: info
    character(*), parameter :: func = 'test_ylm_rotation'

    l_max = max_l_prepared
    n_ylm = (l_max+1)**2
    ! Generate test vector
    allocate(ylm(n_ylm), stat=info)
    call check_allocation(info, 'ylm', func)
    call ylm_real(rvec, l_max, ylm)
    ! Generate reference result
    allocate(refylm(n_ylm), stat=info)
    call check_allocation(info, 'refylm', func)
    rotrvec = matmul(rotmat, rvec)
    call ylm_real(rotrvec, l_max, refylm)

    ! Generate rotation
    allocate(ylmrotmats(-l_max:l_max, -l_max:l_max, 0:l_max), stat=info)
    call check_allocation(info, 'ylmrotmats', func)
    ylmrotmats = 0.d0
    do L = 0, l_max
       call generate_ylm_rotation(L, rotmat, l_max, ylmrotmats(:,:,L))
    end do

    ! Apply rotation
    allocate(resylm(n_ylm), stat=info)
    call check_allocation(info, 'resylm', func)
    resylm = 0.d0
    do L = 0, l_max
       resylm(L**2+1:(L+1)**2) &
       & = matmul(ylmrotmats(-L:L, -L:L, L), ylm(L**2+1:(L+1)**2))
    end do

    ! --- Compare result

    if (any(abs(resylm - refylm) > 1d-10)) then
       call debug_array(0, rotmat, 'rot: ', 'F10.6')
       call debug_array(0, reshape(rvec, (/3, 1/)), 'vec: ', 'F10.6')
       do L = 0, l_max
          do M = -L, L
             i_lm = L**2 + L + M
             write(0,"(5X,2I3,3F10.6)") L, M, &
             & ylm(i_lm), refylm(i_lm), resylm(i_lm)
          end do
       end do
       call aims_stop('Test failed', func)
    end if
    deallocate(ylm, refylm, resylm)
    deallocate(ylmrotmats)

  end subroutine test_ylm_rotation
  !******

!-------------------------------------------------------------------------------
!****s* cartesian_ylm/real_cartesian_coefficient
!  NAME
!    real_cartesian_coefficient
!  SYNOPSIS

real*8 elemental function real_cartesian_coefficient(m, lx, ly, lz)

!  PURPOSE
!    Calculate the Condon-Shortley normalized coefficients that represents
!    the contribution of the cartesian function
!           lx     ly     lz
!          x   *  y   *  z
!    to the spherical harmonics function
!           l
!          r  *  Ylm
!    where         ________________
!               _ / 2     2     2
!          r =   V x  +  y  +  z
!    
!          l = lx + ly + lz
!    
!    To obtain the real valued spherical harmonics, we can use the following
!    linear combinations of complex valued spherical harmonics Y^l_m:
!    
!              /  ___     m      |m| 
!              | V 2  (-1)  Im[ Y   ]  ,  m < 0
!              |                  l
!              |
!              |                 0
!       Ylm = <                 Y      ,  m = 0
!              |                 l
!              |
!              |  ___     m      m 
!              | V 2  (-1)  Re[ Y   ]  ,  m < 0
!              |                 l
!              \
!    
!    which means that only the real or complex valued part are non-zero
!    (or none of them).
!    
!    For more theoretical details please confer
!    [1] Schlegel, H. Bernhard and Frisch, Michael J., 
!        "Transformation Between Cartesian and Pure Spherical Harmonic Gaussians",
!        International Journal of Quantum Chemistry, 54, 83-87 (1995)
!
!  USES
   use constants, only: pi4
   implicit none
!  ARGUMENTS
   integer, intent(in) :: m, lx, ly, lz
!  INPUTS
!   o m -- quantum number m for spherical harmonics
!   o lx, ly, lz -- exponents of the cartesian coordinates x,y, and z
!                   this routine assumes that quantum number l == lx + ly + lz
!  OUTPUT
!   o real_cartesian_coefficient -- see description of function
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   interface
      real*8 elemental function binomial_coefficient(i,j)
         integer, intent(in) :: i,j
      end function
      real*8 elemental function factorial(n)
         integer, intent(in) :: n
      end function
   end interface

   real*8, parameter :: ZERO=0.d0, SQRT2=sqrt(2.d0)

   integer :: abs_m, l, j2, j, sign_i, sign_k, i, k
   real*8 :: sum_i, sum_k

   abs_m = abs(m)
   j2 = int(lx+ly-abs_m)      ! j2 = 2*j
   if ( (modulo(j2, 2) .eq. 1) & ! see comment in [1], after Eq. 15
         .or. ( (m.lt.0) .and. (modulo(ly, 2) .eq. 0) ) &
         .or. ( (m.ge.0) .and. (modulo(ly, 2) .eq. 1) ) ) then
      real_cartesian_coefficient = ZERO
   else
      l = lx + ly + lz
      j = j2 / 2
      ! initial signs for alternating sums:
      sign_i = 1 ! == (-1)**i,  with i=0
      ! instead of 
      !     (-1)**((abs_m-lx+2*k)/2),  with k=0
      ! as in [1], Eq. 15, we use
      sign_k = (-1)**(ly/2-j) ! since 2*j = lx + ly - abs_m

      sum_i = ZERO
      do i = 0, (l-abs_m)/2
         sum_i = sum_i + binomial_coefficient(l,i)*binomial_coefficient(i,j) &
                           * factorial(2*l-2*i)/factorial(l-abs_m-2*i) * sign_i
         sign_i = -sign_i
      enddo ! i

      sum_k = ZERO
      do k = 0, j
         sum_k = sum_k + sign_k * &
                   binomial_coefficient(j,k)*binomial_coefficient(abs_m,lx-2*k)
         sign_k = -sign_k
      enddo

      ! The normalization prefactor is the one of [1], Eq. 7
      real_cartesian_coefficient = sqrt( (2*l+1) * factorial(l-abs_m) / &
                                        ( pi4    * factorial(l+abs_m) ) ) / &
                                    (2**l*factorial(l)) * sum_i * sum_k
      ! apply prefactor for linear combination, see comment above
      ! the prefactor (-1)**m is already included
      if (abs_m .ne. 0) then
         real_cartesian_coefficient = real_cartesian_coefficient * SQRT2
      endif
   endif
end function real_cartesian_coefficient


!-------------------------------------------------------------------------------
!****s* cartesian_ylm/aims_real_cartesian_coefficient
!  NAME
!    aims_real_cartesian_coefficient
!  SYNOPSIS

real*8 elemental function aims_real_cartesian_coefficient(m, lx, ly, lz)

!  PURPOSE
!    Calculate the "FHI-aims convention" normalized coefficients that represents
!    the contribution of the cartesian function
!           lx     ly     lz
!          x   *  y   *  z
!    to the spherical harmonics function
!           l
!          r  *  Ylm
!    where         ________________
!               _ / 2     2     2
!          r =   V x  +  y  +  z
!    
!          l = lx + ly + lz
!    
!    The difference between the FHI-aims and the Condon-Shortley convention
!    is that spherical harmonics with m > 0 get a prefactor of
!              m
!          (-1)
!    
!    The author disproves of this convention, but as long as ONE
!    convention is used consistently, results remain unchanged.
!  USES
   implicit none
!  ARGUMENTS
   integer, intent(in) :: m, lx, ly, lz
!  INPUTS
!   o m -- quantum number m for spherical harmonics
!   o lx, ly, lz -- exponents of the cartesian coordinates x,y, and z
!                   this routine assumes that quantum number l == lx + ly + lz
!  OUTPUT
!   o aims_real_cartesian_coefficient -- see description of function
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   aims_real_cartesian_coefficient = real_cartesian_coefficient(m, lx, ly, lz)
   if (m .gt. 0) then
      aims_real_cartesian_coefficient = aims_real_cartesian_coefficient*(-1)**m
   endif
end function aims_real_cartesian_coefficient


end module cartesian_ylm
