! module cartesian_ylm encapsulates everything for expanding ylm into
! cartesian gaussians
!
! R.Gehrke (2006)
!

module cartesian_ylm

  use cluster
 
  implicit none

  private


  integer, dimension(:,:,:), allocatable :: l ! exponent l_x, l_y, l_z in cartesians
  real*8, dimension(:,:), allocatable :: coeff ! coefficients of expansion of ylms into cartesians
  integer, dimension(:), allocatable :: n_tree ! number of all possible cartesians for a given l
  integer, dimension(:,:), allocatable :: tree ! describes how cartesians for l are calculated from the cartesians for l-1
  integer, dimension(:,:), allocatable :: which_cartesian ! maps coefficients to the corresponding cartesians
  integer, dimension(:,:,:), allocatable :: which_cartesian_gradient ! maps coefficients to the corresponding cartesians for gradient terms
  integer, dimension(:,:,:), allocatable :: which_cartesian_hessian ! maps coefficients to the corresponding cartesians for hessian terms
  integer, dimension(:), allocatable :: n_cartesian ! number of all cartesians with non-zero coefficients
  integer, dimension(:,:), allocatable, public :: index_lm
  real*8, private :: pi
  integer, public :: n_max_cartesian
  real*8, private :: pi4_inv
  integer, public :: l_wave_max
  real*8, private :: coeff_thr

  parameter (l_wave_max = 6)
  parameter (coeff_thr = 1d-8)

  public :: initialize_cartesian_ylm, cleanup_cartesian_ylm, evaluate_cartesians, & 
       tab_ylm_cartesian
  
  contains

    subroutine initialize_cartesian_ylm()

      ! local variables
      integer :: m_abs
      real*8 :: j_in_coefficients
      real*8 :: norm_all
      real*8 :: norm_all_l
      real*8 :: norm_all_lm
      real*8 :: exponent
      real*8 :: exponent_half
      real*8 :: temp_coeff
      real*8 :: sum_one
      real*8 :: sum_two
      integer, dimension(3) :: target_l

      ! functions
      real*8 :: binomial_coefficient
      real*8 :: factorial
      
      ! counter
      integer :: i_coords
      integer :: i_coords_2
      integer :: i_coords_3
      integer :: i_tree
      integer :: i_tree_2
      integer :: i_coeff
      integer :: i_i
      integer :: i_j
      integer :: i_k
      integer :: i_l
      integer :: i_m
      integer :: i_index
      integer :: i_counter

      pi = (4.0D+0)*ATAN(1.0D+0)
      pi4_inv = 1. / (4 * pi)
      ! n_max_cartesians is the number of combinations with repetition binomial_coefficient(n + k - 1    k)
      n_max_cartesian = binomial_coefficient(3 + l_wave_max - 1, l_wave_max)
      if (.not.allocated(l)) then
         allocate(l(3, n_max_cartesian, 0:l_wave_max))
      end if
      if (.not.allocated(coeff)) then
         allocate(coeff(n_max_cartesian, (l_wave_max+1) ** 2))
      end if
      if (.not.allocated(index_lm)) then
         allocate(index_lm(-l_wave_max:l_wave_max, 0:l_wave_max))
      end if
      if (.not.allocated(n_tree)) then
         allocate(n_tree(0:l_wave_max))
      end if
      if (.not.allocated(tree)) then
         allocate(tree(n_max_cartesian, 0:l_wave_max))
      end if
      if (.not.allocated(which_cartesian)) then
         allocate(which_cartesian(n_max_cartesian, (l_wave_max+1) ** 2))
      end if
      if (.not.allocated(which_cartesian_gradient)) then
         allocate(which_cartesian_gradient(3, n_max_cartesian, (l_wave_max+1) ** 2))
      end if
      if (.not.allocated(which_cartesian_hessian)) then
         allocate(which_cartesian_hessian(6, n_max_cartesian, (l_wave_max+1) ** 2))
      end if
      if (.not.allocated(n_cartesian)) then
         allocate(n_cartesian((l_wave_max+1) ** 2))
      end if
      i_index = 0
      do i_l = 0, l_wave_max, 1
         do i_m = -i_l, i_l
            i_index = i_index + 1
            index_lm(i_m, i_l) = i_index
         enddo
      enddo
      ! initialize lookup-tables tree, n_tree, l(:)
      n_tree(0) = 1
      l(:, 1, 0) = 0
      if (l_wave_max .gt. 0) then
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
         do i_l = 2, l_wave_max, 1
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
      norm_all = 1. / sqrt(8*pi)
      do i_l = 0, l_wave_max, 1
         
         norm_all_l = sqrt(2.d0 * (2.d0 * i_l + 1.d0)) / (2**(i_l) * factorial(i_l)) 
         do i_m = - i_l, i_l, 1
            
            m_abs = abs(i_m)
            norm_all_lm = sqrt(factorial(i_l - m_abs) / factorial(i_l + m_abs))
            
            i_coeff = 0
            do i_tree = 1, n_tree(i_l), 1
               
               sum_one = 0.d0
               sum_two = 0.d0
               
               j_in_coefficients = dble(l(1, i_tree, i_l) + l(2, i_tree, i_l) - m_abs) / 2
               if (abs(j_in_coefficients - int(j_in_coefficients)) .lt. 1d-8) then
                  
                  if (j_in_coefficients .ge. 0) then
                     i_j = int(j_in_coefficients + 1d-8)
                  else
                     i_j = int(j_in_coefficients - 1d-8)
                  end if
                  do i_i = 0, (i_l - m_abs) / 2, 1
                     sum_one = sum_one + binomial_coefficient(i_l, i_i) * binomial_coefficient(i_i, i_j) * &
                          (-1) ** i_i * factorial(2*i_l - 2*i_i) / factorial(i_l - m_abs - 2 * i_i)
                  end do
                  ! VB: RALF: DID I FIX THIS CORRECTLY, REPLACING j_in_coefficients by i_j ???
                  do i_k = 0, i_j, 1
                     exponent = m_abs - l(1, i_tree, i_l) + 2 * i_k
                     exponent_half = exponent / 2
                     if ((abs(exponent_half - int(exponent_half)) .gt. 1d-8) .and. (i_m .lt. 0)) then
                        sum_two = sum_two + binomial_coefficient(i_j, i_k) * &
                             binomial_coefficient(m_abs, l(1, i_tree, i_l) - 2 * i_k) * &
                             (-1) ** int((exponent - 1) / 2)
                     end if
                     if ((abs(exponent_half - int(exponent_half)) .lt. 1d-8) .and. (i_m .ge. 0)) then
                        sum_two = sum_two + binomial_coefficient(i_j, i_k) * &
                             binomial_coefficient(m_abs, l(1, i_tree, i_l) - 2 * i_k) * &
                             (-1) ** int(exponent / 2)
                     end if
                  end do
                  temp_coeff = norm_all * norm_all_l * norm_all_lm * sum_one * sum_two
                  if (abs(temp_coeff) .gt. coeff_thr) then
                     i_coeff = i_coeff + 1
                     coeff(i_coeff, index_lm(i_m, i_l)) = temp_coeff
                     ! make assignment in which_cartesian
                     which_cartesian(i_coeff, index_lm(i_m, i_l)) = i_tree
                     if (i_m .ne. 0) then
                        coeff(i_coeff, index_lm(i_m, i_l)) = coeff(i_coeff, index_lm(i_m, i_l)) * sqrt(2.d0)
                     end if
                     if (i_m .gt. 0) then
                        coeff(i_coeff, index_lm(i_m, i_l)) = coeff(i_coeff, index_lm(i_m, i_l)) * ((-1) ** i_m)
                     end if

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

                  end if
               end if
            end do
            n_cartesian(index_lm(i_m, i_l)) = i_coeff
         end do
      end do

    end subroutine initialize_cartesian_ylm
    
    subroutine cleanup_cartesian_ylm()
      
      if (allocated(l)) then
         deallocate(l)
      end if
      if (allocated(coeff)) then
         deallocate(coeff)
      end if
      if (allocated(index_lm)) then
         deallocate(index_lm)
      end if
      if (allocated(n_tree)) then
         deallocate(n_tree)
      end if
      if (allocated(tree)) then
         deallocate(tree)
      end if
      if (allocated(which_cartesian)) then
         deallocate(which_cartesian)
      end if
      if (allocated(which_cartesian_gradient)) then
         deallocate(which_cartesian_gradient)
      end if
      if (allocated(which_cartesian_hessian)) then
         deallocate(which_cartesian_hessian)
      end if
      if (allocated(n_cartesian)) then
         deallocate(n_cartesian)
      end if

    end subroutine cleanup_cartesian_ylm

    subroutine evaluate_cartesians(dir_tab, cartesians)
      
      ! imported variables

      ! input
      real*8, dimension(3, n_atoms), intent(in) :: dir_tab

      ! output
      real*8, dimension(n_max_cartesian, 0:l_wave_max, n_atoms), intent(out) :: cartesians

      ! local variables
      integer :: current_atom

      ! counter
      integer :: i_coords
      integer :: i_l
      integer :: i_tree
      integer :: i_counter
      integer :: i_atom

      ! evaluate all cartesian gaussians with lookup-tables
      do i_atom = 1, n_atoms, 1

         cartesians(1, 0, i_atom) = 1.d0
         
         do i_coords = 1, 3, 1
            cartesians(i_coords, 1, i_atom) = dir_tab(i_coords, i_atom)
         end do
         do i_l = 2, l_wave_max, 1
            i_counter = 1
            do i_tree = 1, n_tree(i_l), 1
               cartesians(i_tree, i_l, i_atom) = cartesians(i_counter, i_l - 1, i_atom) * dir_tab(tree(i_tree, i_l), i_atom)
               if (tree(i_tree, i_l) .eq. 3) then
                  i_counter = i_counter + 1
               end if
            end do
         end do
      end do
      
    end subroutine evaluate_cartesians
    
    subroutine tab_ylm_cartesian(cartesians, l_max, ylm_tab)

      ! evaluates ylm_functions assuming that cartesians are already evaluated !!!
      ! -> subroutine evaluate_cartesians()

      ! imported variables

      ! input 
      real*8, dimension(n_max_cartesian, 0:l_wave_max, n_atoms), intent(in) :: cartesians
      integer, intent(in) :: l_max

      ! output
      real*8, dimension((l_max)*(l_max+1), n_atoms), intent(out) :: ylm_tab

      ! local variables
      integer :: lm_index
      integer :: current_atom

      ! counter
      integer :: i_l
      integer :: i_m
      integer :: i_coeff
      integer :: i_atom

      ! evaluate spherical harmonic with lookup-tables
      do i_atom = 1, n_atoms, 1

         ylm_tab(:, i_atom) = 0.d0
         ylm_tab(1, i_atom) = sqrt(pi4_inv)
         do i_l = 1, l_max, 1
            do i_m = -i_l, i_l, 1

               lm_index = index_lm(i_m, i_l)
               do i_coeff = 1, n_cartesian(lm_index), 1
                  ylm_tab(lm_index, i_atom) = ylm_tab(lm_index, i_atom) + &
                       coeff(i_coeff, lm_index) * cartesians(which_cartesian(i_coeff, lm_index), i_l, i_atom)
!                  ylm_tab(lm_index, i_atom) = ylm_tab(lm_index, i_atom) + cartesians(which_cartesian(i_coeff, lm_index), i_l, i_atom)
!                  write (6,*) i_atom, i_l, i_m, ylm_tab(lm_index, i_atom) 
               end do
            end do
         end do
      end do

    end subroutine tab_ylm_cartesian

end module cartesian_ylm
