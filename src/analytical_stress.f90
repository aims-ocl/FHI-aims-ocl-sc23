!****h* FHI-aims/analytical_stress
!  NAME
!   analytical_stress 
!  SYNOPSIS
module analytical_stress

  !  PURPOSE
  !    This module contains all the routines for 
  !    calculating the analytical stress
  !  USES
  implicit none

  ! CC: Moved here from physics, so that module is accessible without circular references to physics
  real*8,dimension(1:3,1:3)             :: analytical_stress_tensor

  !##################################################################################
  ! ### Variables for sum_up_whole_potential:
  !! Dim = n_atoms
  real*8,dimension(:),allocatable       :: AS_v_at_f
  real*8,dimension(:),allocatable       :: AS_v_at_n
  real*8,dimension(:),allocatable       :: AS_v_at_r
  real*8,dimension(:),allocatable       :: AS_v_at_k
  ! Respective temp. vars
  real*8                                :: AS_v_at_f_temp
  real*8                                :: AS_v_at_n_temp
  real*8                                :: AS_v_at_r_temp
  !real*8                                :: AS_v_at_k_temp = delta_v_hartree_aux

  ! Respective strain derivatives 
  real*8,dimension(:,:,:),allocatable   :: AS_dde_v_at_f
  real*8,dimension(:,:,:),allocatable   :: AS_dde_v_at_n
  real*8,dimension(:,:,:),allocatable   :: AS_dde_v_at_r
  real*8,dimension(:,:,:),allocatable   :: AS_dde_v_at_k
  ! Respective temp. vars
  real*8,dimension(1:3,1:3)             :: AS_dde_v_at_f_temp
  real*8,dimension(1:3,1:3)             :: AS_dde_v_at_n_temp
  real*8,dimension(1:3,1:3)             :: AS_dde_v_at_r_temp
  real*8,dimension(1:3,1:3)             :: AS_dde_v_at_k_temp
  
  ! Dim = n_full_points
  real*8,dimension(:),allocatable       :: AS_rho_f   ! Should be = rho_free_superpos --> Shoudl be replaced accordingly
  real*8,dimension(:),allocatable       :: AS_v_rho_f ! Should be = v_hartree_free --> Shoudl be replaced accordingly
  real*8,dimension(:),allocatable       :: AS_v_rho_n
  real*8,dimension(:),allocatable       :: AS_v_rho_r
  real*8,dimension(:),allocatable       :: AS_v_rho_k
  ! Respective temp. vars
  real*8                                :: AS_rho_f_temp
  real*8                                :: AS_v_rho_f_temp
  real*8                                :: AS_v_rho_n_temp
  real*8                                :: AS_v_rho_r_temp
  real*8                                :: AS_v_rho_k_temp
  ! Respective strain derivatives 
  real*8,dimension(:,:,:),allocatable   :: AS_dde_v_rho_f
  real*8,dimension(:,:,:),allocatable   :: AS_dde_v_rho_n
  real*8,dimension(:,:,:),allocatable   :: AS_dde_v_rho_r
  real*8,dimension(:,:,:),allocatable   :: AS_dde_v_rho_k
  real*8,dimension(:,:,:),allocatable,target   :: AS_dde_potential ! = AS_dde_v_rho_f + AS_dde_v_rho_n + AS_dde_v_rho_r + AS_dde_v_rho_k - AS_dde_v_rho_f_a
  ! Respective temp. vars
  real*8,dimension(1:3,1:3)             :: AS_dde_v_rho_f_temp
  real*8,dimension(1:3,1:3)             :: AS_dde_v_rho_n_temp
  real*8,dimension(1:3,1:3)             :: AS_dde_v_rho_r_temp
  real*8,dimension(1:3,1:3)             :: AS_dde_v_rho_k_temp

  ! Derivative of the density
  real*8,dimension(:,:,:),allocatable   :: AS_dde_rho_free
  real*8,dimension(:,:,:),allocatable   :: AS_dde_rho_mp  
  ! Respective temp. vars
  real*8,dimension(1:3,1:3)             :: AS_dde_rho_free_temp
  real*8,dimension(1:3,1:3)             :: AS_dde_rho_mp_temp 

  ! Average strain derivative of potential
  real*8                                :: AS_average_v_hartree_free
  real*8                                :: AS_average_delta_v_hartree_real
  real*8,dimension(1:3,1:3)             :: AS_dde_v_rho_f_a
  real*8,dimension(1:3,1:3)             :: AS_dde_v_rho_mp_a

  !! Matrices used in hartree potential reciprocal:
  real*8,dimension(:,:,:),allocatable      :: AS_B_matrix
  complex*16,dimension(:,:,:),allocatable  :: AS_Lambda_matrix
  complex*16,dimension(:,:,:),allocatable  :: AS_hartree_coeff_stress

  !! Temporary storage in sum_up_whole_potential_p1.f90
  real*8                                :: AS_delta_v_hartree_aux_save
  real*8,dimension(1:3)                 :: AS_v_hartree_gradient_temp_save

  !! Individual integrals computed in sum_up_whole_pot....
  ! 0.5 * Z_a * (d/de v_(R_a) )
  real*8,dimension(1:3,1:3)             ::   NEW_AS_at_stress
  ! d/de ( - 0.5 \int rho_mp * v )
  real*8,dimension(1:3,1:3)             ::   NEW_AS_rho_stress
  ! Correction for on-site integrals
  ! d/de -0.5 \int sum_{i.eq.j} \rho_f_i v_f_j : Should be zero, but is not
  real*8,dimension(1:3,1:3)             ::   AS_rho_f_f_self_correction 
  ! ~ full_rho contrib:
  ! ( 0.5 \int rho * [d/de v] )
  real*8,dimension(1:3,1:3)             ::   AS_pu_rho_dde_fnrka

  !##################################################################################
  ! Pulay part - aka  Orbital dependent quantities
  real*8,dimension(1:3,1:3)             ::   AS_pulay_stress
  real*8,dimension(1:9)                 ::   AS_pulay_stress_local
  real*8,dimension(:,:,:),allocatable   ::   AS_strain_deriv_wave
  real*8,dimension(:,:),allocatable     ::   AS_strain_deriv_wave_shell
  real*8,dimension(:,:,:,:),allocatable ::   AS_jac_pot_kin_times_psi
  real*8,dimension(2)                   ::   AS_en_xc_minus_pot_xc                                      
  real*8,dimension(1:9)                 ::   AS_dde_potential_local
  ! Here come the corrections (non-rel. case only)
  real*8,dimension(:,:,:),allocatable   ::   AS_Hessian_times_psi
  real*8,dimension(:,:,:),allocatable   ::   AS_Hessian_times_psi_shell
  real*8,dimension(:,:),allocatable     ::   AS_T_times_psi         
  real*8,dimension(:,:,:),allocatable   ::   AS_T_times_psi_shell   
  real*8,dimension(1:3,1:3)             ::   AS_dde_T_sa_at_di_orb         
  real*8,dimension(1:9)                 ::   AS_dde_T_sa_at_di_orb_local   
  real*8,dimension(:,:),allocatable     ::   AS_corr_shell
  ! Respective arrays for efficient indexing of on-site terms
  integer                               ::   AS_n_on_site
  integer,dimension(:,:),allocatable    ::   AS_i_basis_on_site
  ! GGA part - involves terms with derivative of xc-functional with respect to density gradient
  real*8,dimension(:,:,:,:),allocatable ::   AS_hessian_times_xc_deriv_gga
  ! Meta GGA part - involves terms with derivative of xc-functional with respect to kinetic density
  real*8,dimension(:,:,:,:),allocatable ::   AS_hessian_times_xc_deriv_mgga
  ! Relativistic case - strain derivative of kinetic energy wave
  real*8,dimension(:,:,:),allocatable   ::   AS_strain_deriv_kinetic_wave

  !##################################################################################
  ! vdW correction (Hirschfeld)
  real*8,dimension(1:3,1:3)             ::   AS_vdw_stress

  !##################################################################################
  ! Hartree-Fock part
  real*8,dimension(1:9)                 ::   AS_EXX_stress_local
  real*8,dimension(1:3,1:3)             ::   AS_EXX_stress

  !##################################################################################
  !! *** Debug part:
  ! Debug routines for fixing multipole moments / Pulay DMe: 
  !! logical                                    :: AS_flag_write_pulay_dm   = .false.
  !! logical                                    :: AS_flag_read_pulay_dm    = .false.
  !! logical                                    :: AS_flag_write_mm_moments = .false.
  !! logical                                    :: AS_flag_read_mm_moments  = .false.
  !! real*8,dimension(:,:),allocatable          :: AS_mm_save
  !! real*8,dimension(:,:,:),allocatable        :: AS_mm_splines_save

  contains

  !##################################################################################
  ! Routines used in sum_up_whole_potential:
  !
  !! allocate_and_or_zero_all_AS_arrays
  subroutine allocate_and_or_zero_all_AS_arrays(n_points)
    use dimensions
    implicit none

    integer :: n_points

    !! Dim = n_atoms
    if (.not.allocated(AS_v_at_f))  allocate(AS_v_at_f(1:n_atoms))
    if (.not.allocated(AS_v_at_n))  allocate(AS_v_at_n(1:n_atoms))
    if (.not.allocated(AS_v_at_r))  allocate(AS_v_at_r(1:n_atoms))
    if (.not.allocated(AS_v_at_k))  allocate(AS_v_at_k(1:n_atoms))
    AS_v_at_f(:)                   = 0.0d0
    AS_v_at_n(:)                   = 0.0d0
    AS_v_at_r(:)                   = 0.0d0
    AS_v_at_k(:)                   = 0.0d0

    ! Respective strain derivatives 
    if (.not.allocated(AS_dde_v_at_f))  allocate(AS_dde_v_at_f(1:3,1:3,1:n_atoms))
    if (.not.allocated(AS_dde_v_at_n))  allocate(AS_dde_v_at_n(1:3,1:3,1:n_atoms))
    if (.not.allocated(AS_dde_v_at_r))  allocate(AS_dde_v_at_r(1:3,1:3,1:n_atoms))
    if (.not.allocated(AS_dde_v_at_k))  allocate(AS_dde_v_at_k(1:3,1:3,1:n_atoms))
    AS_dde_v_at_f(:,:,:)           = 0.0d0
    AS_dde_v_at_n(:,:,:)           = 0.0d0
    AS_dde_v_at_r(:,:,:)           = 0.0d0
    AS_dde_v_at_k(:,:,:)           = 0.0d0

    ! Dim = n_full_points
    ! AS_v_rho_f = v_hartree_free
    if (.not.allocated(AS_rho_f))    allocate(AS_rho_f(1:n_points))
    if (.not.allocated(AS_v_rho_f))  allocate(AS_v_rho_f(1:n_points))
    if (.not.allocated(AS_v_rho_n))  allocate(AS_v_rho_n(1:n_points))
    if (.not.allocated(AS_v_rho_r))  allocate(AS_v_rho_r(1:n_points))
    if (.not.allocated(AS_v_rho_k))  allocate(AS_v_rho_k(1:n_points))
    AS_rho_f(:)                  = 0.0d0
    AS_v_rho_f(:)                  = 0.0d0
    AS_v_rho_n(:)                  = 0.0d0
    AS_v_rho_r(:)                  = 0.0d0
    AS_v_rho_k(:)                  = 0.0d0

    ! Respective strain derivatives 
    if (.not.allocated(AS_dde_v_rho_f))   allocate(AS_dde_v_rho_f(1:3,1:3,1:n_points))
    if (.not.allocated(AS_dde_v_rho_n))   allocate(AS_dde_v_rho_n(1:3,1:3,1:n_points))
    if (.not.allocated(AS_dde_v_rho_r))   allocate(AS_dde_v_rho_r(1:3,1:3,1:n_points))
    if (.not.allocated(AS_dde_v_rho_k))   allocate(AS_dde_v_rho_k(1:3,1:3,1:n_points))
    if (allocated(AS_dde_potential))    deallocate(AS_dde_potential)
    if (.not.allocated(AS_dde_potential)) allocate(AS_dde_potential(1:3,1:3,1:n_points))
    AS_dde_v_rho_f(:,:,:)         = 0.0d0
    AS_dde_v_rho_n(:,:,:)         = 0.0d0
    AS_dde_v_rho_r(:,:,:)         = 0.0d0
    AS_dde_v_rho_k(:,:,:)         = 0.0d0
    AS_dde_potential(:,:,:)       = 0.0d0

    ! Derivative of the density
    if (.not.allocated(AS_dde_rho_free))  allocate(AS_dde_rho_free(1:3,1:3,1:n_points))
    if (.not.allocated(AS_dde_rho_mp  ))  allocate(AS_dde_rho_mp  (1:3,1:3,1:n_points))
    AS_dde_rho_free(:,:,:)        = 0.0d0
    AS_dde_rho_mp  (:,:,:)        = 0.0d0

    ! Average potentials  
    AS_average_v_hartree_free      = 0.0d0
    AS_average_delta_v_hartree_real= 0.0d0

    ! Average strain derivative of potential
    AS_dde_v_rho_f_a(:,:)         = 0.0d0
    AS_dde_v_rho_mp_a(:,:)        = 0.0d0

    !! N.B. : AS_B_matrix & AS_Lambda_matrix are alloctaed in hartee_potential_recip
    !! Matrices used in hartree potential reciprocal:
    !! real*8,dimension(:,:,:),allocatable      :: AS_B_matrix
    !! complex*16,dimension(:,:,:),allocatable  :: AS_Lambda_matrix

    ! 0.5 * Z_a * (d/de v(R_a) )
    NEW_AS_at_stress(:,:)         = 0.0d0
    ! d/de ( - 0.5 \int rho_mp * v )
    NEW_AS_rho_stress             = 0.0d0
    ! On-site correction
    AS_rho_f_f_self_correction(:,:)= 0.0d0
    ! ~ full_rho contrib:
    ! ( 0.5 \int rho * [d/de v] )
    AS_pu_rho_dde_fnrka(:,:)      = 0.0d0

  end subroutine allocate_and_or_zero_all_AS_arrays

  !! deallocate_all_AS_arrays
  subroutine deallocate_all_AS_arrays()
    use dimensions
    implicit none

    if (allocated(AS_v_at_f))               deallocate(AS_v_at_f)
    if (allocated(AS_v_at_n))               deallocate(AS_v_at_n)
    if (allocated(AS_v_at_r))               deallocate(AS_v_at_r)
    if (allocated(AS_v_at_k))               deallocate(AS_v_at_k)
    if (allocated(AS_dde_v_at_f))           deallocate(AS_dde_v_at_f)
    if (allocated(AS_dde_v_at_n))           deallocate(AS_dde_v_at_n)
    if (allocated(AS_dde_v_at_r))           deallocate(AS_dde_v_at_r)
    if (allocated(AS_dde_v_at_k))           deallocate(AS_dde_v_at_k)
    if (allocated(AS_rho_f))                deallocate(AS_rho_f  )
    if (allocated(AS_v_rho_f))              deallocate(AS_v_rho_f)
    if (allocated(AS_v_rho_n))              deallocate(AS_v_rho_n)
    if (allocated(AS_v_rho_r))              deallocate(AS_v_rho_r)
    if (allocated(AS_v_rho_k))              deallocate(AS_v_rho_k)
    if (allocated(AS_dde_v_rho_f))          deallocate(AS_dde_v_rho_f)
    if (allocated(AS_dde_v_rho_n))          deallocate(AS_dde_v_rho_n)
    if (allocated(AS_dde_v_rho_r))          deallocate(AS_dde_v_rho_r)
    if (allocated(AS_dde_v_rho_k))          deallocate(AS_dde_v_rho_k)
    if (allocated(AS_dde_rho_free))         deallocate(AS_dde_rho_free)
    if (allocated(AS_dde_rho_mp  ))         deallocate(AS_dde_rho_mp  )

    if (allocated(AS_B_matrix  ))           deallocate(AS_B_matrix )
    if (allocated(AS_Lambda_matrix))        deallocate(AS_Lambda_matrix)
    if (allocated(AS_hartree_coeff_stress)) deallocate(AS_hartree_coeff_stress)

  end subroutine deallocate_all_AS_arrays


  ! ###################################################################################
  ! Pulay // orbital routines
  !
  !! AS_eval_jac_pot_kin_times_psi
  ! Construct 
  !           [ delta_lm (en_xc-pot_xc) + (d/de_lm v_h) ] * |i>                          (atomic zora)
  !      or   [ delta_lm (en_xc-pot_xc) + (d/de_lm v_h) ] * |i>  + ( d/dr_l d/dr_m |i> ) (non-relativistic)
  subroutine AS_eval_jac_pot_kin_times_psi( &
    hessian_index, n_compute, en_xc_minus_pot_xc, wave, dde_hartree_pot, new_right_side, relativistic, hessian )
    implicit none

    integer,intent(in)                      :: hessian_index, n_compute
    real*8,intent(in)                       :: en_xc_minus_pot_xc
    real*8,dimension(n_compute),intent(in)  :: wave
    real*8,intent(in)                       :: dde_hartree_pot   
    real*8,dimension(n_compute),intent(out) :: new_right_side      
    logical,intent(in)                      :: relativistic
    real*8,dimension(n_compute),optional    :: hessian

    ! Jacobian
    if ( ( hessian_index .eq. 1) .or. ( hessian_index .eq. 4) .or. ( hessian_index .eq. 6) ) then
      new_right_side(1:n_compute)   = en_xc_minus_pot_xc * wave(1:n_compute)
    else
      new_right_side(1:n_compute)   = 0.0d0
    end if
    ! rho * d/de v
    if (relativistic) then
      new_right_side(1:n_compute) = new_right_side(1:n_compute) + ( dde_hartree_pot * wave(1:n_compute)) 
    else
      new_right_side(1:n_compute) = new_right_side(1:n_compute) + ( dde_hartree_pot * wave(1:n_compute)) + hessian(1:n_compute)
    end if

  end subroutine AS_eval_jac_pot_kin_times_psi

  !! AS_add_gga_jac
  ! Add to LDA Jacobian
  !     - 4 (d/dgrad(rho)^2 f_xc) * grad(rho) * grad(|i>)
  ! This term arises from en_pot_xc = integral(rho * v_xc) where in the GGA case v_xc has an additional term including two gradients.
  ! Here, we add this additional term to the LDA term (calculated by the subroutine above)
  ! and we use integration by parts to shift one gradient from v_xc to rho (wave function).
  ! There is a scalar product between grad(rho) and grad(|i>).
  ! The variable xc_gradient_deriv includes grad(rho).
  subroutine AS_add_gga_jac( &
    new_right_side, n_compute, xc_gradient_deriv_ispin_ipoint, gradient_basis_wave_ipoint )
    implicit none

    real*8,dimension(n_compute),  intent(inout) :: new_right_side
    integer,                      intent(in)    :: n_compute
    real*8,dimension(3),          intent(in)    :: xc_gradient_deriv_ispin_ipoint
    real*8,dimension(n_compute,3),intent(in)    :: gradient_basis_wave_ipoint

    ! local
    integer   :: i_coord

    do i_coord = 1, 3, 1
      new_right_side(1:n_compute) = new_right_side(1:n_compute) &
        - 4 * gradient_basis_wave_ipoint(1:n_compute,i_coord) * xc_gradient_deriv_ispin_ipoint(i_coord)
    end do

  end subroutine AS_add_gga_jac

  !! AS_construct_hessian_op_times_psi
  ! Construct
  ! 	( d/dr_l d/dr_m |i> ) - delta_lm/2 [ ( d^2/dx^2 + d^2/dy^2 + d^2/dz^2 ) |i> ]
  subroutine AS_construct_hessian_op_times_psi( &
    hessian_index, n_compute, hessian, T_times_psi, new_right_side )
    implicit none

    integer,intent(in)                      :: hessian_index
    integer,intent(in)                      :: n_compute
    real*8,dimension(n_compute),intent(in)  :: hessian
    real*8,dimension(n_compute),intent(in)  :: T_times_psi 
    real*8,dimension(n_compute),intent(out) :: new_right_side      

    if ( ( hessian_index .eq. 1) .or. ( hessian_index .eq. 4) .or. ( hessian_index .eq. 6) ) then 
      new_right_side(1:n_compute) = hessian(1:n_compute) + T_times_psi   
    else
      new_right_side(1:n_compute) = hessian(1:n_compute)
    end if

  end subroutine AS_construct_hessian_op_times_psi        

  !! AS_construct_T_times_psi
  ! Construct 
  ! 	( d^2/dx^2 + d^2/dy^2 + d^2/dz^2 ) |i>
  subroutine AS_construct_T_times_psi( &
    n_compute, H_times_psi, local_potential_parts, wave, new_right_side )
    implicit none

    integer,intent(in)                      :: n_compute
    real*8,dimension(n_compute),intent(in)  :: H_times_psi
    real*8                     ,intent(in)  :: local_potential_parts
    real*8,dimension(n_compute),intent(in)  :: wave         
    real*8,dimension(n_compute),intent(out) :: new_right_side      

    new_right_side(1:n_compute) = H_times_psi(1:n_compute) - ( local_potential_parts*wave(1:n_compute) )

  end subroutine AS_construct_T_times_psi             

  !! AS_eval_strain_deriv_wave
  ! Construct
  !   d/de_lm |i> = (d/dr_l |i>) * (r-R_m)
  subroutine AS_eval_strain_deriv_wave(       &
               n_compute_c, n_compute_atoms,  &
               dist_tab,                      &
               dir_tab,                       &
               i_basis,                       &
               i_basis_fns_inv,               &
               atom_index_inv,                &
               gradient_basis_wave,           &
               strain_index,                  &
               AS_strain_der                  )
    use pbc_lists
    use dimensions, only: n_basis_fns, n_centers
    use basis, only: basis_fn
    implicit none

    integer                                                :: n_compute_c
    integer                                                :: n_compute_atoms
    real*8,dimension(                   1:n_compute_atoms) :: dist_tab          ! Need dist_tab since dir_tab already normalized
    real*8,dimension(1:3,               1:n_compute_atoms) :: dir_tab
    integer,dimension(1:n_compute_c)                       :: i_basis
    integer,dimension(n_basis_fns,n_centers)               :: i_basis_fns_inv
    integer,dimension(n_centers)                           :: atom_index_inv
    real*8,dimension(1:n_compute_c,     1:3)               :: gradient_basis_wave
    integer                                                :: strain_index
    real*8,dimension(1:n_compute_c                       ) :: AS_strain_der

    ! local
    integer :: current_basis, i_basis_fn_1, atom_compute
    integer :: i_compute
    integer,dimension(1:9) :: l_index,m_index
 
    ! Map hessian-index -> grad, dir_tab
    l_index(1) = 1   
    l_index(2) = 1   
    l_index(3) = 1   
    l_index(4) = 2   
    l_index(5) = 2   
    l_index(6) = 3   
    l_index(7) = 2   
    l_index(8) = 3   
    l_index(9) = 3   
    m_index(1) = 1
    m_index(2) = 2
    m_index(3) = 3
    m_index(4) = 2
    m_index(5) = 3
    m_index(6) = 3
    m_index(7) = 1
    m_index(8) = 1
    m_index(9) = 2

    ! reset 
    AS_strain_der(:) = 0.0d0
     
    ! basis by basis
    do i_compute = 1, n_compute_c, 1
      current_basis = i_basis(i_compute)
      i_basis_fn_1  = i_basis_fns_inv(basis_fn(Cbasis_to_basis(current_basis)), Cbasis_to_center(current_basis))

      ! NB: if i_basis_fn_1 <= 0, AS_dde_basis_wave = 0 since gradient_basis_wave = 0
      ! If i_basis_fn_1<=0, then atom_compute=0 --> not allowed!
      if ( i_basis_fn_1.gt.0 ) then
        atom_compute             = atom_index_inv(Cbasis_to_center(current_basis))
        AS_strain_der(i_compute) = &
           gradient_basis_wave(i_compute,l_index(strain_index)) &
           * dir_tab(m_index(strain_index),atom_compute) &
           * dist_tab(atom_compute)
      end if
    end do

  end subroutine AS_eval_strain_deriv_wave                

  !! AS_eval_hessian_wave_times_xc_gradient_deriv_gga
  ! Construct
  !   (d/dgrad(rho)^2 f_xc) * grad(rho) * d/de_lm grad(|i>) = (d/dgrad(rho)^2 f_xc) * grad(rho) * (d/dr_l grad(|i>)) * (r-R_m)
  ! there is a scalar product between grad(rho) and (d/dr_l grad(|i>)).
  ! The variable xc_gradient_deriv includes grad(rho).
  subroutine AS_eval_hessian_wave_times_xc_gradient_deriv_gga( &
                               n_compute,                      &
                               n_compute_atoms,                &
                               dist_tab,                       &
                               dir_tab,                        &
                               i_basis,                        &
                               i_basis_fns_inv,                &
                               atom_index_inv,                 &
                               hessian_1,                      &
                               hessian_2,                      &
                               hessian_3,                      &
                               xc_gradient_deriv_ispin_ipoint, &
                               strain_index,                   &
                               AS_gga                          )
    use pbc_lists
    use dimensions, only: n_centers, n_basis_fns
    use basis, only: basis_fn
    implicit none

    integer,                                     intent(in)  :: n_compute
    integer,                                     intent(in)  :: n_compute_atoms
    real*8, dimension(    1:n_compute_atoms),    intent(in)  :: dist_tab        ! Need dist_tab since dir_tab already normalized
    real*8, dimension(1:3,1:n_compute_atoms),    intent(in)  :: dir_tab
    integer,dimension(1:n_compute),              intent(in)  :: i_basis
    integer,dimension(1:n_basis_fns,1:n_centers),intent(in)  :: i_basis_fns_inv
    integer,dimension(1:n_centers),              intent(in)  :: atom_index_inv
    real*8, dimension(1:n_compute),              intent(in)  :: hessian_1
    real*8, dimension(1:n_compute),              intent(in)  :: hessian_2
    real*8, dimension(1:n_compute),              intent(in)  :: hessian_3
    real*8, dimension(1:3),                      intent(in)  :: xc_gradient_deriv_ispin_ipoint
    integer,                                     intent(in)  :: strain_index
    real*8, dimension(1:n_compute),              intent(out) :: AS_gga

    ! local
    integer                 :: current_basis
    integer                 :: i_basis_fn_1
    integer                 :: atom_compute
    integer                 :: i_compute
    integer, dimension(1:9) :: m_index

    ! Map strain_index -> dir_tab
    m_index(1) = 1
    m_index(2) = 2
    m_index(3) = 3
    m_index(4) = 2
    m_index(5) = 3
    m_index(6) = 3
    m_index(7) = 1
    m_index(8) = 1
    m_index(9) = 2

    ! reset 
    AS_gga(:) = 0.0d0

    ! basis by basis
    do i_compute = 1, n_compute, 1
      current_basis = i_basis(i_compute)
      i_basis_fn_1  = i_basis_fns_inv(basis_fn(Cbasis_to_basis(current_basis)), Cbasis_to_center(current_basis))

      ! NB: If i_basis_fn_1<=0, AS_GGA=0 since hessian=0
      ! If i_basis_fn_1<=0, then atom_compute=0 --> not allowed!
      if ( i_basis_fn_1.gt.0 ) then
        atom_compute      = atom_index_inv(Cbasis_to_center(current_basis))
        AS_gga(i_compute) = (   xc_gradient_deriv_ispin_ipoint(1) * hessian_1(i_compute)   &
                              + xc_gradient_deriv_ispin_ipoint(2) * hessian_2(i_compute)   &
                              + xc_gradient_deriv_ispin_ipoint(3) * hessian_3(i_compute) ) &
                            * dir_tab(m_index(strain_index),atom_compute) * dist_tab(atom_compute)
      end if
    end do

  end subroutine AS_eval_hessian_wave_times_xc_gradient_deriv_gga

  !! AS_eval_hessian_wave_times_xc_tau_deriv_mgga
  ! Construct
  !   (d/dgrad(tau) f_xc) * d/de_lm grad(|i>) = (d/dgrad(tau) f_xc) * (d/dr_l grad(|i>)) * (r-R_m)
  subroutine AS_eval_hessian_wave_times_xc_tau_deriv_mgga(&
                               n_compute,                      &
                               n_compute_atoms,                &
                               dist_tab,                       &
                               dir_tab,                        &
                               i_basis,                        &
                               i_basis_fns_inv,                &
                               atom_index_inv,                 &
                               hessian_1,                      &
                               hessian_2,                      &
                               hessian_3,                      &
                               xc_tau_deriv_ispin_ipoint,      &
                               strain_index,                   &
                               AS_mgga                          )

    use dimensions, only: n_centers, n_basis_fns
    use pbc_lists
    implicit none

    integer,                                     intent(in)  :: n_compute
    integer,                                     intent(in)  :: n_compute_atoms
    real*8, dimension(    1:n_compute_atoms),    intent(in)  :: dist_tab        ! Need dist_tab since dir_tab already normalized
    real*8, dimension(1:3,1:n_compute_atoms),    intent(in)  :: dir_tab
    integer,dimension(1:n_compute),              intent(in)  :: i_basis
    integer,dimension(1:n_basis_fns,1:n_centers),intent(in)  :: i_basis_fns_inv
    integer,dimension(1:n_centers),              intent(in)  :: atom_index_inv
    real*8, dimension(1:n_compute),              intent(in)  :: hessian_1
    real*8, dimension(1:n_compute),              intent(in)  :: hessian_2
    real*8, dimension(1:n_compute),              intent(in)  :: hessian_3
    real*8,                                      intent(in)  :: xc_tau_deriv_ispin_ipoint
    integer,                                     intent(in)  :: strain_index
    real*8, dimension(1:n_compute),              intent(out) :: AS_mgga

    ! Local variable to convert xc_tau_deriv in to a format fitting the gga implementation
    real*8, dimension(1:3)                                   :: xc_gradient_deriv_ispin_ipoint_holder

    ! Set up dummy
    xc_gradient_deriv_ispin_ipoint_holder = xc_tau_deriv_ispin_ipoint

    ! Make necessary call to the already implemented method for xc_gradient_deriv
    ! There is probably a way to do this easier, as the outcome is just a scalar product.
    ! For now, however, this will do me as it works. Only necessary for the diagonal elements otherwise
    ! we get double-counting. AJL

    if (strain_index.eq.1.or. &
        strain_index.eq.4.or. &
        strain_index.eq.6) &    
    call AS_eval_hessian_wave_times_xc_gradient_deriv_gga(                      &
                                n_compute,                                      &
                                n_compute_atoms,                                &
                                dist_tab,                                       &
                                dir_tab,                                        &
                                i_basis,                                        &
                                i_basis_fns_inv,                                &
                                atom_index_inv,                                 &
                                hessian_1,                                      &
                                hessian_2,                                      &
                                hessian_3,                                      &
                                xc_gradient_deriv_ispin_ipoint_holder,          &
                                strain_index,                                   &
                                AS_mgga        )

  end subroutine AS_eval_hessian_wave_times_xc_tau_deriv_mgga

  !! AS_compute_i_on_site_basis
  ! Save on-site basis pairs
  subroutine AS_compute_i_on_site_basis( n_compute_a, i_basis, n_on_site, i_on_site_A, i_on_site_B )
   use pbc_lists
   implicit none

   integer :: n_compute_a
   integer :: i_basis(n_compute_a)
   integer :: n_on_site
   integer :: i_on_site_A(n_compute_a**2)
   integer :: i_on_site_B(n_compute_a**2)
   integer :: i_compute_a, i_compute_c

   ! Remember on-site terms once for all
   n_on_site      = 0
   i_on_site_A(:) = 0
   i_on_site_B(:) = 0
   do i_compute_a=1,n_compute_a
     do i_compute_c=1,n_compute_a
       if ( Cbasis_to_center(i_basis(i_compute_a)) .eq. Cbasis_to_center(i_basis(i_compute_c)) ) then
         n_on_site = n_on_site + 1
         i_on_site_A(n_on_site) = i_compute_a
         i_on_site_B(n_on_site) = i_compute_c
       end if
     end do
   end do

  end subroutine AS_compute_i_on_site_basis

  ! ###################################################################################
  ! Sync, output and miscellaneous routines
  !
  !! AS_map_hess_to_symm
  ! Map the 6 independent coordinates of the Hessian matrix to the full 9 coordinates.
  function AS_map_hess_to_symm(i_coord)
   implicit none

   integer :: AS_map_hess_to_symm,i_coord
 
   if (i_coord.le.6) then
     AS_map_hess_to_symm = i_coord
   elseif (i_coord.eq.7) then
     AS_map_hess_to_symm = 2
   elseif (i_coord.eq.8) then
     AS_map_hess_to_symm = 3
   elseif (i_coord.eq.9) then
     AS_map_hess_to_symm = 5
   end if

  end function

  !! AS_symmetrize_stress
  ! If the full analyitical stress tensor was calculated, i.e. AS_components=9,
  ! then this routine can be used to symmetrize the tensor.
  subroutine AS_symmetrize_stress(stress_tensor)
    implicit none
    real*8,dimension(1:3,1:3),intent(inout) :: stress_tensor

    ! local
    integer :: i_coord_1, i_coord_2

    do i_coord_1 = 1, 3, 1
      do i_coord_2 = i_coord_1+1, 3, 1
        stress_tensor(i_coord_1,i_coord_2) = (stress_tensor(i_coord_1,i_coord_2) + stress_tensor(i_coord_2,i_coord_1))/2
        stress_tensor(i_coord_2,i_coord_1) =  stress_tensor(i_coord_1,i_coord_2)
      end do
    end do

  end subroutine AS_symmetrize_stress

  !! AS_sync_analytical_stress
  ! Sync variables from forces_densmat (orbital dependent contributions).
  ! Everything else is synced in sum_up_whole_potential or in the routine below for the EXX part.
  subroutine AS_sync_analytical_stress(pulay_stress_on)
    use dimensions, only: use_AS_Jac_in_pulay
    use runtime_choices, only: flag_rel, REL_none
    use synchronize_mpi
    implicit none

    logical :: pulay_stress_on

    if (.not.use_AS_Jac_in_pulay) then
      call sync_matrix( AS_pu_rho_dde_fnrka, 3, 3 )
    end if
    if (pulay_stress_on) then
      call sync_matrix( AS_pulay_stress, 3, 3 )
      if (flag_rel.eq.REL_none) then
        call sync_matrix( AS_dde_T_sa_at_di_orb , 3, 3 )
      end if
    end if

  end subroutine AS_sync_analytical_stress

  !! AS_sync_EXX_stress
  ! Sync Hartree-Fock part of stress from calculate_fock_matrix_p0.
  ! Map from vector AS_EXX_stress_local to matrix AS_EXX_stress.
  subroutine AS_sync_EXX_stress()
    use synchronize_mpi_basic
    use runtime_choices, only: hybrid_coeff, AS_components
    implicit none

    integer :: i_coord_1
    integer :: i_coord_2
    integer :: i_counter

    ! Synchronize Hartree-Fock stress
    call sync_vector(AS_EXX_stress_local, 9)

    ! Map from one index (1 to 6 or 9) to 2 indices (1 to 3, 1 to 3)
    ! 1 2 3
    ! 7 4 5
    ! 8 9 6
    ! We do not symmetrize here if AS_components=6. This is only allowed if we take the total analytical stress.
    AS_EXX_stress(:,:) = 0.0d0
    i_counter = 0
    do i_coord_1 = 1, 3, 1
      do i_coord_2 = i_coord_1, 3, 1
        i_counter = i_counter + 1
        AS_EXX_stress( i_coord_1 , i_coord_2 ) = AS_EXX_stress_local( i_counter )
      end do
    end do
    ! Set the lower triangle for AS_components=9.
    if (AS_components .eq. 9) then
      AS_EXX_stress(2,1) = AS_EXX_stress_local(7)
      AS_EXX_stress(3,1) = AS_EXX_stress_local(8)
      AS_EXX_stress(3,2) = AS_EXX_stress_local(9)
    end if

  end subroutine AS_sync_EXX_stress

  !! get_total_analytical_stress
  ! Sum up all contribution of the analytical stress tensor, print it and do convergence check
  subroutine get_total_analytical_stress(pulay_on, en_xc_minus_en_pot_xc, en_exx, &
                                         previous_stress_loc, diff_stress_loc, vdw_flag)
    use dimensions
    use physics     
    use runtime_choices
    use numerical_stress
    use sym_base, only: symmetrize_stress_spg, destroy_symmetry_arrays
    use spglib_symmetry, only: write_symm_info, out_symm_mats,&
      & destroy_symmetry, destroy_symmats
    implicit none
 
    logical,                 intent(in)    :: pulay_on
    real*8,                  intent(in)    :: en_xc_minus_en_pot_xc
    real*8,                  intent(in)    :: en_exx
    real*8, dimension(3, 3), intent(inout) :: previous_stress_loc
    real*8,                  intent(out)   :: diff_stress_loc
    logical,                 intent(in)    :: vdw_flag

    ! local
    real*8, dimension(3, 3) :: orbital_dependent_stress
    real*8, dimension(3, 3) :: mp_elec_stress
    real*8, dimension(3, 3) :: mp_nucl_stress
    integer,dimension(9)    :: l_index,m_index
    real*8                  :: pressure
    real*8, dimension(3, 3) :: AS_total_stress
    real*8, dimension(3, 3) :: AS_output_stress
    real*8, dimension(3, 3) :: diff_stress_matrix
    character*400           :: info_str
    integer                 :: i_coord, i_coord_2
    integer                 :: counter

    ! Map hessian-index for output
    !  1 4 5
    !  7 2 6
    !  8 9 3
    l_index(1) = 1
    l_index(2) = 2
    l_index(3) = 3
    l_index(4) = 1
    l_index(5) = 1
    l_index(6) = 2
    l_index(7) = 2
    l_index(8) = 3
    l_index(9) = 3
    m_index(1) = 1
    m_index(2) = 2
    m_index(3) = 3
    m_index(4) = 2
    m_index(5) = 3
    m_index(6) = 3
    m_index(7) = 1
    m_index(8) = 1
    m_index(9) = 2
 
    orbital_dependent_stress(:,:)        = 0.0d0
    mp_elec_stress(:,:)                  = 0.0d0
    mp_nucl_stress(:,:)                  = 0.0d0

    ! Combine the different ingredients...
    ! NB --> Note below
    do i_coord = 1, 3, 1
      do i_coord_2 = 1, 3 ,1

        mp_nucl_stress( i_coord , i_coord_2 ) = mp_nucl_stress( i_coord , i_coord_2 ) - &
           NEW_AS_at_stress( i_coord , i_coord_2 )

        mp_elec_stress( i_coord , i_coord_2 ) = mp_elec_stress( i_coord , i_coord_2 ) + &
           ( NEW_AS_rho_stress(i_coord , i_coord_2 ) - AS_rho_f_f_self_correction( i_coord , i_coord_2 ) )

        if (pulay_on) then
          if (flag_rel .ne. REL_atomic_zora) then
            orbital_dependent_stress( i_coord , i_coord_2 ) = orbital_dependent_stress( i_coord , i_coord_2 ) + &
                   ( AS_pulay_stress( i_coord , i_coord_2 ) &
                   + AS_dde_T_sa_at_di_orb(  i_coord , i_coord_2 ) )
          else
            orbital_dependent_stress( i_coord , i_coord_2 ) = orbital_dependent_stress( i_coord , i_coord_2 ) + &
                   ( AS_pulay_stress( i_coord , i_coord_2 ) ) 
          end if
          if ( .not. use_AS_Jac_in_pulay ) then
            orbital_dependent_stress( i_coord , i_coord_2 ) = orbital_dependent_stress( i_coord , i_coord_2 ) + &
                   ( AS_pu_rho_dde_fnrka(i_coord,i_coord_2) ) 
            if ( i_coord .eq. i_coord_2 ) then
              orbital_dependent_stress( i_coord , i_coord_2 ) = orbital_dependent_stress( i_coord , i_coord_2 ) + &
                     ( en_xc_minus_en_pot_xc )
            end if
          end if !flag_rel
        end if !pulay_on
      end do
    end do

    ! Sum up all contributions to the stress:
    ! Generally AS_components=6, ie, we compute only the upper right of the stress tensor of the orbital dependent contributions, since the stress
    ! tensor is symmetric. However, its individual components are not symmetric: Combining the lower left components of mp_elec_stress with the upper
    ! right components of orbital_dependent_stress gives wrong off-diagonal elements, whenever symmetry is broken.
    AS_total_stress(:,:) = 0.0d0
    do i_coord = 1, 3, 1
      do i_coord_2 = 1, 3, 1

        ! Symmetrize here for AS_components=6
        if ( (AS_components.ne.9) .and. (i_coord_2.lt.i_coord) ) then
          AS_total_stress( i_coord , i_coord_2 ) = AS_total_stress( i_coord_2 , i_coord )
          cycle
        end if

        AS_total_stress( i_coord , i_coord_2 ) = AS_total_stress( i_coord , i_coord_2 ) &
         + ( orbital_dependent_stress( i_coord , i_coord_2 ) + mp_elec_stress( i_coord , i_coord_2 ) + &
         &   mp_nucl_stress( i_coord , i_coord_2 ) )

        ! EXX stress is only calculated when the pulay contributions are calculated
        if ( use_hartree_fock .and. pulay_on ) then
          AS_total_stress( i_coord , i_coord_2 ) = &
             AS_total_stress( i_coord , i_coord_2 ) &
             - hybrid_coeff * AS_EXX_stress( i_coord , i_coord_2 )
          if ( i_coord .eq. i_coord_2 ) then
            AS_total_stress( i_coord , i_coord_2 ) = &
               AS_total_stress( i_coord , i_coord_2 ) &
               - 2.d0 * hybrid_coeff * en_exx
          end if
        end if

        ! Only add vdW term if it has been calculated
        if (( &
            use_vdw_correction_hirshfeld .or. use_mbd_std .or. use_mbd_dev .or. use_libmbd &
        ) .and. (.not. vdw_flag)) then
          AS_total_stress( i_coord , i_coord_2 ) = &
             AS_total_stress( i_coord , i_coord_2 ) &
             + AS_vdw_stress( i_coord , i_coord_2 )
        end if
      end do
    end do

    analytical_stress_tensor(1:3,1:3) = AS_total_stress(1:3,1:3)/cell_volume

    call localorb_info("",use_unit,'(A)')
    if (AS_components.eq.9) then
      write(info_str,'(A,8(18X,A))') "Analytical stress tensor components [eV]         xx","yy","zz","xy","xz","yz","yx","zx","zy"
      call localorb_info(info_str,use_unit,'(2X,A)')
      write(info_str,'(4A)') "---------------------------------------------", &
         "-----------------------------------------------------------------", &
         "-----------------------------------------------------------------", &
         "----------------------------------------"
      call localorb_info(info_str,use_unit,'(2X,A)')
    else
      write(info_str,'(A,5(18X,A))') "Analytical stress tensor components [eV]         xx","yy","zz","xy","xz","yz"
      call localorb_info(info_str,use_unit,'(2X,A)')
      write(info_str,'(3A)') "---------------------------------------------", &
       "-------------------------------------------------------------------", &
       "-------------------------------------------"
      call localorb_info(info_str,use_unit,'(2X,A)')
    end if
    write(info_str,'(2X,A,9E20.10)') "Nuclear Hellmann-Feynman      : ", &
       ( -1.0*NEW_AS_at_stress(l_index(i_coord),m_index(i_coord)) * hartree, &
         i_coord=1,AS_components,1 )
    call localorb_info(info_str,use_unit,'(2X,A)')
    write(info_str,'(2X,A,9E20.10)') "Multipole Hellmann-Feynman    : ", &
       ( NEW_AS_rho_stress(l_index(i_coord),m_index(i_coord)) * hartree, &
         i_coord=1,AS_components,1 )
    call localorb_info(info_str,use_unit,'(2X,A)')
    write(info_str,'(2X,A,9E20.10)') "On-site Multipole corrections : ", &
       ( -1.0 * AS_rho_f_f_self_correction(l_index(i_coord),m_index(i_coord)) &
              * hartree, i_coord=1,AS_components,1 )
    call localorb_info(info_str,use_unit,'(2X,A)')
    write(info_str,'(2X,A,9E20.10)') "Strain deriv. of the orbitals : ", &
       ( AS_pulay_stress(l_index(i_coord),m_index(i_coord)) * hartree, &
         i_coord=1,AS_components,1 )
    call localorb_info(info_str,use_unit,'(2X,A)')
    if (flag_rel .eq. REL_none) then
      write(info_str,'(2X,A,9E20.10)') "On-site kinetic E corrections : ", &
         ( AS_dde_T_sa_at_di_orb(l_index(i_coord),m_index(i_coord)) * hartree, &
           i_coord=1,AS_components,1 )
      call localorb_info(info_str,use_unit,'(2X,A)')
    end if
    if ( .not. use_AS_Jac_in_pulay ) then
      write(info_str,'(2X,A,9E20.10)') "Multipole pot. derivative     : ", &
         ( AS_pu_rho_dde_fnrka(l_index(i_coord),m_index(i_coord)) * hartree, &
           i_coord=1,AS_components,1 )
      call localorb_info(info_str,use_unit,'(2XA)')
      write(info_str,'(2X,A,9E20.10)') "Jacobian contribution         : ", &
      ( en_xc_minus_en_pot_xc * hartree, i_coord=1,3,1 ), &
      ( 0.0d0, i_coord=4,AS_components,1 )
      call localorb_info(info_str,use_unit,'(2X,A)')
    end if
    if ( use_hartree_fock .and. pulay_on ) then
      write(info_str,'(2X,A,9E20.10)') "Exact exchange derivative     : ", &
         ( - hybrid_coeff * AS_EXX_stress(l_index(i_coord),m_index(i_coord)) * hartree, &
           i_coord=1,AS_components,1 )
      call localorb_info(info_str,use_unit,'(2X,A)')
      write(info_str,'(2X,A,9E20.10)') "Jacobian EXX contribution     : ", &
         ( - 2 * hybrid_coeff * en_exx * hartree, i_coord=1,3,1 ), &
         ( 0.0d0, i_coord=4,AS_components,1 )
      call localorb_info(info_str,use_unit,'(2X,A)')
    end if
    ! Only print vdW term if it has been calculated, i.e. we have finished the scf cycle and are doing post-processing
    if (( &
        use_vdw_correction_hirshfeld .or. use_mbd_std .or. use_mbd_dev .or. use_libmbd &
    ) .and. (.not. vdw_flag)) then
      write(info_str,'(2X,A,9E20.10)') "vdW correction   : ", &
         ( AS_vdw_stress(l_index(i_coord),m_index(i_coord)) * hartree, &
           i_coord=1,AS_components,1 )
      call localorb_info(info_str,use_unit,'(2X,A)')
    end if
    if (AS_components.eq.9) then
      write(info_str,'(4A)')  "---------------------------------------------", &
         "-----------------------------------------------------------------", &
         "-----------------------------------------------------------------", &
         "----------------------------------------"
    else
      write(info_str,'(3A)')  "---------------------------------------------", &
       "-------------------------------------------------------------------", &
       "-------------------------------------------"
    end if
    call localorb_info(info_str,use_unit,'(2X,A)')
    write(info_str,'(A,9E20.10)') "Sum of all contributions        : ", &
       ( AS_total_stress(l_index(i_coord),m_index(i_coord)) * hartree, &
         i_coord=1,AS_components,1 )
    call localorb_info(info_str,use_unit,'(2X,A)')                                     
    call localorb_info("",use_unit,'(A)')

    !! Here comes the full output in eV/AA^3 -- see numerical_stress_tensor.f90 and predict_new_geometry.f90
    ! Calculate Pressure: strictly speaking this only makes sense in the hydrostatic case
    ! where the stress tensor is a multiple of unit.
    ! Here we calculate the normed trace of the stress tensor, which corresponds
    ! to the arithmetic mean of the diagonal elements (times minus 1).
    pressure = (  analytical_stress_tensor(1,1) &
                + analytical_stress_tensor(2,2) &
                + analytical_stress_tensor(3,3)  )/(-3.0d0)

    AS_output_stress(1:3,1:3) = analytical_stress_tensor(1:3,1:3)

    ! Print the unsymmetrized analytical stress tensor if it was calculated.
    if (AS_components.eq.9) then
      call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
      if (vdw_flag) then
        call localorb_info("|            Analytical stress tensor w/o vdW correction            |",use_unit,'(2X,A)')
      else
        call localorb_info("|                     Analytical stress tensor                      |",use_unit,'(2X,A)')
      end if
      call localorb_info("|                  Cartesian components [eV/A**3]                   |",use_unit,'(2X,A)')
      call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
      call localorb_info("|                x                y                z                |",use_unit,'(2X,A)')
      call localorb_info("|                                                                   |",use_unit,'(2X,A)')

      write(info_str,'(2X,A1,A,3(2X,F15.8),11X,A1)') "|", "  x  ", &
                        (AS_output_stress(1,i_coord)*hartree / (bohr**3),i_coord=1,3,1), "|"
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(2X,A1,A,3(2X,F15.8),11X,A1)') "|", "  y  ", &
                        (AS_output_stress(2,i_coord)*hartree / (bohr**3),i_coord=1,3,1), "|"
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(2X,A1,A,3(2X,F15.8),11X,A1)') "|", "  z  ", &
                        (AS_output_stress(3,i_coord)*hartree / (bohr**3),i_coord=1,3,1), "|"
      call localorb_info(info_str,use_unit,'(A)')

      call localorb_info("|                                                                   |",use_unit,'(2X,A)')
      call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
      call localorb_info("",use_unit,'(A)')

      ! Symmetrize output of stress
      call AS_symmetrize_stress(AS_output_stress)
    end if

    ! spg symmetrization
    if (use_symmetry_reduced_spg) then
       call symmetrize_stress_spg(AS_output_stress)
       call symmetrize_stress_spg(analytical_stress_tensor)
    endif
    if (use_spglib .and. use_symmetric_forces .and..not.use_symmetry_reduced_spg)then
       call destroy_symmetry()
       call destroy_symmats()
       !    call destroy_sym_maps()
       call destroy_symmetry_arrays()
       call write_symm_info()
       call out_symm_mats()
       call symmetrize_stress_spg(AS_output_stress)
       call symmetrize_stress_spg(analytical_stress_tensor)
    endif

    ! Print the symmetrized analytical stress tensor.
    call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
    if (vdw_flag) then
      call localorb_info("|     Analytical stress tensor w/o vdW correction - Symmetrized     |",use_unit,'(2X,A)')
    else
      call localorb_info("|              Analytical stress tensor - Symmetrized               |",use_unit,'(2X,A)')
    end if
    call localorb_info("|                  Cartesian components [eV/A**3]                   |",use_unit,'(2X,A)')
    call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
    call localorb_info("|                x                y                z                |",use_unit,'(2X,A)')
    call localorb_info("|                                                                   |",use_unit,'(2X,A)')

    write(info_str,'(2X,A1,A,3(2X,F15.8),11X,A1)') "|", "  x  ", &
                      (AS_output_stress(1,i_coord)*hartree / (bohr**3),i_coord=1,3,1), "|"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A1,A,3(2X,F15.8),11X,A1)') "|", "  y  ", &
                      (AS_output_stress(2,i_coord)*hartree / (bohr**3),i_coord=1,3,1), "|"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A1,A,3(2X,F15.8),11X,A1)') "|", "  z  ", &
                      (AS_output_stress(3,i_coord)*hartree / (bohr**3),i_coord=1,3,1), "|"
    call localorb_info(info_str,use_unit,'(A)')

    call localorb_info("|                                                                   |",use_unit,'(2X,A)')
    write(info_str,'(2X,A1,A,2X,F15.8,2X,A,26X,A1)') "|", "  Pressure:", pressure*hartree / (bohr**3), " [eV/A**3] " ,"|"
    call localorb_info(info_str,use_unit,'(A)')
    call localorb_info("|                                                                   |",use_unit,'(2X,A)')
    call localorb_info("+-------------------------------------------------------------------+",use_unit,'(2X,A)')
    call check_pressure(analytical_stress_tensor)
    call localorb_info("",use_unit,'(A)')

    ! Convergence check:
    ! Square the difference between present and previous stress entry.
    do i_coord = 1, 3, 1
      do i_coord_2 = 1, 3, 1
        diff_stress_matrix(i_coord, i_coord_2) = &
           (analytical_stress_tensor(i_coord, i_coord_2) &
            - previous_stress_loc(i_coord, i_coord_2))**2
      end do
    end do

    ! Get the maximum difference, take square root and convert difference to eV/AA^3.
    diff_stress_loc = sqrt(maxval(diff_stress_matrix)) * hartree / (bohr**3)
    ! Store present stress tensor as previous stress for next iteration.
    previous_stress_loc(1:3,1:3) = analytical_stress_tensor(1:3,1:3)
 
  end subroutine get_total_analytical_stress 
 

  !##################################################################################
  !! DEBUG READ/WRITE ROUTINES: As mentioned above, this routines allow to read
  !! write multipole moments and Pulay Density Matrices from/to disk
  subroutine AS_write_mm_moments(mm, mm_dim_1, mm_dim_2)
    implicit none
    integer,intent(in)                                 :: mm_dim_1,mm_dim_2
    real*8,dimension(1:mm_dim_1,1:mm_dim_2),intent(in) :: mm
    integer                                            :: i_dim_1,i_dim_2

    open(111,FILE="multipol_moments_to_file.dat")
    do i_dim_1=1,mm_dim_1
      do i_dim_2=1,mm_dim_2
        write(111,'(I5,I5,EN30.15)') i_dim_1, i_dim_2, mm(i_dim_1,i_dim_2)
      end do
    end do
    close(111)
  end subroutine AS_write_mm_moments

  subroutine AS_read_mm_moments(mm, mm_dim_1, mm_dim_2)
    use localorb_io, only: use_unit
    implicit none
    integer,intent(in)                                  :: mm_dim_1,mm_dim_2
    real*8,dimension(1:mm_dim_1,1:mm_dim_2),intent(out) :: mm
    integer                                             :: i_dim_1,i_dim_2
    integer                                             :: i_dim_1_temp,i_dim_2_temp

    open(111,FILE="multipol_moments_from_file.dat")
    do i_dim_1=1,mm_dim_1
      do i_dim_2=1,mm_dim_2
        read(111,'(I5,I5,EN30.15)') i_dim_1_temp, i_dim_2_temp, mm(i_dim_1,i_dim_2)
        if ( (i_dim_1_temp .ne. i_dim_1) .or. (i_dim_2_temp .ne. i_dim_2) )then
          write(use_unit,*) " *** CONSISTENCY ERROR"
          stop
        end if
      end do
    end do
    close(111)

    ! Check
    call AS_write_mm_moments(mm,mm_dim_1, mm_dim_2)

  end subroutine AS_read_mm_moments

  subroutine AS_write_mm_splines(mm, mm_dim_1, mm_dim_2, mm_dim_3, atom )
    implicit none
    integer,intent(in)                                            :: mm_dim_1,mm_dim_2,mm_dim_3, atom
    real*8,dimension(1:mm_dim_1,1:mm_dim_2,1:mm_dim_3),intent(in) :: mm
    integer                                                       :: i_dim_1,i_dim_2,i_dim_3
    character*200 :: filename

    write(filename,'(A,I6.6,A)')   "multipol_splines_to_file.at_",atom,".dat"
    open(111,FILE=filename)
    do i_dim_1=1,mm_dim_1
      do i_dim_2=1,mm_dim_2
        do i_dim_3=1,mm_dim_3
          write(111,'(I5,I5,I5,EN30.15)') i_dim_1, i_dim_2, i_dim_3, mm(i_dim_1,i_dim_2,i_dim_3)
        end do
      end do
    end do
    close(111)
  end subroutine AS_write_mm_splines

  subroutine AS_read_mm_splines(mm, mm_dim_1, mm_dim_2, mm_dim_3, atom)
    use localorb_io, only: use_unit
    implicit none
    integer,intent(in)                                             :: mm_dim_1,mm_dim_2,mm_dim_3,atom
    real*8,dimension(1:mm_dim_1,1:mm_dim_2,1:mm_dim_3),intent(out) :: mm
    integer                                                        :: i_dim_1,i_dim_2,i_dim_3
    integer                                                        :: i_dim_1_temp,i_dim_2_temp,i_dim_3_temp

    character*200 :: filename

    write(filename,'(A,I6.6,A)')   "multipol_splines_from_file.at_",atom,".dat"
    open(111,FILE=filename)
    do i_dim_1=1,mm_dim_1
      do i_dim_2=1,mm_dim_2
        do i_dim_3=1,mm_dim_3
          read(111,'(I5,I5,I5,EN30.15)') i_dim_1_temp, i_dim_2_temp, i_dim_3_temp, mm(i_dim_1,i_dim_2,i_dim_3)
          if ( (i_dim_1_temp .ne. i_dim_1) .or. (i_dim_2_temp .ne. i_dim_2) .or. (i_dim_3_temp .ne. i_dim_3) )then
            write(use_unit,*) " *** CONSISTENCY ERROR SPLINES"
            stop
          end if
        end do
      end do
    end do
    close(111)
    call AS_write_mm_splines(mm, mm_dim_1, mm_dim_2, mm_dim_3, atom)
  end subroutine AS_read_mm_splines

  subroutine AS_write_pulay_dm( Matrix, dim1, dim2, step)
    implicit none
    integer,intent(in)                   :: dim1,dim2,step
    real*8,dimension(dim1,dim2),intent(in) :: Matrix
    ! local
    integer :: i_dim1, i_dim2

    if (step .eq. 1) then
      open (333,FILE="Pulay_dm_i_term_1.to_file.dat",STATUS='replace',ACTION='write')
    else
      open (333,FILE="Pulay_dm_i_term_2.to_file.dat",STATUS='replace',ACTION='write')
    end if
 
    do i_dim1=1,dim1,1
      do i_dim2=1,dim2,1
        write(333,'(2I8,1EN30.20)') i_dim1,i_dim2, Matrix(i_dim1,i_dim2)
      end do
    end do

    close(333)

  end subroutine AS_write_pulay_dm

  subroutine AS_read_pulay_dm( Matrix, dim1, dim2, step)
    use localorb_io, only: use_unit
    implicit none
    integer,intent(in)                      :: dim1,dim2,step
    real*8,dimension(dim1,dim2),intent(out) :: Matrix
    ! local
    integer :: i_dim1, i_dim2,i_dim1_tmp, i_dim2_tmp

    if (step .eq. 1) then
      open (333,FILE="Pulay_dm_i_term_1.from_file.dat",STATUS='old',ACTION='read')
    else
      open (333,FILE="Pulay_dm_i_term_2.from_file.dat",STATUS='old',ACTION='read')
    end if
 
    do i_dim1=1,dim1,1
      do i_dim2=1,dim2,1
        read(333,'(2I8,1EN30.20)') i_dim1_tmp,i_dim2_tmp, Matrix(i_dim1,i_dim2)
        if ( (i_dim1_tmp .ne. i_dim1) .or. (i_dim2_tmp .ne. i_dim2) ) then
          write(use_unit,*) "READ ERROR *******"
          stop
        end if
      end do
    end do


    close(333)
  end subroutine AS_read_pulay_dm


end module analytical_stress
