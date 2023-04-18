!****h* FHI-aims/hartree_fock_p0
!  NAME
!    hartree_fock_p0
!  SYNOPSIS

module hartree_fock_p0

!  PURPOSE  
!  this module contains the declaration of all the variables, arrays that
!  are commonly used in periodic or cluster Hartree-Fock calculations using RI_method LVL_fast.
!  It contains the subroutines for allocating and deallocating these variables.  
!
!  Subroutines:
!  * allocate_hartree_fock_p0
!  * clean_up_hartree_fock_p0
!
!  USES

      implicit none

!  INPUT
!  o none
!  OUTPUT
!  o none

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
!
!  SOURCE

!   Variales
!   n_homo_max -- the maximal HOMO level of the two spin channels

      real*8, dimension(:,:,:,:), allocatable :: hf_exchange_matr_real
      complex*16, dimension(:,:,:,:), allocatable :: hf_exchange_matr_complex
      
!      Fo LC-wPBEh wich also needs those arrays for the short range exchange
      real*8, dimension(:,:,:,:), allocatable :: hf_exchange_matr_real_SR
      complex*16, dimension(:,:,:,:), allocatable :: hf_exchange_matr_complex_SR
      
      real*8, dimension(:,:,:,:), allocatable, private :: prev_exmat
      complex*16, dimension(:,:,:,:), allocatable, private :: prev_exmat_complex
      real*8, dimension(:,:,:,:,:), allocatable, private :: prev_exmat_diff
      complex*16, dimension(:,:,:,:,:), allocatable, private :: prev_exmat_diff_complex
      real*8, dimension(:,:,:,:,:), allocatable, private :: prev_exmat_error
      complex*16, dimension(:,:,:,:,:), allocatable, private :: prev_exmat_error_complex

contains
  subroutine allocate_hartree_fock_p0
!
!  allocate the arrays needed in Hartree-Fock calculation
!  and set the disk storage file names 
!
    use dimensions
    use prodbas
    use runtime_choices
    use mpi_tasks, only: check_allocation
    use scalapack_wrapper, only: mxld,mxcol
    implicit none

    integer :: info
    character(*), parameter :: func = 'allocate_hartree_fock_p0'

    if(use_scalapack) then

         ! Allocations for Scalapack version
         ! Please note that the number of dimensions is the same as in the Lapack version,
         ! but (n_basis,n_basis,n_k_points_task) is replaced by (mxld, mxcol, 1)

         if(real_eigenvectors)then
            if(.not.allocated(hf_exchange_matr_real)) then
               allocate(hf_exchange_matr_real(mxld,mxcol,1,n_spin))
            endif
         else
            if(.not.allocated(hf_exchange_matr_complex)) then
               allocate(hf_exchange_matr_complex(mxld,mxcol,1,n_spin))
            endif
         endif
         
         if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
            if(real_eigenvectors)then
               if(.not.allocated(hf_exchange_matr_real_SR)) then
                  allocate(hf_exchange_matr_real_SR(mxld,mxcol,1,n_spin))
               endif
            else
               if(.not.allocated(hf_exchange_matr_complex_SR)) then
                  allocate(hf_exchange_matr_complex_SR(mxld,mxcol,1,n_spin))
               endif
            endif
         endif

         if(real_eigenvectors)then
            if (.not.allocated(prev_exmat)) then
               allocate(prev_exmat(mxld, mxcol, 1, n_spin), stat=info)
               call check_allocation(info, 'prev_exmat', func)
            endif
            prev_exmat = 0d0
            if(mixer.eq.MIX_PULAY) then
               if(.not.allocated(prev_exmat_diff)) then
                  allocate(prev_exmat_diff(mxld, mxcol, 1, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_diff', func)
               endif
               prev_exmat_diff = 0.d0
               if(.not.allocated(prev_exmat_error)) then
                  allocate(prev_exmat_error(mxld, mxcol, 1, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_error', func)
               endif
               prev_exmat_error = 0.d0
            endif
            if(mixer.eq.MIX_BROYDEN) then
               if(.not.allocated(prev_exmat_diff)) then
                  allocate(prev_exmat_diff(mxld, mxcol, 1, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_diff', func)
               endif
               prev_exmat_diff = 0.d0
               if(.not.allocated(prev_exmat_error)) then
                  allocate(prev_exmat_error(mxld, mxcol, 1, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_error', func)
               endif
               prev_exmat_error = 0.d0
            endif
         else
            if (.not.allocated(prev_exmat_complex)) then
               allocate(prev_exmat_complex(mxld, mxcol, 1, n_spin), stat=info)
               call check_allocation(info, 'prev_exmat_complex', func)
            endif
            prev_exmat_complex = (0d0, 0d0)
            if(mixer.eq.MIX_PULAY) then
               if(.not.allocated(prev_exmat_diff_complex)) then
                  allocate(prev_exmat_diff_complex(mxld, mxcol, 1, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_diff_complex', func)
               endif
               prev_exmat_diff_complex = (0.d0, 0d0)
               if(.not.allocated(prev_exmat_error_complex)) then
                  allocate(prev_exmat_error_complex(mxld, mxcol, 1, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_error_complex', func)
               endif
               prev_exmat_error_complex = (0.d0, 0d0)
            endif
            if(mixer.eq.MIX_BROYDEN) then
               if(.not.allocated(prev_exmat_diff_complex)) then
                  allocate(prev_exmat_diff_complex(mxld, mxcol, 1, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_diff_complex', func)
               endif
               prev_exmat_diff_complex = (0.d0, 0d0)
               if(.not.allocated(prev_exmat_error_complex)) then
                  allocate(prev_exmat_error_complex(mxld, mxcol, 1, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_error_complex', func)
               endif
               prev_exmat_error_complex = (0.d0, 0d0)
            endif
         endif

    else

         ! Allocations for Lapack version

         if(real_eigenvectors)then
            if(.not.allocated(hf_exchange_matr_real)) then
               allocate(hf_exchange_matr_real(n_basis,n_basis,n_k_points_task,n_spin))
            endif
         else
            if(.not.allocated(hf_exchange_matr_complex)) then
               allocate(hf_exchange_matr_complex(n_basis,n_basis,n_k_points_task,n_spin))
            endif
         endif
         
         if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
            if(real_eigenvectors)then
               if(.not.allocated(hf_exchange_matr_real_SR)) then
                  allocate(hf_exchange_matr_real_SR(n_basis,n_basis,n_k_points_task,n_spin))
               endif
            else
               if(.not.allocated(hf_exchange_matr_complex_SR)) then
                  allocate(hf_exchange_matr_complex_SR(n_basis,n_basis,n_k_points_task,n_spin))
               endif
            endif
         endif

         if(real_eigenvectors)then
            if (.not.allocated(prev_exmat)) then
               allocate(prev_exmat(n_basis, n_basis, n_k_points_task, n_spin), stat=info)
               call check_allocation(info, 'prev_exmat', func)
            endif
            prev_exmat = 0d0
            if(mixer.eq.MIX_PULAY) then
               if(.not.allocated(prev_exmat_diff)) then
                  allocate(prev_exmat_diff(n_basis, n_basis, n_k_points_task, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_diff', func)
               endif
               prev_exmat_diff = 0.d0
               if(.not.allocated(prev_exmat_error)) then
                  allocate(prev_exmat_error(n_basis, n_basis, n_k_points_task, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_error', func)
               endif
               prev_exmat_error = 0.d0
            endif
            if(mixer.eq.MIX_BROYDEN) then
               if(.not.allocated(prev_exmat_diff)) then
                  allocate(prev_exmat_diff(n_basis, n_basis, n_k_points_task, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_diff', func)
               endif
               prev_exmat_diff = 0.d0
               if(.not.allocated(prev_exmat_error)) then
                  allocate(prev_exmat_error(n_basis, n_basis, n_k_points_task, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_error', func)
               endif
               prev_exmat_error = 0.d0
            endif
         else
            if (.not.allocated(prev_exmat_complex)) then
               allocate(prev_exmat_complex(n_basis, n_basis, n_k_points_task, n_spin), stat=info)
               call check_allocation(info, 'prev_exmat_complex', func)
            endif
            prev_exmat_complex = (0d0, 0d0)
            if(mixer.eq.MIX_PULAY) then
               if(.not.allocated(prev_exmat_diff_complex)) then
                  allocate(prev_exmat_diff_complex(n_basis, n_basis, n_k_points_task, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_diff_complex', func)
               endif
               prev_exmat_diff_complex = (0.d0, 0d0)
               if(.not.allocated(prev_exmat_error_complex)) then
                  allocate(prev_exmat_error_complex(n_basis, n_basis, n_k_points_task, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_error_complex', func)
               endif
               prev_exmat_error_complex = (0.d0, 0d0)
            endif
            if(mixer.eq.MIX_BROYDEN) then
               if(.not.allocated(prev_exmat_diff_complex)) then
                  allocate(prev_exmat_diff_complex(n_basis, n_basis, n_k_points_task, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_diff_complex', func)
               endif
               prev_exmat_diff_complex = (0.d0, 0d0)
               if(.not.allocated(prev_exmat_error_complex)) then
                  allocate(prev_exmat_error_complex(n_basis, n_basis, n_k_points_task, n_spin, n_max_pulay), stat=info)
                  call check_allocation(info, 'prev_exmat_error_complex', func)
               endif
               prev_exmat_error_complex = (0.d0, 0d0)
            endif
         endif

    endif

  end subroutine allocate_hartree_fock_p0

!   deallocate the arrays
        subroutine cleanup_hartree_fock_p0 ()

         if (allocated(hf_exchange_matr_real)) then
           deallocate(hf_exchange_matr_real)
         endif
         if (allocated(hf_exchange_matr_complex)) then
           deallocate(hf_exchange_matr_complex)
         endif
         
         if (allocated(hf_exchange_matr_real_SR)) then
           deallocate(hf_exchange_matr_real_SR)
         endif
         if (allocated(hf_exchange_matr_complex_SR)) then
           deallocate(hf_exchange_matr_complex_SR)
         endif

         if(allocated(prev_exmat))then
            deallocate(prev_exmat)
         endif
         if(allocated(prev_exmat_complex))then
            deallocate(prev_exmat_complex)
         endif
         if(allocated(prev_exmat_diff))then
            deallocate(prev_exmat_diff)
         endif
         if(allocated(prev_exmat_diff_complex))then
            deallocate(prev_exmat_diff_complex)
         endif
         if(allocated(prev_exmat_error))then
            deallocate(prev_exmat_error)
         endif
         if(allocated(prev_exmat_error_complex))then
            deallocate(prev_exmat_error_complex)
         endif

       end  subroutine cleanup_hartree_fock_p0

!****s* FHI-aims/exchange_matrix_mixing_p0
!  NAME
!   exchange_matrix_mixing_p0
!  SYNOPSIS

subroutine exchange_matrix_mixing_p0(number_of_loops)

  !  PURPOSE
  !
  !    Mixing exchange matrix for self-consistent Hartree-Fock or hybrid
  !    functional periodic calculations.
  !
  !  USES

  use mpi_tasks
  use dimensions
  use runtime_choices
  use hartree_fock
  use mixing, only: pulay_saved_iter_denmat, mixing_factor
  use scalapack_wrapper, only: mxld,mxcol

  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: number_of_loops

  !   INPUTS 
  !     o number_of_loops -- the current self-consistent iteration
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

  real*8, dimension(:,:,:,:), allocatable :: delta_exmat ! output - input  (error)
  real*8, dimension(:,:,:,:), allocatable :: exmat_diff  ! mixed - output (update)
  complex*16, dimension(:,:,:,:), allocatable :: delta_exmat_complex
  complex*16, dimension(:,:,:,:), allocatable :: exmat_diff_complex
  integer i_store
  integer i_basis_1, i_basis_2, i_spin, i_k_point, i_k, i_count
  character*150 :: info_str
  integer :: info
  character(*), parameter :: func = 'density_matrix_mixing'

  if(use_scalapack) then
     if(real_eigenvectors)then
        allocate(delta_exmat(mxld, mxcol, 1, n_spin), stat=info)
        call check_allocation(info, 'delta_exmat', func)
        allocate(exmat_diff(mxld, mxcol, 1, n_spin), stat=info)
        call check_allocation(info, 'exmat_diff', func)
     else
        allocate(delta_exmat_complex(mxld, mxcol, 1, n_spin), stat=info)
        call check_allocation(info, 'delta_exmat_complex', func)
        allocate(exmat_diff_complex(mxld, mxcol, 1, n_spin), stat=info)
        call check_allocation(info, 'exmat_diff_complex', func)
     endif
  else
     if(real_eigenvectors)then
        allocate(delta_exmat(n_basis, n_basis, n_k_points_task, n_spin), stat=info)
        call check_allocation(info, 'delta_exmat', func)
        allocate(exmat_diff(n_basis, n_basis, n_k_points_task, n_spin), stat=info)
        call check_allocation(info, 'exmat_diff', func)
     else
        allocate(delta_exmat_complex(n_basis, n_basis, n_k_points_task, n_spin), stat=info)
        call check_allocation(info, 'delta_exmat_complex', func)
        allocate(exmat_diff_complex(n_basis, n_basis, n_k_points_task, n_spin), stat=info)
        call check_allocation(info, 'exmat_diff_complex', func)
     endif
  endif

  ! Please note that below there is no need to make a difference between Scalapack and Lapack
  ! since ONLY array assignments without explicit bounds are used!
  ! Please stay with this convention when making changes in this routine !!!

  if (number_of_loops <= 1) then
     ! On first entry, set D^0_in to D^0_out, because this is still better
     ! than using zero.
     if(real_eigenvectors)then
        prev_exmat = hf_exchange_matr_real
     else
        prev_exmat_complex = hf_exchange_matr_complex
     endif
  end if

  ! delta_exmat := density_matrix - prev_exmat
  if(real_eigenvectors)then
     delta_exmat =  hf_exchange_matr_real - prev_exmat
  else
     delta_exmat_complex = hf_exchange_matr_complex - prev_exmat_complex
  endif

  !    linear mixing
  if(mixer.eq.MIX_LINEAR) then

     if(real_eigenvectors)then
        do i_spin = 1, n_spin
           exmat_diff(:,:,:,i_spin) = linear_mix_param(i_spin) * delta_exmat(:,:,:,i_spin)
        enddo
     else
        do i_spin = 1, n_spin
           exmat_diff_complex(:,:,:,i_spin) = linear_mix_param(i_spin) * delta_exmat_complex(:,:,:,i_spin)
        enddo
     endif

  else if (mixer.eq.MIX_PULAY .and. number_of_loops.le.1) then

     ! I am not completely sure if this should better read linear_mix_param to
     ! match the corresponding charge density Pulay mixer.  But first, that
     ! variable is not guaranteed to be initialized, and second, delta_exmat
     ! is zero anyway.
     if(real_eigenvectors)then
        do i_spin = 1, n_spin
           exmat_diff(:,:,:,i_spin) = charge_mix_param(i_spin) * delta_exmat(:,:,:,i_spin)
        enddo
     else
        do i_spin = 1, n_spin
           exmat_diff_complex(:,:,:,i_spin) = charge_mix_param(i_spin) * delta_exmat_complex(:,:,:,i_spin)
        enddo
     endif

  else if (mixer.eq.MIX_BROYDEN .and. number_of_loops.le.1) then

     ! I am not completely sure if this should better read linear_mix_param to
     ! match the corresponding charge density Broyden mixer.  But first, that
     ! variable is not guaranteed to be initialized, and second, delta_exmat
     ! is zero anyway.
     if(real_eigenvectors)then
        do i_spin = 1, n_spin
           exmat_diff(:,:,:,i_spin) = charge_mix_param(i_spin) * delta_exmat(:,:,:,i_spin)
        enddo
     else
        do i_spin = 1, n_spin
           exmat_diff_complex(:,:,:,i_spin) = charge_mix_param(i_spin) * delta_exmat_complex(:,:,:,i_spin)
        enddo
     endif

  else if (mixer.eq.MIX_PULAY) then

     !   Update the density matrix difference

     if(real_eigenvectors)then
        do i_spin = 1, n_spin
           exmat_diff(:,:,:,i_spin) = charge_mix_param(i_spin) * delta_exmat(:,:,:,i_spin)
        enddo
        
        if(pulay_saved_iter_denmat .gt. 0) then
           
           do i_spin = 1, n_spin
              exmat_diff(:,:,:,i_spin) = exmat_diff(:,:,:,i_spin) + &
                   mixing_factor(1) * &
                   ( prev_exmat_diff(:,:,:,i_spin,1) + &
                   charge_mix_param(i_spin) * ( &
                   delta_exmat(:,:,:,i_spin) - prev_exmat_error(:,:,:,i_spin,1) ) &
                   )
           enddo
        endif
        
        do i_store = 2, pulay_saved_iter_denmat, 1
           do i_spin = 1, n_spin
              exmat_diff(:,:,:,i_spin) = exmat_diff(:,:,:,i_spin) + &
                   mixing_factor(i_store) * &
                   ( prev_exmat_diff(:,:,:,i_spin,i_store) + &
                   charge_mix_param(i_spin) * ( &
                   prev_exmat_error(:,:,:,i_spin,i_store-1) - &
                   prev_exmat_error(:,:,:,i_spin,i_store) ) &
                   )
           enddo
        enddo
     else
        do i_spin = 1, n_spin
           exmat_diff_complex(:,:,:,i_spin) = charge_mix_param(i_spin) * delta_exmat_complex(:,:,:,i_spin)
        enddo

        if(pulay_saved_iter_denmat .gt. 0) then

           do i_spin = 1, n_spin
              exmat_diff_complex(:,:,:,i_spin) = exmat_diff_complex(:,:,:,i_spin) + &
                   mixing_factor(1) * &
                   ( prev_exmat_diff_complex(:,:,:,i_spin,1) + &
                   charge_mix_param(i_spin) * ( &
                   delta_exmat_complex(:,:,:,i_spin) - prev_exmat_error_complex(:,:,:,i_spin,1) ) &
                   )
           enddo
        endif

        do i_store = 2, pulay_saved_iter_denmat, 1
           do i_spin = 1, n_spin
              exmat_diff_complex(:,:,:,i_spin) = exmat_diff_complex(:,:,:,i_spin) + &
                   mixing_factor(i_store) * &
                   ( prev_exmat_diff_complex(:,:,:,i_spin,i_store) + &
                   charge_mix_param(i_spin) * ( &
                   prev_exmat_error_complex(:,:,:,i_spin,i_store-1) - &
                   prev_exmat_error_complex(:,:,:,i_spin,i_store) ) &
                   )
           enddo
        enddo
     endif

  else if (mixer.eq.MIX_BROYDEN) then

     !   Update the density matrix difference

     if(real_eigenvectors)then
        do i_spin = 1, n_spin
           exmat_diff(:,:,:,i_spin) = charge_mix_param(i_spin) * delta_exmat(:,:,:,i_spin)
        enddo
        
        if(pulay_saved_iter_denmat .gt. 0) then
           
           do i_spin = 1, n_spin
              exmat_diff(:,:,:,i_spin) = exmat_diff(:,:,:,i_spin) + &
                   mixing_factor(1) * &
                   ( prev_exmat_diff(:,:,:,i_spin,1) + &
                   charge_mix_param(i_spin) * ( &
                   delta_exmat(:,:,:,i_spin) - prev_exmat_error(:,:,:,i_spin,1) ) &
                   )
           enddo
        endif
        
        do i_store = 2, pulay_saved_iter_denmat, 1
           do i_spin = 1, n_spin
              exmat_diff(:,:,:,i_spin) = exmat_diff(:,:,:,i_spin) + &
                   mixing_factor(i_store) * &
                   ( prev_exmat_diff(:,:,:,i_spin,i_store) + &
                   charge_mix_param(i_spin) * ( &
                   prev_exmat_error(:,:,:,i_spin,i_store-1) - &
                   prev_exmat_error(:,:,:,i_spin,i_store) ) &
                   )
           enddo
        enddo
     else
        do i_spin = 1, n_spin
           exmat_diff_complex(:,:,:,i_spin) = charge_mix_param(i_spin) * delta_exmat_complex(:,:,:,i_spin)
        enddo

        if(pulay_saved_iter_denmat .gt. 0) then

           do i_spin = 1, n_spin
              exmat_diff_complex(:,:,:,i_spin) = exmat_diff_complex(:,:,:,i_spin) + &
                   mixing_factor(1) * &
                   ( prev_exmat_diff_complex(:,:,:,i_spin,1) + &
                   charge_mix_param(i_spin) * ( &
                   delta_exmat_complex(:,:,:,i_spin) - prev_exmat_error_complex(:,:,:,i_spin,1) ) &
                   )
           enddo
        endif

        do i_store = 2, pulay_saved_iter_denmat, 1
           do i_spin = 1, n_spin
              exmat_diff_complex(:,:,:,i_spin) = exmat_diff_complex(:,:,:,i_spin) + &
                   mixing_factor(i_store) * &
                   ( prev_exmat_diff_complex(:,:,:,i_spin,i_store) + &
                   charge_mix_param(i_spin) * ( &
                   prev_exmat_error_complex(:,:,:,i_spin,i_store-1) - &
                   prev_exmat_error_complex(:,:,:,i_spin,i_store) ) &
                   )
           enddo
        enddo
     endif
     
  end if

  ! Store result for next cycle
  if (mixer .eq. MIX_PULAY) then

     if(real_eigenvectors)then
     ! Shift old Pulay stores
        do i_store = pulay_saved_iter_denmat, 2, -1
           
           prev_exmat_error(:,:,:,:,i_store) = &
                prev_exmat_error(:,:,:,:,i_store-1)
           
           prev_exmat_diff(:,:,:,:,i_store) = &
                prev_exmat_diff(:,:,:,:,i_store-1)
           
        enddo

     ! Save this Pulay iteration
        prev_exmat_error(:,:,:,:,1) = delta_exmat(:,:,:,:)
        if (number_of_loops <= 1) then
           do i_spin = 1, n_spin
              prev_exmat_diff(:,:,:,i_spin,1) = delta_exmat(:,:,:,i_spin) * charge_mix_param(i_spin)
           enddo
        else
           prev_exmat_diff(:,:,:,:,1) = exmat_diff(:,:,:,:)
        end if
     else
        ! Shift old Pulay stores
        do i_store = pulay_saved_iter_denmat, 2, -1

           prev_exmat_error_complex(:,:,:,:,i_store) = &
                prev_exmat_error_complex(:,:,:,:,i_store-1)

           prev_exmat_diff_complex(:,:,:,:,i_store) = &
                prev_exmat_diff_complex(:,:,:,:,i_store-1)

        enddo

     ! Save this Pulay iteration
        prev_exmat_error_complex(:,:,:,:,1) = delta_exmat_complex(:,:,:,:)
        if (number_of_loops <= 1) then
           do i_spin = 1, n_spin
              prev_exmat_diff_complex(:,:,:,i_spin,1) = &
                   delta_exmat_complex(:,:,:,i_spin) * charge_mix_param(i_spin)
           enddo
        else
           prev_exmat_diff_complex(:,:,:,:,1) = exmat_diff_complex(:,:,:,:)
        end if
     endif

  endif

  ! Store result for next cycle
  if (mixer .eq. MIX_BROYDEN) then

     if(real_eigenvectors)then
     ! Shift old Pulay stores
        do i_store = pulay_saved_iter_denmat, 2, -1
           
           prev_exmat_error(:,:,:,:,i_store) = &
                prev_exmat_error(:,:,:,:,i_store-1)
           
           prev_exmat_diff(:,:,:,:,i_store) = &
                prev_exmat_diff(:,:,:,:,i_store-1)
           
        enddo

     ! Save this Pulay iteration
        prev_exmat_error(:,:,:,:,1) = delta_exmat(:,:,:,:)
        if (number_of_loops <= 1) then
           do i_spin = 1, n_spin
              prev_exmat_diff(:,:,:,i_spin,1) = delta_exmat(:,:,:,i_spin) * charge_mix_param(i_spin)
           enddo
        else
           prev_exmat_diff(:,:,:,:,1) = exmat_diff(:,:,:,:)
        end if
     else
        ! Shift old Pulay stores
        do i_store = pulay_saved_iter_denmat, 2, -1

           prev_exmat_error_complex(:,:,:,:,i_store) = &
                prev_exmat_error_complex(:,:,:,:,i_store-1)

           prev_exmat_diff_complex(:,:,:,:,i_store) = &
                prev_exmat_diff_complex(:,:,:,:,i_store-1)

        enddo

     ! Save this Pulay iteration
        prev_exmat_error_complex(:,:,:,:,1) = delta_exmat_complex(:,:,:,:)
        if (number_of_loops <= 1) then
           do i_spin = 1, n_spin
              prev_exmat_diff_complex(:,:,:,i_spin,1) = &
                   delta_exmat_complex(:,:,:,i_spin) * charge_mix_param(i_spin)
           enddo
        else
           prev_exmat_diff_complex(:,:,:,:,1) = exmat_diff_complex(:,:,:,:)
        end if
     endif

  endif

  ! Update exchange matrix
  if(real_eigenvectors)then
     prev_exmat = prev_exmat + exmat_diff
     hf_exchange_matr_real = prev_exmat
  else
     prev_exmat_complex = prev_exmat_complex + exmat_diff_complex
     hf_exchange_matr_complex = prev_exmat_complex
  endif

  if(allocated(delta_exmat)) deallocate(delta_exmat)
  if(allocated(exmat_diff)) deallocate(exmat_diff)
  if(allocated(delta_exmat_complex)) deallocate(delta_exmat_complex)
  if(allocated(exmat_diff_complex)) deallocate(exmat_diff_complex)

end subroutine exchange_matrix_mixing_p0

end module hartree_fock_p0
!******
