!****h* FHI-aims/hartree_fock
!  NAME
!    hartree_fock
!  SYNOPSIS

      module hartree_fock

!  PURPOSE
!  this module contains the declaration of all the variables, arrays that
!  are commonly used in Hartree-Fock calculations, and all other features based
!  on Hartree-Fock (i.e., GW, MP2, hybrid functional calculations)
!  It contains the subroutines for allocating and deallocating these variables.  
!
!  Subroutines:
!  * allocate_hartree_fock
!  * clean_up_hartree_fock
!
!  USES

      use sparse_tensor, only: sp_ten
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
!   n_homo  -- the HOMO level for each spin
!   ovlp_3fn -- real array, the three-center overlap integral (O integral) over
!          two NAO basis functions and one auxiliary basis function. This is the
!          central quantity of the entire formalism of post-Hartree-Fock calculation.
!          RI-SVS: O_SVS = (ij\nu)
!          RI-V:   O_V   = (ij|\nu).
!          At the end of initialize_hartree_fock(), it is overwritten by
!          RI-SVS: M_SVS = O_SVS * S^(-1) * V^(1/2)
!          RI-V:   M_V   = O_V * V^(-1/2)
!          In both cases, the exchange matrix is then
!            X_ii' = (M * M^T)_iji'j' D_i'j' with the density matrix D.
!          Size: (n_basis_pairs, n_loc_prodbas)
!   O_2bs1HF  --   real array, the 3-point integrals over auxiliary (product) basis, 
!          a regular NAO basis, and a single-partilce (KS/HF) orbital. This is used
!          for Hatree-Fock  calculated version 1.
!          RI-SVS: B = M_SVS c
!          RI-V:   B = M_V c = V^{-1/2} (\nu|in\sigma)
!          In both cases, the exchange matrix is then X = B^T * B
!          Size: (n_loc_prodbas, n_basis, n_homo_max, n_spin)
!   coulomb_matr -- First, the integral of the coulomb interaction between two auxiliary
!          basis functions.  Later on (i.e. most of the time), it contains either
!          S^(-1)*V^(1/2)*S^(-1) (for RI-SVS) or V^(-1/2) (for RI-V).
!   ovlp_prodbas --  the overlap matrix of the auxiliary basis functions 
!   fock_matr -- the exchange matrix here, used for Hartree-Fock  and hybrid 
!          functional calculations
!   screx_matr -- screened exchange matrix, and in case of COHSEX matrix, it contains
!          the coulomb-hole plus screened exchange matrix
!   prev_denmat -- the density matrices from previous iteractions, used here for 
!          pulay mixing

      integer ::  n_homo_max
      integer ::  n_loc_states
      integer, dimension(:), allocatable ::  n_homo

      real*8, dimension(:,:), allocatable :: ovlp_3fn
!     Used for LC-wPBEh for SR-Part
      real*8, dimension(:,:), allocatable :: ovlp_3fn_SR
      type(sp_ten) :: coeff_3fn_ten
      real*8, dimension(:, :), allocatable :: coulomb_matr_lvl
      real*8, dimension(:,:,:,:), allocatable :: O_2bs1HF
      real*8, dimension(:,:,:), allocatable :: fock_matr
!     Used for LC-wPBEh for SR-Part
      real*8, dimension(:,:,:), allocatable :: fock_matr_SR
      real*8, dimension(:,:,:), allocatable :: screx_matr

      integer :: n_nonredundant_denmat, n_local_denmat
      real*8, dimension(:,:), allocatable :: prev_denmat
      real*8, dimension(:,:,:), allocatable :: prev_denmat_diff
      real*8, dimension(:,:,:), allocatable :: prev_denmat_error

!     periodic Hartree Fock in kspace
      integer ::  n_cells_task
      integer, dimension(:,:), allocatable :: kq_point_list
      integer, dimension(:,:), allocatable :: kpq_point_list
      !BLACS parallel
      complex*16, dimension(:,:,:,:,:), allocatable :: lvl_tricoeff_mod_r, lvl_tricoeff_mod_c
      complex*16, dimension(:,:,:), allocatable :: coulomb_matr_blacs
      complex*16, dimension(:,:,:), allocatable :: coulomb_cut_blacs
      !unparallelized
      complex*16, dimension(:,:,:,:), allocatable :: lvl_tricoeff_recip1
      real*8, dimension(:,:,:), allocatable :: lvl_tricoeff_recip2
      complex*16, dimension(:,:,:), allocatable :: coulomb_matr_recip


      contains
        subroutine allocate_hartree_fock &
              ( )
!
!  allocate the arrays needed in Hartree-Fock calculation
!
         use dimensions
         use prodbas
         use runtime_choices
         use mpi_tasks
         use pbc_lists, only: n_cells
         implicit none

         integer :: info
         character(*), parameter :: func = 'allocate_hartree_fock'

         if (.not.allocated(n_homo)) then
            allocate (n_homo(n_spin))
         endif
         if (use_screx) then
            n_loc_states = (n_states-1)/n_tasks + 1
            if (allocated(O_2bs1HF)) deallocate(O_2bs1HF)
            allocate (O_2bs1HF(n_basbas, n_basis, n_loc_states, n_spin), &
            &         stat=info)
            call check_allocation(info, 'O_2bs1HF', func)
         endif

         if (.not.allocated(fock_matr)) then
           allocate (fock_matr(n_basis, n_basis,n_spin))
         endif
         if(use_lc_wpbeh .and. hybrid_coeff .ne. 0.0d0) then
           if (.not.allocated(fock_matr_SR)) then
           allocate (fock_matr_SR(n_basis, n_basis,n_spin))
         endif
         endif
         if(use_screx) then
           if (.not.allocated(screx_matr)) then
             allocate (screx_matr(n_basis, n_basis,n_spin))
           endif
         endif
         if(use_density_matrix_hf) then
            n_nonredundant_denmat = n_basis*(n_basis+1)/2
            n_local_denmat = n_nonredundant_denmat / n_tasks + 1   ! Save value
            if (.not.allocated(prev_denmat)) then
               allocate(prev_denmat(n_local_denmat, n_spin), stat=info)
               call check_allocation(info, 'prev_denmat', func)
            endif
            prev_denmat(:,:) = 0.d0
            if(mixer.eq.MIX_PULAY) then
               if(.not.allocated(prev_denmat_diff)) then
                  allocate(prev_denmat_diff(n_local_denmat, &
                  &                         n_max_pulay,n_spin), stat=info)
                  call check_allocation(info, 'prev_denmat_diff', func)
               endif
               prev_denmat_diff(:,:,:) = 0.d0
               if(.not.allocated(prev_denmat_error)) then
                  allocate(prev_denmat_error(n_local_denmat, &
                  &                         n_max_pulay,n_spin), stat=info)
                  call check_allocation(info, 'prev_denmat_error', func)
               endif
               prev_denmat_error(:,:,:) = 0.d0
            endif
            if(mixer.eq.MIX_BROYDEN) then
               if(.not.allocated(prev_denmat_diff)) then
                  allocate(prev_denmat_diff(n_local_denmat, &
                  &                         n_max_pulay,n_spin), stat=info)
                  call check_allocation(info, 'prev_denmat_diff', func)
               endif
               prev_denmat_diff(:,:,:) = 0.d0
               if(.not.allocated(prev_denmat_error)) then
                  allocate(prev_denmat_error(n_local_denmat, &
                  &                         n_max_pulay,n_spin), stat=info)
                  call check_allocation(info, 'prev_denmat_error', func)
               endif
               prev_denmat_error(:,:,:) = 0.d0
            endif
         endif
         if(use_hf_kspace) then
         ! estimate the largest number of real-space unit cells locally
           if(mod(n_cells, n_tasks).eq.0) then
              n_cells_task = n_cells/n_tasks
           else
              n_cells_task = n_cells/n_tasks+1
           endif

           allocate(lvl_tricoeff_recip1(max_n_basbas_sp,n_basis,n_basis,n_k_points_task),stat=info) 
           call check_allocation(info, 'lvl_tricoeff_recip1', func, max_n_basbas_sp,n_basis,n_basis,n_k_points_task)
           allocate(lvl_tricoeff_recip2(max_n_basbas_sp,n_basis,n_basis),stat=info) 
           call check_allocation(info, 'lvl_tricoeff_recip2', func, max_n_basbas_sp,n_basis,n_basis)
           allocate(kq_point_list(n_k_points,n_k_points),stat=info) 
           call check_allocation(info, 'kq_point_list', func, n_k_points,n_k_points)
           allocate(kpq_point_list(n_k_points,n_k_points),stat=info) 
           call check_allocation(info, 'kpq_point_list', func, n_k_points,n_k_points)
           allocate(coulomb_matr_recip(n_basbas,n_basbas,n_k_points_task),stat=info) 
           call check_allocation(info, 'coulomb_matr_recip', func, n_basbas,n_basbas,n_k_points_task)

         endif


        end  subroutine allocate_hartree_fock

!   deallocate the arrays
        subroutine cleanup_hartree_fock ()

         use dimensions, only: use_hf_kspace
         use sparse_tensor, only: dealloc_sp_ten
         implicit none

         if (allocated(prev_denmat)) then
           deallocate(prev_denmat)
         endif
         if (allocated(prev_denmat_diff)) then
           deallocate(prev_denmat_diff)
         endif
         if (allocated(prev_denmat_error)) then
           deallocate(prev_denmat_error)
         endif
         if (allocated(O_2bs1HF)) then
           deallocate(O_2bs1HF)
         endif
         if (allocated(fock_matr)) then
           deallocate(fock_matr)
         endif
         if (allocated(fock_matr_SR)) then
           deallocate(fock_matr_SR)
         endif
         if (allocated(screx_matr)) then
           deallocate(screx_matr)
         endif
         call dealloc_sp_ten(coeff_3fn_ten)
         if (allocated(coulomb_matr_lvl)) deallocate(coulomb_matr_lvl)
         !IYZ: decide it outside the subroutine
         !if(.not.use_rpa_ene)then
         !  if (allocated(ovlp_3fn)) deallocate(ovlp_3fn)
         !  if (allocated(ovlp_3fn_SR)) deallocate(ovlp_3fn_SR)
         !endif
         if (allocated(ovlp_3fn)) deallocate(ovlp_3fn)
         if (allocated(ovlp_3fn_SR)) deallocate(ovlp_3fn_SR)

         if(allocated(lvl_tricoeff_recip1)) then
           deallocate(lvl_tricoeff_recip1)
         endif
         if(allocated(lvl_tricoeff_recip2)) then
           deallocate(lvl_tricoeff_recip2)
         endif
         if(allocated(kq_point_list)) then
           deallocate(kq_point_list)
         endif
         if(allocated(kpq_point_list)) then
           deallocate(kpq_point_list)
         endif
         if(allocated(coulomb_matr_recip)) then
           deallocate(coulomb_matr_recip)
         endif
        end  subroutine cleanup_hartree_fock
        
!      BSSE       
        subroutine allocate_hartree_fock_bsse &
              ( )
              
         use dimensions
         use prodbas
         use runtime_choices
         use mpi_tasks, only: check_allocation, n_tasks
         implicit none

         integer :: info
         character(*), parameter :: func = 'allocate_hartree_fock_bsse'
         

         if (allocated(O_2bs1HF)) O_2bs1HF(:,:,:,:)=0.d0 
         if (allocated(screx_matr)) screx_matr(:,:,:)=0.d0

          if(use_density_matrix_hf) then
            n_nonredundant_denmat = n_basis*(n_basis+1)/2
            n_local_denmat = n_nonredundant_denmat / n_tasks + 1   ! Save value
            if (.not.allocated(prev_denmat)) then
               allocate(prev_denmat(n_local_denmat, n_spin), stat=info)
               call check_allocation(info, 'prev_denmat', func)
            endif
            prev_denmat(:,:) = 0.d0
            if(mixer.eq.MIX_PULAY) then
               if(.not.allocated(prev_denmat_diff)) then
                  allocate(prev_denmat_diff(n_local_denmat, &
                  &                         n_max_pulay,n_spin), stat=info)
                  call check_allocation(info, 'prev_denmat_diff', func)
               endif
               prev_denmat_diff(:,:,:) = 0.d0
               if(.not.allocated(prev_denmat_error)) then
                  allocate(prev_denmat_error(n_local_denmat, &
                  &                         n_max_pulay,n_spin), stat=info)
                  call check_allocation(info, 'prev_denmat_error', func)
               endif
               prev_denmat_error(:,:,:) = 0.d0
            endif
            if(mixer.eq.MIX_BROYDEN) then
               if(.not.allocated(prev_denmat_diff)) then
                  allocate(prev_denmat_diff(n_local_denmat, &
                  &                         n_max_pulay,n_spin), stat=info)
                  call check_allocation(info, 'prev_denmat_diff', func)
               endif
               prev_denmat_diff(:,:,:) = 0.d0
               if(.not.allocated(prev_denmat_error)) then
                  allocate(prev_denmat_error(n_local_denmat, &
                  &                         n_max_pulay,n_spin), stat=info)
                  call check_allocation(info, 'prev_denmat_error', func)
               endif
               prev_denmat_error(:,:,:) = 0.d0
            endif
         endif
         end  subroutine allocate_hartree_fock_bsse        

      end module hartree_fock
!******
