!****s* FHI-aims/evaluate_ex_and_xc_matr_kspace
!  NAME
!   evaluate_ex_and_xc_matr_kspace
!  SYNOPSIS

      subroutine evaluate_ex_and_xc_matr_kspace &
           ( n_recip_points, n_recip_points_task, n_low_state, n_high_state, &
             KS_eigenvector_irk, KS_eigenvector_complex_irk, &
             xc_realspace, &
             exact_x_kspace, xc_kspace &
           )

!  PURPOSE
!  Subroutine "evaluate_ex_and_xc_matr_kspace" evaluates the 
!  the GW self-energy for a periodic system.
!

! USES
      use dimensions
      use runtime_choices
      use hartree_fock
      use hartree_fock_p0
      use mpi_tasks
      use synchronize_mpi_basic
      use localorb_io, only: use_unit
      use pbc_lists, only: inv_irk_point_mapping

      implicit none

! ARGUMENTS 

      integer :: n_recip_points
      integer :: n_recip_points_task
      integer :: n_low_state
      integer :: n_high_state
      integer :: n_freq
      integer :: n_full_freq
      real*8  :: KS_eigenvector_irk(n_basis,n_states,n_spin,n_recip_points_task)
      complex*16  :: KS_eigenvector_complex_irk(n_basis,n_states,n_spin,n_recip_points_task)
      real*8  :: xc_realspace(n_hamiltonian_matrix_size,n_spin)

!     output
      real*8  :: exact_x_kspace(n_low_state:n_high_state,n_spin,n_recip_points_task)
      real*8  :: xc_kspace(n_low_state:n_high_state,n_spin,n_recip_points_task)

! INPUTS
! o  n_recip_points -- number of k points (not necessarily on the regular k mesh)
! o  n_recip_points_task -- number of k points per CPU core 
! o  n_low_state  -- integer number,
!            the lowest KS/HF eigenstate taken into account in the GW self-energy calculations
! o  n_high_state -- integer number,
!            the highest KS/HF eigenstate taken into account
! o  KS_eigenvector_irk -- real array,
!            the eigenvector of the single-particle calculation
! o  KS_eigenvector_complex_irk -- complex array,
!            the complex eigenvector of the single-particle calculation,
!            used when "real_eigenvectors == .false."
! o  xc_realspace -- the XC contribution to the DFT Hamiltonian in realspace
!           
!
! OUTPUT
! o  exact_x_kspace -- the exact-exchange term for each single-particle level
! o  xc_kspace -- the DFT-XC term for each single-particle level
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

!  local variables

!     auxiliary matrices for Level 3 Blas matrix multiplications

      integer :: info
      integer :: n_workxc

      real*8, dimension(:,:,:), allocatable :: work_xc
      real*8, dimension(:,:), allocatable :: xc_matr
      real*8, dimension(:,:), allocatable :: xc_matr_w
      real*8, dimension(:), allocatable :: xc_real_tmp
      real*8, dimension(:,:,:), allocatable :: hf_exchange_tmp
      complex*16, dimension(:,:), allocatable :: xc_matr_complex
      complex*16, dimension(:,:), allocatable :: xc_matr_w_complex
      complex*16, dimension(:), allocatable :: xc_complex_tmp
      complex*16, dimension(:,:,:), allocatable :: hf_exchange_complex_tmp


!     external functions
      real*8  :: ddot
      complex*16  :: zdotc

!     counters
      integer :: i_state
      integer :: i_spin
      integer :: i_k_point
      integer :: i_k_point_local
      integer :: i_recip_point
      integer :: i_recip_point_local
      integer :: i_basis_1, i_basis_2
      integer :: i_index
      integer :: id_send, id_recv

!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2A)')"----------------------------------------", &
                  "--------------------------------------"
        write(use_unit,'(2X,2A)') &
         "Start to calculate the exact-exchange and DFT xc contribution", &
         " for each single-particle level ..."
      endif

      if(packed_matrix_format==PM_none)then
          allocate(work_xc(1:n_centers_basis_I, 1:n_centers_basis_I, n_spin),stat=info)
          call check_allocation(info, 'work_xc                      ')
          n_workxc=n_centers_basis_I
      else
          allocate(work_xc(1, 1, 1),stat=info)
          call check_allocation(info, 'work_xc                      ')
          n_workxc=1
      end if

      if(real_eigenvectors) then
         allocate(xc_matr_w(n_basis*(n_basis+1)/2, n_spin),stat=info)
         call check_allocation(info, 'xc_matr_w             ')
         allocate(xc_matr(n_basis, n_basis),stat=info)
         call check_allocation(info, 'xc_matr               ')
         allocate(xc_real_tmp(n_basis),stat=info)
         call check_allocation(info, 'xc_real_tmp                ')
         allocate(hf_exchange_tmp(n_basis,n_basis,n_spin),stat=info)
         call check_allocation(info, 'hf_exchange_tmp            ')
         !
         ! Require dummy allocations of corresponding complex arrays to avoid
         ! triggering -check bounds error for unallocated array in subroutine call.
         ! The dummy arrays should never be used.
         allocate(xc_matr_w_complex(1,1),stat=info)
         call check_allocation(info, 'dummy xc_matr_w       ')
         allocate(xc_matr_complex(1, 1),stat=info)
         call check_allocation(info, 'dummy xc_matr         ')
         allocate(xc_complex_tmp(1),stat=info)
         call check_allocation(info, 'dummy xc_real_tmp          ')
         allocate(hf_exchange_complex_tmp(1,1,1),stat=info)
         call check_allocation(info, 'dummy hf_exchange_complex_tmp')
      else
         allocate(xc_matr_w_complex(n_basis*(n_basis+1)/2, n_spin),stat=info)
         call check_allocation(info, 'xc_matr_w             ')
         allocate(xc_matr_complex(n_basis, n_basis),stat=info)
         call check_allocation(info, 'xc_matr               ')
         allocate(xc_complex_tmp(n_basis),stat=info)
         call check_allocation(info, 'xc_real_tmp                ')
         allocate(hf_exchange_complex_tmp(n_basis,n_basis,n_spin),stat=info)
         call check_allocation(info, 'hf_exchange_complex_tmp    ')
         !
         ! Require dummy allocations of corresponding real arrays to avoid
         ! triggering -check bounds error for unallocated array in subroutine call.
         ! The dummy arrays should never be used.
         allocate(xc_matr_w(1, 1),stat=info)
         call check_allocation(info, 'dummy xc_matr_w       ')
         allocate(xc_matr(1, 1),stat=info)
         call check_allocation(info, 'dummy xc_matr         ')
         allocate(xc_real_tmp(1),stat=info)
         call check_allocation(info, 'dummy xc_real_tmp          ')
         allocate(hf_exchange_tmp(1,1,1),stat=info)
         call check_allocation(info, 'dummy hf_exchange_tmp      ')
      endif

      do i_recip_point = 1, n_recip_points, 1
        if(use_inv_symmetry) then
         i_k_point = inv_irk_point_mapping(i_recip_point)
        else
         i_k_point = i_recip_point
        endif
         i_k_point_local = (i_k_point-1)/n_tasks + 1
      
         id_send=mod(i_k_point, n_tasks)
         id_recv=mod(i_recip_point, n_tasks)
         if(myid.eq.id_send) then
           if(real_eigenvectors) then
              hf_exchange_tmp(:,:,:) = hf_exchange_matr_real(:,:,i_k_point_local,:)
           else
              hf_exchange_complex_tmp(:,:,:) = hf_exchange_matr_complex(:,:,i_k_point_local,:)
           endif
         endif

         if(id_send .ne. id_recv) then
           if(real_eigenvectors) then
              if(myid .eq. id_send) then
                  call send_real_vector(hf_exchange_tmp, n_basis*n_basis*n_spin, id_recv)
              elseif(myid .eq. id_recv) then
                  call receive_real_vector(hf_exchange_tmp, n_basis*n_basis*n_spin, id_send)
              endif
            else
              if(myid .eq. id_send) then
                  call send_complex_vector(hf_exchange_complex_tmp, n_basis*n_basis*n_spin, id_recv)
              elseif(myid .eq. id_recv) then
                  call receive_complex_vector(hf_exchange_complex_tmp, n_basis*n_basis*n_spin, id_send)
              endif
            endif
         endif
         
         if(myid .ne. mod(i_recip_point, n_tasks)) cycle
           i_recip_point_local = (i_recip_point-1)/n_tasks + 1

!  construct the XC matrix for each k point
          call construct_xc_matr_kspace &
          ( xc_realspace, xc_matr_w, xc_matr_w_complex, i_k_point, n_workxc, work_xc)

          do i_spin = 1, n_spin, 1

            if(real_eigenvectors) then

               i_index=0
               do i_basis_1 = 1, n_basis, 1
                  do i_basis_2 = 1, i_basis_1, 1
                    i_index = i_index + 1
                    xc_matr(i_basis_2,i_basis_1) = xc_matr_w(i_index,i_spin)
!                    xc_matr(i_basis_1,i_basis_2) = xc_matr_w(i_index,i_spin)
                  enddo
               enddo

! get the XC contribution for each state
               do i_state = n_low_state, n_high_state, 1
                 call dsymv( 'U', n_basis, 1.d0, xc_matr, n_basis, &
                             KS_eigenvector_irk(:,i_state,i_spin,i_recip_point_local), 1, &
                             0.d0, xc_real_tmp, 1)
                 xc_kspace(i_state,i_spin,i_recip_point_local) = &
                 ddot(n_basis, KS_eigenvector_irk(:,i_state,i_spin,i_recip_point_local), 1, &
                      xc_real_tmp, 1)
               enddo

! get the exact-exchange contribution for each state
               do i_state = n_low_state, n_high_state, 1
                 call dsymv( 'U', n_basis, 1.d0, &
                             hf_exchange_tmp(:,:,i_spin), n_basis, &
                             KS_eigenvector_irk(:,i_state,i_spin,i_recip_point_local), 1, &
                             0.d0, xc_real_tmp, 1)
                 exact_x_kspace(i_state,i_spin,i_recip_point_local) = &
                 - ddot(n_basis, KS_eigenvector_irk(:,i_state,i_spin,i_recip_point_local), 1, &
                        xc_real_tmp, 1)
               enddo

!               do i_basis_1=1, n_basis, 1
!                 do i_basis_2=1, n_basis, 1
!                   write(use_unit,'(2I4,f16.8)') i_basis_2, i_basis_1, xc_matr(i_basis_2, i_basis_1)
!                 enddo
!               enddo
            else

               i_index=0
               do i_basis_1 = 1, n_basis, 1
                  do i_basis_2 = 1, i_basis_1, 1
                    i_index = i_index + 1
                    xc_matr_complex(i_basis_2,i_basis_1) = xc_matr_w_complex(i_index,i_spin)
!                    xc_matr_complex(i_basis_1,i_basis_2) = conjg(xc_matr_w_complex(i_index,i_spin))
                  enddo
               enddo

               do i_state = n_low_state, n_high_state, 1
                 call zhemv( 'U', n_basis, (1.d0,0.d0), xc_matr_complex, n_basis, &
                             KS_eigenvector_complex_irk(:,i_state,i_spin,i_recip_point_local), 1, &
                             (0.d0,0.d0), xc_complex_tmp, 1)

!                 call zgemv( 'N', n_basis, n_basis, (1.d0,0.d0), xc_matr_complex, n_basis, &
!                             KS_eigenvector_complex_irk(:,i_state,i_spin,i_recip_point_local), 1, &
!                             (0.d0,0.d0), xc_complex_tmp, 1)

                 xc_kspace(i_state,i_spin,i_recip_point_local) = &
                 real(zdotc(n_basis, KS_eigenvector_complex_irk(:,i_state,i_spin,i_recip_point_local), &
                      1,  xc_complex_tmp, 1))

!                 write(use_unit,'(2I4,2f16.8)') i_recip_point, i_state, xc_kspace(i_state,i_spin,i_recip_point_local), &
!                       xc_kspace(i_state,i_spin,i_recip_point_local)*27.2113845
               enddo

               do i_state = n_low_state, n_high_state, 1
                 call zhemv( 'U', n_basis,  (1.d0,0.d0), &
                             hf_exchange_complex_tmp(:,:,i_spin), n_basis, &
                             KS_eigenvector_complex_irk(:,i_state,i_spin,i_recip_point_local), 1, &
                             (0.d0,0.d0), xc_complex_tmp, 1)
!                 call zgemv( 'N', n_basis, n_basis,  (1.d0,0.d0), &
!                             hf_exchange_matr_complex(:,:,i_recip_point_local,i_spin), n_basis, &
!                             KS_eigenvector_complex_irk(:,i_state,i_spin,i_recip_point_local), 1, &
!                             (0.d0,0.d0), xc_complex_tmp, 1)

                 exact_x_kspace(i_state,i_spin,i_recip_point_local) = &
                 -real(zdotc(n_basis, KS_eigenvector_complex_irk(:,i_state,i_spin,i_recip_point_local), &
                      1,  xc_complex_tmp, 1))

!                 write(use_unit,'(2I4,2f16.8)') i_recip_point, i_state, exact_x_kspace(i_state,i_spin,i_recip_point_local), &
!                       exact_x_kspace(i_state,i_spin,i_recip_point_local)*27.2113845
               enddo

            endif

! end loop over i_spin
          enddo

! end loop over i_recip_point
      enddo
      if(allocated(work_xc)) then
        deallocate(work_xc)
      endif
      if(allocated(xc_matr_w)) then
        deallocate(xc_matr_w)
      endif
      if(allocated(xc_matr)) then
        deallocate(xc_matr)
      endif
      if(allocated(xc_matr_w_complex)) then
        deallocate(xc_matr_w_complex)
      endif
      if(allocated(xc_matr_complex)) then
        deallocate(xc_matr_complex)
      endif
      if(allocated(xc_real_tmp)) then
        deallocate(xc_real_tmp)
      endif
      if(allocated(xc_complex_tmp)) then
        deallocate(xc_complex_tmp)
      endif
      if(allocated(hf_exchange_tmp)) then
        deallocate(hf_exchange_tmp)
      endif
      if(allocated(hf_exchange_complex_tmp)) then
        deallocate(hf_exchange_complex_tmp)
      endif

      end subroutine evaluate_ex_and_xc_matr_kspace
