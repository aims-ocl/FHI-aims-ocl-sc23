!****s* FHI-aims/evaluate_exx_matr_kspace
!  NAME
!   evaluate_exx_matr_kspace
!  SYNOPSIS

      subroutine evaluate_exx_matr_kspace &
           ( n_recip_points, n_recip_points_task, n_low_state, n_high_state, &
             KS_eigenvector_irk, KS_eigenvector_complex_irk, &
             exact_x_kspace &
           )

!  PURPOSE
!  Subroutine "evaluate_exx_matr_kspace" evaluates the 
!  the GW self-energy for a periodic system.
!

! USES
      use dimensions
      use runtime_choices
      use hartree_fock
      use hartree_fock_p0
      use scalapack_wrapper
      use mpi_tasks
      use synchronize_mpi_basic
      use localorb_io, only: use_unit
      use pbc_lists, only: inv_irk_point_mapping
      use aims_memory_tracking

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

!     output
      real*8  :: exact_x_kspace(n_low_state:n_high_state,n_spin,n_recip_points_task)

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
!           
!
! OUTPUT
! o  exact_x_kspace -- the exact-exchange term for each single-particle level
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

      integer, dimension(:), allocatable :: kpoint_task
      integer, dimension(:), allocatable :: id_kpoint
!      real*8, dimension(:,:,:), allocatable :: work_xc
      real*8, dimension(:), allocatable :: exx_real_aux
      real*8, dimension(:,:,:), allocatable :: hf_exchange_tmp
      real*8, dimension(:,:), allocatable :: tmp
      real*8, dimension(:), allocatable :: exx_global
      complex*16, dimension(:), allocatable :: exx_complex_aux
      complex*16, dimension(:,:,:), allocatable :: hf_exchange_complex_tmp
      complex*16, dimension(:,:), allocatable :: tmp_complex


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
      integer :: i_task

!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2A)')"----------------------------------------", &
                  "--------------------------------------"
        write(use_unit,'(2X,2A)') &
         "Start to calculate the exact-exchange", &
         " for each single-particle level ..."
      endif

      if(use_scalapack .and. n_recip_points .le. n_tasks) then
           if(real_eigenvectors) then
              call aims_allocate(tmp,mxld,mxcol,'tmp')
           else
              call aims_allocate(tmp_complex,mxld,mxcol,'tmp_complex')
           endif         
           call aims_allocate(exx_global,n_states,'exx_global')
      else
              if(real_eigenvectors) then
                 allocate(exx_real_aux(n_basis),stat=info)
                 call check_allocation(info, 'exx_real_aux                ')
                 allocate(hf_exchange_tmp(n_basis,n_basis,n_spin),stat=info)
                 call check_allocation(info, 'hf_exchange_tmp            ')
                 allocate(exx_complex_aux(1),stat=info)
                 call check_allocation(info, 'dummy exx_real_aux          ')
                 allocate(hf_exchange_complex_tmp(1,1,1),stat=info)
                 call check_allocation(info, 'dummy hf_exchange_complex_tmp')
              else
                 allocate(exx_complex_aux(n_basis),stat=info)
                 call check_allocation(info, 'exx_complex_aux                ')
                 allocate(hf_exchange_complex_tmp(n_basis,n_basis,n_spin),stat=info)
                 call check_allocation(info, 'hf_exchange_complex_tmp    ')
                 allocate(exx_real_aux(1),stat=info)
                 call check_allocation(info, 'dummy exx_real_aux          ')
                 allocate(hf_exchange_tmp(1,1,1),stat=info)
                 call check_allocation(info, 'dummy hf_exchange_tmp      ')
              endif
      endif


!      write(use_unit,*) "myid, my_k_point", myid, my_k_point

!      exact_x_kspace(:,:,:) = 0.d0

      if(use_scalapack .and. n_recip_points .le. n_tasks) then

          call aims_allocate(kpoint_task, n_tasks, 'kpoint_task')
          call aims_allocate(id_kpoint, n_recip_points, 'id_kpoint')
!     determine the first MPI task that a give k point lives on 
 
          kpoint_task(:) = 0
          do i_task = 1, n_tasks, 1
            if(myid .eq. i_task -1) then
              kpoint_task (i_task) = my_k_point
            endif
          enddo 
          call sync_integer_vector(kpoint_task, n_tasks)

          do i_recip_point = 1, n_recip_points, 1
            if(use_inv_symmetry) then
              i_k_point = inv_irk_point_mapping(i_recip_point)
            else
              i_k_point = i_recip_point
            endif

            do i_task = 1, n_tasks, 1
              if(i_k_point .eq. kpoint_task(i_task)) then
                id_kpoint(i_recip_point) = i_task-1   
                if(i_recip_point .eq. i_task - 1) then
                  exit
                endif
              endif
            enddo
!          write(use_unit,*) "i_recip_point, id_kpoint", i_recip_point, id_kpoint(i_recip_point)
          enddo

          do i_spin = 1, n_spin, 1
            if (real_eigenvectors) then
              call pdgemm('N','N',n_basis,n_states,n_basis,1.d0, hf_exchange_matr_real(1,1,1,i_spin),1,1,sc_desc, &
                                     eigenvec(1,1,i_spin),1,1,sc_desc,0.d0,tmp,1,1,sc_desc)
            else
              call pzgemm('N','N',n_basis,n_states,n_basis,(1.d0,0.d0),hf_exchange_matr_complex(1,1,1,i_spin),1,1,sc_desc, & 
                      eigenvec_complex(1,1,i_spin),1,1,sc_desc,(0.d0,0.d0),tmp_complex,1,1,sc_desc)
            endif

            exx_global(:) = 0.d0
            do i_state = 1, n_states, 1
              if(l_col(i_state) == 0) cycle
              do i_basis_1 = 1, n_basis
                if(l_row(i_basis_1) == 0) cycle
                if(real_eigenvectors)then
                   exx_global(i_state) = exx_global(i_state) - &
                        eigenvec(l_row(i_basis_1),l_col(i_state),i_spin) * tmp (l_row(i_basis_1),l_col(i_state))
                else
                   exx_global(i_state) = exx_global(i_state) - &
                        dble(conjg(eigenvec_complex(l_row(i_basis_1),l_col(i_state),i_spin)) * tmp_complex (l_row(i_basis_1),l_col(i_state)))
!                   write(use_unit,'(2I4,4f18.6)') i_basis_1, i_state, eigenvec_complex(l_row(i_basis_1),l_col(i_state),i_spin), &
!                            hf_exchange_matr_complex(i_basis_1,i_state,1,i_spin)
                endif
              enddo
!              write(use_unit,'(2I4,f18.6)') my_k_point, i_state, exx_global(i_state,i_spin)
            enddo    

            call sync_vector(exx_global,n_states,my_scalapack_comm_all)

            do i_recip_point = n_recip_points, 1, -1
!          do i_recip_point = 1, n_recip_points, 1
                if(use_inv_symmetry) then
                  i_k_point = inv_irk_point_mapping(i_recip_point)
                else
                  i_k_point = i_recip_point
                endif
                i_k_point_local = (i_k_point-1)/n_tasks + 1
                i_recip_point_local = (i_recip_point-1)/n_tasks + 1

                id_send=mod(id_kpoint(i_recip_point),n_tasks)
                id_recv=mod(i_recip_point, n_tasks)

!                write(use_unit,'(A,5I4)') "myid, i_k_point, i_recip_point, id_send, id_recv", myid, i_k_point, i_recip_point, id_send, id_recv
                if((id_send .eq. id_recv).and. (myid.eq.id_recv)) then
                    exact_x_kspace(n_low_state:n_high_state,i_spin,i_recip_point_local) = exx_global(n_low_state:n_high_state) 
                else
                   if(myid .eq. id_send) then
                        call send_real_vector(exx_global(n_low_state), (n_high_state-n_low_state+1), id_recv)
                   elseif(myid .eq. id_recv) then
                        call receive_real_vector(exact_x_kspace(n_low_state,i_spin,i_recip_point_local), (n_high_state-n_low_state+1), id_send)
                   endif
                endif
!                if(myid.eq. id_recv) then
!                    exact_x_kspace(n_low_state:n_high_state,i_spin,i_recip_point_local) = -exx_global(n_low_state:n_high_state,i_spin) 
!                endif
! end loop over i_recip_point
             enddo
! end loop over i_spin
          enddo

      else

!       do i_recip_point = 1, n_recip_points, 1
           do i_recip_point = n_recip_points, 1, -1
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
!                    if(i_recip_point.eq.1) then
!                      do i_basis_1 = 1, n_basis, 1
!                        write(use_unit,'(I4,2f18.6)')i_basis_1, KS_eigenvector_complex_irk(i_basis_1, 20, 1, i_recip_point_local)
!                      enddo
!                    endif
                    
                    do i_spin = 1, n_spin, 1

                      if(real_eigenvectors) then

! get the exact-exchange contribution for each state
                         do i_state = n_low_state, n_high_state, 1
                            call dsymv( 'U', n_basis, 1.d0, &
                                       hf_exchange_tmp(:,:,i_spin), n_basis, &
                                       KS_eigenvector_irk(:,i_state,i_spin,i_recip_point_local), 1, &
                                       0.d0, exx_real_aux, 1)
                             exact_x_kspace(i_state,i_spin,i_recip_point_local) = &
                          -  ddot(n_basis, KS_eigenvector_irk(:,i_state,i_spin,i_recip_point_local), 1, &
                                  exx_real_aux, 1)

                         enddo

                       else

                          do i_state = n_low_state, n_high_state, 1
                             call zhemv( 'U', n_basis,  (1.d0,0.d0), &
                                         hf_exchange_complex_tmp(:,:,i_spin), n_basis, &
                                         KS_eigenvector_complex_irk(:,i_state,i_spin,i_recip_point_local), 1, &
                                        (0.d0,0.d0), exx_complex_aux, 1)

                              exact_x_kspace(i_state,i_spin,i_recip_point_local) = &
                                 -real(zdotc(n_basis, KS_eigenvector_complex_irk(:,i_state,i_spin,i_recip_point_local), &
                                   1,  exx_complex_aux, 1))
                              
!                             do i_basis_1 = 1, n_basis
!                               write(use_unit,'(2I4,2f18.6)') i_basis_1, i_state, KS_eigenvector_complex_irk(i_basis_1,i_state,i_spin,i_recip_point_local)
!                             enddo
!            write(use_unit,'(2I4,f18.6)') i_k_point, i_state, exact_x_kspace(i_state,i_spin,i_recip_point_local)
 
                          enddo

                       endif

! end loop over i_spin
                     enddo

! end loop over i_recip_point
                enddo
!    end of use_scalapack
      endif

      if(allocated(kpoint_task))then
         call aims_deallocate(kpoint_task, 'kpoint_task')
      endif
      if(allocated(id_kpoint))then
        call aims_deallocate(id_kpoint, 'id_kpoint')
      endif
      if(allocated(exx_real_aux)) then
        deallocate(exx_real_aux)
      endif
      if(allocated(exx_complex_aux)) then
        deallocate(exx_complex_aux)
      endif
      if(allocated(hf_exchange_tmp)) then
        deallocate(hf_exchange_tmp)
      endif
      if(allocated(hf_exchange_complex_tmp)) then
        deallocate(hf_exchange_complex_tmp)
      endif
      if(allocated(tmp)) then
        deallocate(tmp)
      endif
      if(allocated(tmp_complex)) then
        deallocate(tmp_complex)
      endif

      end subroutine evaluate_exx_matr_kspace
