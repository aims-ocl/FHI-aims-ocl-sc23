!****s* FHI-aims/orthonormalize_eigenvectors
!  NAME
!   orthonormalize_eigenvectors
!  SYNOPSIS
subroutine orthonormalize_eigenvectors( )

!  PURPOSE
! modified Gram-Schmidt orthonormalization
!
!  USES
  use mpi_tasks
  use localorb_io
  use dimensions
  use physics
  use inner_product
  use aims_memory_tracking, only : aims_allocate, aims_deallocate

  implicit none

!  ARGUMENTS
  real*8, dimension(:), allocatable :: overlap_matrix_w
  complex*16, dimension(:), allocatable :: overlap_matrix_w_complex
  real*8, dimension(:,:),allocatable :: work_ovl

  real*8 :: r_ii, r_ij
  integer :: i_state, i_state_2, i_spin, i_k_point, i_k
  integer :: i, n, info
  real*8, allocatable :: ovlp_full(:,:), ovlp_ev(:,:), ev_ovlp_ev(:,:)
  complex*16, allocatable :: ovlp_full_complex(:,:), ovlp_ev_complex(:,:), ev_ovlp_ev_complex(:,:)

  character*100 :: info_str
!  INPUTS
!  none
!  OUTPUTS
!  orthonormalized eigenvectors
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals:
!    FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject
!   to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!  This routine uses the Loewdin transform for orthogonalization:
!
!  Given a matrix of vectors W, calculate M = W^T S W and
!  the Cholesky factorisation of M = LL^T (or M = U^T U).
!  Set Z = W L^-T, then Z^T S Z = I.



  write(info_str,'(2X,A)') "Orthonormalizing eigenvectors"
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  if (real_eigenvectors) then
     call aims_allocate( overlap_matrix_w, n_basis*(n_basis+1)/2,                 "overlap_matrix_w" )
     ! dummy complex matrix only
     call aims_allocate( overlap_matrix_w_complex, 1,                     "overlap_matrix_w_complex" )
  else
     call aims_allocate( overlap_matrix_w_complex, n_basis*(n_basis+1)/2, "overlap_matrix_w_complex" )
     ! dummy real matrix only
     call aims_allocate( overlap_matrix_w, 1,                                     "overlap_matrix_w" )
  end if

  if (n_periodic == 0) then

     call aims_allocate( ovlp_full, n_basis, n_basis,                                    "ovlp_full" )
     call aims_allocate( ovlp_ev, n_basis, n_states,                                       "ovlp_ev" )
     call aims_allocate( ev_ovlp_ev, n_states, n_states,                                "ev_ovlp_ev" )

     ! Store full overlap matrix in ovlp_full

     n = 0
     do i = 1, n_basis
        ovlp_full(1:i,i) = overlap_matrix(n+1:n+i)
        ovlp_full(i,1:i-1) = overlap_matrix(n+1:n+i-1)
        n = n+i
     enddo

     do i_spin = 1, n_spin, 1

        ! ovlp_ev = ovlp_full * KS_eigenvector
        call dgemm('N','N', n_basis,n_states,n_basis, 1.d0, &
                ovlp_full,ubound(ovlp_full,1), &
                KS_eigenvector(1,1,i_spin,1),ubound(KS_eigenvector,1), &
                0.d0, ovlp_ev,ubound(ovlp_ev,1))


        ! ev_ovlp_ev = KS_eigenvector**T * ovlp_ev
        call dgemm('T','N', n_states,n_states,n_basis, 1.d0, &
                KS_eigenvector(1,1,i_spin,1),ubound(KS_eigenvector,1), &
                ovlp_ev,ubound(ovlp_ev,1), &
                0.d0, ev_ovlp_ev,ubound(ev_ovlp_ev,1))

        ! Cholesky decomposition of ev_ovlp_ev
        ! The matrix should never be singular because it is an overlap matrix,
        ! the test of info is for safety only.
        call dpotrf('U',n_states,ev_ovlp_ev,ubound(ev_ovlp_ev,1),info)
        if(info/=0) then
           print *,'dpotrf info = ',info
           stop
        endif

        ! Invert upper Cholesky factor, set lower half to 0
        ! The matrix should never be singular, the test of info is for safety only

        call dtrtri('U','N',n_states,ev_ovlp_ev,ubound(ev_ovlp_ev,1),info)
        if(info/=0) then
           print *,'dtrtri info = ',info
           stop
        endif
        do i = 1, n_states-1
           ev_ovlp_ev(i+1:n_states,i) = 0
        enddo

        ! Orthogonalized eigenvectors = KS_eigenvector * U**-1
        call dgemm('N','N', n_basis,n_states,n_states, 1.d0, &
                KS_eigenvector(1,1,i_spin,1),ubound(KS_eigenvector,1), &
                ev_ovlp_ev,ubound(ev_ovlp_ev,1), &
                0.d0, ovlp_ev,ubound(ovlp_ev,1))
        KS_eigenvector(1:n_basis,1:n_states,i_spin,1) = ovlp_ev(1:n_basis,1:n_states)

     end do

     call aims_deallocate( ovlp_full,   "ovlp_full" )
     call aims_deallocate( ovlp_ev,       "ovlp_ev" )
     call aims_deallocate( ev_ovlp_ev, "ev_ovlp_ev" )

  else ! periodic system

     i_k = 0
     if (real_eigenvectors) then

        call aims_allocate(ovlp_full, n_basis, n_basis,     "ovlp_full" )
        call aims_allocate(ovlp_ev, n_basis, n_states,        "ovlp_ev" )
        call aims_allocate(ev_ovlp_ev, n_states, n_states, "ev_ovlp_ev" )

        do i_k_point = 1, n_k_points, 1

           if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then

              i_k = i_k + 1

              call construct_overlap( overlap_matrix, overlap_matrix_w, &
                   overlap_matrix_w_complex, i_k_point )

              ! Store full overlap matrix in ovlp_full

              n = 0
              do i = 1, n_basis
                 ovlp_full(1:i,i) = overlap_matrix_w(n+1:n+i)
                 ovlp_full(i,1:i-1) = overlap_matrix_w(n+1:n+i-1)
                 n = n+i
              enddo

              do i_spin = 1, n_spin, 1

                 ! ovlp_ev = ovlp_full * KS_eigenvector
                 call dgemm('N','N', n_basis,n_states,n_basis, 1.d0, &
                         ovlp_full,ubound(ovlp_full,1), &
                         KS_eigenvector(1,1,i_spin,i_k),&
                         ubound(KS_eigenvector,1), &
                         0.d0, ovlp_ev,ubound(ovlp_ev,1))

                 ! ev_ovlp_ev = KS_eigenvector**T * ovlp_ev
                 call dgemm('T','N', n_states,n_states,n_basis, 1.d0, &
                         KS_eigenvector(1,1,i_spin,i_k), &
                         ubound(KS_eigenvector,1), &
                         ovlp_ev,ubound(ovlp_ev,1), &
                         0.d0, ev_ovlp_ev,ubound(ev_ovlp_ev,1))


                 ! Cholesky decomposition of ev_ovlp_ev
                 ! The matrix should never be singular, the test of info is for safety only

                 call dpotrf('U',n_states,ev_ovlp_ev,ubound(ev_ovlp_ev,1),info)
                 if(info/=0) then
                    print *,'dpotrf info = ',info
                    stop
                 endif

                 ! Invert upper Cholesky factor, set lower half to 0
                 ! The matrix should never be singular, the test of info is for safety only

                 call dtrtri('U','N',n_states,ev_ovlp_ev,ubound(ev_ovlp_ev,1),info)
                 if(info/=0) then
                    print *,'dtrtri info = ',info
                    stop
                 endif
                 do i = 1, n_states-1
                    ev_ovlp_ev(i+1:n_states,i) = 0
                 enddo

                 ! Orthogonalized eigenvectors = KS_eigenvector * U**-1
                 call dgemm('N','N', n_basis,n_states,n_states, 1.d0, &
                         KS_eigenvector(1,1,i_spin,i_k), &
                         ubound(KS_eigenvector,1), &
                         ev_ovlp_ev,ubound(ev_ovlp_ev,1), &
                         0.d0, ovlp_ev,ubound(ovlp_ev,1))

                 KS_eigenvector(1:n_basis,1:n_states,i_spin,i_k) = &
                    ovlp_ev(1:n_basis,1:n_states)

              end do
           end if

        end do

        call aims_deallocate( ovlp_full,   "ovlp_full" )
        call aims_deallocate( ovlp_ev,       "ovlp_ev" )
        call aims_deallocate( ev_ovlp_ev, "ev_ovlp_ev" )

     else ! complex eigenvectors

        call aims_allocate( ovlp_full_complex, n_basis, n_basis,     "ovlp_full_complex" )
        call aims_allocate( ovlp_ev_complex, n_basis, n_states,        "ovlp_ev_complex" )
        call aims_allocate( ev_ovlp_ev_complex, n_states, n_states, "ev_ovlp_ev_complex" )

        if(packed_matrix_format==PM_none)then
           call aims_allocate(work_ovl, n_centers_basis_I, n_centers_basis_I, "work_ovl" )
        else
           ! dummy allocation for the normal packed matrix case
           call aims_allocate(work_ovl, 1, 1,                                 "work_ovl" )
        end if

        do i_k_point = 1, n_k_points, 1

           if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then

              i_k = i_k + 1

              call construct_overlap( overlap_matrix, overlap_matrix_w, &
                   overlap_matrix_w_complex, i_k_point, work_ovl )

              ! Store full overlap matrix in ovlp_full_complex

              n = 0
              do i = 1, n_basis
                 ovlp_full_complex(1:i,i) = overlap_matrix_w_complex(n+1:n+i)
                 ovlp_full_complex(i,1:i-1) = conjg(overlap_matrix_w_complex(n+1:n+i-1))
                 n = n+i
              enddo

              do i_spin = 1, n_spin, 1

                 ! ovlp_ev_complex = ovlp_full_complex * KS_eigenvector_complex

                 call zgemm('N','N', n_basis,n_states,n_basis, (1.d0,0.d0), &
                         ovlp_full_complex,ubound(ovlp_full_complex,1), &
                         KS_eigenvector_complex(1,1,i_spin,i_k),ubound(KS_eigenvector_complex,1), &
                         (0.d0,0.d0), ovlp_ev_complex,ubound(ovlp_ev_complex,1))

                 ! ev_ovlp_ev_complex = KS_eigenvector_complex**H * ovlp_ev_complex

                 call zgemm('C','N', n_states,n_states,n_basis, (1.d0,0.d0), &
                          KS_eigenvector_complex(1,1,i_spin,i_k),ubound(KS_eigenvector_complex,1), &
                          ovlp_ev_complex,ubound(ovlp_ev_complex,1), &
                          (0.d0,0.d0), ev_ovlp_ev_complex,ubound(ev_ovlp_ev_complex,1))

                 ! Cholesky decomposition of ev_ovlp_ev_complex
                 ! The matrix should never be singular, the test of info is for safety only

                 call zpotrf('U',n_states,ev_ovlp_ev_complex,ubound(ev_ovlp_ev_complex,1),info)
                 if(info/=0) then
                    print *,'zpotrf info = ',info
                    stop
                 endif

                 ! Invert upper Cholesky factor, set lower half to 0
                 ! The matrix should never be singular, the test of info is for safety only

                 call ztrtri('U','N',n_states,ev_ovlp_ev_complex,ubound(ev_ovlp_ev_complex,1),info)
                 if(info/=0) then
                    print *,'ztrtri info = ',info
                    stop
                 endif
                 do i = 1, n_states-1
                    ev_ovlp_ev_complex(i+1:n_states,i) = 0
                 enddo

                 ! Orthogonalized eigenvectors = KS_eigenvector_complex * U**-1

                 call zgemm('N','N', n_basis,n_states,n_states, (1.d0,0.d0), &
                          KS_eigenvector_complex(1,1,i_spin,i_k),ubound(KS_eigenvector_complex,1), &
                          ev_ovlp_ev_complex,ubound(ev_ovlp_ev_complex,1), &
                          (0.d0,0.d0), ovlp_ev_complex,ubound(ovlp_ev_complex,1))
                 KS_eigenvector_complex(1:n_basis,1:n_states,i_spin,i_k) = ovlp_ev_complex(1:n_basis,1:n_states)

              end do

           end if

        end do

        call aims_deallocate( work_ovl,                     "work_ovl" )
        call aims_deallocate( ovlp_full_complex,   "ovlp_full_complex" )
        call aims_deallocate( ovlp_ev_complex,       "ovlp_ev_complex" )
        call aims_deallocate( ev_ovlp_ev_complex, "ev_ovlp_ev_complex" )

     end if ! if (real_eigenvectors)

  end if ! if (n_periodic == 0)

  ! Allocatable array that are tracked
  if (allocated(overlap_matrix_w))         call aims_deallocate( overlap_matrix_w,                 "overlap_matrix_w" )
  if (allocated(overlap_matrix_w_complex)) call aims_deallocate( overlap_matrix_w_complex, "overlap_matrix_w_complex" )
  if (allocated(work_ovl))                 call aims_deallocate( work_ovl,                                 "work_ovl" )
  if (allocated(ovlp_full))                call aims_deallocate( ovlp_full,                               "ovlp_full" )
  if (allocated(ovlp_ev))                  call aims_deallocate( ovlp_ev,                                   "ovlp_ev" )
  if (allocated(ev_ovlp_ev))               call aims_deallocate( ev_ovlp_ev,                             "ev_ovlp_ev" )
  if (allocated(ovlp_full_complex))        call aims_deallocate( ovlp_full_complex,               "ovlp_full_complex" )
  if (allocated(ovlp_ev_complex))          call aims_deallocate( ovlp_ev_complex,                   "ovlp_ev_complex" )
  if (allocated(ev_ovlp_ev_complex))       call aims_deallocate( ev_ovlp_ev_complex,             "ev_ovlp_ev_complex" )

end subroutine orthonormalize_eigenvectors
!******
