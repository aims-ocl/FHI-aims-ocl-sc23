!****s* FHI-aims/evaluate_invs_times_sqrtv
!  NAME
!   evaluate_invs_times_sqrtv
!  SYNOPSIS
!
      subroutine evaluate_invs_times_sqrtv &
                  (ovlp_prodbas,coulomb_matr)

!  PURPOSE
!   calculate the inverse of overlap matrix S times square root of the Coulomb matrix 
!   V, namely S^(-1)*V^(1/2).
!   This is the lapack version, where the the matrix inversion is performed only
!   at one processor.
!
!  USES

       use dimensions
       use runtime_choices
       use prodbas
       use mpi_tasks
       use synchronize_mpi
       use localorb_io, only: localorb_info, use_unit, OL_norm
       implicit none

!  ARGUMENTS
      real*8  ovlp_prodbas(n_basbas,n_loc_prodbas)
      real*8  coulomb_matr(n_basbas,n_loc_prodbas)
!
!  INPUTS
!  o ovlp_prodbas  :  real array, overlap matrix of the auxiliary (product) basis
!  o coulomb_matr  :  real array, bare coulomb matrix within the auxiliary basis
!  OUTPUTS
!  o ovlp_prodbas  :  real array, the inverse of the overlap matrix
!  o coulomb_matr  :  real array, the inverse of the overlap matrix time the
!                     squre root of the bare coulomb matrix within the auxiliary basis
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

!  local variable

      real*8, dimension(:,:), allocatable ::  temp_coulomb_matr(:,:)
      real*8, dimension(:), allocatable ::    ovlp_eigenvalues
      real*8, dimension(:,:), allocatable ::  ovlp_transform
      real*8, dimension(:,:), allocatable ::  invs_times_sqrtv(:,:)
      real*8  ev_sqrt

!     working array

      integer   n_nonsingular

!  counter
      integer i_basis_1
      integer i_basis_2
      integer i_basis_3

      integer :: info
      character(*), parameter :: func = 'temp_coulomb_matr'

!  start to work

      call localorb_info('Multiplying S^-1 x V^0.5 (lapack)', &
      &                  use_unit, '(2X,A)', OL_norm)

       if(myid.eq.0) then
         if(.not.allocated(temp_coulomb_matr)) then
            allocate(temp_coulomb_matr(n_basbas,n_basbas), stat=info)
            call check_allocation(info, 'temp_coulomb_matr', func)
         endif
         temp_coulomb_matr(:,:) = 0.d0
       endif

       call gather_auxmat(temp_coulomb_matr, coulomb_matr, n_basbas)

!  Averaging the coulomb matrix on processer 0.
      if(myid.eq.0) then

          do i_basis_1 =1, n_basbas
            do i_basis_2 = i_basis_1 , n_basbas
               temp_coulomb_matr(i_basis_1,i_basis_2) = &
               0.5d0* ( temp_coulomb_matr(i_basis_1,i_basis_2) + &
                        temp_coulomb_matr(i_basis_2,i_basis_1) )

               temp_coulomb_matr(i_basis_2,i_basis_1) = &
               temp_coulomb_matr(i_basis_1,i_basis_2)

            enddo
          enddo
       endif

!  The diagonalization is only performed on the first processer.
      if(myid.eq.0) then

          allocate(ovlp_eigenvalues(n_basbas), stat=info)
          call check_allocation(info, 'ovlp_eigenvalues', func)
          allocate(ovlp_transform(n_basbas,n_basbas), stat=info)
          call check_allocation(info, 'ovlp_transform', func)
          ovlp_eigenvalues=0.d0
          ovlp_transform=0.d0

          call diagonalize_auxmat_lapack &
          (  n_basbas, temp_coulomb_matr, safe_minimum, &
             0.d0, &
             n_nonsingular, ovlp_eigenvalues, ovlp_transform, &
             "Coulomb" &
           )

!          write(use_unit,'(2X,A,E10.4)') &
!             "Lowest eigenvalue of the Coulomb matrix:", &
!               ovlp_eigenvalues(1)
!          write(use_unit,'(2X,A,I8,A,I8,A)') "Using ", n_nonsingular,  &
!                   "  eigenvalues out of rank ", &
!                    n_basbas, "  Coulomb matrix (auxiliary basis)."
!          write(use_unit,*)

!           if (ovlp_eigenvalues(1).lt. 1.d-5) then
!             if(myid.eq.0) then
!               write(use_unit,'(2X,A)') &
!              "Be careful! The Coulomb matrix may be ill-conditioned."
!               write(use_unit,*)
!             endif
!           endif

           do i_basis_1 = 1, n_basbas, 1
!             if(myid.eq.0) then
!              write(use_unit,*)i_basis_1, ovlp_eigenvalues(i_basis_1)
!              endif
             if(i_basis_1.le.n_nonsingular) then 
               ev_sqrt = sqrt(ovlp_eigenvalues(i_basis_1))
               do i_basis_2 = 1, n_basbas, 1
!              if(myid.eq.0) then
!                if (abs(ovlp_transform(i_basis_2, i_basis_1)) .gt. 0.1) then
!                 write(use_unit,'(5I5,f20.10)')i_basis_1, i_basis_2, basbas_atom(i_basis_2), &
!                       basbas_l(i_basis_2), basbas_m(i_basis_2), &
!                      ovlp_transform(i_basis_2, i_basis_1)
!                endif
!              endif

                   ovlp_transform(i_basis_2, i_basis_1) = &
                     ovlp_transform(i_basis_2,i_basis_1)*ev_sqrt
               enddo
             else
                   ovlp_transform(:, i_basis_1) = 0.d0
             endif
            enddo


!            temp_coulomb_matr = 0.d0
!            call dgemm('N', 'T', n_basbas, n_basbas, &
!                       n_nonsingular, 1.0d0, &
!                       ovlp_transform(:,1:n_nonsingular), &
!                       n_basbas, &
!                       ovlp_transform(:,1:n_nonsingular), &
!                       n_basbas, 0.d0, &
!                       temp_coulomb_matr, n_basbas &
!                      )

!  end of if myid.eq.0
        endif
!            if(myid.eq.0) then
!            write(use_unit,*) "inverted coulomb"
!            do i_basis_1 = 1, n_basbas, 1
!             do i_basis_2= 1,  i_basis_1, 1
!               write(use_unit,'(A,2I6,f20.10)')"inv", i_basis_1, i_basis_2,
!     +         temp_coulomb_matr(i_basis_1,i_basis_2)
!             enddo
!             enddo
!            endif

        call scatter_auxmat(ovlp_transform, coulomb_matr, n_basbas)

        call power_auxmat_lapack(ovlp_prodbas, -1.d0, "overlap")

       allocate(invs_times_sqrtv(n_basbas,n_loc_prodbas))
!       call check_allocation(i_index, 'invs_times_sqrtv              ')


       call auxiliary_matrix_multi &
               ( n_basbas,n_loc_prodbas,map_prodbas, &
                 ovlp_prodbas, coulomb_matr, invs_times_sqrtv)

        coulomb_matr (:,:) = invs_times_sqrtv(:,:)
       
!       if(myid.eq.0) then
!         do i_basis_1 = 1, n_basbas, 1
!         do i_basis_2 = 1, n_loc_prodbas, 1
!          write(use_unit,'(2I6,f20.10)') i_basis_1,i_basis_2, &
!                         coulomb_matr(i_basis_1,i_basis_2)
!         enddo
!         enddo
!       endif

!  deallocate arrays

       if(allocated(ovlp_transform)) then
         deallocate(ovlp_transform)
       endif
       if(allocated(ovlp_eigenvalues)) then
         deallocate(ovlp_eigenvalues)
       endif
       if(allocated(temp_coulomb_matr)) then
         deallocate(temp_coulomb_matr)
       endif
       if(allocated(invs_times_sqrtv)) then
         deallocate(invs_times_sqrtv)
       endif
       end subroutine evaluate_invs_times_sqrtv

!******
