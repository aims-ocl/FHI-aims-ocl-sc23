!****s* FHI-aims/power_auxmat_lapack_supercell
!  NAME
!   power_auxmat_lapack_supercell
!  SYNOPSIS
!
      subroutine power_auxmat_lapack_supercell &
                  (auxmat, power, name)

!  PURPOSE
!   Calculate some power (e.g. inverse square root V^-(1/2))
!   of a square matrix within the global auxiliary basis.
!   This is done by diagonalization.  Eigenmodes corresponding to
!   eigenvalues smaller then prodbas_threshold are set to zero.
!   This is the lapack version
!
!  USES

       use dimensions
       use runtime_choices
       use prodbas
       use mpi_tasks
       use synchronize_mpi
       use localorb_io, only: use_unit
       implicit none

!  ARGUMENTS
      real*8, intent(INOUT) :: auxmat(n_basbas_supercell,n_loc_prodbas_supercell)
      real*8, intent(IN)    :: power
      character*(*), intent(IN) :: name
!
!  INPUTS
!    o auxmat -- real array, e.g. bare coulomb matrix within the auxiliary basis
!    o power -- real number, e.g. -0.5d0
!    o name -- name of matrix (for output only)
!  OUTPUTS
!    o auxmat -- real array, e.g. the inverse square root of the bare coulomb matrix
!                within the auxiliary basis
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

      real*8, dimension(:,:), allocatable ::  temp_auxmat
      real*8, dimension(:), allocatable ::    eigenvalues
      real*8, dimension(:,:), allocatable ::  transform
      real*8  :: ev_sqrt, threshold
      integer, dimension(:), allocatable :: stride
      integer :: mpierr

!     working array

      integer   n_nonsingular

!  counter
      integer i_basis_1
      integer i_basis_2
      integer i_basis_3
      integer i_task

!  start to work

       if(myid.eq.0) then
         if(.not.allocated(temp_auxmat)) then
           allocate(temp_auxmat(n_basbas_supercell,n_basbas_supercell))
         endif
         temp_auxmat(:,:) = 0.d0
       endif

!       call gather_auxmat_supercell(temp_auxmat, auxmat, n_basbas_supercell)

       if(use_mpi)then
          if(myid.eq.0)then
             allocate(stride(n_tasks))
             stride(1)=0
             do i_task = 2, n_tasks
                stride(i_task) = stride(i_task-1)+n_prodbas_per_proc(i_task-1)*n_basbas_supercell
             enddo
          endif
          call mpi_gatherv(auxmat,n_loc_prodbas_supercell*n_basbas_supercell,&
               MPI_DOUBLE_PRECISION,temp_auxmat,n_prodbas_per_proc*n_basbas_supercell,&
               stride,MPI_DOUBLE_PRECISION,0,mpi_comm_global,mpierr)
       else
          temp_auxmat = auxmat
       endif


!  Averaging the coulomb matrix on processer 0.
      if(myid.eq.0) then

          do i_basis_1 =1, n_basbas_supercell
            do i_basis_2 = i_basis_1 , n_basbas_supercell
               temp_auxmat(i_basis_1,i_basis_2) = &
               0.5d0* ( temp_auxmat(i_basis_1,i_basis_2) + &
                        temp_auxmat(i_basis_2,i_basis_1) )

               temp_auxmat(i_basis_2,i_basis_1) = &
               temp_auxmat(i_basis_1,i_basis_2)

            enddo
          enddo
       endif

!  The diagonalization is only performed on the first processer.
      if(myid.eq.0) then

          allocate(eigenvalues(n_basbas_supercell))
          allocate(transform(n_basbas_supercell,n_basbas_supercell))
          eigenvalues=0.d0
          transform=0.d0

          if (power < 0.d0) then
             threshold = prodbas_threshold
          else
             threshold = 0.d0   ! No threshold needed.
          end if

          call diagonalize_auxmat_lapack &
          (  n_basbas_supercell, temp_auxmat, safe_minimum, &
             threshold, &
             n_nonsingular, eigenvalues, transform, &
             name &
           )

          if (name /= '') then
             write(use_unit,'(2X,3A,E10.4)') &
             "Lowest eigenvalue of the ", trim(name), " matrix:", &
             eigenvalues(1)
             write(use_unit,'(2X,A,I8,A,I8,3A)') "Using ", n_nonsingular,  &
             "  eigenvalues out of rank ", &
             n_basbas_supercell, "  ", trim(name), " matrix (auxiliary basis)."
             write(use_unit,*)
          end if

           if (eigenvalues(1).lt. 1.d-5 .and. power < 0.d0) then
             if(myid.eq.0) then
               write(use_unit,'(2X,3A)') &
              "Be careful! The ", trim(name), " matrix may be ill-conditioned."
               write(use_unit,*)
             endif
           endif

           do i_basis_1 = 1, n_nonsingular, 1
!             if(myid.eq.0) then
!              write(use_unit,*)i_basis_1, eigenvalues(i_basis_1)
!              endif
              ev_sqrt = sqrt(eigenvalues(i_basis_1))
              do i_basis_2 = 1, n_basbas_supercell, 1
                 transform(i_basis_2, i_basis_1) = &
                   transform(i_basis_2,i_basis_1) * ev_sqrt**power
              enddo
            enddo


            ! BL: not necessary, is done by dgemm
            !temp_auxmat = 0.d0
            call dgemm('N', 'T', n_basbas_supercell, n_basbas_supercell, &
                    n_nonsingular, 1.0d0, &
                    transform(:,1:n_nonsingular), &
                    n_basbas_supercell, &
                    transform(:,1:n_nonsingular), &
                    n_basbas_supercell, 0.d0, &
                    temp_auxmat, n_basbas_supercell &
                   )

        endif   ! myid.eq.0

!  Distribute the inverse coulomb matrix over different threads.
!        call scatter_auxmat_supercell(temp_auxmat, auxmat, n_basbas_supercell)

        if(use_mpi)then
           call mpi_scatterv(temp_auxmat,n_prodbas_per_proc*n_basbas_supercell,stride,&
                MPI_DOUBLE_PRECISION,auxmat,n_loc_prodbas_supercell*n_basbas_supercell,&
                MPI_DOUBLE_PRECISION,0,mpi_comm_global,mpierr)
        else
           auxmat = temp_auxmat
        endif


!  deallocate arrays

       if(allocated(temp_auxmat)) then
         deallocate(temp_auxmat)
       endif
       if(allocated(eigenvalues)) then
         deallocate(eigenvalues)
       endif
       if(allocated(transform)) then
         deallocate(transform)
       endif
       if(allocated(stride)) then
         deallocate(stride)
       endif

     end subroutine power_auxmat_lapack_supercell

!******
