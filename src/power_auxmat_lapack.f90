!****s* FHI-aims/power_auxmat_lapack
!  NAME
!   power_auxmat_lapack
!  SYNOPSIS
!
      subroutine power_auxmat_lapack &
                  (auxmat, power, name)

!  PURPOSE
!   Calculate some power (e.g. inverse square root V^-(1/2))
!   of a square matrix within the auxiliary basis.
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
       use localorb_io, only: localorb_info
       implicit none

!  ARGUMENTS
      real*8, intent(INOUT) :: auxmat(n_basbas,n_loc_prodbas)
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

!     working array

      integer   n_nonsingular
      real*8 :: diff, thisdiff
      character*150 :: info_str

!  counter
      integer i_basis_1
      integer i_basis_2
      integer i_basis_3

      integer :: info
      character(*), parameter :: func = 'power_auxmat_lapack'

!  start to work

       if(myid.eq.0) then
         if(.not.allocated(temp_auxmat)) then
            allocate(temp_auxmat(n_basbas,n_basbas), stat=info)
            call check_allocation(info, 'temp_auxmat', func)
         endif
       else
         ! need dummy allocation! This array will never be touched unless myid = 0
            allocate(temp_auxmat(1,1), stat=info)
       endif
       temp_auxmat(:,:) = 0.d0

       call gather_auxmat(temp_auxmat, auxmat, n_basbas)

!  Averaging the coulomb matrix on processer 0.
      if(myid.eq.0) then
         diff = 0.d0
          do i_basis_1 =1, n_basbas
            do i_basis_2 = i_basis_1 , n_basbas
               thisdiff = abs(temp_auxmat(i_basis_1,i_basis_2) - &
                              temp_auxmat(i_basis_2,i_basis_1))
               diff = max(diff, thisdiff)
               temp_auxmat(i_basis_1,i_basis_2) = &
               0.5d0* ( temp_auxmat(i_basis_1,i_basis_2) + &
                        temp_auxmat(i_basis_2,i_basis_1) )

               temp_auxmat(i_basis_2,i_basis_1) = &
               temp_auxmat(i_basis_1,i_basis_2)

            enddo
          enddo
          if (name /= '') then
             write(info_str, &
             &     "(2X,'Difference of ',A,' matrix to its transposed:',ES12.4)") &
             & trim(name), diff
             call localorb_info(info_str)
          end if
       endif

!  The diagonalization is only performed on the first processer.
      if(myid.eq.0) then

         if (power < 0.d0) then
            threshold = prodbas_threshold
         else
            if (prodbas_threshold > 0.0) then
               threshold = prodbas_threshold
            else
               ! This is still a significant threshold as it removes negative EVs.
               threshold = 0.d0
            end if
         end if

         call power_genmat_lapack(n_basbas, temp_auxmat, power, &
         &                        safe_minimum, threshold, name)

        endif   ! myid.eq.0

!  Distribute the inverse coulomb matrix over different threads.
        call scatter_auxmat(temp_auxmat, auxmat, n_basbas)

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
       end subroutine power_auxmat_lapack

!******
