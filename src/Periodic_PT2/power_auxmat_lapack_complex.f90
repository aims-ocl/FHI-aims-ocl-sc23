!****s* FHI-aims/power_auxmat_lapack_complex
!  NAME
!   power_auxmat_lapack_complex
!  SYNOPSIS
!
      subroutine power_auxmat_lapack_complex &
                  (auxmat, power, name)

!  PURPOSE
!   Calculate some power (e.g. inverse square root V^-(1/2))
!   of a complex square matrix within the auxiliary basis.
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
       implicit none

!  ARGUMENTS
      complex*16, intent(INOUT) :: auxmat(n_basbas,n_basbas)
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

      !complex*16, dimension(:,:), allocatable ::  temp_auxmat
      real*8, dimension(:), allocatable ::    eigenvalues
      complex*16, dimension(:,:), allocatable ::  transform
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
      character(*), parameter :: func = 'power_auxmat_lapack_complex'

!  start to work

! To purge the coulomb matrix "auxmat" to be the exact Hermite matrix
      !diff = 0.d0
      !do i_basis_1 =1, n_basbas
      !  do i_basis_2 = i_basis_1 , n_basbas
      !     thisdiff = abs(auxmat(i_basis_1,i_basis_2) - &
      !                    conjg(auxmat(i_basis_2,i_basis_1)))
      !     diff = max(diff, thisdiff)
      !     auxmat(i_basis_1,i_basis_2) = &
      !     0.5d0* ( auxmat(i_basis_1,i_basis_2) + &
      !              conjg(auxmat(i_basis_2,i_basis_1) ))

      !     auxmat(i_basis_2,i_basis_1) = &
      !     conjg(auxmat(i_basis_1,i_basis_2))

      !  enddo
      !enddo
      !if (name /= '') then
      !   write(info_str, &
      !   &     "(2X,'Difference of ',A,' matrix to its transposed:',ES12.4)") &
      !   & trim(name), diff
      !   call localorb_info(info_str)
      !end if

      if (power < 0.d0) then
         threshold = prodbas_threshold
      else
         !if (prodbas_threshold > 0.0) then
         !   threshold = prodbas_threshold
         !else
            ! This is still a significant threshold as it removes negative EVs.
            threshold = 0.d0
         !end if
      end if

      call power_genmat_lapack_complex(n_basbas, auxmat, &
         & power, safe_minimum, threshold, name)

!  deallocate arrays

       if(allocated(eigenvalues)) then
         deallocate(eigenvalues)
       endif
       if(allocated(transform)) then
         deallocate(transform)
       endif
       end subroutine power_auxmat_lapack_complex

!******
