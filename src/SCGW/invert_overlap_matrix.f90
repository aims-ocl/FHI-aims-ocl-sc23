!****s* FHI-aims/power_auxmat_lapack
!  NAME
!   power_auxmat_lapack
!  SYNOPSIS
!
      subroutine invert_overlap_matrix (auxmat)
!                  (auxmat, power, name)

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
       use physics
       implicit none

!  ARGUMENTS
     real*8 auxmat(n_basis, n_basis)
!      real*8, intent(INOUT) :: auxmat(n_basbas,n_loc_prodbas)
      real*8 :: power
!      character*(*) :: name
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
      integer i_index

      power = -1.0
 
      i_index = 0
      do i_basis_2 = 1, n_basis, 1
         do i_basis_1 = 1, i_basis_2, 1
            i_index = i_index+1
            auxmat(i_basis_1,i_basis_2) = overlap_matrix(i_index)
            auxmat(i_basis_2,i_basis_1) = overlap_matrix(i_index)
         enddo
      enddo
      
!  start to work
!  Averaging the coulomb matrix on processer 0.
!      if(myid.eq.0) then
!!$         if (power < 0.d0) then
            threshold =  1.d-5
!!$         else
         call power_genmat_lapack(n_basis, auxmat, power, &
         &                        safe_minimum, threshold, '')

!       endif   ! myid.eq.0
!  Distribute the inverse coulomb matrix over different threads.
!        call scatter_auxmat(temp_auxmat, auxmat, n_basis)

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
!       end subroutine power_auxmat_lapack
        end subroutine invert_overlap_matrix
!******
!****s* FHI-aims/power_auxmat_lapack
!  NAME
!   power_auxmat_lapack
!  SYNOPSIS
!
      subroutine invert_overlap_matrix_2 (auxmat)
!                  (auxmat, power, name)

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
       use physics
       implicit none

!  ARGUMENTS
     real*8 auxmat(n_basis, n_basis)
!      real*8, intent(INOUT) :: auxmat(n_basbas,n_loc_prodbas)
      real*8 :: power
!      character*(*) :: name
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
      integer i_index

      power = 0.5
 
      i_index = 0
      do i_basis_2 = 1, n_basis, 1
         do i_basis_1 = 1, i_basis_2, 1
            i_index = i_index+1
            auxmat(i_basis_1,i_basis_2) = overlap_matrix(i_index)
            auxmat(i_basis_2,i_basis_1) = overlap_matrix(i_index)
         enddo
      enddo
      
!  start to work
!  Averaging the coulomb matrix on processer 0.
!      if(myid.eq.0) then
!!$         if (power < 0.d0) then
            threshold =  1.d-5
!!$         else
         call power_genmat_lapack(n_basis, auxmat, power, &
         &                        safe_minimum, threshold, '')

!       endif   ! myid.eq.0
!  Distribute the inverse coulomb matrix over different threads.
!        call scatter_auxmat(temp_auxmat, auxmat, n_basis)

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
!       end subroutine power_auxmat_lapack
        end subroutine invert_overlap_matrix_2
!******
