!****s* FHI-aims/diagonalize_auxmat_scalapack_supercell
!  NAME
!   diagonalize_auxmat_scalapack_supercell
!  SYNOPSIS

      subroutine diagonalize_auxmat_scalapack_supercell &
      ( dist_auxmat, threshold, &
        n_nonsingular, sc_desc_am, &
        eigenvalues, auxmat_transform, &
        name &
      )

!  PURPOSE
!  to diagonalize the Coulomb matrix using ScaLapack
!
!  USES

      use dimensions, only: n_basbas_supercell
      use mpi_tasks
      use prodbas
      use scalapack_wrapper

      implicit none


! ARGUMENTS
      integer sc_desc_am(dlen_)
      character*(*), intent(IN) :: name
      real*8 dist_auxmat(sc_max_row,sc_max_col)
      real*8, intent(IN) :: threshold

      integer n_nonsingular
      real*8, dimension(n_basbas_supercell) :: eigenvalues
      real*8 auxmat_transform(sc_max_row,sc_max_col)

! INPUTS 
! o  sc_desc_am -- integer array, the descriptor for scalapack matrix manipulation
! o  dist_auxmat -- the distributed auxiliary matrix to be diagonalized using
!           scalapack
! o  threshold -- the cutoff threshold for the eigenvalues of the concerned
!         overlap or Coulomb matrix, only those above this threshold will be included.
! o  name -- name of matrix (for output only); Use empty string for no output.
! OUTPUT
! o  n_nonsingular -- real number,  number of nonsingular eigenvales larger than 
!          a certain (positive) cutoff threshold.
! o  eigenvalues -- Eigenvalues of the auxiliary matrix
! o  auxmat_transform -- Non-singular eigenvectors of the Coulomb matrix
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

!     n_found: total number of (overlap matrix) eigenvalues found
!     eigenvalues : The desired eigenvalues, including those to be discarded
!     info : Status value to indicate potential errors

      integer n_found
      integer nb_prod
      integer info

      integer :: lwork
      integer :: liwork2

      integer :: i_row, i_col
      integer :: trilwmin, lwormtr, lworsrt
      integer :: np0, nq0
      integer :: iarow, iacol

      real*8, dimension(:), allocatable :: work
      integer, dimension(:), allocatable :: iwork

      character*100 :: info_str

!  external function
      integer, external :: numroc, indxg2p

!     counters

      integer i_eval

      if (threshold < 0.d0) then
         write(info_str, "('** ',A)") &
         & 'diagonalize_auxmat_scalapack: threshold must be non-negative'
         call aims_stop(info_str)
      end if

!  begin work

!      if(myid.eq.0.and.flag_info_output.eq.0) then
!       write (use_unit,'(2X,2A)') &
!        "Diagonalizing the Coulomb matrix and ", &
!        "checking for singularities with ScaLapack ."
!       write (use_unit,*)
!      endif

      nb_prod = sc_nb_aux
!      np0 = NUMROC( MAX(n_basbas_supercell,nb_prod,2), nb_prod, 0, 0, nprow_aux )
!      nq0 = NUMROC( MAX(n_basbas_supercell,nb_prod,2), nb_prod, 0, 0, npcol_aux )

      iarow = INDXG2P( 0, nb_prod, sc_myprow_aux, 0, sc_nprow_aux )
      iacol = INDXG2P( 0, nb_prod, sc_mypcol_aux, 0, sc_npcol_aux )

      np0 = NUMROC( n_basbas_supercell, nb_prod, sc_myprow_aux, iarow, sc_nprow_aux )
      nq0 = NUMROC( n_basbas_supercell, nb_prod, sc_myprow_aux, iacol, sc_npcol_aux )

!      np0 = sc_max_row
!      nq0 = sc_max_col

      TRILWMIN = 3*n_basbas_supercell + MAX( nb_prod*( NP0+1 ), 3*nb_prod )
!      lwormtr = MAX(( nb_prod*(nb_prod-1))/2, &
!                    ( np0 + nq0)*nb_prod + 2*nb_prod*nb_prod)
      lwormtr = MAX(( nb_prod*(nb_prod-1))/2, &
                    ( np0 + nq0)*nb_prod) + nb_prod*nb_prod
      lworsrt = max(n_basbas_supercell, np0*(nb_prod+nq0))
      lwork = MAX( 1+6*n_basbas_supercell+2*NP0*NQ0, TRILWMIN, lwormtr,lworsrt ) + &
                  2*n_basbas_supercell
      allocate(work(lwork),stat=i_row)
      call check_allocation(i_row, 'work                          ')
	


      liwork2 = MAX( 7*n_basbas_supercell + 8*sc_npcol_aux + 2, n_basbas_supercell + &
                   2*nb_prod + 2*sc_npcol_aux)
!       liwork2= 7*n_basbas_supercell + 8*npcol_aux +2
      allocate(iwork(liwork2),stat=i_row)
      call check_allocation(i_row, 'iwork                         ')

! then, solve the problem
      dist_auxmat(:,:) = - dist_auxmat(:,:)
      call PDSYEVD('V', 'U', n_basbas_supercell, dist_auxmat, 1, 1, sc_desc_am, &
                   eigenvalues, auxmat_transform, 1, 1, sc_desc_am, work, &
                   lwork, iwork, liwork2, info)
      if(info /= 0) call scalapack_err_exit(info,"PDSYEVD")

      eigenvalues(1:n_basbas_supercell) = -eigenvalues(1:n_basbas_supercell)
      do i_row = 1, n_basbas_supercell, 1
         if(eigenvalues(i_row) < threshold) exit
!          write(use_unit,*) i_row, eigenvalues(i_row)
      enddo

      n_nonsingular = i_row -1

       if (allocated(work)) then
         deallocate(work)
       endif
       if (allocated(iwork)) then
         deallocate(iwork)
       endif


!  that's all folks

     end subroutine diagonalize_auxmat_scalapack_supercell
!-----------------------------------------------------------------------
!******
