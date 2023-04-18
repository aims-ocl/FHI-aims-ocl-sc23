!****s* FHI-aims/power_auxmat_scalapack_supercell
!  NAME
!   power_auxmat_scalapack_supercell
!  SYNOPSIS

      subroutine power_auxmat_scalapack_supercell &
                  (auxmat, power, name)

!  PURPOSE
!   Calculate some power (e.g. inverse square root V^-(1/2))
!   of a square matrix within the auxiliary basis.
!   This is done by diagonalization.  Eigenmodes corresponding to
!   eigenvalues smaller then prodbas_threshold are set to zero.
!   This is the scalapack version
!
!  USES
       use dimensions
       use runtime_choices
       use prodbas
       use mpi_tasks
       use synchronize_mpi
       use scalapack_wrapper
       use localorb_io
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

      real*8, dimension(:,:), allocatable ::  dist_auxmat
      real*8, dimension(:), allocatable ::    eigenvalues
      real*8, dimension(:,:), allocatable ::  transform
      real*8  :: ev_sqrt, threshold

      integer, dimension(:), allocatable :: l_row_prod
      integer, dimension(:), allocatable :: l_col_prod

!     working array

      integer   n_nonsingular

      integer   id_send
      integer   id_recv
      integer   tag
      integer   my_status(MPI_STATUS_SIZE)
      integer   sc_desc_am(dlen_)
      integer   info
      integer   mpierr
      integer   nb_prod

!  external function
      integer, external ::   numroc

!  counter
      integer i_basis_1
      integer i_basis_2
      integer i_basis_3
      integer i_index
      integer i_task, i_task_1
      integer lr, lc

!  start to work

      if(myid.eq.0 .and. name /= '') then
       write (use_unit,'(2X,4A)') &
        "Diagonalizing the ", trim(name), " matrix and ", &
        "checking for singularities with ScaLapack ."
       write (use_unit,*)
      endif


      nb_prod = 64

!      call descinit( sc_desc_am, n_basbas_supercell, n_basbas_supercell, mb_prod, nb_prod, &
!                      rsrc, csrc, &
!                      my_blacs_ctxt, MAX(1,max_row), info )

      call descinit( sc_desc_am, n_basbas_supercell, n_basbas_supercell, sc_mb_aux, sc_nb_aux, &
                      0, 0, &
                      sc_my_blacs_ctxt_aux, MAX(1,sc_max_row), info )

      if(.not.allocated(l_row_prod)) then
        allocate(l_row_prod(n_basbas_supercell))
      endif
      if(.not.allocated(l_col_prod)) then
        allocate(l_col_prod(n_basbas_supercell))
      endif

      l_row_prod(:) = 0
      l_col_prod(:) = 0

! ATTENTION: The following code assumes rsrc==0 and csrc==0 !!!!
! copied from "scalapack_wrapper.f90"

      lr = 0 ! local row counter
      lc = 0 ! local column counter

      do i_basis_1 = 1, n_basbas_supercell

        if( MOD((i_basis_1-1)/sc_mb_aux,sc_nprow_aux) .eq. sc_myprow_aux) then
! row i is on local processor
           lr = lr+1
           l_row_prod(i_basis_1) = lr
        endif

        if( MOD((i_basis_1-1)/sc_nb_aux,sc_npcol_aux) .eq. sc_mypcol_aux) then
! column i is on local processor
           lc = lc+1
           l_col_prod(i_basis_1) = lc
        endif

      enddo

!  redistribute the Coulomb interaction matrix

      if(.not.allocated(dist_auxmat)) then
         allocate(dist_auxmat(sc_max_row,sc_max_col))
         dist_auxmat(:,:) = auxmat(:,1:sc_max_col)
      endif

      if (.not.allocated(transform)) then
        allocate(transform(sc_max_row,sc_max_col))
        transform (:,:) = 0.d0
      endif
! averaging the upper and lower triganles
      call pdtran( n_basbas_supercell, n_basbas_supercell, 1.0d0, dist_auxmat, &
                   1, 1, sc_desc_am, 0.0d0, transform, &
                   1, 1, sc_desc_am )
       dist_auxmat(:,:) = 0.5d0* (dist_auxmat(:,:) + &
                                        transform(:,:) )
! diagonalizing the matrix

       if (power < 0.d0) then
          threshold = prodbas_threshold
       else
          threshold = 0.d0   ! No threshold needed.
       end if

       allocate(eigenvalues(n_basbas_supercell))
       eigenvalues(:) = 0.d0
       transform(:,:) = 0.d0
       call diagonalize_auxmat_scalapack_supercell &
          (  dist_auxmat, threshold, &
             n_nonsingular, sc_desc_am, &
             eigenvalues, transform, &
             name &
           )

!       write(use_unit,*) "n_nonsingular", n_nonsingular
!       write(use_unit,*) eigenvalues(:)

      if(myid.eq.0 .and. name /= '') then
        write(use_unit,'(2X,3A,E10.4)') &
          "Lowest eigenvalue of the ", trim(name), " matrix:", &
             eigenvalues(n_nonsingular)
        write(use_unit,'(2X,A,I8,A,I8,3A)') "Using ", n_nonsingular,  &
               "   eigenvalues out of rank ", &
                n_basbas_supercell, "   ", trim(name), " matrix (auxiliary basis)."
        write(use_unit,*)
      endif

      if (myid.eq.0 .and. &
         eigenvalues(n_nonsingular).lt. 1.d-5 .and. power < 0.d0) then
          write(use_unit,'(2X,3A)') &
            "Be careful! The ", trim(name), " matrix may be ill-conditioned."
             write(use_unit,*)
      endif

!      do i_basis_1 = 1, n_basbas_supercell, 1
!         if(myid.eq.0) then
!           write(use_unit,*)i_basis_1, eigenvalues(i_basis_1)
!         endif
!      enddo

      do i_basis_1 = 1, n_nonsingular, 1
!         if(myid.eq.0) then
!           write(use_unit,*)i_basis_1, eigenvalues(i_basis_1)
!         endif
         lc = l_col_prod(i_basis_1)
         if(lc .gt. 0) then
            ev_sqrt = sqrt(eigenvalues(i_basis_1))
            transform(:, lc) = &
               transform(:, lc) * ev_sqrt**power
         endif
       enddo


       dist_auxmat = 0.d0
       call pdgemm('N', 'T', n_basbas_supercell, n_basbas_supercell, &
                   n_nonsingular, 1.0d0, &
                   transform(1,1), &
                   1, 1, sc_desc_am, &
                   transform(1,1), &
                   1, 1, sc_desc_am, 0.d0, &
                   dist_auxmat, &
                   1, 1, sc_desc_am &
                  )

       auxmat(:,1:sc_max_col)=dist_auxmat(:,:)

!  deallocate arrays

       if(allocated(l_row_prod)) then
         deallocate(l_row_prod)
       endif
       if(allocated(l_col_prod)) then
         deallocate(l_col_prod)
       endif
       if(allocated(eigenvalues)) then
         deallocate(eigenvalues)
       endif
       if(allocated(transform)) then
         deallocate(transform)
       endif
       if(allocated(dist_auxmat)) then
         deallocate(dist_auxmat)
       endif
     end subroutine power_auxmat_scalapack_supercell

!******
