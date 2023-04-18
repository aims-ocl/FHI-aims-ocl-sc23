!****s* FHI-aims/evaluate_rpa_integrand
!  NAME
!   evaluate_rpa_integrand
!  SYNOPSIS

subroutine evaluate_rpa_integrand &
           (polar_freq,rpa_c_integrand)

!  PURPOSE
!  for a given polarisability and Coulomb matrix, evaluate the integrand 
!  of RPA correlation energy. This is given by
!
!  E_c_integrand =  -ln(det(1-chi0*v)) + tr(chi0*v) }

! USES
      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use scalapack_wrapper
      use scalapack_utils, only : sclpck_loc_ind
      use localorb_io
      use runtime_choices
      implicit none

! ARGUMENTS 
      real*8  :: polar_freq(n_basbas, n_loc_prodbas)
      real*8  :: rpa_c_integrand

! INPUTS
! o  polar_freq -- real array,
!            the RPA polarizability 
!
! OUTPUT
! o  rpa_c_integrand  -- real number,
!            the integrand for RPA correlation energy (-Tr(ln(1-v*chi_0)+Tr(v*chi_0))
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

      real*8  det_v_times_polar
      real*8  trace_v_times_polar

      real*8, dimension(:,:), allocatable :: v_times_polar

!     working array
      integer :: ipiv(n_basbas)
      integer :: info

!     timing

      character*50  filename

!     counters

! for MPI
      integer :: id_send
      integer :: id_recv
      integer :: tag
      integer :: mpierr
      integer :: my_status(MPI_STATUS_SIZE)

! for scalapack
      integer, dimension(:), allocatable :: ipiv_scal(:)

      integer :: i_state
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_basis_3
      integer :: i_index
      integer :: i_spin
      integer :: i_task
      integer :: lr, lc

!     begin work

!   invert the bare Coulomb matrix for calculating the screened Coulomb interaction .

       if(use_scalapack) then

         if(.not.allocated(v_times_polar)) then
             allocate(v_times_polar(MAX(1,max_row_2d),max_col_2d),stat=i_index)
             call check_allocation(i_index, 'v_times_polar                 ')
         endif

        endif

        rpa_c_integrand = 0.d0 

        if(use_scalapack) then
! scalapack version
! parameters are defined in module prodbas.f90  

! Copy polar_freq to v_times_polar changing the matrix distribution from 1D to 2D:

             call PDGEMR2D(n_basbas, n_basbas, &
                           polar_freq, 1, 1, aux_sc_desc, &
                           v_times_polar, 1, 1, aux_sc_desc_2d, &
                           my_blacs_ctxt_aux)

             trace_v_times_polar = 0.d0
             do i_basis_1 = 1, n_basbas
                lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
                lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
                if(lr>0 .and. lc>0) then
                   trace_v_times_polar = trace_v_times_polar + &
                         v_times_polar(lr, lc)
                   v_times_polar(lr,lc)  = &
                         v_times_polar(lr,lc) - 1.d0
!  Hubbard correction, need to be improved
!  E_c_integrand =  -ln(det(1-1/2*chi0*v)) + tr(chi0*v) }
                endif
             enddo
             v_times_polar(:,:) = - v_times_polar(:,:)

             call pdpotrf( 'U', n_basbas, v_times_polar, 1, 1, aux_sc_desc_2d, info)

             if (info.ne.0) then
               write(use_unit,*) " * Failure of Cholesky decomposition! "
               write(use_unit,*) " * Error info = ", info
               call aims_stop
             endif

             det_v_times_polar = 1.d0
             do i_basis_1 = 1, n_basbas, 1
                lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
                lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
                if(lr>0 .and. lc>0) then
                   det_v_times_polar = det_v_times_polar *  &
                         v_times_polar(lr, lc)**2
                endif
             enddo

             rpa_c_integrand = log (abs(det_v_times_polar)) + &
                                    trace_v_times_polar
! scalapack version
        else
          if(myid.eq.0) then
             if(.not.allocated(v_times_polar)) then
                allocate(v_times_polar(n_basbas,n_basbas),stat=i_index)
                call check_allocation(i_index, 'v_times_polar                 ')
             endif
             v_times_polar(:,:) = 0.d0
          endif

          id_recv = 0
          do i_task = 1, n_tasks

            id_send = i_task - 1
            tag= id_send
            if(id_recv .ne. id_send) then
              if( myid .eq. id_send ) then
                call MPI_SEND (polar_freq, &
                     n_basbas*n_loc_prodbas, &
                     MPI_DOUBLE_PRECISION, id_recv, &
                     tag, mpi_comm_global,mpierr)
              endif

              if(myid .eq. id_recv) then
                call MPI_RECV (polar_freq, &
                     n_basbas*n_loc_prodbas, &
                     MPI_DOUBLE_PRECISION, id_send, &
                     tag, mpi_comm_global,my_status,mpierr)
              endif
            endif

            if(myid.eq.0) then
              ! IGOR debug
              !write(filename,'(A,I1,A)') 'chi_',i_task,'.out'
              !open(UNIT=666,FILE=filename,STATUS='NEW',FORM='FORMATTED')
              !write(666,*) polar_freq(:,:)
              !close(666)
              do i_basis_1 = 1, n_loc_prodbas

               i_index=map_prodbas(i_basis_1,i_task)

               if(i_index.gt.0) then
                  v_times_polar(:,i_index) = &
                  polar_freq(:,i_basis_1)
               endif

              enddo
            endif
! end of loop over i_task
         enddo

        if(myid.eq.0) then 

           ! IGOR debug
           !open(UNIT=666,FILE='chi_d.out',STATUS='NEW',FORM='FORMATTED')
           !write(666,*) v_times_polar(:,:)
           !close(666)

           trace_v_times_polar = 0.d0
           do i_basis_1 = 1, n_basbas

             trace_v_times_polar = trace_v_times_polar + &
               v_times_polar(i_basis_1, i_basis_1)

!          write(use_unit,'(I6,2f18.8)') i_basis_1, v_times_polar(i_basis_1,i_basis_1), &
!            trace_v_times_polar

             v_times_polar(i_basis_1,i_basis_1) = &
               v_times_polar(i_basis_1,i_basis_1) - 1.d0

           enddo

           v_times_polar(:,:) = - v_times_polar(:,:)


!    Cholesky factorization of matrix 1-v*chi
!        call dpotf2('U', n_basbas, v_times_polar, n_basbas, info)

!    LU factorization of matrix 1-v*chi
           call dgetrf( n_basbas, n_basbas, v_times_polar, &
                        n_basbas, ipiv, info )

           if (info.ne.0) then
             write(use_unit,*) " * Failure of LU decomposition! "
             write(use_unit,*) " * Error info = ", info
             stop
           endif

           det_v_times_polar = 1.d0
           do i_basis_1=1, n_basbas
             det_v_times_polar = det_v_times_polar * &
               v_times_polar(i_basis_1, i_basis_1)
          
!            write(use_unit,*) i_basis_1, v_times_polar(i_basis_1,i_basis_1)
           enddo


           if(det_v_times_polar.lt.0) then
             write(use_unit,*) " * Detminant of V_TIMES_POLAR is negetive ! "

!          write(use_unit,*) " * STOP! "
!             stop
           endif
           rpa_c_integrand = log (abs(det_v_times_polar)) + &
                                 trace_v_times_polar
! end of if myid == 0
         endif

! end of if use scalapack
       endif

       call sync_real_number(rpa_c_integrand)     

      if(allocated(ipiv_scal)) then
         deallocate(ipiv_scal)
      endif
      if (allocated (v_times_polar)) then
        deallocate (v_times_polar)
      endif

      end subroutine evaluate_rpa_integrand
!---------------------------------------------------------------------
!******
!****s* FHI-aims/evaluate_rpa_integrand
!  NAME
!   evaluate_rpa_integrand
!  SYNOPSIS

subroutine evaluate_rpa_integrand_2 &
           (polar_freq,rpa_c_integrand)

!  PURPOSE
!  for a given polarisability and Coulomb matrix, evaluate the integrand 
!  of RPA correlation energy. This is given by
!
!  E_c_integrand =  -ln(det(1-chi0*v)) + tr(chi0*v) }

! USES
      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use scalapack_wrapper
      use scalapack_utils, only : sclpck_loc_ind
      use localorb_io
      implicit none

! ARGUMENTS 
      real*8  :: polar_freq(max_row_2d, max_col_2d)
      real*8  :: rpa_c_integrand

! INPUTS
! o  polar_freq -- real array,
!            the RPA polarizability 
!
! OUTPUT
! o  rpa_c_integrand  -- real number,
!            the integrand for RPA correlation energy (-Tr(ln(1-v*chi_0)+Tr(v*chi_0))
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

      real*8  det_v_times_polar
      real*8  trace_v_times_polar

      real*8, dimension(:,:), allocatable :: v_times_polar

!     working array
      integer :: ipiv(n_basbas)
      integer :: info

!     timing

      character*50  filename

!     counters

! for MPI
      integer :: id_send
      integer :: id_recv
      integer :: tag
      integer :: mpierr
      integer :: my_status(MPI_STATUS_SIZE)

! for scalapack
      integer, dimension(:), allocatable :: ipiv_scal(:)

      integer :: i_state
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_basis_3
      integer :: i_index
      integer :: i_spin
      integer :: i_task
      integer :: lr, lc

!     begin work

!   invert the bare Coulomb matrix for calculating the screened Coulomb interaction .

      allocate(v_times_polar(max_row_2d,max_col_2d),stat=i_index)
      call check_allocation(i_index, 'v_times_polar                 ')

      rpa_c_integrand = 0.d0 

! Please note: This 2D version only works with scalapack

! Copy polar_freq to v_times_polar changing the matrix distribution from 1D to 2D:

      v_times_polar(:,:) = polar_freq(:,:)

      trace_v_times_polar = 0.d0
      do i_basis_1 = 1, n_basbas
         lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
         lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
         if(lr>0 .and. lc>0) then
            trace_v_times_polar = trace_v_times_polar + &
                  v_times_polar(lr, lc)
            v_times_polar(lr,lc)  = &
                  v_times_polar(lr,lc) - 1.d0
!  Hubbard correction, need to be improved
!  E_c_integrand =  -ln(det(1-1/2*chi0*v)) + tr(chi0*v) }
         endif
      enddo
      v_times_polar(:,:) = - v_times_polar(:,:)

      call pdpotrf( 'U', n_basbas, v_times_polar, 1, 1, aux_sc_desc_2d, info)

      if (info.ne.0) then
        write(use_unit,*) " * Failure of Cholesky decomposition! "
        write(use_unit,*) " * Error info = ", info
        call aims_stop
      endif

      det_v_times_polar = 1.d0
      do i_basis_1 = 1, n_basbas, 1
         lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
         lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
         if(lr>0 .and. lc>0) then
            det_v_times_polar = det_v_times_polar *  &
                  v_times_polar(lr, lc)**2
         endif
      enddo

      rpa_c_integrand = log (abs(det_v_times_polar)) + &
                             trace_v_times_polar

      call sync_real_number(rpa_c_integrand)     

      if(allocated(ipiv_scal)) then
         deallocate(ipiv_scal)
      endif
      if (allocated (v_times_polar)) then
        deallocate (v_times_polar)
      endif

      end subroutine evaluate_rpa_integrand_2
!---------------------------------------------------------------------
!******
