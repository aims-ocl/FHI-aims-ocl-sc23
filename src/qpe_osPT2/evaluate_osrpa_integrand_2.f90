!****s* FHI-aims/evaluate_osrpa_integrand
!  NAME
!   evaluate_osrpa_integrand
!  SYNOPSIS

subroutine evaluate_osrpa_integrand_2 &
           (polar_freq,rpa_c_integrand_os, &
            rpa_c_integrand_ss, &
            rpa_c_integrand_total,c_osrpa_integrand)

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
      use runtime_choices, only: prodbas_threshold, safe_minimum, sc_check, n_lambda_osrpa, c_osrpa_order, c_osrpa_threshold
      use localorb_io
      implicit none

! ARGUMENTS
      real*8  :: polar_freq(max_row_2d, max_col_2d, n_spin)
      real*8  :: rpa_c_integrand_os, rpa_c_integrand_ss, rpa_c_integrand_total
      real*8  :: c_osrpa_integrand(c_osrpa_order)

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
      real*8  delta_trace(2)
      real*8  i_term(2)

      real*8, dimension(:,:), allocatable :: v_times_polar
      real*8, dimension(:), allocatable :: eigenvalues
      real*8, dimension(:), allocatable :: lambda_grid
      real*8, dimension(:), allocatable :: lambda_weight
      real*8 special_radius
      integer n_lambda
      integer i_lambda, i_order

      character :: spin_name(2)
      real*8 :: rpa_c_integrand_spin
      real*8 :: c_osrpa_integrand_spin(c_osrpa_order)


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
      integer :: i_spin, j_spin
      integer :: i_task
      integer :: lr, lc

!     begin work

!   invert the bare Coulomb matrix for calculating the screened Coulomb interaction .

     if(.not.allocated(v_times_polar)) then
        allocate(v_times_polar(max_row_2d,max_col_2d),stat=i_index)
        call check_allocation(i_index, 'v_times_polar                 ')
     endif
     v_times_polar(:,:) = 0.d0

     rpa_c_integrand_os    = 0.0d0
     rpa_c_integrand_ss    = 0.0d0
     rpa_c_integrand_total = 0.0d0

     if (n_spin .eq. 1) then
       ! IGOR debug
       !open(UNIT=666,FILE='chi_c.out',STATUS='NEW',FORM='FORMATTED')
       !write(666,*) v_times_polar(:,:,1)
       !close(666)

       ! check if the system is stronly correlated or not.
       !temp_v_times_polar_a = v_times_polar(:,:,1)
       !if (sc_check(1)) then
       !  call diagonalize_auxmat_lapack(n_basbas, temp_v_times_polar_a, &
       !  & safe_minimum, -1.0d10, n_nonsingular, eigenvalues, &
       !  polar_transform, " ")
       !  special_radius=maxval(abs(eigenvalues))
       !  temp_v_times_polar_a = v_times_polar(:,:,1)
       !  write(use_unit,'(2X,A,f16.8,A)') &
       !      "Special radius of non-interacting response matrix in each spin channel = ", &
       !      special_radius,"."
       !  sc_check = .false.
       !endif

       !=================================
       ! for the dRPA calculations
       !=================================
       trace_v_times_polar = 0.d0
       v_times_polar(:,:) = polar_freq(:,:,1)*2.0d0
       do i_basis_1 = 1, n_basbas
         lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
         lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
         if (lr>0 .and. lc>0) then
             trace_v_times_polar = trace_v_times_polar + &
               v_times_polar(lr, lc)

             v_times_polar(lr,lc) = &
               v_times_polar(lr,lc) - 1.d0
         endif
       enddo

       v_times_polar(:,:) = - v_times_polar(:,:)

       !    LU factorization of matrix 1-v*chi
       !call dgetrf( n_basbas, n_basbas, temp_v_times_polar_a(:,:), &
       !             n_basbas, ipiv, info )
       call pdpotrf( 'U', n_basbas, v_times_polar, 1, 1, aux_sc_desc_2d, info)

       if (info.ne.0) then
         write(use_unit,*) " * Failure of LU decomposition! "
         write(use_unit,*) " * Error info = ", info
         stop
       endif

       det_v_times_polar = 1.d0
       do i_basis_1=1, n_basbas, 1
         lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
         lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
         if(lr>0 .and. lc>0) then
           det_v_times_polar = det_v_times_polar * &
             v_times_polar(lr, lc)**2.0d0
         endif
       enddo
       !write(use_unit,*) "Igor debug, det chi = ", det_v_times_polar

       if(det_v_times_polar.lt.0) then
         write(use_unit,*) " * Detminant of V_TIMES_POLAR is negetive ! "

!      write(use_unit,*) " * STOP! "
!         stop
       endif

       rpa_c_integrand_total = log (abs(det_v_times_polar)) + &
                             trace_v_times_polar

       call sync_real_number(rpa_c_integrand_total)
       !if (myid.eq.0) write(use_unit,*) "Igor debug ", rpa_c_integrand_total, trace_v_times_polar

       !======================================
       ! Now for os-RPA and ss-RPA
       !======================================

       trace_v_times_polar = 0.d0
       v_times_polar(:,:) = polar_freq(:,:,1)
       do i_basis_1 = 1, n_basbas
         lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
         lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
         if (lr>0 .and. lc>0) then
             trace_v_times_polar = trace_v_times_polar + &
               v_times_polar(lr, lc)*2.0d0

             v_times_polar(lr,lc) = &
               v_times_polar(lr,lc) - 1.d0
         endif
       end do
       v_times_polar(:,:) = - v_times_polar(:,:)
       !    Cholesky factorization of matrix 1-v*chi
       !        call dpotf2('U', n_basbas, v_times_polar, n_basbas, info)

       !    LU factorization of matrix 1-v*chi
       !call dgetrf( n_basbas, n_basbas, v_times_polar(:,:,1), &
       !             n_basbas, ipiv, info )
       call pdpotrf( 'U', n_basbas, v_times_polar, 1, 1, aux_sc_desc_2d, info)

       if (info.ne.0) then
         write(use_unit,*) " * Failure of LU decomposition! "
         write(use_unit,*) " * Error info = ", info
         stop
       endif

       det_v_times_polar = 1.d0
       do i_basis_1=1, n_basbas
         lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
         lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
         if(lr>0 .and. lc>0) then
           det_v_times_polar = det_v_times_polar * &
             v_times_polar(lr,lc)**2.0d0
         endif
       enddo
       !write(use_unit,*) "Igor debug, det chi = ", det_v_times_polar

       if(det_v_times_polar.lt.0) then
         write(use_unit,*) " * Detminant of V_TIMES_POLAR is negetive ! "

!      write(use_unit,*) " * STOP! "
!         stop
       endif
       rpa_c_integrand_os = 2.0*log (abs(det_v_times_polar)) + &
                             trace_v_times_polar

       call sync_real_number(rpa_c_integrand_os)
       ! for close-shell cases, ss-rpa always equals to os-rpa at the first order
       rpa_c_integrand_ss = rpa_c_integrand_os

     else ! n_spin !=1

        n_lambda = n_lambda_osrpa
        if(.not.allocated(lambda_grid)) then
           allocate(lambda_grid(n_lambda),stat=i_index)
           call check_allocation(i_index, 'n_lambda_grid               ')
        endif
        if(.not.allocated(lambda_weight)) then
           allocate(lambda_weight(n_lambda),stat=i_index)
           call check_allocation(i_index, 'n_lambda_weight             ')
        endif
        call gauleg(0d0,1.d0,lambda_grid,lambda_weight,n_lambda)

        !if(.not.allocated(temp_v_times_polar_a)) then
        !   allocate(temp_v_times_polar_a(max_row_2d,max_col_2d),stat=i_index)
        !   call check_allocation(i_index, 'temp_v_times_polar          ')
        !endif
        !temp_v_times_polar_a(:,:) = 0.0d0
        !if(.not.allocated(temp_v_times_polar_b)) then
        !   allocate(temp_v_times_polar_b(max_row_2d,max_col_2d),stat=i_index)
        !   call check_allocation(i_index, 'temp_v_times_polar          ')
        !endif
        !temp_v_times_polar_b(:,:) = 0.0d0
        if(.not.allocated(eigenvalues)) then
           allocate(eigenvalues(n_basbas),stat=i_index)
           call check_allocation(i_index, 'eigenvalues                 ')
        endif
        eigenvalues(:) = 0.0d0

       !=================================
       ! for the dRPA calculations
       !=================================
       v_times_polar(:,:) = polar_freq(:,:,1) + polar_freq(:,:,2)

       trace_v_times_polar = 0.d0
       do i_basis_1 = 1, n_basbas
         lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
         lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
         if (lr>0 .and. lc>0) then
             trace_v_times_polar = trace_v_times_polar + &
               v_times_polar(lr, lc)

             v_times_polar(lr,lc) = &
               v_times_polar(lr,lc) - 1.d0
         endif
       enddo

       v_times_polar(:,:) = - v_times_polar(:,:)
       !    Cholesky factorization of matrix 1-v*chi
       !        call dpotf2('U', n_basbas, v_times_polar, n_basbas, info)

       !    LU factorization of matrix 1-v*chi
       call pdpotrf( 'U', n_basbas, v_times_polar, 1, 1, aux_sc_desc_2d, info)

       if (info.ne.0) then
         write(use_unit,*) " * Failure of LU decomposition! "
         write(use_unit,*) " * Error info = ", info
         stop
       endif

       det_v_times_polar = 1.d0
       do i_basis_1=1, n_basbas
         lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
         lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
         if(lr>0 .and. lc>0) then
           det_v_times_polar = det_v_times_polar * &
             v_times_polar(lr,lc)**2.0d0
         endif
       enddo
       !write(use_unit,*) "Igor debug, det chi = ", det_v_times_polar

       if(det_v_times_polar.lt.0) then
         write(use_unit,*) " * Detminant of V_TIMES_POLAR is negetive ! "

!      write(use_unit,*) " * STOP! "
!         stop
       endif
       rpa_c_integrand_total = log (abs(det_v_times_polar)) + &
                             trace_v_times_polar

       call sync_real_number(rpa_c_integrand_total)

       !=================================
       ! for the ss-RPA calculations
       !=================================
       rpa_c_integrand_ss = trace_v_times_polar
       do i_spin = 1, n_spin
         v_times_polar = polar_freq(:,:,i_spin)

         do i_basis_1 = 1, n_basbas
            lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
            lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
            if(lr>0 .and. lc>0) then
              v_times_polar(lr, lc) = &
                v_times_polar(lr, lc) - 1.d0
            endif
         enddo

         v_times_polar(:,:) = - v_times_polar(:,:)
         !    Cholesky factorization of matrix 1-v*chi
         !        call dpotf2('U', n_basbas, v_times_polar, n_basbas, info)

         !    LU factorization of matrix 1-v*chi
         !call dgetrf( n_basbas, n_basbas, temp_v_times_polar_a(:,:), &
         !             n_basbas, ipiv, info )
         call pdpotrf( 'U', n_basbas, v_times_polar, 1, 1, aux_sc_desc_2d, info)

         if (info.ne.0) then
           write(use_unit,*) " * Failure of LU decomposition! "
           write(use_unit,*) " * Error info = ", info
           stop
         endif

         det_v_times_polar = 1.d0
         do i_basis_1=1, n_basbas
            lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
            lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
            if(lr>0 .and. lc>0) then
              det_v_times_polar = det_v_times_polar * &
                v_times_polar(lr, lc)**2.0d0
            endif
         enddo
         !write(use_unit,*) "Igor debug, det chi = ", det_v_times_polar

         if(det_v_times_polar.lt.0) then
           write(use_unit,*) " * Detminant of V_TIMES_POLAR is negetive ! "

!        write(use_unit,*) " * STOP! "
!           stop
         endif
         rpa_c_integrand_ss = rpa_c_integrand_ss + log (abs(det_v_times_polar))
       end do

       call sync_real_number(rpa_c_integrand_ss)

       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! IGOR debug
       !open(UNIT=666,FILE='chi_a.out',STATUS='NEW',FORM='FORMATTED')
       !write(666,*) v_times_polar(:,:,1)
       !close(666)
       !open(UNIT=888,FILE='chi_b.out',STATUS='NEW',FORM='FORMATTED')
       !write(888,*) v_times_polar(:,:,2)
       !close(888)
       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       !=================================
       ! for the os-RPA calculations
       !=================================
       delta_trace = 0.0d0
       i_term = 1.0d0

       do i_spin = 1, n_spin
         rpa_c_integrand_spin = 0.0d0
         !c_osrpa_integrand_spin = 0.0d0
         if (i_spin.eq.1) then
             j_spin=2
         else
             j_spin=1
         end if
         call evaluate_osrpa_response_2(i_spin,j_spin, spin_name,&
           polar_freq,eigenvalues, &
           n_lambda, lambda_grid,lambda_weight,&
           rpa_c_integrand_spin,c_osrpa_integrand_spin)

         if (myid.eq.0) write(use_unit,'(2X,A,X,A,A,f19.8,A)') &
             "osRPA correlation (",spin_name(1),"-spin) :", &
             rpa_c_integrand_spin," Ha."
         rpa_c_integrand_os = rpa_c_integrand_os + rpa_c_integrand_spin
         !do i_order = 1, c_osrpa_order
         !  c_osrpa_integrand(i_order) =  c_osrpa_integrand(i_order) + &
         !      c_osrpa_integrand_spin(i_order)
         !end do
       enddo ! i_spin
      end if ! if (n_spin == 2)

      !call sync_real_number(rpa_c_integrand_os)
      !call sync_real_number(rpa_c_integrand_ss)
      !call sync_real_number(rpa_c_integrand_total)

      if (allocated (v_times_polar)) then
        deallocate (v_times_polar)
      endif
      if (allocated (eigenvalues)) then
        deallocate (eigenvalues)
      endif
      if (allocated (lambda_grid)) then
        deallocate (lambda_grid)
      endif
      if (allocated (lambda_weight)) then
        deallocate (lambda_weight)
      endif

end subroutine evaluate_osrpa_integrand_2

subroutine evaluate_osrpa_response_2(i_spin,j_spin, spin_name,&
      v_times_polar, eigenvalues, &
      n_lambda, lambda_grid,lambda_weight,&
      rpa_c_integrand_spin,c_osrpa_integrand_spin)

    use dimensions, only: n_basbas, n_spin
    use runtime_choices, only: prodbas_threshold, safe_minimum, sc_check, c_osrpa_order, c_osrpa_threshold
    use scalapack_utils, only : sclpck_loc_ind
    use localorb_io
    use prodbas
    use mpi_tasks
    use synchronize_mpi
    implicit none
    integer :: i_spin, j_spin
    character :: spin_name(2)
    real*8  :: v_times_polar(max_row_2d,max_col_2d,n_spin)
    real*8  :: eigenvalues(n_basbas)
    integer :: n_lambda
    real*8  :: lambda_grid(n_lambda)
    real*8  :: lambda_weight(n_lambda)
    real*8  :: rpa_c_integrand_spin
    real*8  :: c_osrpa_integrand_spin(c_osrpa_order)
    integer :: n_nonsingular

    ! local variables
    real*8  :: special_radius
    real*8  :: delta_trace
    real*8  :: i_term
    integer :: i_basis_1, i_basis_2, i_order
    integer :: i_lambda
    integer :: i_index
    integer :: lr, lc

    real*8, dimension(:,:), allocatable :: temp_v_times_polar_a
    real*8, dimension(:,:), allocatable :: temp_v_times_polar_b

    if(.not.allocated(temp_v_times_polar_a)) then
       allocate(temp_v_times_polar_a(max_row_2d,max_col_2d),stat=i_index)
       call check_allocation(i_index, 'temp_v_times_polar_a           ')
    endif
    if(.not.allocated(temp_v_times_polar_b)) then
       allocate(temp_v_times_polar_b(max_row_2d,max_col_2d),stat=i_index)
       call check_allocation(i_index, 'temp_v_times_polar_b           ')
    endif

    if (i_spin .eq. 1 .and. j_spin.eq.2) then
        spin_name(1)='alpha'
        spin_name(2)='beta'
    else if (i_spin .eq. 2 .and. j_spin.eq.1) then
        spin_name(2)='alpha'
        spin_name(1)='beta'
    else if (i_spin .eq. 1 .and. j_spin.eq.1) then
        spin_name(2)='alpha'
        spin_name(1)='alpha'
    else if (i_spin .eq. 2 .and. j_spin.eq.2) then
        spin_name(2)='beta'
        spin_name(1)='beta'
    end if

    temp_v_times_polar_a = v_times_polar(:,:,j_spin)

    delta_trace = 0.0d0
    special_radius = 0.0d0

    if (sc_check(j_spin)) then
      !call diagonalize_auxmat_scalapack_real(n_basbas, temp_v_times_polar_a, &
      !& safe_minimum, -1.0d10, n_nonsingular, eigenvalues, &
      !polar_transform, " ")
      !special_radius=maxval(abs(eigenvalues))
      call diagonalize_auxmat_scalapack_real(n_basbas,temp_v_times_polar_a, &
                  prodbas_threshold, n_nonsingular, eigenvalues)
      special_radius=maxval(abs(eigenvalues))
      temp_v_times_polar_a = v_times_polar(:,:,j_spin)
    else
      special_radius = 0.0d0
    endif

    if (special_radius .lt. c_osrpa_threshold) then
        i_term = 1.0d0
        if (sc_check(j_spin)) then
            if (myid.eq.0) then
              write(use_unit,'(2X,A,X,A,A,f4.2,A,f9.6,A)') &
                  "Special radius of non-interacting response matrix in the", &
                  trim(spin_name(2)),"-spin channel < ",c_osrpa_threshold,&
                  " (",special_radius,")."
              write(use_unit,'(2X,A,X,A,X,A,A)') &
                  "Strong correlation is small in the OS-type", &
                  "particle-hole summation in the",trim(spin_name(1)),&
                  "-spin channel."
              write(use_unit,'(2X,A)') "Invoke standard perturbative solution."
            endif
            sc_check(j_spin)=.false.
        endif
        !write(use_unit,*) "igor debug dgemm", i_term, n_basbas
        !call dgemm &
        !( 'N','N',n_basbas,n_basbas,n_basbas,1.0d0, temp_v_times_polar_a,&
        !  n_basbas, v_times_polar(:,:,i_spin), n_basbas, 0.0d0, &
        !  temp_v_times_polar_b, n_basbas &
        !)
        call pdgemm('N','N',n_basbas,n_basbas,n_basbas, &
        &      1.0d0, temp_v_times_polar_a,     1,1, aux_sc_desc_2d, &
        &             v_times_polar(:,:,i_spin),1,1, aux_sc_desc_2d, &
        &      0.0d0, temp_v_times_polar_b,     1,1, aux_sc_desc_2d)
        !write(use_unit,*)  'debug 3', temp_v_times_polar_b(:, 1)
        ! compute the trace
        delta_trace = 0.0d0
        do i_basis_1 = 1, n_basbas
            lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
            lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
            if(lr>0 .and. lc>0) then
              delta_trace = delta_trace + &
                  temp_v_times_polar_b(lr, lc)
            endif
        enddo
        delta_trace = -1.0d0/(i_term+1.0d0)*delta_trace
        !write(use_unit,*) "igor debug", i_spin, int(i_term), delta_trace
        call sync_real_number(delta_trace)

        rpa_c_integrand_spin = rpa_c_integrand_spin + delta_trace
        !write(use_unit,*) "igor debug", i_spin, int(i_term), rpa_c_integrand_spin

        osRPA_LOOP: do while (abs(delta_trace) .gt. 1.0d-09 .and. &
                              int(i_term) .lt. 300 &
                              )
          i_term = i_term + 1.0d0
          !write(use_unit,*) "igor debug dgemm", i_term(i_spin), n_basbas
          !call dgemm &
          !( 'N','N',n_basbas,n_basbas,n_basbas,1.0d0, v_times_polar(:,:,j_spin),&
          !  n_basbas, temp_v_times_polar_b, n_basbas, 0.0d0, &
          !  temp_v_times_polar_a, n_basbas &
          !)
          call pdgemm('N','N',n_basbas,n_basbas,n_basbas, &
          &      1.0d0, v_times_polar(:,:,j_spin),1,1, aux_sc_desc_2d, &
          &             temp_v_times_polar_b,     1,1, aux_sc_desc_2d, &
          &      0.0d0, temp_v_times_polar_a,     1,1, aux_sc_desc_2d)
          ! compute the trace
          delta_trace = 0.0d0
          do i_basis_1 = 1, n_basbas
            lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
            lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
            if(lr>0 .and. lc>0) then
              delta_trace = delta_trace + &
                  temp_v_times_polar_a(lr, lc)
            endif
          enddo

          temp_v_times_polar_b = temp_v_times_polar_a

          delta_trace = -1.0d0/(i_term+1.0d0)*delta_trace
          !write(use_unit,*) "igor debug", i_spin, int(i_term), delta_trace

          call sync_real_number(delta_trace)

          rpa_c_integrand_spin = rpa_c_integrand_spin + delta_trace
          !write(use_unit,*) "igor debug", i_spin, int(i_term), rpa_c_integrand_spin
        end do osRPA_LOOP
    else
        if (myid.eq.0) then
          if (special_radius.ge.1.0d0) then
              write(use_unit,'(2X,A,X,A,A,f9.6,A)') &
                  "Special radius of non-interacting response matrix in the", &
                  trim(spin_name(2)),"-spin channel >= 1.0 (",special_radius,")."
              write(use_unit,'(2X,A,X,A,X,A,A)') &
                  "Strong correlation breaks the perturbative OS-type", &
                  "particle-hole summation in the",trim(spin_name(1)),&
                  "-spin channel."
          else
              write(use_unit,'(2X,A,X,A,A,f4.2,A,f9.6,A)') &
                  "Special radius of non-interacting response matrix in the", &
                  trim(spin_name(2)),"-spin channel >= ",c_osrpa_threshold,&
                  " (",special_radius,")."
              write(use_unit,'(2X,A,X,A,X,A,A)') &
                  "Strong correlation is large in the OS-type", &
                  "particle-hole summation in the",trim(spin_name(1)),&
                  "-spin channel."
          endif
          write(use_unit,'(2X,A)') &
              "Invoke renormalization to properly consider the strong correlation"
        endif

        rpa_c_integrand_spin = 0.d0
        c_osrpa_integrand_spin = 0.d0

        do i_basis_1 = 1, n_basbas
          lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
          lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
          if(lr>0 .and. lc>0) then
            rpa_c_integrand_spin = rpa_c_integrand_spin + &
              v_times_polar(lr,lc,i_spin)
          endif
        enddo
        !write(use_unit,*) "igor debug 1", rpa_c_integrand_spin

        temp_v_times_polar_a = v_times_polar(:,:,j_spin)

        temp_v_times_polar_b(:,:) = 0.0d0
        do i_lambda = 1, n_lambda
            temp_v_times_polar_a =  -lambda_grid(i_lambda)*v_times_polar(:,:,j_spin)
            do i_basis_1 = 1, n_basbas
              lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
              lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
              if(lr>0 .and. lc>0) then
                temp_v_times_polar_a(lr,lc) = &
                    1.0d0 + temp_v_times_polar_a(lr,lc)
              endif
            enddo
            !call  power_auxmat_lapack(temp_v_times_polar_b,-1.0d0," ")
            !write(use_unit,*) "igor debug inverse"
            !MARK IGOR
            call power_auxmat_scalapack_2d(temp_v_times_polar_a, &
                                           -1.0d0, &
                                           '')
            temp_v_times_polar_b(:,:) = temp_v_times_polar_b(:,:) + &
                lambda_weight(i_lambda)*temp_v_times_polar_a
        enddo

        call pdgemm('N','N',n_basbas,n_basbas,n_basbas, &
        &      1.0d0, temp_v_times_polar_b,     1,1, aux_sc_desc_2d, &
        &             v_times_polar(:,:,i_spin),1,1, aux_sc_desc_2d, &
        &      0.0d0, temp_v_times_polar_a,     1,1, aux_sc_desc_2d)
        do i_basis_1 = 1, n_basbas
            lr = sclpck_loc_ind(i_basis_1, myprow_aux_2d, nprow_aux_2d, nb_aux_2d, 0)
            lc = sclpck_loc_ind(i_basis_1, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
            if(lr>0 .and. lc>0) then
              rpa_c_integrand_spin = rpa_c_integrand_spin - &
                  temp_v_times_polar_a(lr, lc)
            endif
        enddo
        call sync_real_number(rpa_c_integrand_spin)
    end if
    if (allocated (temp_v_times_polar_a)) then
      deallocate (temp_v_times_polar_a)
    endif
    if (allocated (temp_v_times_polar_b)) then
      deallocate (temp_v_times_polar_b)
    endif
end subroutine evaluate_osrpa_response_2
