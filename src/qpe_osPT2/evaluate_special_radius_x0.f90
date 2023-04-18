!****s* FHI-aims/evaluate_special_radius_x0
!  NAME
!   evaluate_special_radius_x0
!  SYNOPSIS

subroutine evaluate_special_radius_x0 &
           (polar_freq,special_radius)

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
      use runtime_choices, only: prodbas_threshold, safe_minimum
      use localorb_io
      implicit none

! ARGUMENTS 
      real*8  :: polar_freq(n_basbas, n_loc_prodbas, n_spin)
      real*8  :: rpa_c_integrand, rpa_c_integrand_total
      real*8  :: special_radius(2)

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

      real*8, dimension(:,:,:), allocatable :: v_times_polar
      real*8, dimension(:,:), allocatable :: temp_v_times_polar_a
      real*8, dimension(:,:), allocatable :: temp_v_times_polar_b
      real*8, dimension(:,:), allocatable :: polar_transform
      real*8, dimension(:), allocatable :: eigenvalues

      real*8, dimension(:,:,:), allocatable :: polar_freq_remote
      integer n_loc_prodbas_remote

      real*8 rpa_c_integrand_spin
      real*8 n_nonsingular
      

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
      integer :: i_state
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_basis_3
      integer :: i_index
      integer :: i_spin, j_spin
      integer :: i_task
      integer :: lr, lc
      character :: spin_name(2)

!     begin work

    special_radius = 0.0d0

    if(myid.eq.0) then
       if(.not.allocated(v_times_polar)) then
          allocate(v_times_polar(n_basbas,n_basbas,n_spin),stat=i_index)
          call check_allocation(i_index, 'v_times_polar                 ')
       endif
       v_times_polar(:,:,:) = 0.d0
       if(.not.allocated(polar_transform)) then
          allocate(polar_transform(n_basbas,n_basbas),stat=i_index)
          call check_allocation(i_index, 'polar_transform             ')
       endif
       polar_transform(:,:) = 0.0d0
       if(.not.allocated(eigenvalues)) then
          allocate(eigenvalues(n_basbas),stat=i_index)
          call check_allocation(i_index, 'eigenvalues                 ')
       endif
       eigenvalues(:) = 0.0d0
    endif

    id_recv = 0
    do i_task = 1, n_tasks

      id_send = i_task - 1
      tag= id_send
      if(id_send.ne.id_recv) then
        n_loc_prodbas_remote = COUNT(map_prodbas(:,i_task)>0)
        if (myid.eq.id_send) then
            allocate(polar_freq_remote(n_basbas,n_loc_prodbas_remote,n_spin),stat=i_index)
            polar_freq_remote(:,:,:) = polar_freq(:,:,:)
            call MPI_SEND (polar_freq_remote, &
                 n_basbas*n_loc_prodbas_remote*n_spin, &
                 MPI_DOUBLE_PRECISION, id_recv, &
                 tag, mpi_comm_global,mpierr)
            deallocate(polar_freq_remote)
        else if (myid.eq.id_recv) then
           allocate(polar_freq_remote(n_basbas,n_loc_prodbas_remote,n_spin),stat=i_index)
           polar_freq_remote = 0.0d0
           call MPI_RECV (polar_freq_remote, &
                n_basbas*n_loc_prodbas_remote*n_spin, &
                MPI_DOUBLE_PRECISION, id_send, &
                tag, mpi_comm_global,my_status,mpierr)

           do i_basis_1 = 1, n_loc_prodbas_remote
             i_index=map_prodbas(i_basis_1,i_task)
             do i_spin = 1, n_spin
               if(i_index.gt.0) then
                  v_times_polar(:,i_index,i_spin) = &
                  polar_freq_remote(:,i_basis_1,i_spin)
               endif
             enddo
           enddo
           deallocate(polar_freq_remote)
        end if
      else ! if(id_recv .ne. id_send)
        if (myid.eq.0) then
           do i_basis_1 = 1, n_loc_prodbas
            i_index=map_prodbas(i_basis_1,i_task)
            do i_spin = 1, n_spin
               v_times_polar(:,i_index,i_spin) = &
               polar_freq(:,i_basis_1,i_spin)
            enddo
           enddo
        end if
      endif
! end of loop over i_task
     enddo

     if(myid.eq.0) then 
        if (n_spin .eq. 1) then
          ! IGOR debug
          !open(UNIT=666,FILE='chi_c.out',STATUS='NEW',FORM='FORMATTED')
          !write(666,*) v_times_polar(:,:,1)
          !close(666)

          ! check if the system is stronly correlated or not.
          !temp_v_times_polar_a = v_times_polar(:,:,1)
          call diagonalize_auxmat_lapack(n_basbas, v_times_polar(:,:,1), &
          & safe_minimum, -1.0d10, n_nonsingular, eigenvalues, &
          polar_transform, " ")
          special_radius=maxval(abs(eigenvalues))
          !temp_v_times_polar_a = v_times_polar(:,:,1)
        else ! n_spin !=1

          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! IGOR debug
          !open(UNIT=666,FILE='chi_a.out',STATUS='NEW',FORM='FORMATTED')
          !write(666,*) v_times_polar(:,:,1)
          !close(666)
          !open(UNIT=888,FILE='chi_b.out',STATUS='NEW',FORM='FORMATTED')
          !write(888,*) v_times_polar(:,:,2)
          !close(888)
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          do i_spin = 1, n_spin 
            if (i_spin .eq. 1) then 
                j_spin = 2 
                spin_name(1)='alpha'
                spin_name(2)='beta'
            else
                j_spin = 1 
                spin_name(2)='alpha'
                spin_name(1)='beta'
            end if

            !temp_v_times_polar_a = v_times_polar(:,:,j_spin)

            call diagonalize_auxmat_lapack(n_basbas, v_times_polar(:,:,j_spin), &
            & safe_minimum, -1.0d10, n_nonsingular, eigenvalues, &
            polar_transform, " ")
            special_radius(j_spin)=maxval(abs(eigenvalues))
          end do
        endif
        write(use_unit,'(2X,A,2f16.8,A)') &
            "Special radius of non-interacting response matrix in each spin channel = ", &
            special_radius(1),special_radius(2),"." 
     endif

      call MPI_BCAST(special_radius,2,MPI_DOUBLE_PRECISION,0,mpi_comm_global,mpierr)

     if (allocated (v_times_polar)) then
       deallocate (v_times_polar)
     endif
     if (allocated (polar_transform)) then
       deallocate (polar_transform)
     endif
     if (allocated (eigenvalues)) then
       deallocate (eigenvalues)
     endif

      end subroutine evaluate_special_radius_x0
