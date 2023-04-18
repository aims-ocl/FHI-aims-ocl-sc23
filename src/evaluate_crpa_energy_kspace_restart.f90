!****s* FHI-aims/evaluate_crpa_energy_kspace_restart
!  NAME
!   evaluate_crpa_energy_kspace
!  SYNOPSIS

      subroutine evaluate_crpa_energy_kspace_restart &
           (n_low_state, occ_numbers, n_full_freq, &
            omega_full, womega_full, &
            KS_eigenvalue, KS_eigenvector, &
            KS_eigenvector_complex, &
            rpa_c_energy &
           )

!  PURPOSE
!  Subroutine evaluate_crpa_energy_kspace_restart evaluates the correlation
!  energy at the RPA level using the adiabatic connection fluctuation
!  dissipation theorem.
!
!  E_RPA = 1/2pi \int dw { ln(det(1-v_times_polar)) + tr(v_times_polar) }

! USES
      use dimensions
      use prodbas
      use pbc_lists
      use hartree_fock
      use constants
      use mpi_tasks
      use synchronize_mpi
      use timing
      ! Igor for the restarting utility
      use runtime_choices
      use restart_rpa
      use evaluate_polarisability_kspace_mod
      use localorb_io, only: use_unit
      implicit none

! ARGUMENTS 

      integer :: n_full_freq
      integer :: n_low_state
      integer :: n_high_state
!      integer :: n_lumo(n_spin)
!      integer :: n_homo(n_spin)

      real*8  :: occ_numbers(n_states,n_spin,n_k_points)
      real*8  :: omega_full(n_full_freq)
      real*8  :: womega_full(n_full_freq)
      real*8  :: KS_eigenvalue(n_states,n_spin,n_k_points)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin,n_k_points_task)
      complex*16  :: KS_eigenvector_complex(n_basis,n_states,n_spin,n_k_points_task)

!     output
      real*8  :: rpa_c_energy

! INPUTS
! o  n_full_freq -- integer number,
!            the number of frequency points for the screened Coulomb interaction W
! o  n_low_state  -- integer number,
!            the lowest KS/HF eigenstate taken into account in the polarisability calcua            ltions
! o  n_high_state -- integer number,
!            the highest KS/HF eigenstate. In the present case, n_high_state >= n_homo
!            should be fine. 
! o  n_electrons -- real number
!            the total number of electrons in the system
! o  occ_numbers -- real 2-dimentianal array of length (n_states, n_spin)
!            the occupation number of the electrons for each eigenstate and each spin
! o  omega_full(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  womega_full(n_freq) -- real array
!            the weigth of the Gauss-Legendre frequency grid for the screened Coulomb 
!            in teraction
! o  chemical_potential -- real number, the chemical potential of the system
! o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
! o  KS_eigenvector -- real array,
!            the eigenvector of the single-particle calculation
! o  KS_eigenvector_complex -- complex array,
!            the complex eigenvector of the single-particle calculation,
!            used when "real_eigenvectors == .false."
!           
!
! OUTPUT
! o  rpa_c_energy -- real number, the calculated RPA correlation energy
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

      complex*16  det_v_times_polar
      complex*16  trace_v_times_polar
      complex*16  trace_v(n_irk_points), trace_polar(n_irk_points), trace_v_polar_output(n_irk_points)
      real*8  rpa_c_integrand
      real*8, dimension(:), allocatable :: rpa_c_kgrid

      ! for restart
      real*8  :: crpa_each_thread, current_crpa_original
      
!      integer :: n_k_points_original
!      real*8, dimension(:,:), allocatable :: k_point_list_old

!    local timing
      real*8  temp_time_rpa
      real*8  temp_clock_time_rpa

!     auxiliary matrices for Level 3 Blas matrix multiplications
!     n_first : the first orbital which is NOT fully occupied

      integer :: n_first(n_spin)
      integer :: ipiv(n_basbas)
      integer :: info

      complex*16, dimension(:,:,:), allocatable :: polar_kspace
      complex*16, dimension(:,:), allocatable :: v_times_polar
      complex*16, dimension(:,:), allocatable :: coulomb_tmp

!     timing

!     parameters of the fitting tails
!       real*8  s1, s2, omega_1, omega_2
!       real*8  alpha, beta, a, b

      character*50  filename

!     counters

      integer :: i_state
      integer :: i_freq
      integer :: i_spin
      integer :: i_index
      integer :: i_k_point
      integer :: i_k_point_local
      integer :: i_irk_point
      integer :: i_irk_point_local
      integer :: i_prodbas_1
      integer :: i_prodbas_2
      integer :: i_task
      integer :: id_recv
      integer :: id_send

!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"-----------------------------------------------------------------------"
        write(use_unit,'(2X,A)') &
              "Start to calculate the periodic RPA correlation energy  ... "
      endif

      if(flag_KS_eigenfunc_conjg) then
         KS_eigenvector_complex = conjg(KS_eigenvector_complex)
      endif
!     determine the highest occupied orbital level
!     such complication occurs when some of the orbitals are
!     either not fully occupied or not fully empty
      n_first(:) = 1
      do i_k_point = 1, n_k_points, 1
        do i_spin = 1, n_spin
          do i_state = 1, n_states
           if (abs(occ_numbers(i_state,i_spin,i_k_point)-dble(2/n_spin)) &
                           .lt.1.d-8) then
             n_first(i_spin)= i_state + 1
           endif
          enddo
          if(n_first(i_spin) .gt. n_states) then
           n_first(i_spin) = n_states
          endif
        enddo
      enddo

      allocate(polar_kspace(n_basbas, n_basbas, n_irk_points_task),stat=i_index)
      call check_allocation(i_index, 'polar_kspace                    ')

      allocate(v_times_polar(n_basbas, n_basbas),stat=i_index)
      call check_allocation(i_index, 'v_times_poalr                   ')

      allocate(coulomb_tmp(n_basbas, n_basbas),stat=info)
      call check_allocation(info, 'coulomb_tmp                     ')

! work array

      allocate(rpa_c_kgrid(n_irk_points),stat=i_index)

!      do i_k_point = 1, n_k_points, 1
!         do i_prodbas_1 = 1, n_basbas, 1
!           write(use_unit,'(3I4,4f18.8)') i_k_point, i_prodbas_1, basbas_l(i_prodbas_1), &
!                 coulomb_matr_recip(i_prodbas_1,i_prodbas_1,i_k_point) , &
!!                multipole_basbas_fn(basbas_fn(i_prodbas_1))
!         enddo
!      enddo

      time_rpa_corr = 0.d0
      clock_time_rpa_corr = 0.d0
      time_polar = 0.d0
      clock_time_polar = 0.d0

      rpa_c_energy = 0.d0 
      rpa_c_kgrid = 0.d0 
      trace_v =(0.d0,0.d0)
      trace_polar =(0.d0,0.d0)
      trace_v_polar_output =(0.d0,0.d0)
      current_crpa_original = 0.d0

      n_low_state=max(1,n_low_state)
      if (flag_frozen_core_postSCF) then ! count the frozen core states, Igor
          call count_frozen_core_states(n_low_state)
      endif
      if (myid .eq. 0) then
          write(use_unit,'(2X,A,I12)') &
              'The first valence state in the periodic RPA calculation :', &
              n_low_state
      endif

      
      if (restart_rpa_read) then ! reading the restart info., Igor
          call read_restart_rpa_info()
          current_crpa_original = current_crpa
      endif

      do i_freq = rpa_freq_start+1, n_full_freq, 1

        crpa_each_thread = 0.d0

        rpa_c_integrand = 0.d0 

        if(myid.eq.0) then
           write(use_unit,'(A,I4,f20.6)') " | i_freq ", i_freq, & 
                   omega_full(i_freq)
        endif

         call get_timestamps(temp_time_rpa, temp_clock_time_rpa )

!    evaluate the polarisability at frequency point i_freq
!         write(use_unit,*) "occ_numbers now", occ_numbers
         call  evaluate_polarisability_kspace &
              (n_low_state,  omega_full(i_freq), KS_eigenvalue, KS_eigenvector, &
                 KS_eigenvector_complex, &
                 occ_numbers, polar_kspace )


         call get_timestamps(rtime, clock_rtime)
         time_polar = time_polar + rtime - temp_time_rpa
         clock_time_polar = clock_time_polar + clock_rtime - temp_clock_time_rpa

!          polar_kspace(:,:,:) = polar_kspace(:,:,:) * 2.d0/dble(n_spin)
! XR: we don't need this in peridoc calculations. Note that this is different from the cluster case.

! n_k_points_original is the number of k points in the original mesh
          do i_k_point = 1, n_k_points, 1

             if(.not. irk_point_included(i_k_point) ) cycle

             i_k_point_local = (i_k_point-1)/n_tasks + 1
             i_irk_point = irk_point_mapping(i_k_point)

             i_task = mod(i_k_point, n_tasks)
             id_recv = mod(i_irk_point, n_tasks)
             if(myid.eq.i_task) then
               coulomb_tmp(:,:) = coulomb_matr_recip(:,:,i_k_point_local)
             endif
             if(i_task.ne.id_recv) then
               if(myid.eq.i_task) then
                 call send_complex_vector(coulomb_tmp,n_basbas*n_basbas,id_recv)
               elseif (myid.eq.id_recv) then
                 call receive_complex_vector(coulomb_tmp,n_basbas*n_basbas,i_task)
               endif
             endif

            if(myid .ne. mod(i_irk_point, n_tasks)) cycle
            i_irk_point_local = (i_irk_point-1)/n_tasks + 1

            if(i_freq == 1) then
              do i_prodbas_1 = 1, n_basbas, 1
                 trace_v(i_irk_point) = trace_v(i_irk_point) + coulomb_tmp(i_prodbas_1,i_prodbas_1)
                 trace_polar(i_irk_point) = trace_polar(i_irk_point) + &
                     polar_kspace(i_prodbas_1,i_prodbas_1,i_irk_point_local)
              enddo
            endif

!         if(i_freq .eq. 1) then
!              do i_prodbas_1 = 1, n_basbas, 1
!                  write(use_unit,'(3I4,4f18.8)')i_k_point, i_prodbas_1, basbas_l(i_prodbas_1), &
!                      coulomb_matr_recip(1,i_prodbas_1,i_k_point) , &
!                      polar_kspace(1, i_prodbas_1,i_irk_point_local)
!               enddo
!          endif
!  Multiply \chi_0 with bare coulomb matrix v.
            call zgemm('N', 'N', n_basbas, n_basbas, n_basbas, (1.d0,0.d0), &
                   coulomb_tmp, n_basbas, &
                   polar_kspace(1,1,i_irk_point_local), n_basbas, (0.d0, 0.d0), &
                   v_times_polar, n_basbas)
            
!            v_times_polar(:,:) = (0.d0,0.d0)
!            do i_prodbas_1 = 1, n_basbas, 1
!              do i_prodbas_2 = 1, n_basbas, 1
!                 do i_index = 1, n_basbas, 1
!                    v_times_polar(i_prodbas_2,i_prodbas_1) = &
!                       v_times_polar(i_prodbas_2,i_prodbas_1) + &
!                       coulomb_tmp (i_index,i_prodbas_2)* &
!                       polar_kspace(i_index,i_prodbas_1,i_irk_point_local)
!                 enddo
!              enddo
!            enddo

            trace_v_times_polar = (0.d0,0.d0)
            do i_prodbas_1 = 1, n_basbas, 1
              trace_v_times_polar = trace_v_times_polar + &
               v_times_polar(i_prodbas_1, i_prodbas_1)
!              if(i_freq .eq. 1) then
!               write(use_unit,'(3I4,4f18.8)') i_k_point, i_prodbas_1, basbas_l(i_prodbas_1), v_times_polar(i_prodbas_1, i_prodbas_1), trace_v_times_polar
!              endif

               v_times_polar(i_prodbas_1, i_prodbas_1) = v_times_polar(i_prodbas_1, i_prodbas_1) - (1.d0,0.d0)
            enddo
            if(i_freq.eq.1) then
              trace_v_polar_output (i_irk_point) = trace_v_times_polar
            endif

            v_times_polar(:,:) = - v_times_polar(:,:)

            call  zgetrf (n_basbas,n_basbas,v_times_polar,n_basbas,ipiv,info)

            det_v_times_polar = (1.d0,0.d0)
            do i_prodbas_1 = 1, n_basbas, 1
             det_v_times_polar = det_v_times_polar * v_times_polar(i_prodbas_1,i_prodbas_1)
!              if(i_freq .eq. 1) then
!               write(use_unit,'(I4,4f20.8)') i_prodbas_1, v_times_polar(i_prodbas_1,i_prodbas_1), det_v_times_polar
!              endif
            enddo

            rpa_c_integrand = log (abs(real(det_v_times_polar))) + &
                                 real(trace_v_times_polar)

            crpa_each_thread = crpa_each_thread + &
                         rpa_c_integrand * womega_full(i_freq) * &
                         irk_weight(i_irk_point)


            rpa_c_energy = rpa_c_energy + &
                         rpa_c_integrand * womega_full(i_freq) * &
                         irk_weight(i_irk_point)


            rpa_c_kgrid(i_irk_point) = rpa_c_kgrid(i_irk_point) + &
                         rpa_c_integrand * womega_full(i_freq) 

!           if(i_k_point.eq.1 .and. i_freq .eq. 1) then
!            write(use_unit,'(2I4,6f20.8)') i_k_point, i_freq, rpa_c_integrand, rpa_c_energy, &
!                    det_v_times_polar, trace_v_times_polar
!           endif
! end of loop over i_k_point            
          enddo

!      write(81,'(I6,3f20.8)')i_freq, omega_full(i_freq), &
!                  womega_full(i_freq), rpa_c_integrand

       call get_timestamps(temp_time_rpa, temp_clock_time_rpa )
       time_rpa_corr = time_rpa_corr +  temp_time_rpa - rtime
       clock_time_rpa_corr = clock_time_rpa_corr + temp_clock_time_rpa - clock_rtime

       rpa_freq_start = rpa_freq_start + 1
       call sync_real_number(crpa_each_thread)
       crpa_each_thread=crpa_each_thread/2.d0/pi
       current_crpa = current_crpa + crpa_each_thread
       call write_restart_rpa_info()
       if(myid .eq. 0) then
           write(use_unit,'(2X,I8,I8,f21.10)') n_full_freq, rpa_freq_start, current_crpa
       endif

! end of loop over i_freq
      enddo
      call sync_timing(time_polar)
      call sync_timing(time_rpa_corr)
      call sync_vector(rpa_c_kgrid,n_irk_points)
      call sync_vector_complex(trace_v,n_irk_points)
      call sync_vector_complex(trace_polar,n_irk_points)

      if(myid.eq.0) then
        do i_k_point = 1, n_k_points, 1
            if(.not. irk_point_included(i_k_point) ) cycle

            i_irk_point = irk_point_mapping(i_k_point)
            write(use_unit,'(2I6,4f18.8)') i_k_point, i_irk_point, k_point_list(i_k_point,:), &
                  rpa_c_kgrid(i_irk_point)
            write(use_unit,'(6f18.8)') trace_v(i_irk_point), trace_polar(i_irk_point), trace_v_polar_output(i_irk_point)
        enddo
      endif

      call sync_real_number(rpa_c_energy)     
      rpa_c_energy=rpa_c_energy/2.d0/pi + current_crpa_original
!      rpa_c_energy= - rpa_c_energy/2.d0

! Delta term correction
! This deals with the case when there are frational occupations and 
! the intra-orbtial excitations should be taken into account. The contribution
! is a delta function at zero frequency and a finite contribution when integrated
! out over frequency axis.

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"----------------------------------------------------", &
                  "-------------------------"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            " RPA correlation energy :", rpa_c_energy, "Ha,", &
             rpa_c_energy*hartree, "eV"
!        write(use_unit,*)"----------------------------------------------------", &
!                 "-------------------------"
        write(use_unit,*)
      endif

      if (allocated (polar_kspace)) then
        deallocate (polar_kspace)
      endif
      if (allocated (v_times_polar)) then
        deallocate (v_times_polar)
      endif
      if (allocated (coulomb_tmp)) then
        deallocate (coulomb_tmp)
      endif

      return

      end subroutine evaluate_crpa_energy_kspace_restart

