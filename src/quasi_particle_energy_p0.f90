!****s* FHI-aims/quasi_particle_energy_p0
!  NAME
!   quasi_particle_energy_p0
!  SYNOPSIS

      subroutine quasi_particle_energy_p0 &
           (anacon_type, n_max_par, &
            n_low_state, n_high_state, &
            n_freq, omega, &
            sigma_par_p0, occ_numbers, KS_eigenvalue, &
            chemical_potential_spin, &
            exact_x_kspace, xc_kspace, &
            qp_energy )

!  PURPOSE
!  Subroutine quasi_partilce_energy_p0 calculate the quasi particle energies
!  using precomputed self-energies for a periodic system.  
!
!  USES

      use dimensions
      use runtime_choices
      use pbc_lists
      use mpi_tasks
      use synchronize_mpi_basic
      use localorb_io
      use constants

      implicit none

!  ARGUMENTS

      integer :: anacon_type 
      integer :: n_max_par
      integer :: n_low_state, n_high_state
      integer :: n_freq

      real*8 :: omega(n_freq)
      real*8 :: occ_numbers(n_states,n_spin,n_k_points)
      real*8 :: chemical_potential_spin(n_spin)

      real*8, dimension(n_states,n_spin,n_k_points) :: KS_eigenvalue
      real*8, dimension(n_low_state:n_high_state,n_spin,n_irk_points_task) :: xc_kspace
      real*8, dimension(n_low_state:n_high_state,n_spin,n_irk_points_task) :: exact_x_kspace
      complex*16 :: sigma_par_p0(n_max_par,n_low_state:n_high_state,n_spin,n_irk_points_task)

      real*8 :: qp_energy(n_low_state:n_high_state,n_spin,n_irk_points_task)

!  INPUTS
!  o  anacon_type -- integer number, if 0, the two-pole fitting for analytical
!          continuation; if 1, using Pade approximation for ana. cont.        
!  o  n_max_par -- the number of parameters used for analytical continuation 
!          For anacon_type = 0, recommended n_max_par is  4 or 5. If 4, this will 
!          be the normal two-pole fitting, else if 5, it will be two-pole plus a 
!          (small) constant number
!          For anacon_type = 1, recommended n_max_par is the half of n_freq  
!  o  n_low_state  -- integer number,
!          the lowest KS/HF eigenstate for self-energy correction
!  o  n_high_state -- integer number,
!          the highest KS/HF eigenstate for self-energy correction
!  o  n_freq -- integer number, the number of frequency points for the GW self-energy
!  o  omega(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the self-energy
!  o  occ_numbers -- occupation numbers of single-particle energy levels
!  o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
!  o  chemical_potential_spin -- real number, the chemical potential of the system
!  o  xc_kspace -- the DFT xc term for the single-particle orbital energy
!  o  exact_x_kspace -- the exact-exchage contribution to the GW single-particle orbital energy
!  o  sigma_par_p0 -- complex array, the fitting parameters from analytical continuation
! OUTPUTS
!  o  qp_energy -- real array, the final calculated quasiparticle energy for each
!            concerned state (from n_low_state to n_high_state).

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

      real*8, allocatable :: qp_energy_tmp(:,:)
      real*8, allocatable :: correl_energy(:,:,:)
      real*8, allocatable :: correl_energy_tmp(:,:)
      real*8, allocatable :: exact_x_tmp(:,:)
      real*8, allocatable :: xc_tmp(:,:)
      real*8, allocatable :: qp_vbm_k(:,:)
      real*8, allocatable :: qp_cbm_k(:,:)
      real*8, dimension(n_spin)  :: qp_vbm
      real*8, dimension(n_spin)  :: qp_cbm
      real*8, dimension(n_spin)  :: direct_band_gap
      integer, dimension(n_spin)  :: direct_gap_k_point
      integer, dimension(n_spin)  :: vbm_k_point
      integer, dimension(n_spin)  :: cbm_k_point
!      real*8 :: correl_energy(n_low_state:n_high_state,n_spin,n_k_points_task)

      real*8  e_diff
      real*8  en
      real*8  mu
      real*8  delta_mu

!      real*8 :: qpe_sum

      complex*16 selfe, dselfe
      complex*16 exchange_tmp

      real*8, parameter ::  qp_energy_thr = 1.d-5

      integer :: n_states_count
      integer :: info, mpierr
      character(*), parameter :: func = 'quasi_particle_energy_p0.f90'
      character*150 :: info_str

!     counters


      integer :: i_state
      integer :: i_state_1
      integer :: i_count
      integer :: i_spin
      integer :: i_k_point, i_k_point_local
      integer :: i_k_point_1, i_k_point_local_1
      integer :: i_irk_point, i_irk_point_local
      integer :: i_basis_1, i_basis_2
      integer :: i_task, i_task_1

!   external function


!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,2A)')   "Quasi particle energy calculation ", &
          "for periodic systems starts ..."
      endif

      if(use_hartree_fock .and. (.not.use_screx) &
            .and. .not. use_gw_and_hse) then
        exact_x_kspace = exact_x_kspace * (1.0-hybrid_coeff)
      endif

!  read in DFT exchange-correlation energy for all states
!      allocate(xc_kspace(n_low_state:n_high_state,n_spin,n_k_points_task))
!      call check_allocation(info, 'xc_kspace', func)
!      if(use_split_xc_gw) then
!          xc_kspace(n_low_state:n_high_state,:,:) = &
!          x_KS_array(n_low_state:n_high_state,:,:)+c_KS_array(n_low_state:n_high_state,:,:)
!      else
!          xc_kspace(n_low_state:n_high_state,:,:) = &
!          xc_KS_matr(n_low_state:n_high_state,:,:)
!      endif

!  quasi particle energy calculation

!      if(allocated(qp_energy_tmp)) then
!        deallocate(qp_energy_tmp)
!     endif
      allocate(qp_energy_tmp(n_low_state:n_high_state,n_spin),stat=info)
      call check_allocation(info, 'qp_energy_tmp',func)
      allocate(correl_energy(n_low_state:n_high_state,n_spin,n_irk_points_task),stat=info)
      call check_allocation(info, 'correl_energy',func)

      do i_irk_point = 1, n_irk_points, 1
         if (myid .ne. mod(i_irk_point,n_tasks) ) cycle

         i_k_point = inv_irk_point_mapping(i_irk_point)
         i_k_point_local = (i_k_point-1)/n_tasks + 1
         i_irk_point_local = (i_irk_point-1)/n_tasks + 1

         qp_energy(n_low_state:n_high_state,:,i_irk_point_local)= &
                KS_eigenvalue(n_low_state:n_high_state,:,i_k_point)

         qp_energy_tmp(n_low_state:n_high_state,:)= &
         KS_eigenvalue(n_low_state:n_high_state,:,i_k_point)

         do i_spin = 1, n_spin, 1
           do i_state = n_low_state, n_high_state, 1
   
             e_diff = 1.d-3
             i_count =0
             do while (abs(e_diff).gt.qp_energy_thr)
               i_count = i_count +1
               qp_energy(i_state,i_spin,i_irk_point_local) = &
                     qp_energy_tmp(i_state,i_spin) + 0.5d0* e_diff
               qp_energy_tmp(i_state,i_spin) = qp_energy(i_state,i_spin,i_irk_point_local)
   
               mu =  chemical_potential_spin(i_spin)
   
               en = qp_energy(i_state,i_spin,i_irk_point_local)-mu
   
               call get_real_selfenergy(anacon_type,n_freq,omega, &
                         dcmplx(en,0.d0), n_max_par, &
                         sigma_par_p0(1:n_max_par,i_state,i_spin,i_irk_point_local), selfe)
   
   
               qp_energy(i_state,i_spin,i_irk_point_local) = &
                          KS_eigenvalue(i_state,i_spin,i_k_point) &
                        + real(selfe) &
                        + exact_x_kspace(i_state,i_spin,i_irk_point_local) &
                        - xc_kspace(i_state,i_spin,i_irk_point_local)
   
               e_diff =  qp_energy(i_state,i_spin,i_irk_point_local) &
                        - qp_energy_tmp(i_state,i_spin)
   
!               if(myid.eq.0 .and. i_irk_point_local .eq.1) then
!                   write(use_unit,'(I4,4f18.6)') i_state, KS_eigenvalue(i_state,i_spin,i_k_point), &
!                            exact_x_kspace(i_state,i_spin,i_irk_point_local), &
!                            xc_kspace(i_state,i_spin,i_irk_point_local), &
!                            real(selfe)
!               endif
              if(i_count .gt. 100) then
                   write(use_unit,'(2X,3A,I4,A,I4 )') &
                  " * Error: QUASI_PARTILCE_ENERGY: self-consistent", &
                  " quasiparticle solution can not be found for", &
                  " i_state = ",  i_state, "  i_spin = ", i_spin
   
               exit
              endif
   
! end of do while
            enddo
   
            correl_energy(i_state,i_spin,i_irk_point_local) = real(selfe)
!            write(use_unit,'(4I4,f16.8)') i_irk_point, i_k_point, i_spin, i_state, &
!                      correl_energy(i_state,i_spin,i_irk_point_local)
! end of do i_state
            enddo
! end of do i_spin
           enddo
! end of loop over i_irk_point
       enddo

! add the information of other k points    
!      do i_k_point = 1, n_k_points, 1
!         i_k_point_local = (i_k_point-1)/n_tasks + 1
!
!         if(.not.irk_point_included(i_k_point)) then
!            i_irk_point = irk_point_mapping(i_k_point)
!            i_k_point_1 = inv_irk_point_mapping(i_irk_point)
!            
!            if(myid.eq. mod(i_k_point_1, n_tasks)) then
!              i_k_point_local_1 = (i_k_point_1-1)/n_tasks + 1
!              correl_energy(:,:,i_k_point_local) = correl_energy(:,:,i_k_point_local_1)
!            endif
!
!         endif
!      enddo

      if(myid.eq.0) then

        write(use_unit,*)
        write(use_unit,'(2A)')"--------------------------------------", &
         "------------------------------------------------------"
        write(use_unit,'(15X,A)')"GW quasi-particle energy levels"
        write(use_unit,*)
        write(use_unit,'(15X,A,I6)')"| # of k_points in 1st BZ: ", n_k_points
        write(use_unit,*)
        write(use_unit,'(15X,A)')"e_qp = e_gs + e_x^ex - e_xc^gs + e_c^nloc"
        write(use_unit,*)
        write(use_unit, '(2X, A, 5X, A,8X, A, 8X,A,8X,A,8X,A, 8X, A)') &
             "state", "occ_num", "e_gs", "e_x^ex", "e_xc^gs", &
              "e_c^nloc", "e_qp"
        write(use_unit,'(2A)')"-----------------------------------------", &
         "-------------------------------------------------------------"
! end of if myid.eq.0
      endif
      call mpi_barrier(mpi_comm_global,mpierr)

      n_states_count = n_high_state - n_low_state + 1
      allocate(exact_x_tmp(n_states_count,n_spin),stat=info)
      call check_allocation(info, 'exact_x_tmp',func)
      allocate(xc_tmp(n_states_count,n_spin),stat=info)
      call check_allocation(info, 'xc_tmp',func)
      allocate(correl_energy_tmp(n_states_count,n_spin),stat=info)
      call check_allocation(info, 'correl_energy_tmp',func)

      do i_k_point =1, min(n_k_points,out_k_points_eigenvalues), 1

        i_k_point_local = (i_k_point-1)/n_tasks+1
        i_irk_point = irk_point_mapping(i_k_point)
        i_irk_point_local = (i_irk_point-1)/n_tasks+1
   
        i_task = mod(i_irk_point,n_tasks)
        if(myid.eq.i_task) then
          exact_x_tmp(1:n_states_count,:) = exact_x_kspace(n_low_state:n_high_state,:,i_irk_point_local)
          xc_tmp(1:n_states_count,:) = xc_kspace(n_low_state:n_high_state,:,i_irk_point_local)
          correl_energy_tmp(1:n_states_count,:) = correl_energy(n_low_state:n_high_state,:,i_irk_point_local)
          qp_energy_tmp(n_low_state:n_high_state,:) = qp_energy(n_low_state:n_high_state,:,i_irk_point_local)
        endif

        if(i_task.gt.0) then
          if(myid.eq.i_task) then
            call send_real_vector(exact_x_tmp,n_states_count*n_spin, 0)
            call send_real_vector(xc_tmp,n_states_count*n_spin, 0)
            call send_real_vector(correl_energy_tmp,n_states_count*n_spin, 0)
            call send_real_vector(qp_energy_tmp,n_states_count*n_spin, 0)
          elseif(myid.eq.0) then
            call receive_real_vector(exact_x_tmp,n_states_count*n_spin, i_task)
            call receive_real_vector(xc_tmp,n_states_count*n_spin, i_task)
            call receive_real_vector(correl_energy_tmp,n_states_count*n_spin, i_task)
            call receive_real_vector(qp_energy_tmp,n_states_count*n_spin, i_task)
          endif
        endif

!        if(i_task_1.gt.0) then
!          if(myid.eq.i_task_1) then
!            call send_real_vector(qp_energy_tmp,n_states_count*n_spin, 0)
!          elseif(myid.eq.0) then
!            call receive_real_vector(qp_energy_tmp,n_states_count*n_spin, i_task_1)
!          endif
!        endif

        if(myid.eq.0) then
          do i_spin = 1, n_spin
            if(n_spin.eq.2.and.i_spin.eq.1) then
             write(use_unit,'(35X, A)') "Spin Up"
            endif

           if(n_spin.eq.2.and.i_spin.eq.2) then
             write(use_unit,'(35X, A)') "Spin Down"
           endif

           write(use_unit,'(2X,A,I4,A,3f16.4)')"K_point ", i_k_point, " : ",  k_point_list(i_k_point,1:3)
           write(use_unit,'(2A)')"-----------------------------------------", &
           "----------------------------------------------------------"
            do i_state = n_low_state, n_high_state, 1
               i_state_1 = i_state - n_low_state + 1
               write(use_unit, '(2X, I6, 2X, F8.4, 5F14.4)') &
                   i_state, occ_numbers(i_state,i_spin,i_k_point), &
                   KS_eigenvalue(i_state,i_spin,i_k_point)*hartree, &
                   exact_x_tmp(i_state_1,i_spin)*hartree, &
                   xc_tmp(i_state_1,i_spin)*hartree, &
                   correl_energy_tmp(i_state_1,i_spin)*hartree, &
                   qp_energy_tmp(i_state,i_spin)*hartree
            enddo

 
            if(n_spin.eq.2.and.i_spin.eq.1) then
             write(use_unit,'(2A)')"-----------------------------------------", &
              "-------------------------------------------------------"
            endif
! end of i_spin
           enddo
           write(use_unit,'(2A)')"------------------------------------------", &
            "---------------------------------------------------------"
           write(use_unit,*)
! end  if (myid.eq. 0)
        endif
! end of loop i_k_point
      enddo

      allocate(qp_vbm_k(n_k_points,n_spin),stat=info)
      call check_allocation(info, 'qp_vbm',func)
      allocate(qp_cbm_k(n_k_points,n_spin),stat=info)
      call check_allocation(info, 'qp_cbm',func)

      do i_spin = 1, n_spin
        qp_vbm_k(:, i_spin) = 0.d0
        qp_cbm_k(:, i_spin) = 0.d0
        do i_k_point =1, n_k_points, 1

           i_irk_point = irk_point_mapping(i_k_point)

           if(myid .ne. mod(i_irk_point,n_tasks) ) cycle
           i_k_point_local = (i_k_point-1)/n_tasks+1
           i_irk_point_local = (i_irk_point-1)/n_tasks+1

           qp_vbm_k(i_k_point, i_spin) = -1.e9
           qp_cbm_k(i_k_point, i_spin) = 1.e9

           do i_state = n_low_state, n_high_state, 1
              if(occ_numbers(i_state,i_spin,i_k_point) .gt. 1.d-6) then
                 qp_vbm_k(i_k_point,i_spin) = qp_energy(i_state,i_spin,i_irk_point_local)
              endif
           enddo

           do i_state = n_high_state, n_low_state, -1
              if(occ_numbers(i_state,i_spin,i_k_point) .lt. 1.d0) then
                 qp_cbm_k(i_k_point,i_spin) = qp_energy(i_state,i_spin,i_irk_point_local)
              endif
           enddo
!          write(use_unit,*)i_k_point, qp_vbm_k(i_k_point,i_spin), qp_cbm_k(i_k_point,i_spin)
          
!  end of loop over i_k_point
        enddo
        call sync_vector(qp_vbm_k(1,i_spin),n_k_points)
        call sync_vector(qp_cbm_k(1,i_spin),n_k_points)

!          write(use_unit,*)"now",  qp_vbm_k(:,i_spin), qp_cbm_k(:,i_spin)
        qp_vbm(i_spin) = qp_vbm_k(1, i_spin)
        qp_cbm(i_spin) = qp_cbm_k(1, i_spin)
        direct_band_gap(i_spin) = qp_cbm(i_spin) - qp_vbm(i_spin)
        vbm_k_point(i_spin) = 1
        cbm_k_point(i_spin) = 1
        direct_gap_k_point(i_spin) = 1
!        write(use_unit,*) "qp_vbm, qp_cbm :", qp_vbm, qp_cbm
        do i_k_point = 2, n_k_points, 1
           if (qp_vbm(i_spin) .lt. qp_vbm_k(i_k_point, i_spin)) then
             qp_vbm(i_spin) = qp_vbm_k(i_k_point, i_spin) 
             vbm_k_point(i_spin) = i_k_point
           endif
           if (qp_cbm(i_spin) .gt. qp_cbm_k(i_k_point, i_spin)) then
             qp_cbm(i_spin) = qp_cbm_k(i_k_point, i_spin) 
             cbm_k_point(i_spin) = i_k_point
           endif

           if( direct_band_gap(i_spin) .gt. &
               (qp_cbm_k(i_k_point, i_spin) - qp_vbm_k(i_k_point, i_spin)) )  then
                direct_band_gap(i_spin) = qp_cbm_k(i_k_point, i_spin) - qp_vbm_k(i_k_point, i_spin)
                direct_gap_k_point(i_spin) = i_k_point
           endif
        enddo
!  end of loop over i_spin
      enddo

      call mpi_barrier(mpi_comm_global,mpierr)
      if(n_spin .eq. 1) then
          write(info_str,'(2X,A,f15.8,A)') &
             "Valence band maximum (VBM) from the GW calculation is ", &
             qp_vbm(1)*hartree, " eV (relative to internal zero)"
          call localorb_info(info_str)
          write(info_str,'(10X,A,3f8.4)') "at the k point ", k_point_list(vbm_k_point(1),1:3)
          call localorb_info(info_str)
          if(myid.eq.0) write(use_unit,*) 

          write(info_str,'(2X,A,f15.8,A)') &
             "Conduction band mininum (CBM) from the GW calculation is ", &
             qp_cbm(1)*hartree,  " eV (relative to internal zero)" 
          call localorb_info(info_str)
          write(info_str,'(10X,A,3f8.4)') "at the k point ", k_point_list(cbm_k_point,1:3)
          call localorb_info(info_str)
          if(myid.eq.0) write(use_unit,*) 

          write(info_str,'(2X,A,F15.8,A, 3f8.4)') &
             "Direct band gap from the GW calculation is ", &
             direct_band_gap(1)*hartree,  " eV at the k point ", k_point_list(direct_gap_k_point,1:3)  
          call localorb_info(info_str)
          if(myid.eq.0) write(use_unit,*) 
          write(info_str,'(2X,A,F15.8,A)') &
             "Indirect band gap from the GW calculation is ", &
             (qp_cbm(1)-qp_vbm(1))*hartree,  " eV " 
          call localorb_info(info_str)
      else
          write(info_str,'(2X,A,f15.8,A)') &
          "Spin-up valence band maximum (VBM) from the GW calculation is ", &
           qp_vbm(1)*hartree, " eV (relative to internal zero)"
          call localorb_info(info_str)
          write(info_str,'(10X,A,3f8.4)') "at the k point ", k_point_list(vbm_k_point(1),1:3)
          call localorb_info(info_str)
          if(myid.eq.0) write(use_unit,*) 

          write(info_str,'(2X,A,f15.8,A)') &
             "Spin-up conduction band mininum (CBM) from the GW calculation is ", &
             qp_cbm(1)*hartree,  " eV (relative to internal zero)" 
          call localorb_info(info_str)
          write(info_str,'(10X,A,3f8.4)') "at the k point ", k_point_list(cbm_k_point(1),1:3)
          call localorb_info(info_str)
          if(myid.eq.0) write(use_unit,*) 

          write(info_str,'(2X,A,F15.8,A, 3f8.4)') &
             "Spin-up direct band gap from the GW calculation is ", &
             direct_band_gap(1)*hartree,  " eV at the k point ", k_point_list(direct_gap_k_point(1),1:3)  
          call localorb_info(info_str)
          if(myid.eq.0) write(use_unit,*) 
          write(info_str,'(2X,A,F15.8,A)') &
             "Spin-up indirect band gap from the GW calculation is ", &
             (qp_cbm(1)-qp_vbm(1))*hartree,  " eV " 
          call localorb_info(info_str)
          if(myid.eq.0) write(use_unit,*) 

          write(info_str,'(2X,A,f15.8,A)') &
           "Spin-Down valence band maximum (VBM) from the GW calculation is ", &
            qp_vbm(2)*hartree, " eV (relative to internal zero)"
          call localorb_info(info_str)
          write(info_str,'(10X,A,3f8.4)') "at the k point ", k_point_list(vbm_k_point(2),1:3)
          call localorb_info(info_str)
          write(use_unit,*) 

          write(info_str,'(2X,A,f15.8,A)') &
             "Spin-down conduction band mininum (CBM) from the GW calculation is ", &
             qp_cbm(1)*hartree,  " eV (relative to internal zero)" 
          call localorb_info(info_str)
          write(info_str,'(10X,A,3f8.4)') "at the k point ", k_point_list(cbm_k_point(2),1:3)
          call localorb_info(info_str)
          if(myid.eq.0) write(use_unit,*) 

          write(info_str,'(2X,A,F15.8,A, 3f8.4)') &
             "Spin-down direct band gap from the GW calculation is ", &
             direct_band_gap(2)*hartree,  " eV at the k point ", k_point_list(direct_gap_k_point(2),1:3)  
          call localorb_info(info_str)
          if(myid.eq.0) write(use_unit,*) 

          write(info_str,'(2X,A,F15.8,A)') &
             "Spin-down indirect band gap from the GW calculation is ", &
             (qp_cbm(2)-qp_vbm(2))*hartree,  " eV " 
          call localorb_info(info_str)
      endif


      

!      call output_real_eigenfunctions(qp_energy,KS_eigenvector,occ_numbers)
      deallocate (qp_energy_tmp)
      deallocate (correl_energy)
      deallocate (qp_vbm_k)
      deallocate (qp_cbm_k)
      deallocate (exact_x_tmp)
      deallocate (xc_tmp)
      deallocate (correl_energy_tmp)

      return
      end subroutine quasi_particle_energy_p0


!---------------------------------------------------------------------
!******
