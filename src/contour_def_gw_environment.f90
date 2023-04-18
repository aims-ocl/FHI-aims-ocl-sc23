!****h* FHI-aims/contour_def_gw
!  NAME
!    Routines for contour deformation in single-shot GW
!  SYNOPSIS

module contour_def_gw_environment
   use constants,                  only: pi,&
                                         hartree
   use contour_def_gw_types,       only: cd_environment_type,&
                                         qp_adjust_type
   use dimensions,                 only: use_ev_scgw,&
                                         use_ev_scgw0
   use localorb_io,                only: localorb_info,&
                                         use_unit
   use mpi_tasks,                  only: myid,&
                                         aims_stop,&
                                         aims_warn
   use output_handling,            only: open_file, close_file
   use synchronize_mpi,            only: sync_timing,&
                                         bcast_real,&
                                         bcast_integer,&
                                         bcast_logical

   implicit none

   private

   public :: basic_init_contour_def_env, init_contour_def_env, init_contour_sc_loop,&
             update_contour_sc_loop, write_sc_info, contour_def_integrate,&
             contour_def_add_residues, sync_contour_def_timings, print_restart_information,&
             read_restart_information, adjust_KS_eigenvalues, adjust_qp_energies 
 
contains
 
! **************************************************************************************************
!> brief initializes and checks the environment for the contour deformation --iterative solution
!  o gw_cd  --  contour deformation environment
!  o n_low_state  -- lowest KS/HF eigenstate for self-energy correction
!  o n_high_state -- highest KS/HF eigenstate for self-energy correction
!  o n_spin -- number of spins
! **************************************************************************************************
   subroutine basic_init_contour_def_env(gw_cd,n_low_state,n_high_state,n_states,&
                                         n_spin,setup_qp_iteration,n_homo)
      
      type(cd_environment_type)                          :: gw_cd
      integer, intent(in)                                :: n_low_state, n_high_state
      integer, intent(in)                                :: n_states
      integer, intent(in)                                :: n_spin
      logical, intent(in)                                :: setup_qp_iteration
      integer, dimension(:), intent(in), optional        :: n_homo
 
      integer                                            :: i_level
   
      if(setup_qp_iteration.and..not.(present(n_homo))) then
        call aims_stop('CD initialization failed')
      endif

      !*** check CD range
      if(any(gw_cd%contour_def_start < n_low_state)) then
         gw_cd%contour_def_start(:) = n_low_state
      endif
 
      if(any(gw_cd%contour_def_end > n_high_state)) then
         gw_cd%contour_def_end(:) = n_high_state
      endif

      if(.not.allocated(gw_cd%num_levels)) then
         allocate(gw_cd%num_levels(n_spin))
      endif
   
      !*** set-up QP iteration
      if(setup_qp_iteration) then 
        call init_qp_iteration_env(gw_cd,n_spin,n_states,n_homo)
      else
        gw_cd%num_levels(:) = n_high_state
        if(.not.allocated(gw_cd%corrected_levels)) then
           allocate(gw_cd%corrected_levels(n_high_state,n_spin))
        endif
        do i_level = 1, n_high_state
           gw_cd%corrected_levels(i_level,1:n_spin) = i_level
        enddo
        if(.not.allocated(gw_cd%index_cd_levels)) then
          allocate(gw_cd%index_cd_levels(n_states,n_spin))
        endif
        gw_cd%index_cd_levels(:,:) = 0
        do i_level = 1, n_states
           gw_cd%index_cd_levels(i_level,1:n_spin) = i_level
        enddo
      endif
      gw_cd%num_levels_tot=MAXVAL(gw_cd%num_levels) 
   
      !*** allocate self-energy structures 
      if(.not.allocated(gw_cd%self_energy)) then
        allocate(gw_cd%self_energy)
        if(.not.allocated(gw_cd%self_energy%re)) then
         allocate (gw_cd%self_energy%re(gw_cd%num_levels_tot,n_spin)) 
        endif
        if(.not.allocated(gw_cd%self_energy%im)) then
         allocate (gw_cd%self_energy%im(gw_cd%num_levels_tot,n_spin)) 
        endif
        if(.not.allocated(gw_cd%self_energy%complx)) then
         allocate (gw_cd%self_energy%complx(gw_cd%num_levels_tot,n_spin)) 
        endif
      endif
      gw_cd%self_energy%re(:,:) = 0.0d0
      gw_cd%self_energy%im(:,:) = 0.0d0
      gw_cd%self_energy%complx(:,:) = (0.0d0, 0.0d0)

      !***allocate array for Hedin-shifts
      if(.not.allocated(gw_cd%shift)) then
         allocate(gw_cd%shift(gw_cd%num_levels_tot,n_spin))
      endif
      gw_cd%shift(:,:) = 0.0d0

   end subroutine basic_init_contour_def_env

! **************************************************************************************************
!> brief initialize QP iteration
!  o gw_cd  --  contour deformation environment
!  o n_states -- number of states (occupied and virtual)
!  o n_spin -- number of spins
!  o n_homo -- HOMO index for each spin
! **************************************************************************************************
   subroutine init_qp_iteration_env(gw_cd,n_spin,n_states,n_homo)

      type(cd_environment_type)                          :: gw_cd
      integer, intent(in)                                :: n_spin
      integer, intent(in)                                :: n_states
      integer, dimension(:), intent(in)                  :: n_homo

      integer                                            :: i_level, i_spin, diff
      integer                                            :: level_index, tmp_nr_level
      integer                                            :: start_occ_level, tmp_total, max_level
      integer                                            :: n_occ, n_virt, max_virt


      !***initialize self-consistent environment
      if(.not.allocated(gw_cd%sc_env)) then
        allocate(gw_cd%sc_env(n_spin))
      endif

      do i_spin = 1, n_spin
         gw_cd%num_levels(i_spin) =  gw_cd%contour_def_end(i_spin) &
                                     - gw_cd%contour_def_start(i_spin) + 1
      enddo
    
      !***set-up num_levels for sc-evGW0 or sc-evGW 
      if(gw_cd%self_consistent) then
        do i_spin = 1, n_spin
          n_occ = gw_cd%sc_env(i_spin)%n_occ
          n_virt = gw_cd%sc_env(i_spin)%n_virt
          if(n_homo(i_spin)-n_occ+1.le.0) then
            start_occ_level = 1
          else
            start_occ_level = n_homo(i_spin) - n_occ + 1
          endif
          if(gw_cd%contour_def_start(i_spin) < start_occ_level) then
            level_index = gw_cd%contour_def_start(i_spin)
          else
            level_index = start_occ_level
          endif
          max_level = MAX(n_homo(i_spin)+n_virt,gw_cd%contour_def_end(i_spin))
          tmp_total = 0  
          do i_level = 1, max_level 
             tmp_total =  tmp_total + 1
             if(level_index+1 > gw_cd%contour_def_end(i_spin).and.&
                level_index+1 < start_occ_level) then
                level_index = start_occ_level
             elseif(level_index+1 > n_homo(i_spin)+n_virt.and.&
                    level_index+1 < gw_cd%contour_def_start(i_spin)) then
                level_index = gw_cd%contour_def_start(i_spin)
             else
                level_index = level_index + 1 
             endif
             if(level_index > max_level) exit
          enddo
          gw_cd%num_levels(i_spin) =  tmp_total
        enddo
      endif
 
      !*** set-up array that contains indices of corrected levels
      tmp_nr_level = MAXVAL(gw_cd%num_levels) 
      if(.not.allocated(gw_cd%corrected_levels)) then
         allocate(gw_cd%corrected_levels(tmp_nr_level,n_spin))
      endif

      if(gw_cd%self_consistent) then
        do i_spin = 1, n_spin
          n_occ = gw_cd%sc_env(i_spin)%n_occ
          if(n_homo(i_spin)-n_occ+1.le.0) then
            start_occ_level = 1
          else
            start_occ_level = n_homo(i_spin) - n_occ + 1
          endif
          if(gw_cd%contour_def_start(i_spin) < start_occ_level) then
            level_index = gw_cd%contour_def_start(i_spin)
          else
            level_index = start_occ_level
          endif
          do i_level = 1, gw_cd%num_levels(i_spin)
             gw_cd%corrected_levels(i_level, i_spin)  = level_index    
             if(level_index+1 > gw_cd%contour_def_end(i_spin).and.&
                level_index+1 < start_occ_level) then
                level_index = start_occ_level
             elseif(level_index+1 > n_homo(i_spin)+n_virt.and.&
                    level_index+1 < gw_cd%contour_def_start(i_spin)) then
                level_index = gw_cd%contour_def_start(i_spin)
             else
                level_index = level_index + 1 
             endif 
          enddo
        enddo
      else
        do i_spin = 1, n_spin
           level_index = gw_cd%contour_def_start(i_spin)
           do i_level =  1, gw_cd%num_levels(i_spin)
              gw_cd%corrected_levels(i_level,i_spin) = level_index
              level_index = level_index + 1
           enddo
        enddo
      endif

      !*** set-up index array 
      if(.not.allocated(gw_cd%index_cd_levels)) then
         allocate(gw_cd%index_cd_levels(n_states,n_spin))
      endif
      gw_cd%index_cd_levels(:,:) = 0

      do i_spin = 1, n_spin
         do i_level = 1, gw_cd%num_levels(i_spin)
            level_index = gw_cd%corrected_levels(i_level,i_spin)
            gw_cd%index_cd_levels(level_index,i_spin) = i_level          
         enddo
      enddo

      !***initialize sc qp convergence
      if(allocated(gw_cd%sc_env)) then
        do i_spin = 1, n_spin
            allocate(gw_cd%sc_env(i_spin)%qp_converged(gw_cd%num_levels(i_spin)))
            allocate(gw_cd%sc_env(i_spin)%converged_last_step(gw_cd%num_levels(i_spin)))
            gw_cd%sc_env(i_spin)%qp_converged(:) = .false.
            gw_cd%sc_env(i_spin)%converged_last_step(:) = .false.
            gw_cd%sc_env(i_spin)%last_step = .false.
        enddo
      endif
      
   end subroutine init_qp_iteration_env

! **************************************************************************************************
!> brief initializes and checks the environment for the contour deformation solving the quasi-
!>       particle equations iteratively
!  o gw_cd  --  contour deformation environment
!  o real_omega --  real frequency, i.e. we have Sigma_istate(real_omega)
!  o KS_eigenvalue -- KS/HF eigenvalues of the single-particle calculation
!  o chemical_potential -- the chemical potential of the system
!  o i_state -- i-th KS/HF state 
!  o i_spin -- ith spin 
!  o n_low_state  -- lowest KS/HF eigenstate for self-energy correction
!  o n_high_state -- highest KS/HF eigenstate for self-energy correction
!  o n_states -- number of KS/HF eigenstates
!  o n_spin -- number of spins
!  o n_homo -- the HOMO level, i.e., the number of occupied state 
! **************************************************************************************************
   subroutine init_contour_def_env(gw_cd, real_omega, KS_eigenvalue, chemical_potential, i_state, &
                                   i_spin, n_states, n_homo)
      
      type(cd_environment_type)                          :: gw_cd
      real(kind=8), intent(in)                           :: real_omega
      real(kind=8), dimension(:,:), &
        intent(in)                                       :: KS_eigenvalue
      real(kind=8), dimension(:), intent(in)             :: chemical_potential
      integer, intent(in)                                :: i_state, i_spin
      integer, intent(in)                                :: n_states
      integer, dimension(:), intent(in)                  :: n_homo
 
      integer                                            :: i_ener, j_state, &
                                                            num_residues, index_contour_def
 
      index_contour_def  = gw_cd%index_cd_levels(i_state,i_spin)
      gw_cd%self_energy%re(index_contour_def,i_spin) = 0.0d0
      gw_cd%self_energy%im(index_contour_def,i_spin) = 0.0d0
      gw_cd%self_energy%complx(index_contour_def,i_spin) = 0.0d0
 
      !*** determine how many residues we need
 
      num_residues = 0
 
      if(real_omega < chemical_potential(i_spin)) then
        do j_state = 1, n_homo(i_spin)
           if(real_omega < KS_eigenvalue(j_state,i_spin)) then
              num_residues = num_residues + 1
           endif
        enddo 
      else
        do j_state = n_homo(i_spin)+1, n_states
           if(real_omega > KS_eigenvalue(j_state,i_spin)) then
              num_residues = num_residues + 1
           endif
        enddo 
      endif
 
      gw_cd%num_residues = num_residues

      ! allocate real frequencies of the residues
      allocate(gw_cd%real_freq(num_residues))
      gw_cd%real_freq =0.0d0
      !
      !! for getting the residue (connected with the MO m)
      allocate(gw_cd%residue_from_freq(num_residues))
 
      ! 
      i_ener = 0
      if(real_omega < chemical_potential(i_spin)) then
        do j_state = 1, n_homo(i_spin)
           if(real_omega < KS_eigenvalue(j_state,i_spin)) then
              i_ener = i_ener+1
              gw_cd%residue_from_freq(i_ener) = j_state
              gw_cd%real_freq(i_ener) = KS_eigenvalue(j_state,i_spin)-real_omega 
           endif
        enddo  
      else
        do j_state = n_homo(i_spin)+1, n_states
           if(real_omega > KS_eigenvalue(j_state,i_spin)) then
              i_ener = i_ener+1
              gw_cd%residue_from_freq(i_ener) = j_state
              gw_cd%real_freq(i_ener) = KS_eigenvalue(j_state,i_spin)-real_omega 
           endif
        enddo  
      endif   

   end subroutine init_contour_def_env

! **************************************************************************************************
!> brief initialize self-consistent GW0 loop
!  o gw_cd  --  contour deformation environment
!  o n_iteration -- maximal number of iterations in sc-ev loops
!  o threshold_sc -- threshold for the sc-evGW0 or sc-evGW convergence; note that it can't be too
!                    tight because the CD would be no longer stable
!  o scgw_converged -- if sc-ev loop is converged
!  o n_spin -- number of spins
! **************************************************************************************************
   subroutine init_contour_sc_loop(gw_cd,n_iteration,threshold_sc,scgw_converged,n_spin)
 
      type(cd_environment_type)                          :: gw_cd
      integer, intent(out)                               :: n_iteration
      real(kind=8), intent(out)                          :: threshold_sc   
      logical, intent(out)                               :: scgw_converged
      integer, intent(in)                                :: n_spin

      integer                                            :: itemp

      itemp =1 

      if(gw_cd%self_consistent) then
        n_iteration = gw_cd%n_iter_sc
      else
        n_iteration = 1
      endif
   
      scgw_converged = .false.
      threshold_sc = 0.005d0/hartree

      if(use_ev_scgw) then
        gw_cd%sctype="ev-scGW"
      elseif(use_ev_scgw0) then
        gw_cd%sctype="ev-scGW0"
      endif

      if (myid.eq.0 .and. n_iteration .gt. 1) then
         write(use_unit,'(2X,A)') TRIM(gw_cd%sctype) // " loop initialization ... "
         write(use_unit,*) " "
         write(use_unit,'(T2,A)')"************************************************************"
         write(use_unit, '(T3,A19,T25,I3)') TRIM(gw_cd%sctype) // " Iteration:", itemp
         write(use_unit,'(T2,A)')"************************************************************"
      endif

   end subroutine init_contour_sc_loop

! **************************************************************************************************
!> brief update sc-evGW0 or sc-evGW loops
!  o gw_cd  --  contour deformation environment
!  o scgw_converged -- if sc-ev loop is converged
!  o i_iter -- iteration step in sc-evGW0 or sc-evGW loops
!  o n_iteration -- maximal number of iterations in sc-ev loops
!  o n_spin -- number of spins
!  o n_states -- number of KS/HF eigenstates
!  o n_homo -- the HOMO level, i.e., the number of occupied state 
!  o threshold_sc -- threshold for the sc-evGW0 or sc-evGW convergence
!  o qp_energy -- QP energies within the sc-ev cycle calculate in step i_iter
!  o KS_eigenvalue_last -- in 1st loop corresponds to the KS/HF eigenvalues
! **************************************************************************************************
   subroutine update_contour_sc_loop(gw_cd,scgw_converged,i_iter,n_iteration,n_spin,n_states,n_homo,&
                                     threshold_sc,qp_energy,KS_eigenvalue_last)

      type(cd_environment_type)                          :: gw_cd
      logical, intent(out)                               :: scgw_converged
      integer, intent(in)                                :: i_iter
      integer, intent(in)                                :: n_iteration
      integer, intent(in)                                :: n_spin
      integer, intent(in)                                :: n_states
      integer, dimension(:), intent(in)                  :: n_homo
      real(kind=8), intent(in)                           :: threshold_sc
      real(kind=8), dimension(:,:), intent(inout)        :: qp_energy
      real(kind=8), dimension(:,:), intent(inout)        :: KS_eigenvalue_last

      integer                                            :: i_spin, i_state, i_level
      integer                                            :: index_contour_def,&
                                                            n_explicit_virt,&
                                                            n_explicit_occ
      logical                                            :: tmp_converged
      real(kind=8)                                       :: diff

      tmp_converged = .true.

      if(gw_cd%self_consistent) then

        do i_spin = 1, n_spin
           do i_level = 1, gw_cd%num_levels(i_spin)
              i_state = gw_cd%corrected_levels(i_level,i_spin)
              diff = (ABS(qp_energy(i_state,i_spin)-KS_eigenvalue_last(i_state,i_spin)))
              gw_cd%sc_env(i_spin)%converged_last_step(i_level) &
                   = gw_cd%sc_env(i_spin)%qp_converged(i_level)
              if(diff < threshold_sc) then
                 gw_cd%sc_env(i_spin)%qp_converged(i_level) = .true.
              endif
           enddo
           !***scissor shift for virtuals
           diff = 0.0d0
           n_explicit_virt = 0
           do i_state = n_homo(i_spin)+1, n_states
              index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
              if(index_contour_def == 0) cycle
              diff = diff + qp_energy(i_state,i_spin)-KS_eigenvalue_last(i_state,i_spin)
              n_explicit_virt = n_explicit_virt + 1
           enddo
           diff = diff/REAL(n_explicit_virt,KIND=8)
           if(n_explicit_virt < n_states-n_homo(i_spin)) then
             do i_state = n_homo(i_spin)+1, n_states
                index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
                if(index_contour_def /= 0) cycle
                qp_energy(i_state,i_spin) = KS_eigenvalue_last(i_state,i_spin) + diff
             enddo
             if(myid.eq.0) then
               write(use_unit,'(T3,A42,T49,F14.4,T64,A2)') "Scissor-shift the rest of " //&
                    "the virtuals by:", diff*hartree, "eV" 
             endif
           endif
           !*** scissor shift for occupied
           diff = 0.0d0
           n_explicit_occ = 0
           do i_state = 1, n_homo(i_spin)
              index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
              if(index_contour_def == 0) cycle
              !do not include core states in the shift
              if(KS_eigenvalue_last(i_state,i_spin)*hartree < -100.0d0) cycle
              diff = diff + qp_energy(i_state,i_spin)-KS_eigenvalue_last(i_state,i_spin)
              n_explicit_occ = n_explicit_occ + 1
           enddo
           diff = diff/REAL(n_explicit_occ,KIND=8)
           if(n_explicit_occ < n_homo(i_spin)) then
             do i_state = 1, n_homo(i_spin)
                index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
                if(index_contour_def /= 0) cycle
                if(KS_eigenvalue_last(i_state,i_spin)*hartree < -100.0d0) then
                  !account for core states
                  qp_energy(i_state,i_spin) = KS_eigenvalue_last(i_state,i_spin) + 10.0d0*diff
                else
                  qp_energy(i_state,i_spin) = KS_eigenvalue_last(i_state,i_spin) + diff
                endif
             enddo
           endif
        enddo
        call initiate_reiteration(gw_cd,n_spin)
        tmp_converged = .true.
        do i_spin = 1, n_spin
           if(.not.all(gw_cd%sc_env(i_spin)%qp_converged)) then 
             tmp_converged = .false.
           endif
        enddo
        if(tmp_converged) then 
          scgw_converged = .true.
        else
           KS_eigenvalue_last(:,:) = qp_energy(1:n_states,1:n_spin)
        endif

      else
        scgw_converged = .true.

      endif

   end subroutine update_contour_sc_loop
 
! **************************************************************************************************
!> brief initialize re-iteration of state of interest when ev-scGW or ev-scGW0 loop is converged.
!        Background is that CD is unstabile when QP-result is very close to any incoming
!        "new" KS-value; thus once converged, these states leave the ev-scGW or ev-scGW0 loop, but
!        might be re-iterated once all other states are converged
!  o gw_cd  --  contour deformation environment
! **************************************************************************************************
  subroutine initiate_reiteration(gw_cd,n_spin)

      type(cd_environment_type)                          :: gw_cd
      integer, intent(in)                                :: n_spin

      integer                                            :: i_state, i_spin,&
                                                            index_contour_def

      do i_spin = 1, n_spin
         if(gw_cd%sc_env(i_spin)%reiterate) then
            if(all(gw_cd%sc_env(i_spin)%qp_converged).and..not.gw_cd%sc_env(i_spin)%last_step) then
               gw_cd%sc_env(i_spin)%last_step = .true.
               do i_state = gw_cd%contour_def_start(i_spin),gw_cd%contour_def_end(i_spin)
                  index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
                  if(gw_cd%sc_env(i_spin)%converged_last_step(index_contour_def).and.&
                     gw_cd%sc_env(i_spin)%qp_converged(index_contour_def)) then
                    gw_cd%sc_env(i_spin)%qp_converged(index_contour_def) = .false.   
                    if(myid.eq.0) write(use_unit,'(T2,A16,I4,2X,A8,I4)') "Re-iterate state", i_state,&
                                  "for spin", i_spin
                  endif
               enddo
            else
               if(gw_cd%sc_env(i_spin)%last_step) then
                 gw_cd%sc_env(i_spin)%qp_converged(:) = .true.
               endif
            endif
         endif
      enddo

  end subroutine initiate_reiteration
! **************************************************************************************************
!> brief initialize self-consistent GW0 loop
!  o gw_cd  --  contour deformation environment
!  o n_iteration -- maximal number of iterations in sc-ev loops
!  o i_iter -- iteration step in sc-evGW0 or sc-evGW loops
! **************************************************************************************************
   subroutine write_sc_info(gw_cd,n_iteration,i_iter,scgw_converged)

      type(cd_environment_type)                          :: gw_cd
      integer, intent(in)                                :: n_iteration
      integer, intent(in)                                :: i_iter
      logical, intent(in)                                :: scgw_converged

      if(.not.scgw_converged) then 
        if(myid == 0.and.i_iter+1 <= n_iteration) then
          write(use_unit,'(T2,A)')"************************************************************"
          write(use_unit, '(T3,A19,T25,I3)') TRIM(gw_cd%sctype) // " Iteration:", i_iter+1
          write(use_unit,'(T2,A)')"************************************************************"  
        endif
      endif

   end subroutine write_sc_info

! **************************************************************************************************
!> brief numerical integration for the contour deformation (iterative solution)
!  o gw_cd  --  contour deformation environment
!  o Wmn_element -- screened Coulomb matrix element in the MO basis, i.e. (mn|W(iomega)|mn)
!  o real_omega -- real frequency
!  o KS_eigenvalue -- KS/HF eigenvalues of the single-particle calculation
!  o i_spin -- ith spin 
!  o i_state -- i-th KS/HF state
!  o j_state -- j-th KS/HF state
!  o i_freq -- i-th frequency for screened Coulomb interaction W or real frequency 
!  o omega_full -- the Gauss-Legendre frequency grid for the screened Coulomb interaction
!  o womega(n_freq) -- the weight of the Gauss-Legendre frequency grid for the self-energy
! **************************************************************************************************
   subroutine contour_def_integrate(gw_cd, Wmn_element,real_omega, KS_eigenvalue, i_spin, i_state,&
                                    j_state, i_freq, omega_full, womega_full, chemical_potential)
      
      type(cd_environment_type)                          :: gw_cd
      complex(kind=8), intent(in)                        :: Wmn_element
      real(kind=8), intent(in)                           :: real_omega
      real(kind=8), dimension(:,:), intent(in)           :: KS_eigenvalue
      integer, intent(in)                                :: i_spin, i_state, j_state, &
                                                            i_freq
      real(kind=8), dimension(:), intent(in)             :: omega_full, womega_full    
      real(kind=8), intent(in)                           :: chemical_potential   

      integer                                            :: index_contour_def
      real(kind=8)                                       :: omega
         
       index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
       if(gw_cd%corrected_levels(index_contour_def,i_spin) == i_state) then
         omega = omega_full(i_freq)  
         call numerical_integrate_cd(gw_cd%self_energy%complx(index_contour_def,i_spin),&
                                     Wmn_element, womega_full(i_freq), omega,&
                                     real_omega, KS_eigenvalue(j_state,i_spin),&
                                     gw_cd%full_cmplx_sigma, gw_cd%eta,&
                                     chemical_potential)
      endif

   end subroutine contour_def_integrate

! **************************************************************************************************
!> brief adding residues for contour deformation (iterative solution)
!  o gw_cd  --  contour deformation environment
!  o Wmn_element -- screened Coulomb matrix element in the MO basis, i.e. (mn|W(e_m-omega)|mn)
!  o real_omega -- real frequency omega, i.e. we have Sigma_istate(real_omega)
!  o chemical_potential -- the chemical potential of the system
!  o i_spin -- ith spin 
!  o i_state -- i-th KS/HF state
!  o j_state -- j-th KS/HF state
!  o i_freq -- i-th real frequency, i.e. (epsilon_m-epsilon_n) 
! **************************************************************************************************
   subroutine contour_def_add_residues(gw_cd, self_energy, Wmn_element, real_omega, chemical_potential,&
                                       i_spin, i_state, j_state, i_freq_real)

      type(cd_environment_type)                          :: gw_cd
      complex(kind=8), dimension(:,:), intent(inout)     :: self_energy
      complex(kind=8), intent(in)                        :: Wmn_element
      real(kind=8), intent(in)                           :: chemical_potential
      real(kind=8), intent(in)                           :: real_omega
      integer, intent(in)                                :: i_spin, i_state, j_state,&
                                                            i_freq_real

      integer                                            :: index_contour_def, &
                                                            m_level_of_real_energy
      real(kind=8)                                       :: omega
     
 
      omega = gw_cd%real_freq(i_freq_real) ! epsilon_m - real_omega
      m_level_of_real_energy = gw_cd%residue_from_freq(i_freq_real) 
      ! only add residue to specific level
      if(m_level_of_real_energy == j_state) then
         index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
         if(chemical_potential-real_omega < omega ) then
            self_energy(index_contour_def,i_spin) = &
                   self_energy(index_contour_def,i_spin) + Wmn_element
         elseif(chemical_potential-real_omega  > omega) then
            self_energy(index_contour_def,i_spin) = &
                   self_energy(index_contour_def,i_spin) - Wmn_element
         endif
      endif
 
   end subroutine contour_def_add_residues

! **************************************************************************************************
!> brief numerical integration for the contour deformation technique 
!  o Wmn_element -- screened Coulomb matrix element in the MO basis, i.e. (mn|W(iomega)|mn)
!  o weight -- integration weight for the respective grid
!  o omega -- frequency
!  o real_energy -- real energy
!  o eigen_val -- KS/HF eigenstate 
! **************************************************************************************************

   subroutine numerical_integrate_cd(self_energy_real_update, Wmn_element, weight, omega, real_energy,&
                                     eigen_val, full_cmplx_sigma, eta, chemical_potential)

      complex(kind=8), intent(inout)                      :: self_energy_real_update
      complex(kind=8), intent(in)                         :: Wmn_element
      real(kind=8), intent(in)                            :: weight, omega, real_energy,&
                                                             eigen_val
      logical, intent(in), optional                       :: full_cmplx_sigma              
      real(kind=8), intent(in), optional                  :: eta, chemical_potential


      real(kind=8)                                        :: my_eta, sgn
      complex(kind=8)                                     :: im_unit


      im_unit = (0.0d0, 1.0d0)

      my_eta = 0.0d0
      sgn = 1.0d0

      if(present(full_cmplx_sigma).and.present(eta)) then
        if(full_cmplx_sigma) then
          my_eta = eta
          sgn = sign(1.0d0,chemical_potential - eigen_val)
        endif
      endif
      
      self_energy_real_update = self_energy_real_update- &
                          ( &
                          0.5d0/pi*weight*Wmn_element/(im_unit*(omega-my_eta*sgn)+real_energy-eigen_val)+ &
                          0.5d0/pi*weight*Wmn_element/(im_unit*(-omega-my_eta*sgn)+real_energy-eigen_val) &
                          )

   end subroutine numerical_integrate_cd

! **************************************************************************************************
!> brief sync cpu times
! **************************************************************************************************
   subroutine sync_contour_def_timings(time_self_energy_cd, time_Wmn_realfreq_cd, &
                                       time_WPQ_realfreq_cd, time_polar_realfreq_cd)

      real(kind=8), intent(inout)                         :: time_self_energy_cd,&
                                                             time_Wmn_realfreq_cd,&
                                                             time_WPQ_realfreq_cd,&
                                                             time_polar_realfreq_cd
   
  
      call sync_timing(time_self_energy_cd)
      call sync_timing(time_Wmn_realfreq_cd)
      call sync_timing(time_WPQ_realfreq_cd)
      call sync_timing(time_polar_realfreq_cd)

   end subroutine sync_contour_def_timings

! **************************************************************************************************
!> print restart information to separate file
! **************************************************************************************************
   subroutine print_restart_information(gw_cd, qp_energy, qp_non_convergence, e_diff, i_count,&
                                        current_state, current_spin) 

     type(cd_environment_type)                            :: gw_cd
     real(kind=8), dimension(:,:), intent(in)             :: qp_energy
     logical, dimension(:,:), intent(in)                  :: qp_non_convergence
     real(kind=8), intent(in)                             :: e_diff
     integer, intent(in)                                  :: i_count
     integer, intent(in)                                  :: current_state
     integer, intent(in)                                  :: current_spin

     character(len=40)                                    :: filename 
     integer                                              :: i_state, i_spin,&
                                                             i_level, unit_number,&
                                                             my_i_count
     real(kind=8)                                         :: my_e_diff
 

     if(trim(adjustl(gw_cd%restart))=='write'.or.trim(adjustl(gw_cd%restart))=='read_and_write') then
       if(myid == 0) then
          call open_file(filename,&
                         file_action = 'write',&
                         front_name='contour_gw_qp_energies',&
                         extension='dat',&
                         unit_number=unit_number)
          do i_spin = 1, current_spin
             do i_level = 1, gw_cd%num_levels(i_spin)
                i_state = gw_cd%corrected_levels(i_level,i_spin)
                if(i_state > current_state.and.i_spin == current_spin) cycle
                if(i_state==current_state.and.i_spin==current_spin) then
                  my_e_diff = e_diff
                  my_i_count = i_count
                else
                  my_e_diff = 0.0d0
                  my_i_count = 0
                  if(qp_non_convergence(i_state,i_spin)) my_i_count = -1
                endif 
                write(unit_number,'(2I6,2F18.12,I6)') i_spin, i_state, qp_energy(i_state,i_spin),&
                                                      my_e_diff, my_i_count
                                                    
             enddo
          enddo
          call close_file(unit_number)
       endif 
     endif 
     
 
   end subroutine print_restart_information

! **************************************************************************************************
!> read restart information from file
! **************************************************************************************************
   subroutine read_restart_information(gw_cd, qp_energy, qp_energy_old,qp_non_convergence,&
                                       e_diff, i_count, i_spin, i_state, already_converged)

     type(cd_environment_type)                            :: gw_cd
     real(kind=8), dimension(:,:), intent(inout)          :: qp_energy
     real(kind=8), dimension(:,:), intent(inout)          :: qp_energy_old
     logical, dimension(:,:), intent(inout)               :: qp_non_convergence
     real(kind=8), intent(inout)                          :: e_diff
     integer, intent(inout)                               :: i_count
     integer, intent(in)                                  :: i_state, i_spin
     logical, intent(out)                                 :: already_converged 
   
     character(len=40)                                    :: filename
     logical                                              :: open_failed
     integer                                              :: unit_number, istat,&
                                                             tmp_i_count, tmp_state,&
                                                             tmp_spin
     integer                                              :: spin_to_restart
     integer                                              :: state_to_restart
     real(kind=8)                                         :: tmp, tmp_e_diff


     already_converged = .false.
     open_failed = .false.

     if(trim(adjustl(gw_cd%restart))=='read'.or.trim(adjustl(gw_cd%restart))=='read_and_write') then
       if(myid == 0) then
          
          call open_file(filename,&
                         front_name='contour_gw_qp_energies',&
                         extension='dat',&
                         file_action='read',&
                         unit_number=unit_number,&
                         open_failed=open_failed)
          if(.not.open_failed) then
            do
              read(unit_number,'(2I6,2F18.12,I6)',iostat=istat) tmp_spin, tmp_state, tmp, &
              tmp_e_diff, tmp_i_count
              if(tmp_spin == i_spin.and.tmp_state == i_state) then
                qp_energy(i_state,i_spin) = tmp
                e_diff = tmp_e_diff
                i_count = tmp_i_count
              endif
              if(istat /= 0) then
                spin_to_restart  = tmp_spin
                state_to_restart = tmp_state
                exit
              endif
            enddo
                                              
            call close_file(unit_number)

          endif
       endif
       call bcast_logical(open_failed,0)
       if(open_failed) then  
         if(myid.eq.0) call aims_warn("* Warning! Restart file does not exist. Not restarting")
       else
         call bcast_real(qp_energy(i_state,i_spin),0)
         call bcast_real(e_diff,0)
         call bcast_integer(i_count,0)
         call bcast_integer(spin_to_restart,0)
         call bcast_integer(state_to_restart,0)
         qp_energy_old(i_state,i_spin) = qp_energy(i_state,i_spin) - e_diff
         if(i_spin < spin_to_restart) already_converged = .true.
         if(i_state < state_to_restart .and. i_spin == spin_to_restart) already_converged = .true.
         if(i_count < 0) qp_non_convergence(i_state,i_spin) = .true.
       endif
     endif 
     
 
   end subroutine read_restart_information

! **************************************************************************************************
!> brief add relativistic correction to KS eigenvalues before GW calculation starts
!  o qp_adjust  --  environment for relativistic correction
!  o KS_eigenvalue -- KS/HF eigenvalues of the single-particle calculation
!  o n_states -- number of states
!  o n_spin -- number of spins 
! **************************************************************************************************
   subroutine adjust_KS_eigenvalues(qp_adjust, KS_eigenvalue, n_states, n_spin)

     type(qp_adjust_type)                                        :: qp_adjust
     real(kind=8), dimension(n_states,n_spin), intent(inout)     :: KS_eigenvalue
     integer, intent(in)                                         :: n_spin, n_states
     character(len=120)                                          :: info_str

     integer                                                     :: i, i_spin


     call localorb_info('',use_unit,'(A)')

     do i = 1 , qp_adjust%n_qps
        do i_spin = 1, n_spin
           write (info_str,'(2X,A,I2,A,I2,A,F14.8,A,F14.8)') &
                "Adjusting KS eigenvalue ", qp_adjust%state(i), " in spin channel  ",&
                i_spin, " to ",KS_eigenvalue(qp_adjust%state(i), i_spin)*hartree,' + ',&
                qp_adjust%values(i)
           call localorb_info(info_str,use_unit,'(A)')
           KS_eigenvalue(qp_adjust%state(i), i_spin) &
                = KS_eigenvalue(qp_adjust%state(i), i_spin) + qp_adjust%values(i)/hartree
        end do
     end do
     
     call localorb_info('',use_unit,'(A)')
     
     
   end subroutine adjust_KS_eigenvalues

! **************************************************************************************************
!> brief add relativistic correction to QP energies
!  o gw_cd  --  contour deformation environment
!  o qp_adjust  --  environment for relativistic correction
!  o qp_energy -- quasi-particle energies
!  o n_states -- number of states
!  o n_spin -- number of spins 
! **************************************************************************************************
   subroutine adjust_qp_energies(gw_cd, qp_adjust, qp_energy, n_states, n_spin)

     type(cd_environment_type)                                   :: gw_cd
     type(qp_adjust_type)                                        :: qp_adjust
     real(kind=8), dimension(:,:), intent(inout)                 :: qp_energy
     integer, intent(in)                                         :: n_spin, n_states

     character(len=120)                                          :: info_str
     integer                                                     :: i, i_spin


     call localorb_info('',use_unit,'(A)')

     do i = 1 , qp_adjust%n_qps
        do i_spin = 1, n_spin
           if(i_spin/=gw_cd%spin_channel(i_spin)) cycle
           if(qp_adjust%state(i) >= gw_cd%contour_def_start(i_spin) .and.&
              qp_adjust%state(i) <= gw_cd%contour_def_end(i_spin)) then
              write (info_str,'(2X,A,I2,A,I2,A,F14.8,A,F14.8)') &
                   "Adjusting QP energy ", qp_adjust%state(i), " in spin channel  ",&
                   i_spin, " to ",qp_energy(qp_adjust%state(i), i_spin)*hartree,' + ',&
                   qp_adjust%values(i)
              call localorb_info(info_str,use_unit,'(A)')
              qp_energy(qp_adjust%state(i), i_spin) &
                   = qp_energy(qp_adjust%state(i), i_spin) + qp_adjust%values(i)/hartree
           endif
        end do
     end do
     
     call localorb_info('',use_unit,'(A)')
     
     
   end subroutine adjust_qp_energies

end module contour_def_gw_environment
