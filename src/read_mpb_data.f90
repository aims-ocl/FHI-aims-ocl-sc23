!****s* FHI-aims/read_mpb_data
!  NAME
!   read_mpb_data
!  AUTHOR
!   Stefan Ringe
!  REFERENCE
!   This routine is part of the MPBE solver for the modeling of 
!   electrolytes with FHI-aims. Large parts of this code are still 
!   experimental, that is why we highly recommend you to contact
!   the main authors of the corresponding paper (below) before starting
!   to calculate anything. We are highly motivated to help and 
!   cooperate with interested FHI-aims users.
!  SEE ALSO
!    Ringe, S., Oberhofer, H., Hille, C., Matera, S., Reuter, K., "Function-Space oriented Solution Scheme for the 
!    Size-Modified Poisson-Boltzmann Equation in Full-Potential DFT", JCTC 12, 4052-4066 (2016).
!    DOI: 10.1021/acs.jctc.6b00435
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2016).
!  SYNOPSIS

subroutine read_mpb_data()

!  PURPOSE
!
!
!  Subroutine read_mpb_data reads in all necessary parameters for
!  the MPBE solver
! 
!  VARIABLES
!
!  For a short description of the MPBE variables check out the modules
!  lpb_solver_utilities and mpb_solver_utilities which define all
!  additional global variables used by the MPBE solver

   use dimensions
   use lpb_solver_utilities
   use mpb_solver_utilities

   use dimensions
   use grids
   use species_data
   use mpi_tasks
   use localorb_io
   use constants
   use read_fixed_grid
   use runtime_choices
   use vdw_correction
   use pseudodata

  implicit none

  integer :: linecount, i_code
  character*500 inputline
  character*100 desc_str
  character*20 solvent_type
  logical :: flag_dielec_func, flag_ions_temp, flag_ions_charge,&
    flag_ions_size, flag_ions_kind, flag_ions_conc, &
    flag_SPE_lmax,flag_constant_lpb_solver,&
    flag_SPE_conv, flag_dynamic_cavity, flag_start_mpb_solver, &
    flag_delta_rho_in_merm,flag_converge_rho_mpb,&
    flag_SPE_cut_and_lmax_ff,flag_ions_mod_alpha,flag_solve_lpbe_only,&
    flag_MERM_in_SPE_solver, flag_reg_method, flag_correct_dielectric_decrement,&
    flag_set_nonelstat_params, flag_dynamic_ions, flag_ions_mod_alpha_anion,flag_forces_mpb_force_off,&
    flag_nonsc_Gnonmf, flag_Gnonmf_FD_delta, flag_MERM_atom_wise
  character*300 :: info_str
  real*8 :: shift_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   //First catch some Aims calculations that do not work, yet with MPBE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (use_forces .and. atomic_MERM .and. use_density_matrix) then
    call aims_stop('MPBE: Forces are under current development and not working for atomic_MERM together with the density matrix update method. Please use either orbital based update method or set atomic_MERM to False.')
  end if
  if (n_periodic.ne.0) then
    call aims_stop('MPBE: MPBE is not working for periodic systems, yet, only for cluster systems without PBCs')    
  end if
  if (use_embedding_potential) then
    call aims_stop('MPBE: MPBE is not working, yet, with embedding.')
  end  if
  if (use_embedding_pp) then
    call aims_stop('MPBE: We did not test, if MPBE is actually working with pseudo-potentials.&
      &In principle it', 'should be possible to make it work. If you nead pseudo-potentials, please contact the authors')
  end if
  if (use_load_balancing) then
    call aims_stop('MPBE: MPBE does not work, yet, with load balancing')
  end if
  if (use_local_index) then
    call aims_stop('MPBE: stopping here, because we did not actually test if MPBE is working with&
      &use_local_index if one needs this please try it out by yourself and remove this stop and contact the authors about the results')
  end if
  if (.not.force_new_functional) then
    call aims_stop('MPBE: stopping here, we are not sure, if the MPBE solver really works with force_new_functional = .false. has to be tested&
      &please contact the authors.')
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   //flags of things to read in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  flag_dielec_func = .False.
  flag_ions_temp = .False.
  flag_ions_charge = .false.
  flag_ions_size = .false.
  flag_ions_kind = .false.
  flag_ions_conc = .false.
  flag_SPE_lmax = .False.
  flag_constant_lpb_solver = .False.
  flag_SPE_conv = .False.
  flag_dynamic_cavity = .False.
  flag_start_mpb_solver = .False.
  flag_delta_rho_in_merm = .False.
  flag_converge_rho_mpb = .False.
  flag_SPE_cut_and_lmax_ff = .False.
  flag_ions_mod_alpha = .False.
  flag_solve_lpbe_only = .False.
  flag_reg_method = .False.
  flag_MERM_in_SPE_solver = .False.
  flag_correct_dielectric_decrement = .False.
  flag_set_nonelstat_params = .False.
  flag_dynamic_ions  = .False.
  flag_ions_mod_alpha_anion = .False.
  flag_forces_mpb_force_off = .False.
  flag_nonsc_Gnonmf = .False.
  flag_Gnonmf_FD_delta = .False.
  flag_MERM_atom_wise = .False.

   if (myid.eq.0) then
      write(use_unit,*)
      write (use_unit,'(2X,A)') &
            "| Reading configuration options for the MPB solver "
   end if
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  //Start reading
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
  linecount = 0
   
!    radial_multiplier_mpb = 10
!    r_radial_mpb_min = 2.2
!    r_radial_mpb_max = 3.0   
   
!    call modify_grid_mpb()
!     !add grid shells at dielectric funcion transition
!     n_radial = n_radial + 30

   lineloop: do

      read(7,'(A)',iostat=i_code) inputline
      if(i_code<0) then
         backspace(7)
        exit lineloop        ! end of file reached
      end if
      if(i_code>0) then
         write (use_unit,'(2X,A)') &
               "Error reading file 'control.in'..."
         stop
      end if

      linecount = linecount + 1

      read(inputline,*,iostat=i_code) desc_str
      if(i_code/=0) cycle lineloop             ! skip empty line

      if (desc_str(1:1).eq.'#') cycle lineloop ! skip comment

!       //SYSTEM SPECIFICATIONS
      
      if (desc_str.eq."dielec_func") then
         read(inputline,*,end=88,err=99) desc_str, dieleckind
	 flag_dielec_func = .True.

	 if (dieleckind.eq.0) then
	    read(inputline,*,end=88,err=99) desc_str, desc_str, epsinf_mpb, rho0_eps, beta_eps
	 else if (dieleckind.eq.1) then
	    read(inputline,*,end=88,err=99) desc_str, desc_str, epsinf_mpb, rhomin_mpb, rhomax_mpb
	 end if

!       //IONS

      else  if (desc_str.eq."ions_temp") then
	 read(inputline,*,end=88,err=99) desc_str,T_mpb
	 flag_ions_temp = .True. 

      else  if (desc_str.eq."ions_conc") then
	 read(inputline,*,end=88,err=99) desc_str,c_mpb
	 flag_ions_conc = .True. 

      else  if (desc_str.eq."ions_charge") then
	 read(inputline,*,end=88,err=99) desc_str,z_mpb
	 flag_ions_charge = .True. 

      else  if (desc_str.eq."ions_size") then
	 read(inputline,*,end=88,err=99) desc_str,a_mpb
	 flag_ions_size = .True. 

      else  if (desc_str.eq."ions_kind") then
	 read(inputline,*,end=88,err=99) desc_str,alphakind
	 flag_ions_kind = .True. 
	 
      else  if (desc_str.eq."ions_mod_alpha") then
	  read(inputline,*,end=88,err=99) desc_str, d_alpha_ion,xi_alpha_ion
	  flag_ions_mod_alpha = .True. 

      else  if (desc_str.eq."ions_mod_alpha_anion") then
	  read(inputline,*,end=88,err=99) desc_str, d_alpha_ion_anion,xi_alpha_ion_anion
	  use_separate_alpha = .True.
	  flag_ions_mod_alpha_anion = .True.  		 

!       //SPE SOLVER
	 
      else  if (desc_str.eq."SPE_lmax") then
	  read(inputline,*,end=88,err=99) desc_str, limit_l_max_SPE
	  flag_SPE_lmax = .True.
      else  if (desc_str.eq."SPE_conv") then
	  read(inputline,*,end=88,err=99) desc_str, tau_merm, eta_merm
	  flag_SPE_conv = .True.  
	  
      else  if (desc_str.eq."SPE_cut_and_lmax_ff") then
	  read(inputline,*,end=88,err=99) desc_str, multipole_radius_SPE_sq, limit_l_max_SPE_ff
	  flag_SPE_cut_and_lmax_ff = .True.  
	  
	  if (multipole_radius_SPE_sq.eq.0.0) then
	    !using standard settings with free atom radius
	    changed_SPE_ff_cutoffs = .False.
	  else
	    !using our own multipole radius given in the above statement
	    changed_SPE_ff_cutoffs = .True.
	  end if
	  
	  
!       //OTHER TECHNICAL SETTINGS
	  
      else  if (desc_str.eq."start_mpb_solver") then
	  read(inputline,*,end=88,err=99) desc_str, mpb_start_crit
	  if (mpb_start_crit.eq.'step') then
	    read(inputline,*,end=88,err=99) desc_str, desc_str, start_at_scf_step
	  else
	    start_at_scf_step = 0
	  end if
	  flag_start_mpb_solver = .True. 
	  
      else  if (desc_str.eq."dynamic_cavity_off") then
	  dynamic_cavity = .false.
	  flag_dynamic_cavity = .True.  

      else  if (desc_str.eq."dynamic_ions_off") then
	  dynamic_ions = .false.
	  flag_dynamic_ions = .True.  

      else  if (desc_str.eq."nonsc_Gnonmf") then
	  evaluate_sc_Gnonmf = .false.
	  flag_nonsc_Gnonmf = .True.  

      else  if (desc_str.eq."Gnonmf_FD_delta") then
	  read(inputline,*,end=88,err=99) desc_str, delta_Gnonmf
	  flag_Gnonmf_FD_delta = .True.  

      else  if (desc_str.eq."delta_rho_in_merm") then
	  delta_rho_in_merm = .true.
	  flag_delta_rho_in_merm = .True.  

      else  if (desc_str.eq."not_converge_rho_mpb") then
	  not_converge_rho_mpb = .True.
	  flag_converge_rho_mpb = .True.  



      else if (desc_str.eq."lpbconstcoeff") then
	  flag_constant_lpb_solver = .True.
	  solve_lpb_with_constant_dielec = .True.
	  write(use_unit,'(2X,A)') '| WARNING! Setting lpbconstcoeff is just for testing purposes'
	  
      else  if (desc_str.eq."solve_lpbe_only") then
	  read(inputline,*,end=88,err=99) desc_str, use_mpbe_free_energy
	  solve_lpbe_only = .True.
	  flag_solve_lpbe_only = .True.  

      else  if (desc_str.eq."reg_method") then
	  read(inputline,*,end=88,err=99) desc_str, reg_method
	  flag_reg_method = .True.  

      else  if (desc_str.eq."MERM_in_SPE_solver") then
	  read(inputline,*,end=88,err=99) desc_str, MERM_in_SPE_solver
	  flag_MERM_in_SPE_solver = .True.  
	  
      else  if (desc_str.eq."MERM_atom_wise") then
	  read(inputline,*,end=88,err=99) desc_str, atomic_MERM
	  flag_MERM_atom_wise = .True.  

      else  if (desc_str.eq."correct_dielectric_decrement") then
      read(inputline,*,end=88,err=99) desc_str, dielectric_decrement_method
	  correct_dielectric_decrement = .True.
	  flag_correct_dielectric_decrement = .True.  

      else  if (desc_str.eq."set_nonelstat_params") then
	   read(inputline,*,end=88,err=99) desc_str, alphaplusgamma_mpb, beta_mpb
	   flag_set_nonelstat_params = .True.

      else  if (desc_str.eq."forces_mpb_force_off") then
           forces_mpb_force_off = .true.
           flag_forces_mpb_force_off = .True.
	   
      else
!         must have reached the end of the current species description.

         backspace(7)
         exit lineloop

      end if

   enddo lineloop


  if (linecount == 0) then
    if (myid.eq.0) then
	write(use_unit,*) "No data for MPB solver. Assuming defaults."
    end if
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   //Setting up default values, if nothing was given as INPUT.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if (.not.flag_set_nonelstat_params) then
    !these are the defaults of the paper from Andreussi and Marzari for neutral molecules. for charged take different params!
    alphaplusgamma_mpb = 50
    beta_mpb = -0.35
  end if

  if (.not.flag_correct_dielectric_decrement) then
    correct_dielectric_decrement = .False.
  end if

  if (.not.flag_nonsc_Gnonmf) then
    evaluate_sc_Gnonmf = .True.
  end if

  if (.not.flag_MERM_in_SPE_solver) then
    !Default: do not update potential in each LPBE step, but run
    !loop inside kerker_mpb routine
    MERM_in_SPE_solver = .True.
  end if
  
  if (.not.flag_MERM_atom_wise) then
    atomic_MERM = .False.
  end if
  
  if (.not.flag_solve_lpbe_only) then
    !Default: use MPBE solver, not pure LPBE solver
    solve_lpbe_only = .False.
    use_mpbe_free_energy = .False.
  end if

  if (.not.flag_reg_method) then
    !Default: regularize MPBE solver with free Hartree potential, not with LPBE solution
    reg_method = 'none' !'vacandconst'
  end if

  if (.not. flag_ions_mod_alpha) then
    d_alpha_ion = 0.5 !-0.07
    xi_alpha_ion = 1.0 !0.4
  end if

  if (.not. flag_Gnonmf_FD_delta) then
    delta_Gnonmf = 1e-8
  end if
  
  if (.not. flag_ions_mod_alpha_anion) then
    use_separate_alpha = .false.
  end if

  if (.not. flag_dielec_func) then
    dieleckind = 1
    epsinf_mpb = 78.36
    rho0_eps = 0.0004
    beta_eps = 1.3
     rhomin_mpb = 0.0001
     rhomax_mpb = 0.005
!      rhomin_mpb = 0.0001
!      rhomax_mpb = 0.0015
  end if
  log_rhomin_mpb = log(rhomin_mpb)
  log_rhomax_mpb = log(rhomax_mpb)
  log_epsinf_mpb = log(epsinf_mpb)

  if (.not. flag_converge_rho_mpb) then
    !not converge rho with solvent potential. if true the scf will be converged in vacuum until conv crit is hit and 
    !after the last step the solvation potential will be calculated (so will not effect rho)
    not_converge_rho_mpb = .False.
  end if

  if (.not.flag_ions_charge) then
    z_mpb = 1
  end if
  if (.not.flag_ions_conc) then
    c_mpb = 1.0
  end if
  if (.not.flag_ions_kind) then
    alphakind = 1
  end if
  if (.not.flag_ions_size) then
    a_mpb = 5.0
  end if
  if (.not.flag_ions_temp) then
    T_mpb = 300.
  end if

  if (.not. flag_dynamic_cavity) then
    dynamic_cavity = .True.
  end if

  if (.not. flag_dynamic_ions) then
    dynamic_ions = .True.
  end if

  if (.not. flag_SPE_lmax) then
    !precondition with ions is default
    limit_l_max_SPE = 0 !maxval(l_hartree) will be automatically set to maxval(l_hartree) later
  end if
  precondition_with_ions = .True.


  if (.not. flag_constant_lpb_solver) then
    solve_lpb_with_constant_dielec = .False.
  end if

  if (.not. flag_SPE_conv) then
    tau_merm = 1.0d-10 !1.0d-8
    eta_merm = .5
  end if

  if (.not.flag_start_mpb_solver) then
    mpb_start_crit = 'step'
    start_at_scf_step = 2
  end if

  if (.not.flag_SPE_cut_and_lmax_ff) then
    limit_l_max_SPE_ff = limit_l_max_SPE
    changed_SPE_ff_cutoffs = .False.
    !cutoff radius will be set to multipole_radius_free + 
    !extra_adding_to_hartree_potential_distance (2Ang)
  end if



  if (.not.flag_delta_rho_in_merm) then
    delta_rho_in_merm = .False.
  else
    write(info_str,'(2X,2A)') '| WARNING: Delta-SPE source term update',&
    'is not supported, yet, therefore turning it off.'
    call localorb_info(info_str) 
    delta_rho_in_merm = .False.
  end if
  
  if (.not.flag_forces_mpb_force_off) then
    !the default setting does not turn off the force calculation for MPB solvation
    forces_mpb_force_off = .False.
  end if
  
  
  !settings for solvation energies
    !solve_lpbe_only = .true. !false.
    !MERM_in_SPE_solver = .true.
    !c_mpb = 0.0
    !d_alpha_ion = 0.5
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    //FINISHED setting defaults. Now modify the input for FHI-aims (units,...) and output
!    warning if any INPUT was not correct or will not work
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    multipole_radius_SPE_sq = multipole_radius_SPE_sq**2
    
   if (evaluate_sc_Gnonmf) then
    write(info_str,'(2X,A)') '| Using self-consistent calculation of dispersion/Pauli-'
    call localorb_info(info_str)  
    write(info_str,'(4X,A)') 'repulsion/cavity creation free energy. Could lead to worse convergence of SCF and'
    call localorb_info(info_str) 
    write(info_str,'(4X,A)') 'is not always needed, see Ringe,Oberhofer,Hille,Matera,Reuter, JCTC, 2016.'
    call localorb_info(info_str) 
   end if
    !per default: correct KS Hamiltonian, i.e. consider deps/dnel and dalpha/dnel != 0
    KS_mpb_corr1=.true.
    KS_mpb_corr2=.true.

    !in the case of cavities that do not depend on nel but on nelfree, these terms are zero
    if (.not. dynamic_cavity) then
      KS_mpb_corr1=.false.
    end if
    if (.not. dynamic_ions) then
      KS_mpb_corr2=.false.
    end if
  
  if (solve_lpbe_only .and. use_separate_alpha) then
	call aims_stop('MPBE: LPBE does not work together with use_separate_alpha, choose either MPBE or only 1 alpha function!') 
  end if 
  if (MERM_in_SPE_solver .and. use_separate_alpha) then
    call aims_stop('MPBE: MERM_in_SPE_solver does not work with use_separate_alpha at the moment')
  end if
  if (use_separate_alpha .and. correct_dielectric_decrement) then
	call aims_stop('MPBE: Dielectric decrement does not work with two separate ion exclusion functions at the moment!') 
  end if  
  if (correct_dielectric_decrement .and. (decrement_kind == 2 .or. decrement_kind == 3)) then
    call aims_stop('MPBE: Dielectric decrement is only working with linear eps(c_ion) function at the moment, use decrement_kind = 1 or 4 instead!') 
  end if
  if (.not. MERM_in_SPE_solver .and. atomic_MERM) then
    write(info_str,'(2X,A)') '| WARNING: If "atomic_MERM = True" works only, if MERM_in_SPE_solver is used, setting atomic_MERM to False.'
   call localorb_info(info_str)
    atomic_MERM = .False.
  end if
  
  if (not_converge_rho_mpb) then
   write(info_str,'(2X,A)') '| WARNING: If "not_converge_rho_mpb = True" no forces will be'
   call localorb_info(info_str)
   write(info_str,'(4X,A)') 'calculated. The SCF is converged in vacuum and the MPBE potential corresponding'
   call localorb_info(info_str)
   write(info_str,'(4X,A)') 'to this electron density is evaluated once.'
   call localorb_info(info_str)
    mpb_start_crit = 'converge'
  end if  
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    //Writing out the parameters to be used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  write(info_str,'(2X,A,A,E14.6,A)') '| The temperature of the pseudoatomic',&
    'charges is set to:'
  call localorb_info(info_str)
  write(info_str,'(4X,A,F15.8,2X,A)') '| T_mpb = ', T_mpb, 'K'
  call localorb_info(info_str)
  write(info_str,'(2X,A,A,E14.6,A)') '| The following parameters are used in',&
    'the MPB method to model a smooth dielectric function:'
  call localorb_info(info_str)
  write(info_str,'(4X,A,F15.8,2X,A)') '| eps_infinity = ', epsinf_mpb
  call localorb_info(info_str)
  if (dieleckind == 0) then
    write(info_str,'(4X,A,F15.8,2X,A)') '| rho_0 = ', rho0_eps, 'e/bohr**3'	
    call localorb_info(info_str)
    write(info_str,'(4X,A,F15.8)') 'beta = ', beta_eps
    call localorb_info(info_str)
  else
    write(info_str,'(4X,A,F15.8,2X,A)') '| n_min = ', rhomin_mpb, 'e/bohr**3'	
    call localorb_info(info_str)
    write(info_str,'(4X,A,F15.8,2X,A)') '| n_max = ', rhomax_mpb, 'e/bohr**3'	
    call localorb_info(info_str)
  end if
  write(info_str,'(2X,A,A,E14.6,A)') '| A z:z electrolyt will be incorporated with'
  call localorb_info(info_str)
  write(info_str,'(4X,A,F15.8,2X,A)') '| concentration: c_mpb = ', c_mpb, 'mol/L' 
  call localorb_info(info_str)
  write(info_str,'(4X,A,I2)') '| charge: z_mpb = ', z_mpb 
  call localorb_info(info_str)
  write(info_str,'(4X,A,F15.8,2X,A)') '| size: a_mpb = ', a_mpb, 'Ang' 
  call localorb_info(info_str)
  if (alphakind==1) then
    write(info_str,'(4X,A,F15.8)') '| d_alpha_ion = ', d_alpha_ion
    call localorb_info(info_str)
    write(info_str,'(4X,A,F15.8)') '| xi_alpha_ion = ', xi_alpha_ion
    call localorb_info(info_str)
  end if
  write(info_str,'(2X,A,A,E14.6,A)') '| The linear PB equation is solved by multipole expansion up to'
  call localorb_info(info_str)
  write(info_str,'(4X,A,I2)') '| Using maximum angular momentum for SPE multipole expansion of l_max = ', limit_l_max_SPE 
  call localorb_info(info_str)
  write(info_str,'(4X,A,I2)') '| and in the far field l_max = ', limit_l_max_SPE_ff 
  call localorb_info(info_str)
  write(info_str,'(2X,A,E14.6)') '| tau_merm = ', tau_merm 
  call localorb_info(info_str)
  if (dynamic_cavity) then
    write(info_str,'(2X,A)') '| The dielectric function will be updated in each SCF step.' 
    call localorb_info(info_str)
  else
    write(info_str,'(2X,A)') '| The dielectric function will be built up by the use of the superposition of free atom electron densities.' 
    call localorb_info(info_str)
  end if
  if (dynamic_ions) then
    write(info_str,'(2X,A)') '| The exclusion function will be updated in each SCF step.' 
    call localorb_info(info_str)
  else
    write(info_str,'(2X,A)') '| The exclusion function will be built up by the use of the superposition of free atom electron densities.' 
    call localorb_info(info_str)
  end if
  if (KS_mpb_corr1) then
    write(info_str,'(2X,A)') '| The dielectric function dependence on the electron density is added to the Kohn-Sham Hamiltonian.'
    call localorb_info(info_str)
  else
    write(info_str,'(2X,A)') '| The dielectric function dependence on the electron density is ignored.' 
    call localorb_info(info_str)
  end if
  if (KS_mpb_corr2) then
    write(info_str,'(2X,A)') '| The exclusion function dependence on the electron density is added to the Kohn-Sham Hamiltonian.'
    call localorb_info(info_str)
  else
    write(info_str,'(2X,A)') '| The exclusion function dependence on the electron density is ignored.' 
    call localorb_info(info_str)
  end if
 

  
  if (forces_mpb_force_off) then
    write(info_str,'(2X,A,A,E14.6,A)') '| Forces (from MPB scheme) will NOT be computed.'
    call localorb_info(info_str)
  else if (.not.forces_mpb_force_off.and.use_forces) then
    write(info_str,'(2X,A,A,E14.6,A)') '| Forces (from MPB scheme) will be computed.'
    call localorb_info(info_str)  
  end if
  
  if (use_separate_alpha) then
    write(info_str,'(2X,A,A,E14.6,A)') '| Two separate alpha-functions will be used, one for the cations, one for the anions.'
    call localorb_info(info_str)
  end if

  !Now we check the parameters and modify them. we calculate also additional parameters...


  epsinf_mpb_inv = 1.0D0 / epsinf_mpb


  if (c_mpb == 0d0) then !we cannot calculate with c = 0 M, have to give it a low value 
    c_mpb = 1d-30	!(below c_cut in kerker_mpb)
    write(info_str,'(2X,2A,E14.6)') '| Giving concentration a finite value due to numerical ',&
    'issues, c_mpb = ',c_mpb 
    call localorb_info(info_str)
    mpbe_no_ion = .true. 
    if (.not.solve_lpbe_only) then
      write(info_str,'(2X,A)') '| Switching from MPBE to LPBE solver which is faster and more stable for c = 0M (no ion) case.'
      call localorb_info(info_str)
    end if
    solve_lpbe_only = .true.
  end if
  
  if (a_mpb.eq.0.0) then
    if (solve_lpbe_only) then
      solve_pbe_only = .false.
      solve_debye_only = .true.
    else
      solve_pbe_only = .true.	
      solve_debye_only = .false.
    end if
  else
    solve_pbe_only = .false.
  end if  
    
  if (solve_lpbe_only) then
    if (solve_debye_only) then
      write(info_str,'(2X,A)') '| Using Debye-LPB-DFT method.'
      call localorb_info(info_str)       
      write(info_str,'(2X,A)') '| Using iterative/relaxation SC-LPBE solver.'
      call localorb_info(info_str)    
    else
      write(info_str,'(2X,A)') '| Using LPB-DFT method.'
      call localorb_info(info_str)       
      write(info_str,'(2X,A)') '| Using iterative/relaxation SC-LPBE solver.'
      call localorb_info(info_str)
    end if
  else if (solve_pbe_only) then
    write(info_str,'(2X,A)') '| Using PB-DFT method.'
    call localorb_info(info_str)      
    write(info_str,'(2X,A)') '| Using Newton solver combined with iterative iterative/relaxation SC-LPBE solver.'
    call localorb_info(info_str)
  else
    write(info_str,'(2X,A)') '| Using MPB-DFT method.'  
    call localorb_info(info_str)      
    write(info_str,'(2X,A)') '| Using Newton solver combined with iterative iterative/relaxation SC-LPBE solver.'
    call localorb_info(info_str)  
  end if
  
  if (MERM_in_SPE_solver) then
    write(info_str,'(2X,A)') '| Using faster SC-LPBE solver which does not require full update of potential in each iteration step.'  
    call localorb_info(info_str)      
  end if

  if (correct_dielectric_decrement) then
    write(info_str,'(2X,A)') '| Applying Dielectric Decrement Method, making the dielectric function dependent on the ionic concentrations.'
    call localorb_info(info_str)  
  end if

  c_mpb = c_mpb * 6.02214129D23 * 1.0D-27 * bohr**3 !concentration in a_0**-3
  kappainf_mpb = SQRT(2d0 * pi4 * c_mpb * z_mpb) !calculate  parameter for MPB

  kappa_div_sqrt_eps = kappainf_mpb / SQRT(epsinf_mpb)


  a_mpb = 1/bohr * a_mpb 
  phi_zero_mpb = 2 * c_mpb * a_mpb**3

  if (phi_zero_mpb.ge.1d0) then
    write(info_str,'(2X,A,F10.5,2A)') '| WARNING: phi_zero_mpb is equal or larger than 1 (',phi_zero_mpb,'). This is physically unreasonable',&
    'and computationally instable.'
    call localorb_info(info_str)
    phi_zero_mpb = 0.99
    a_mpb = (phi_zero_mpb/2d0/c_mpb)**(1d0/3d0)
    write(info_str,'(2X,A,F15.8,A)') '| Changed a_mpb to a_mpb = ',a_mpb*bohr,' Ang to get phi_zero_mpb = 0.99.'
    call localorb_info(info_str)
  end if
  
  kBT_mpb = boltzmann_kB * T_mpb

    if (precondition_with_ions) then
      kappa_debye = kappa_div_sqrt_eps * sqrt(z_mpb/kBT_mpb)
    else
      kappa_debye = 5d0*kappa_div_sqrt_eps * sqrt(z_mpb/kBT_mpb)
    end if
  sqrt_z_mpb_div_kBT_mpb = sqrt(z_mpb/kBT_mpb)

  if (d_alpha_ion.ne.0d0.or.xi_alpha_ion.ne.1d0) then
    rhomin_kappa = log(rhomin_mpb) + (log(rhomin_mpb)-log(rhomax_mpb))*d_alpha_ion
    rhomax_kappa = log(rhomax_mpb) + (log(rhomin_mpb)-log(rhomax_mpb))*d_alpha_ion
    shift_dist = (rhomax_kappa-rhomin_kappa)*(1d0-xi_alpha_ion)/2d0
    rhomin_kappa = exp(rhomin_kappa + shift_dist)
    rhomax_kappa = exp(rhomax_kappa - shift_dist)
  else
    rhomin_kappa = rhomin_mpb !take care that this is NOT negative, and the function is NOT 100% symmetric
    rhomax_kappa = rhomax_mpb !if you change this, you have to correct the kappa gradient if correct_dielectric_d is true
  end if

  if (use_separate_alpha) then
    if (d_alpha_ion_anion.ne.0d0.or.xi_alpha_ion_anion.ne.1d0) then
      rhomin_kappa_anion = log(rhomin_mpb) + (log(rhomin_mpb)-log(rhomax_mpb))*d_alpha_ion_anion
      rhomax_kappa_anion = log(rhomax_mpb) + (log(rhomin_mpb)-log(rhomax_mpb))*d_alpha_ion_anion
      shift_dist = (rhomax_kappa_anion-rhomin_kappa_anion)*(1d0-xi_alpha_ion_anion)/2d0
      rhomin_kappa_anion = exp(rhomin_kappa_anion + shift_dist)
      rhomax_kappa_anion = exp(rhomax_kappa_anion - shift_dist)
    else
      rhomin_kappa_anion = rhomin_mpb !make sure by yourself, that this is NOT negative, and keep in mind that the function is NOT 100% symmetric
      rhomax_kappa_anion = rhomax_mpb !if you change this, you have to correct the kappa gradient if correct_dielectric_d is true
    end if 
  end if
  log_rhomin_kappa = log(rhomin_kappa)
  log_rhomax_kappa = log(rhomax_kappa)
  if (rhomin_kappa_anion > 0.0d0) log_rhomin_kappa_anion = log(rhomin_kappa_anion)
  if (rhomax_kappa_anion > 0.0d0) log_rhomax_kappa_anion = log(rhomax_kappa_anion)

  write(info_str,'(2X,A,F15.8,A)') '| volume fraction of ions, phi_zero_mpb = ', phi_zero_mpb, ' a.u.'
  call localorb_info(info_str)

  return

88 continue
   write(use_unit,'(A)') & 
   "Syntax error reading data for mpb solver from 'control.in' (missing arguments)"
   write(use_unit,'(A)') "in line: '"//trim(inputline)//"'"
   stop

99 continue
   write(use_unit,'(A)') "Syntax error reading data for mpb solver from 'control.in'"
   write(use_unit,'(A)') "in line: '"//trim(inputline)//"'"
   stop

end subroutine read_mpb_data
