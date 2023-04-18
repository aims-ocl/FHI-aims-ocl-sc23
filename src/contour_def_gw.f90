!****h* FHI-aims/contour_def_gw
!  NAME
!    Routines for contour deformation in single-shot GW
!  SYNOPSIS

module contour_def_gw
   use constants,                  only: hartree
   use contour_def_gw_environment, only: init_contour_sc_loop,&
                                         update_contour_sc_loop,&
                                         write_sc_info,&
                                         adjust_KS_eigenvalues,&
                                         adjust_qp_energies
   use contour_def_gw_types,       only: cd_environment_type
   use dimensions,                 only: flag_calc_spectral_func,&
                                         flag_print_self_energy,&
                                         flag_adjust_qp_pre,&
                                         flag_adjust_qp_post,&
                                         use_ev_scgw,&
                                         use_contour_def_gw,&
                                         n_spin
   use gw_para,                    only: state_to_print,& 
                                         print_range,&
                                         z_value
   use localorb_io,                only: localorb_info,&
                                         use_unit
   use mpi_tasks,                  only: myid,&
                                         aims_warn,&
                                         aims_stop
   use evaluate_self_energy_freq,  only: evaluate_Wmn_imag_freq_for_cd
   use quasi_particle_energies,    only: quasi_particle_energy_contour_def,&
                                         quasi_particle_energy_contour_def_zshot,&
                                         get_exchange_terms,&
                                         output_quasi_particle_energies
   use print_self_energies,        only: print_self_energy_analytic,&
                                         print_self_energy_cd
   use spectral_func_cd,           only: calc_spectral_func
 
   implicit none

   private

   public :: contour_def_gw_calculation

contains

! **************************************************************************************************
!> brief perform G0W0, ev-scGW0 or ev-scGW calculation with contour deformation 
!  o gw_cd  --  contour deformation environment
!  o qp_energy -- real array, the final calculated quasiparticle energy for each
!  o qp_non_convergence -- flag whether the qp calculation from the previous calculation with
!                          analytic continuation has converged 
!  o KS_eigenvalue -- real array,
!           the eigenvalues of the single-particle calculation. For DFT calculation,
!           this is the KS eigenvalue, but for HF calculation, this is then the HF
!           eigenvalue
!  o exchange_self_energy -- real array the exact exchange part of the self-energy 
!           of each KS/HF state for each spin channel
!  o xc_KS_matr -- real array, the matrix elements of the exchange correlation 
!           potential witin KS/HF orbitals
!  o x_KS_array -- real array, the matrix elements of the exchange
!           potential witin KS/HF orbitals
!  o c_KS_array -- real array, the matrix elements of the correlation
!           potential witin KS/HF orbitals
!  o ovlp_3KS -- transformed 3-center integrals
!  o Wmn_freq_cd_cmplx --  screened Coulomb matrix element in the MO basis, i.e. (mn|W(iomega)|mn)
!                          dim: (ndim2_o3KS,ndim1_o3KS,n_full_freq,n_spin)
!  o n_freq -- integer number, the number of frequency points for the GW self-energy
!  o n_full_freq -- the number of frequency points for the screened Coulomb interaction W
!  o omega -- the Gauss-Legendre frequency grid for the GW self-energy
!  o omega_full -- the Gauss-Legendre frequency grid for the screened Coulomb interaction
!  o womega_full -- the weight of the Gauss-Legendre frequency grid for the self-energy
!  o chemical_potential -- the chemical potential of the system
!  o n_low_state  -- the lowest KS/HF eigenstate for self-energy correction
!  o n_high_state  -- the highest KS/HF eigenstate for self-energy correction
!  o n_states -- number of KS/HF eigenstates
!  o n_homo -- the HOMO level for each spin channel
! **************************************************************************************************
   subroutine contour_def_gw_calculation(gw_cd,qp_energy,qp_non_convergence,KS_eigenvalue,&
                                         KS_eigenvector,occ_numbers,correl_energy,&
                                         exchange_self_energy,xc_KS_matr,x_KS_array,c_KS_array,&
                                         ovlp_3KS,Wmn_freq_cd_cmplx,n_freq,n_full_freq,omega_grid,&
                                         omega_full_grid,womega_full,chemical_potential,partition_tab,&
                                         ndim1_o3KS,ndim2_o3KS,n_low_state,n_high_state,n_states,&
                                         n_homo,tot_time_polar,tot_time_sigma,l_shell_max)

      type(cd_environment_type)                                :: gw_cd
      real(kind=8), dimension(:,:), intent(inout)              :: qp_energy
      logical, dimension(:,:), intent(inout)                   :: qp_non_convergence
      real(kind=8), dimension(:,:), intent(inout)              :: KS_eigenvalue
      real(kind=8), dimension(:,:,:), intent(in)               :: KS_eigenvector
      real(kind=8), dimension(:,:), intent(in)                 :: occ_numbers
      real(kind=8), dimension(:,:), intent(inout)              :: correl_energy
      real(kind=8), dimension(:,:), intent(in)                 :: exchange_self_energy
      real(kind=8), dimension(:,:,:), intent(in)               :: xc_KS_matr
      real(kind=8), dimension(:,:), intent(in)                 :: x_KS_array      
      real(kind=8), dimension(:,:), intent(in)                 :: c_KS_array      
      real(kind=8), dimension(:,:,:,:), allocatable,&
        intent(inout)                                          :: ovlp_3KS
      complex(kind=8), dimension(:,:,:,:), allocatable,&
        intent(inout)                                          :: Wmn_freq_cd_cmplx
      integer, intent(in)                                      :: n_freq, n_full_freq
      real(kind=8), dimension(:), intent(in)                   :: omega_grid
      real(kind=8), dimension(:), intent(in)                   :: omega_full_grid
      real(kind=8), dimension(:), intent(in)                   :: womega_full
      real(kind=8), dimension(:), intent(in)                   :: chemical_potential   
      real(kind=8), dimension(:), intent(in)                   :: partition_tab
      integer, intent(in)                                      :: ndim1_o3KS, ndim2_o3KS
      integer, intent(in)                                      :: n_low_state, n_high_state
      integer, intent(in)                                      :: n_states
      integer, dimension(:), intent(in)                        :: n_homo
      real(kind=8), intent(inout)                              :: tot_time_polar, tot_time_sigma
      integer, dimension(:), intent(in)                        :: l_shell_max


      logical                                                  :: flag_coh, out_self_energy
      logical                                                  :: scgw_converged 
      integer                                                  :: i_iter, n_iteration
      real(kind=8)                                             :: threshold_sc
      real(kind=8), dimension(:,:), allocatable                :: KS_eigenvalue_last
      real(kind=8), dimension(:,:), allocatable                :: qp_energy_ac
      real(kind=8), dimension(:,:), allocatable                :: my_exchange_self_energy
      real(kind=8), dimension(:,:), allocatable                :: xc_energy
      complex(kind=8), dimension(:,:,:), allocatable           :: self_energy_omega

      flag_coh = .false.
      out_self_energy = .false.
      allocate(KS_eigenvalue_last(n_states,n_spin))
      allocate(self_energy_omega(n_low_state:n_high_state,n_freq,n_spin))
      self_energy_omega(:,:,:) = (0.d0,0.d0)
 
      !*** add relativistic correction to KS eigenvalues 
      if(flag_adjust_qp_pre) then
         call adjust_KS_eigenvalues(gw_cd%qp_adjust, KS_eigenvalue, n_states, n_spin)
      endif

      if(.not.gw_cd%iterative) then
         KS_eigenvalue_last(:,:) = KS_eigenvalue(1:n_states,1:n_spin) 
         call quasi_particle_energy_contour_def_zshot(gw_cd, Wmn_freq_cd_cmplx,qp_energy,correl_energy, &
                                                  KS_eigenvalue,KS_eigenvalue_last,exchange_self_energy,&
                                                  xc_KS_matr,x_KS_array,c_KS_array,occ_numbers,&
                                                  n_low_state,n_high_state,n_homo,n_full_freq, &
                                                  omega_full_grid,womega_full,chemical_potential, &
                                                  ovlp_3KS)
         !**** add relativistic correction to QP energies
         if(flag_adjust_qp_post) then
            call adjust_qp_energies(gw_cd,gw_cd%qp_adjust,qp_energy,n_states,n_spin)
         endif
         !**** output QP energies
         call qp_output(gw_cd, qp_energy,KS_eigenvalue,occ_numbers,&
                        exchange_self_energy,xc_KS_matr,x_KS_array,c_KS_array,correl_energy,&
                        n_low_state,n_high_state,n_states,n_homo,z_value=z_value)
      endif

      allocate(qp_energy_ac(n_low_state:n_high_state,n_spin))
      qp_energy_ac(:,:) = qp_energy
      if(gw_cd%iterative) then
         KS_eigenvalue_last(:,:) = KS_eigenvalue(1:n_states,1:n_spin) 
         !*** initialize GW0 loop
         call init_contour_sc_loop(gw_cd,n_iteration,threshold_sc,scgw_converged,n_spin)
         !*** start GW0 lop
         do i_iter = 1, n_iteration
            if(use_ev_scgw.and.i_iter > 1) then
              Wmn_freq_cd_cmplx = (0.0d0, 0.0d0)
              call evaluate_Wmn_imag_freq_for_cd(gw_cd, Wmn_freq_cd_cmplx, self_energy_omega,&
                                                 n_low_state, n_high_state, ndim2_o3KS, &
                                                 ndim1_o3KS, n_spin, n_homo, occ_numbers,&
                                                 n_freq, n_full_freq, omega_grid, omega_full_grid,&
                                                 womega_full, chemical_potential,&
                                                 KS_eigenvalue_last, KS_eigenvector, ovlp_3KS,&
                                                 out_self_energy, partition_tab, l_shell_max,&
                                                 flag_coh, tot_time_polar, tot_time_sigma)
            endif
            call quasi_particle_energy_contour_def(gw_cd, Wmn_freq_cd_cmplx, qp_energy, correl_energy, &
                                                   KS_eigenvalue, KS_eigenvalue_last, &
                                                   exchange_self_energy, xc_KS_matr, x_KS_array, &
                                                   c_KS_array, occ_numbers, n_low_state,&
                                                   n_high_state, n_homo, n_full_freq, omega_full_grid, &
                                                   womega_full, chemical_potential, ovlp_3KS, &
                                                   qp_non_convergence,i_iter)
            if(any(qp_non_convergence).and.gw_cd%self_consistent) then
              call replace_nonconverged_qps(gw_cd,qp_energy,qp_energy_ac, Wmn_freq_cd_cmplx, correl_energy, &
                                            KS_eigenvalue, KS_eigenvalue_last, exchange_self_energy,&
                                            xc_KS_matr, x_KS_array, c_KS_array, occ_numbers,&
                                            n_low_state, n_high_state, n_homo, n_full_freq, &
                                            omega_full_grid, womega_full, chemical_potential, &
                                            ovlp_3KS, qp_non_convergence, i_iter)
            endif
            call update_contour_sc_loop(gw_cd,scgw_converged,i_iter,n_iteration,n_spin,n_states,n_homo,&
                                        threshold_sc,qp_energy,KS_eigenvalue_last)
            if(gw_cd%self_consistent.and..not.scgw_converged.and.i_iter/=n_iteration) then
               call qp_output(gw_cd, qp_energy,KS_eigenvalue,occ_numbers,exchange_self_energy,&
                              xc_KS_matr,x_KS_array,c_KS_array,correl_energy,n_low_state,&
                              n_high_state,n_states,n_homo,qp_non_convergence=qp_non_convergence)
               call write_sc_info(gw_cd,n_iteration,i_iter,scgw_converged)
            endif
            if(scgw_converged) exit
         enddo
 
         !**** add relativistic correction to QP energies
         if(flag_adjust_qp_post) then
            call adjust_qp_energies(gw_cd,gw_cd%qp_adjust,qp_energy,n_states,n_spin)
         endif
         !**** output QP energies
         call qp_output(gw_cd, qp_energy,KS_eigenvalue,occ_numbers,exchange_self_energy,&
                        xc_KS_matr,x_KS_array,c_KS_array,correl_energy,n_low_state,&
                        n_high_state,n_states,n_homo,qp_non_convergence=qp_non_convergence)


         if(gw_cd%self_consistent.and.myid.eq.0) then
           write(use_unit,'(T2,A)')"************************************************************"
           if(scgw_converged) then
              write(use_unit,'(T3,A18)') TRIM(gw_cd%sctype) // " converged"
           else
              write(use_unit,'(T3,A33)') "* Warning! " // TRIM(gw_cd%sctype) // " NOT converged" 
           endif
           write(use_unit,'(T2,A)')"************************************************************"
           write(use_unit,*) " "
         endif

         if(flag_print_self_energy) then
            call print_self_energy_cd(gw_cd, Wmn_freq_cd_cmplx, KS_eigenvalue, KS_eigenvalue_last,&
                                      occ_numbers, n_full_freq, n_homo, omega_full_grid, womega_full,&
                                      chemical_potential, ovlp_3KS, state_to_print, print_range)
         endif
         if(flag_calc_spectral_func) then
            call calc_spectral_func(gw_cd, Wmn_freq_cd_cmplx, KS_eigenvalue, &
                                    KS_eigenvalue_last, exchange_self_energy,&
                                    xc_KS_matr, x_KS_array, c_KS_array, occ_numbers, &
                                    n_high_state, n_full_freq, n_homo, omega_full_grid, womega_full, &
                                    chemical_potential, ovlp_3KS)
         endif
      endif
      deallocate(KS_eigenvalue_last)
      deallocate(self_energy_omega)
   end subroutine contour_def_gw_calculation
 
! **************************************************************************************************
!> brief if QP energies not converged, rewrite with QPs from analytic continuation or do Z-shot 
!  o gw_cd  --  contour deformation environment
!  o qp_energy -- QP energies from CD
!  o qp_energy_ac -- QP energies from analytic continuation
!  o i_iter -- iteration step in ev-scGW0 or ev-scGW cycle
! **************************************************************************************************
   subroutine replace_nonconverged_qps(gw_cd,qp_energy,qp_energy_ac,Wmn_freq_cd,correl_energy, &
                                       KS_eigenvalue,KS_eigenvalue_last,exchange_self_energy,&
                                       xc_KS_matr,x_KS_array,c_KS_array,occ_numbers,&
                                       n_low_state,n_high_state,n_homo,n_full_freq, &
                                       omega_full,womega_full,chemical_potential, &
                                       ovlp_3KS,qp_non_convergence,i_iter)
     type(cd_environment_type)                                 :: gw_cd
     real(kind=8), dimension(:,:), intent(inout)               :: qp_energy
     real(kind=8), dimension(:,:), intent(in)                  :: qp_energy_ac
     complex(kind=8), dimension(:,:,:,:), intent(in)           :: Wmn_freq_cd
     real(kind=8), dimension(:,:), intent(inout)               :: correl_energy
     real(kind=8), dimension(:,:), intent(in)                  :: KS_eigenvalue
     real(kind=8), dimension(:,:), intent(in)                  :: KS_eigenvalue_last
     real(kind=8), dimension(:,:), intent(in)                  :: exchange_self_energy
     real(kind=8), dimension(:,:,:), intent(in)                :: xc_KS_matr
     real(kind=8), dimension(:,:),  intent(in)                 :: x_KS_array      
     real(kind=8), dimension(:,:), intent(in)                  :: c_KS_array      
     real(kind=8), dimension(:,:)                              :: occ_numbers
     integer, intent(in)                                       :: n_low_state,n_high_state,&
                                                                  n_full_freq
     integer, dimension(:), intent(in)                         :: n_homo
     real(kind=8), dimension(:), intent(in)                    :: omega_full
     real(kind=8), dimension(:), intent(in)                    :: womega_full
     real(kind=8), dimension(:), intent(in)                    :: chemical_potential   
     real(kind=8), dimension(:,:,:,:), intent(in)              :: ovlp_3KS
     logical, dimension(:,:), intent(inout)                    :: qp_non_convergence
     integer, intent(in)                                       :: i_iter

     character(*), parameter :: func = 'replace_nonconverged_qps'  
 
     character(len=120)                                       :: info_str
     integer                                                  :: i_spin, i_state, i_level
   
     if(i_iter == 1) then 
       do i_spin = 1, n_spin
          do i_level = 1, gw_cd%num_levels(i_spin)
             i_state = gw_cd%corrected_levels(i_level,i_spin)
             if(qp_non_convergence(i_state,i_spin)) then
               qp_energy(i_state,i_spin) = qp_energy_ac(i_state,i_spin)
               qp_non_convergence(i_state,i_spin) = .false.
               if(myid.eq.0) then
                  write(use_unit,'(T2,A23,T26,I4,A13,T44,I2)') "| * Non-converged state", i_state, &
                                                     " spin-channel", i_spin
                  write(use_unit,'(T2,A45,T49,F14.4,T64,A2)') "| * Use QP energy from analytic" // & 
                                                              " continuation for state:", &
                                                              qp_energy(i_state,i_spin)*hartree,"eV"
                  write(use_unit,*) ""                  
               endif
             endif
          enddo
       enddo
     endif

     if(gw_cd%try_zshot.and.i_iter > 1) then
       gw_cd%contour_def_offset = 0.002d0
       call quasi_particle_energy_contour_def_zshot(gw_cd, Wmn_freq_cd,qp_energy,correl_energy, &
                                                KS_eigenvalue,KS_eigenvalue_last,&
                                                exchange_self_energy,&
                                                xc_KS_matr,x_KS_array,c_KS_array,occ_numbers,&
                                                n_low_state,n_high_state,n_homo,n_full_freq, &
                                                omega_full,womega_full,chemical_potential, &
                                                ovlp_3KS,qp_non_convergence=qp_non_convergence)
       do i_spin = 1, n_spin
          do i_level = 1, gw_cd%num_levels(i_spin)
             i_state = gw_cd%corrected_levels(i_level,i_spin)
             if(qp_non_convergence(i_state,i_spin)) then
               qp_non_convergence(i_state,i_spin) = .false.
               if(myid.eq.0) then
                  write(use_unit,'(T2,A23,T26,I4,A13,T44,I2)') "| * Non-converged state", i_state, &
                                                     " spin-channel", i_spin
                  write(use_unit,'(T2,A47,T49,F14.4,T64,A2)') "| * Performed Z-shot (linearization)" // & 
                                                              " for state:", &
                                                              qp_energy(i_state,i_spin)*hartree,"eV"
                  write(use_unit,*) ""                  
               endif
             endif
          enddo
       enddo
     endif
   end subroutine replace_nonconverged_qps

! **************************************************************************************************
!> brief if QP energies not converged, rewrite with QPs from analytic continuation or do Z-shot 
!  o gw_cd  --  contour deformation environment
! **************************************************************************************************
   subroutine qp_output(gw_cd, qp_energy,KS_eigenvalue,occ_numbers,exchange_self_energy,xc_KS_matr,&
                        x_KS_array,c_KS_array,correl_energy,n_low_state,n_high_state,n_states,n_homo,&
                        qp_non_convergence,z_value)

      type(cd_environment_type)                                :: gw_cd
      real(kind=8), dimension(:,:), intent(in)                 :: qp_energy
      real(kind=8), dimension(:,:), intent(in)                 :: KS_eigenvalue
      real(kind=8), dimension(:,:), intent(in)                 :: occ_numbers
      real(kind=8), dimension(:,:), intent(in)                 :: exchange_self_energy
      real(kind=8), dimension(:,:,:), intent(in)               :: xc_KS_matr
      real(kind=8), dimension(:,:), intent(in)                 :: x_KS_array      
      real(kind=8), dimension(:,:), intent(in)                 :: c_KS_array      
      real(kind=8), dimension(:,:), intent(in)                 :: correl_energy
      integer, intent(in)                                      :: n_low_state, n_high_state
      integer, intent(in)                                      :: n_states
      integer, dimension(:), intent(in)                        :: n_homo
      logical, dimension(:,:), intent(in), optional            :: qp_non_convergence
      real(kind=8), dimension(:,:), intent(in), optional       :: z_value

      integer                                                  :: i_spin, i_level, i_state
      real(kind=8), dimension(:,:), allocatable                :: my_exchange_self_energy
      real(kind=8), dimension(:,:), allocatable                :: xc_energy

      ! get correct XC contributions
      allocate(my_exchange_self_energy(n_high_state,n_spin))
      allocate(xc_energy(n_states,n_spin))
      call get_exchange_terms(exchange_self_energy,xc_KS_matr,x_KS_array,c_KS_array,&
                              my_exchange_self_energy,xc_energy)

      if(present(z_value).and.present(qp_non_convergence)) then
        call aims_stop('Only Z-shot or iterative procedure possible')
      endif

      if(present(z_value)) then
        do i_spin = 1, n_spin
           if(gw_cd%spin_channel(i_spin) /= i_spin)  cycle
           do i_level = 1, gw_cd%num_levels(i_spin)
              i_state = gw_cd%corrected_levels(i_level,i_spin)

              my_exchange_self_energy(i_state,i_spin) = z_value(i_state,i_spin)&
                                                      *my_exchange_self_energy(i_state,i_spin)
              xc_energy(i_state,i_spin) = z_value(i_state,i_spin)*xc_energy(i_state,i_spin)
           enddo
        enddo
        call output_quasi_particle_energies(qp_energy,KS_eigenvalue,my_exchange_self_energy,&
                                            x_KS_array,c_KS_array,correl_energy,xc_energy,&
                                            occ_numbers,n_low_state,n_high_state,n_spin,n_states,&
                                            n_homo,gw_cd%contour_def_start,gw_cd%contour_def_end,&
                                            gw_cd%spin_channel,z_value=z_value)
      elseif(present(qp_non_convergence)) then
        call output_quasi_particle_energies(qp_energy,KS_eigenvalue,my_exchange_self_energy,&
                                            x_KS_array,c_KS_array,correl_energy,xc_energy,&
                                            occ_numbers,n_low_state,n_high_state,n_spin,n_states,&
                                            n_homo,gw_cd%contour_def_start,gw_cd%contour_def_end,&
                                            gw_cd%spin_channel,qp_non_convergence=qp_non_convergence)
      else 
        if(myid.eq.0) call aims_warn("* Warning! No output of QP energies")
      endif

      deallocate(my_exchange_self_energy)
      deallocate(xc_energy)
   
   end subroutine qp_output
end module
