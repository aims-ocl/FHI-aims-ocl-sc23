      subroutine self_consistent_gw ( )

      use dimensions
      use runtime_choices
      use species_data
      use physics
      use prodbas
      use hartree_fock
      use gw_para
      use constants
      use mpi_tasks
      use synchronize_mpi
      use numerical_utilities 
      use scgw_grid
      use poles_fit
      use gt
      use localorb_io, only: use_unit
   use scgw_allocations 

      implicit none
!     everything is allocated in the module "scgw_allocations"
      real*8       error_on_G (n_spin)
!
! begin work

      if (myid.eq.0) then
        write(use_unit,'(A)')
        write(use_unit,'(A)')"--------------------------------------------"
        write(use_unit,'(10X,A)') "Self-Consistent GW calculation starts ..."
      endif

      if ((maxval(species_z,n_species).gt.22)&
               .and.( pole_dist_type ==0 ))then
        if(myid.eq.0)then
           write(use_unit,*)"****************** WARNING **********************"
           write(use_unit,*)"*** Self-consistent GW is presently working only "
           write(use_unit,*)"*** for light elements (Z<22).                   "
           write(use_unit,*)'*** At your own risk, you can override this      '
           write(use_unit,*)'*** warning by setting the flag                  '
           write(use_unit,*)'*** "dist_pole_type 1".                          '
           write(use_unit,*)'*** This flag chooses a different parametrization'
           write(use_unit,*)'*** of the pole expansion employed in the        '
           write(use_unit,*)'*** computation of the Fourier transforms (FT).  '
           write(use_unit,*)'*** This has shown to work in some cases, but    '
           write(use_unit,*)'*** does not guaratee the usual accuracy for FT. '
           write(use_unit,*)'************************************************ '
           stop
        endif
      endif


      if(use_mpi) then
         call MPI_Bcast(KS_eigenvector, &
              n_basis*n_states*n_spin*n_k_points, &
              MPI_DOUBLE_PRECISION, &
              0, mpi_comm_global, mpierr)
      endif
      n_high_state=n_states     

       !calculate the KS xc pot (in NAO e KS basis)
       if (.not. allocated(xc_matr))then
          allocate(xc_matr(n_basis,n_basis,n_spin))
       endif
       call integrate_xc_matrix &
          ( partition_tab, &
            rho, rho_gradient, &
            kinetic_density, &
            l_shell_max, &
            xc_matr )

        !determine the HOMO level
        if(.not.allocated(n_homo)) then
           allocate(n_homo(n_spin))
           n_homo(:) = 0
        endif
        do i_spin = 1, n_spin
         do i_state = 1, n_states
           if(occ_numbers(i_state,i_spin,1).gt.1.d-6) then
             n_homo(i_spin)=i_state
           endif
         enddo
        enddo

        !check for fractional occupation -- stop if any 
        if (n_spin .eq.1) then
          do i_state = 1, n_homo(1)
            if ( abs(occ_numbers(i_state,1,1)-2.d0).gt. 1.d-4 )then
                 write(use_unit,*)' Error: sc-GW works only with ', &
                            'integer occupation! STOP. '
               stop
            endif
          enddo
        elseif (n_spin.eq.2)then
          do i_spin = 1, n_spin
            do i_state = 1, n_homo(i_spin)
              if ( abs(occ_numbers(i_state,i_spin,1)-1.d0).gt. 1.d-4 )then  
                   write(use_unit,*)' Error: sc-GW works only with ', &
                              'integer occupation! STOP.  '
                 stop
              endif
            enddo
          enddo
        endif

        !set the chemical_potential 
        if( specify_mu )then 
           chemical_potential = scgw_chem_pot/hartree
        else 
          if(.not.(n_spin .gt.1 .and. n_homo(n_spin).eq.0))then
             chemical_potential = (KS_eigenvalue(n_homo(1),1,1)&
              + KS_eigenvalue(n_homo(1)+1,1,1))/2.d0
             if(myid.eq.0) write(use_unit,*) ' ... setting chemical potential'
             if( n_spin .gt. 1)then
               if((chemical_potential .lt. KS_eigenvalue(n_homo(2),2,1)) &
                   .or. (chemical_potential .gt. KS_eigenvalue(n_homo(2)+1,2,1)) )then
                  chemical_potential = (KS_eigenvalue(n_homo(2),2,1)&
                 + KS_eigenvalue(n_homo(2)+1,2,1))/2.d0
               endif     
             endif     
          endif
        endif

        !get the product basis overlap matrix
        if (.not.allocated(ovlp_3KS)) then
           allocate(ovlp_3KS(n_loc_prodbas,n_states,n_high_state,n_spin))
        end if
        if (.not. use_hartree_fock)then
          if (sparse_o3fn) call aims_stop('No RI-LVL with SCGW', func)
          call transform_ovlp3fn(n_high_state,KS_eigenvector, ovlp_3fn, ovlp_3KS)
        endif

        !initialize the (logarithmic) time and frequency grids
        call init_grid()

!--------------------------------------------------------------
!GET THE DFT/HF GREEN FUNCTION
         
        if(myid.eq.0)then
          write(use_unit,*)"   Allocating matrices for self-consistent GW."
          write(use_unit,'(2X,A,f12.3,A)') &
             "| The Green function in frequency takes :", &
              (dble(n_basis*(n_basis))* nomega )*8.d0/1.024d3/1.024d3/1.024d3," Gbs"
          write(use_unit,'(2X,A,f12.3,A)') &
             "| The Green function in time takes :", &
              (dble(n_basis*(n_basis))* (ntau*2+1) )*8.d0/1.024d3/1.024d3/1.024d3," Gbs"
          write(use_unit,'(2X,A,f12.3,A)') &
             "| The polarizability takes :", &
              (dble(n_max_loc_prodbas*(n_basbas))* (ntau+1) )*8.d0/1.024d3/1.024d3/1.024d3," Gbs"
           write(use_unit,*) " " 
        endif 

        if (.not.allocated(green_fn_freq)) then
           allocate( green_fn_freq(n_basis,n_basis,nomega,n_spin) )
        end if
        if (.not.allocated(green_fn_time)) then !this is updated at each iteration
            allocate( green_fn_time(n_basis,n_basis,-ntau:ntau,n_spin) )
        end if
        if (.not.allocated(green_0_t_0)) then
            allocate( green_0_t_0 (n_basis,n_basis,n_spin) )
        end if
        if (.not.allocated(inv_green_fn_freq)) then
           allocate( inv_green_fn_freq(n_basis,n_basis,nomega,n_spin) )
        end if

        call cpu_time (green_fn_timing)
        call  get_green_function_freq (&
          green_fn_freq, chemical_potential,&
          occ_numbers, hamiltonian, overlap_matrix)
        call get_green_function_time ( n_homo, KS_eigenvalue,&
         KS_eigenvector, green_fn_time, chemical_potential, &
         overlap_matrix,occ_numbers)

!---------------------------------------------------------
!         !for calculations using ensemble Green functions 
!         if(write_G) then !write the green function to file  
!            call write_green_fn ()
!         endif 
!
!        get_ensemble_G = .true.
!!         get_ensemble_G = .false.
!         if(get_ensemble_G .and. read_G .and. n_spin == 2) then ! for fractional charge sc-GW (if possible)
!
!           if (.not.allocated(green_fn_freq_ens)) then
!              allocate( green_fn_freq_ens(n_basis,n_basis,nomega,n_spin))
!              green_fn_freq_ens (:,:,:,:) = 0.d0
!           endif
!           if (.not.allocated(green_fn_time_ens)) then
!               allocate( green_fn_time_ens(n_basis,n_basis,-ntau:ntau,n_spin))
!               green_fn_time_ens (:,:,:,:) = 0.d0
!           endif
!           if (.not.allocated(hartree_dft)) then
!               allocate(hartree_dft (n_basis,n_basis))
!               hartree_dft (:,:) = 0.d0
!           endif
! 
!           call get_hartree_pot &
!             (ovlp_3fn, green_fn_time(:,:,0,:), &
!              hartree_dft)
!                      
!           call read_green_fn  ()
! 
!           fraction_p = first_state_included 
!           !the ensemble G is defined as a linear combination of 2 (or more) different Gs
!           ! where one of the two is red from file
!           green_fn_freq (:,:,:,:) = &
!                 (1.d0-fraction_p) * green_fn_freq     (:,:,:,:) + &
!                      +fraction_p  * green_fn_freq_ens (:,:,:,:)
!           green_fn_time (:,:,:,:) = &
!                 (1.d0-fraction_p) * green_fn_time     (:,:,:,:) + &
!                      +fraction_p  * green_fn_time_ens (:,:,:,:)
!
!           if(allocated(green_fn_freq_ens)) deallocate (green_fn_freq_ens) !allocated in read_green_fn.f90 
!           if(allocated(green_fn_time_ens)) deallocate (green_fn_time_ens)
!        endif
!---------------------------------------------------------------

        do i_spin = 1, n_spin
          call invert_green (green_fn_freq(:,:,:,i_spin),&
               nomega,inv_green_fn_freq(:,:,:,i_spin), n_basis)
        enddo
        call cpu_time (green_fn_timing1)        
        green_fn_timing = green_fn_timing1 - green_fn_timing
        green_0_t_0 (:,:,:) = green_fn_time (:,:,0,:)
        if (myid.eq.0) then
         write(use_unit,'(10X,A,F12.3,A)') &
       "| Total time for calculating the Green's function:            ", &
              green_fn_timing, " s"
        endif

!------------------------------------------
      !OUTPUT GENERAL SETTINGS  
      if(myid.eq.0)then
      write(use_unit,'(A,F10.4,A,F10.4,A)') "          | Chemical Potential : ", &
           chemical_potential , ' Ha', chemical_potential*hartree, ' eV'
write(use_unit,*) "         | Parameter for the poles expansion of G and Sigma"
        write(use_unit,*) "         | Number of poles used : ", n_poles
        write(use_unit,'(A,F10.4)') "          | Smallest pole :", pole_min
        write(use_unit,'(A,F10.4)') "          | Largest  pole :", pole_max
      endif

!-------------TEST FOURIER TRANSFORM ON G

        test_FT = .true.
        if(test_FT) then
           call cpu_time(fourier)
           if(.not. allocated (green_fn_time_FT))then
               allocate (green_fn_time_FT (n_basis,n_basis, -ntau:ntau,n_spin))
               green_fn_time_FT(:,:,:,:) = 0.d0
           endif

           do i_spin = 1, n_spin, 1
             call transform_G &
               (green_fn_freq (:,:,:,i_spin) , n_basis, &
               n_basis, green_fn_time_FT (:,:,:,i_spin))
           enddo

           call cpu_time(fourier1)
           name_of_quantity = "G(t)"
           do i_spin = 1, n_spin            
           call check_the_error&
             (green_fn_time(:,:,:,i_spin), &
              green_fn_time_FT(:,:,:,i_spin), &
              n_basis, n_basis,&
              tau, ntau, wtau, &
              name_of_quantity,&
              average_error, max_error,.false.)
            enddo

            fourier1 = fourier1 - fourier
            if (myid.eq.0) then
              write(use_unit,'(10X,A,F12.3,A)') &
       "| Total time for Fourier transform G(t):                      ", &
                  fourier1, " s"
              write(use_unit,*) " "
            endif
           if(allocated (green_fn_time_FT))then
               deallocate (green_fn_time_FT)
           endif
        endif

!           if(.false.)then !old version of FT -- based on tail fit
!               green_fn_time_FT (:,:,:,i_spin) = 0.d0
!               call transform_and_interpolate_green_fn_to_time &
!                    (  green_fn_freq (:,:,:,i_spin),  &
!                       tau,ntau,wtau, &
!                       omega, nomega, &
!                       womega1,  &
!                       green_fn_time_FT (:,:,:,i_spin), &
!                       inv_overlap_matrix , &
!                       omegamax, ovlp_NAO_KS(:,:,i_spin))
!           endif

!-------------------------------------------------------------------
!Loop initialization

      number_reiterations = 0
      scgw_converged      = .false.
      max_reiteration     = scgw_it_limit
      threshold_green     = scgw_scf_accuracy

      if (myid.eq.0) then
        write(use_unit,*) "  Self-Consistent loop initialization ... "
        write(use_unit,*) " "
      endif

      do while(.not. scgw_converged)
         number_reiterations =  number_reiterations +1 
         if (myid.eq.0) then
            write(use_unit,*) " ----------------------------------------------"
            write(use_unit,*) "             Iteration ", number_reiterations
            write(use_unit,*) " ----------------------------------------------"
            write(use_unit,*) " "
         endif

         if(.false.)then 
           call scgw_densmat_analysis ()
         endif 

!------------------------------------------------------------------
      !SET TYPE OF SELF-CONSISTENT CALCULATION 
      if(use_scgw)then 
        update_W = .true. ! if false perform a GW0 calculation (only partially self-consistent)
      elseif(use_scgw0)then     
        update_W = .false.
        if(myid.eq.0 .and. number_reiterations.eq.1)then
          write(use_unit,*) " A PARTIALLY self-consistent",&
                     " GW_0 calculation will be done, ",&
                     "the W will be calculated once",&
                     " but it won't be updated!"
          write(use_unit,*) " "
        endif
      endif

!------------EVALUATE POLARIZABILITY (IN TIME)
      if(number_reiterations .eq.1 .or. (update_W) )then 

        if (.not.allocated(polar_green_time)) then
           allocate(polar_green_time(n_basbas,n_loc_prodbas,0:ntau))
        endif

        call cpu_time (time_chi_0)
        call get_polar ( green_fn_time, tau, &
           ntau, wtau, polar_green_time)
        call cpu_time (time_chi_1)
        time_chi_1 = time_chi_1 - time_chi_0

        if (myid.eq.0) then
           write(use_unit,'(10X,A,F12.3,A)') &
        "| Total time for calculating the Polarizability:              ", &
             time_chi_1, " s"
        endif

!---------------FOURIER TRANSFORM THE POLARIZABILITY (TO FREQUENCY)

      if (.not.allocated(polar_green_freq)) then
         if(n_loc_prodbas .gt. 0)then 
           allocate(polar_green_freq(n_basbas,n_loc_prodbas,nomega))
         else
           ! WPH, 18 October 2018: changed the second dimension from 1 to
           !      n_loc_prodbas.  This is needed even when n_loc_prodbas.eq.0
           !      because otherwise gfortran will throw an array mismatch
           !      error when aux_polar_freq is set equal to polar_green_freq.
           allocate(polar_green_freq(n_basbas,n_loc_prodbas,nomega))
           polar_green_freq(:,:,:) = 0.d0
         endif
      endif
      call cpu_time(fourier)
      call transform_polar (polar_green_time,  &
       n_basbas, n_loc_prodbas, &
       polar_green_freq )
      call cpu_time(fourier1)
      fourier1 = fourier1 - fourier
      if (myid.eq.0) then
        write(use_unit,'(10X,A,F12.3,A)') &
        "| Total time for Fourier transforming the Polarizability:     ", &
            fourier1, " s"
      endif
      if (allocated(polar_green_time)) then
          deallocate(polar_green_time)
      end if

!---------------EVALUATE \int \chi_0(r1,r2,w) v(r2,r1) dr1 dr2 dw -- UNNEEDED
        if (.false.)then 
          integral_vchi = (0.d0,0.d0) 
          do i_freq=1, nomega
            do i_basbas = 1, n_loc_prodbas
              do j_basbas = 1, n_basbas
                if (i_freq==1)then
                  integral_vchi = integral_vchi + &
                  polar_green_freq (j_basbas, i_basbas, i_freq)&
                  * womega1(i_freq) *2.d0/dble(n_spin)
                else 
                  integral_vchi = integral_vchi + &
                  polar_green_freq (j_basbas, i_basbas, i_freq) &
                  *womega1(i_freq)*2.0 *2.d0/dble(n_spin)
                endif
              enddo
            enddo
          enddo 
          if(myid.eq.0 .and. .false. )then
            write(use_unit,*) "  \int \chi_0(r1,r2,w) v(r2,r1) dr1 dr2 dw =  ",&
                        real(integral_vchi)
          endif
        endif 

!---------- EVALUATE THE RPA CORRELATION ENERGY -- for comparison pourpose only
      if(.true.)then
        rpa_c_energy = 0.d0
        if (.not.allocated(aux_polar_freq)) then
             allocate(aux_polar_freq(n_basbas,n_loc_prodbas))
        endif
  
        do i_freq =1, nomega, 1
         aux_polar_freq (:,:) = polar_green_freq(:,:,i_freq)* 2.d0/dble(n_spin)
         call evaluate_rpa_integrand(aux_polar_freq,rpa_c_integrand)
         rpa_c_energy = rpa_c_energy + &
                        rpa_c_integrand * womega1(i_freq)         
        enddo
        rpa_c_energy=rpa_c_energy/2.d0/pi

        if (allocated(aux_polar_freq)) then
           deallocate(aux_polar_freq)
        endif
      endif


!----------EVALUATE THE SCREENED COULOMB INTERACTION (IN FREQUENCY)
    
      if (.not.allocated(screened_coul_int_freq)) then
         if(n_loc_prodbas .gt. 0)then 
           allocate(screened_coul_int_freq(n_basbas,n_loc_prodbas,nomega))
         else
           allocate(screened_coul_int_freq(n_basbas,1,nomega)) !fake allocation. This isn't used anywhere, in practice.
           screened_coul_int_freq (:,:,:) = 0.d0
         endif
      endif         

      call cpu_time(time_W)
      screened_coul_int_freq(:,:,:) = 0.d0
      if (myid.eq.0) then
        write(use_unit,*) " "
        write(use_unit,*) "  --- Evaluating the Screened Coulomb interaction"
      endif

      do i_freq = 1, nomega,1
        call screened_coulomb_interaction     &
              (polar_green_freq(:,:,i_freq),  &
              screened_coul_int_freq(:,:,i_freq))

        ! subtract the bare Coulomb potential from W, Wc = W-V 
        do j_basbas = 1, n_loc_prodbas, 1
          i_index = map_prodbas(j_basbas,myid+1)
          do i_basbas = 1, n_basbas, 1
            if(i_index.eq.i_basbas) then
               screened_coul_int_freq(i_basbas, j_basbas,i_freq) = &
                 screened_coul_int_freq(i_basbas, j_basbas,i_freq) - 1.d0
            endif
          enddo
        enddo
      enddo
      

      call cpu_time(time_W1)
      time_W1 = time_W1 - time_W
      if (myid.eq.0) then
        write(use_unit,'(10X,A,F12.3,A)') &
       "| Total time for calculating the Screened Coulomb interaction:", &
             time_W1, " s"
      endif

      if (allocated(polar_green_freq)) then
         deallocate(polar_green_freq)
      end if

!!----------------TRANSFORM THE SCI TO TIME  

        if (.not.allocated(screened_coul_int_time)) then
           allocate(screened_coul_int_time(n_basbas,n_loc_prodbas,0:ntau))
        endif
        call cpu_time(fourier)
        call transform_W (screened_coul_int_freq,&
           n_basbas, n_loc_prodbas,&
           screened_coul_int_time )
        call cpu_time(fourier1)

        fourier1 = fourier1 - fourier
        if (myid.eq.0) then
            write(use_unit,'(10X,A,F12.3,A)') &
        "| Total time for Fourier transforming the Scr. Coul. Int.:    " &
            ,  fourier1, " s"
            write(use_unit,*) " "
        endif

        if (allocated(screened_coul_int_freq)) then
           deallocate(screened_coul_int_freq)
        endif

      endif !update_W

!------------- EVALUATE THE SELF ENERGY  (IN TIME)

      if(myid.eq.0)then
        write(use_unit,*)"  --- Evaluating the Self Energy "
      endif
      if (.not.allocated(self_energy_time)) then
         allocate(self_energy_time(n_basis,n_basis,-ntau:ntau,n_spin))
      endif

      call cpu_time (time_self_0)
      do i_spin = 1, n_spin
        call get_self_energy ( green_fn_time(:,:,:,i_spin), &
           screened_coul_int_time, &
           self_energy_time (:,:,:,i_spin) )
      enddo
      call cpu_time (time_self_1)
      time_self_1 = time_self_1 - time_self_0

      if (myid.eq.0) then
        write(use_unit,'(10X,A,F12.3,A)') &
      "| Total time for calculating the Self Energy:                 ", &
             time_self_1, " s"
      endif
      if(update_W)then
        if (allocated(screened_coul_int_time)) then
           deallocate(screened_coul_int_time)
        endif
      endif

!----------------FOURIER TRANSFORM THE SELF ENERGY (TO FREQUENCY)
      if(myid.eq.0) then
        write(use_unit,*)&
"         | Fourier transform of the Self Energy from time to frequency "
      endif
      if (.not.allocated(self_energy_freq)) then
         allocate(self_energy_freq(n_basis,n_basis,nomega,n_spin))
      endif

      call cpu_time(fourier)
      do i_spin = 1, n_spin
         call transform_sigma (&
         self_energy_time(:,:,:,i_spin), n_basis, n_basis,  &
         self_energy_freq(:,:,:,i_spin))
      enddo
      call cpu_time(fourier1)
      fourier1 = fourier1 - fourier
      if (myid.eq.0) then
        write(use_unit,'(10X,A,F12.3,A)') &
      "| Total time for Fourier transforming the Self Energy:        ", &
             fourier1, " s"
      endif

      if (allocated(self_energy_time)) then
         deallocate(self_energy_time)
      endif

!---------------------------------------------------------
      if(myid.eq.0)then
        write(use_unit,*)"         | Evaluating the Exact-exchange and Hartree operators"
      endif

      if (.not. allocated(exchange_self_energy0))then
        allocate(exchange_self_energy0(n_basis,n_basis, n_spin))
      endif
      if (.not. allocated(exchange_self_energy))then
         allocate(exchange_self_energy(n_basis,n_basis, n_spin))
      endif
      if (.not. allocated(hartree_pot))then
         allocate(hartree_pot(n_basis,n_basis))
      endif
      if (.not. allocated(hartree_pot0))then
         allocate(hartree_pot0(n_basis,n_basis))
      endif

      call get_hartree_pot &
           (ovlp_3fn, green_fn_time(:,:,0,:), &
            hartree_pot)
      call get_exchange_self_energy_v1 &
           (ovlp_3fn, green_fn_time(:,:,0,:), &
            exchange_self_energy)

      if(number_reiterations.eq.1)then
        hartree_pot0 (:,:)        = hartree_pot (:,:)
        exchange_self_energy0 (:,:,:) = exchange_self_energy (:,:,:)
      endif

!-------------------------------------------------------

!      if(.not. read_G )then
          call total_energy_calculation_v1 (&
                green_fn_time(:,:,0,:),    &
                green_fn_freq,             &
                exchange_self_energy,      &
                exchange_self_energy0,     &
                hartree_pot,           &
                hartree_pot0,          &
                self_energy_freq,          &
                xc_matr,                   &
                rpa_c_energy               &
                )
!      else 
!          call total_energy_calculation_v1 (&
!                green_fn_time(:,:,0,:),    &
!                green_fn_freq,             &
!                exchange_self_energy,      &
!                exchange_self_energy0,     &
!                hartree_pot,           &
!                hartree_dft,          &
!                self_energy_freq,          &
!                xc_matr,                   &
!                rpa_c_energy               &
!                )
!      endif
!
!!-----------SOLVE DYSON'S EQUATION
      if(scgw_print_all_spectrum)then
           if(myid.eq.0)then
              write(use_unit,*) "  --- Calculating the Spectrum"
           endif
         do i_spin =1, n_spin
           call  solve_dyson_equation_re_v2 (         &
                self_energy_freq      (:,:,:,i_spin), &
                xc_matr               (:,:,i_spin),   &
                hartree_pot       (:,:),          &
                hartree_pot0      (:,:),          &
                exchange_self_energy0 (:,:,i_spin),   &
                exchange_self_energy  (:,:,i_spin),   &
                number_reiterations,                  &
                scgw_converged,                       &
                i_spin)

!          call solve_dyson_equation_re_v0 &
!              ( inv_green_fn_freq(:,:,:,i_spin), &
!                self_energy_freq(:,:,:,i_spin), &
!                green_fn_freq(:,:,:,i_spin), &
!                exchange_self_energy(:,:,i_spin), &
!                xc_matr(:,:,i_spin), hartree_pot(:,:,i_spin),&
!                exchange_self_energy0(:,:,i_spin), &
!                new_chem_pot, number_reiterations,&
!                scgw_converged, i_spin)
         enddo
      endif

      call solve_dyson_equation_v1 &
         ( inv_green_fn_freq,      &
           self_energy_freq,       & 
           green_fn_freq,          & 
           exchange_self_energy,   &
           xc_matr,                &
           hartree_pot,        &
           hartree_pot0,       &
           exchange_self_energy0 )

        
!------------------------------------------CHECK THE SPARSITY of the Green's function
! Print out the green's function, in a 
! gnuplot readable format! Good to see what are the
! distribution of non-null matrix elements in G. 
! Sparse matrix multiplication could become an option at a certain point!

      if (.false.)then
        if( number_reiterations.lt.10 ) then
           write(iter,'(A,I1)') "0",number_reiterations
        else
           write(iter,'(I2)') number_reiterations
        endif
        filename = "hys_G_"//iter//".dat"
        open(555,file = filename)
        do i_basis = 1, n_basis, 1
          do j_basis = 1, n_basis, 1    
            g1 = 0
            do i_freq = 1, nomega, 1
              g1 = g1 + abs(green_fn_freq (i_basis,j_basis,i_freq,1))*&
                   womega1(i_freq)
            enddo              
            write(555,'(I4, I4, f9.5)') i_basis, j_basis, g1
          enddo
          write(555,*) " "
        enddo         
        close(555)
      endif

!---------------------TRANSFROM G(iw) -> G(it)
      if (.not. allocated(new_green_fn_time))then
         allocate(new_green_fn_time(n_basis,n_basis,-ntau:ntau,n_spin))
      endif
      call cpu_time(fourier)
      do i_spin = 1, n_spin, 1
           call transform_G &
           (green_fn_freq (:,:,:,i_spin), n_basis, &
           n_basis, new_green_fn_time(:,:,:,i_spin) )  
      enddo
      call cpu_time(fourier1)
      fourier1 = fourier1 - fourier
      if (myid.eq.0) then
        write(use_unit,'(10X,A,F12.3,A)') &
      "| Total time for Fourier transforming the G(t)                ", &
            fourier1, " s"
      endif

      !check deviation from previous iteration 
      name_of_quantity = "G(it)"
      do i_spin = 1, n_spin
        call check_the_error&
         (green_fn_time(:,:,:,i_spin), &
          new_green_fn_time(:,:,:,i_spin), &
          n_basis, n_basis,&
          tau, ntau, wtau, &
          name_of_quantity,&
          error_on_G(i_spin), max_error, .false.)
      enddo

      !linear mixing -- if required
      if(number_reiterations.gt.2)then
        if(myid.eq.0 .and. scgw_mix_param .gt. 1.d-8) &
          write(use_unit,*) '         | Using mixing parameter alpha = ', scgw_mix_param
        green_fn_time (:,:,:,:) = (1-scgw_mix_param) * new_green_fn_time(:,:,:,:)&
        +scgw_mix_param * green_fn_time(:,:,:,:)
      else 
        green_fn_time (:,:,:,:) = new_green_fn_time(:,:,:,:)
      endif

      if (allocated(new_green_fn_time))then
         deallocate(new_green_fn_time)
      endif

!------CALCULATE NUMBER OF PARTICLEs
      if(.not. allocated(aux_overlap_matrix)) then 
         allocate(aux_overlap_matrix(n_basis,n_basis))
      endif
      i_index = 0
      do i_basis=1, n_basis, 1
       do j_basis = 1, i_basis, 1
          i_index = i_index + 1
          aux_overlap_matrix(i_basis,j_basis) =  overlap_matrix(i_index)
          if(.not.i_basis.eq.j_basis)then
            aux_overlap_matrix(j_basis,i_basis) = overlap_matrix(i_index)
          endif
        enddo
      enddo
      n_part = 0
      do i_spin = 1, n_spin
        do i_basis=1, n_basis, 1
          do j_basis = 1, n_basis, 1
             n_part = n_part + green_fn_time (i_basis,j_basis,0,i_spin)*&
             aux_overlap_matrix(i_basis,j_basis)
          enddo
        enddo
      enddo
      n_part =(2.d0/n_spin)*n_part
      if(myid.eq.0)then
         write(use_unit,'(A,A,F12.7)') "          | Number of particle from ",&
             " G   : ",  n_part
         write(use_unit,*) " "
      endif

      if(allocated(aux_overlap_matrix)) then 
         deallocate(aux_overlap_matrix)
      endif

      !CHECK IF SELF CONSISTENCY HAS BEEN REACHED
      if (n_spin ==2)then
        if ((error_on_G(1).lt. threshold_green).and.&
            (error_on_G(2).lt. threshold_green))then
          scgw_converged = .true.
        endif
      elseif(n_spin==1)then  
        if (error_on_G(1).lt. threshold_green)then
          scgw_converged = .true.
        endif
      endif
      if ( number_reiterations.ge.max_reiteration ) scgw_converged = .true.
      enddo ! while(.not. scgw_convergedd), end of the self consistent loop

        !----------CALCULATE THE DENSITY OF STATES as ImG(w) 
        if(.true.)then
             if(myid.eq.0)then
                write(use_unit,*) "  --- Calculating the Spectrum"
             endif
           do i_spin =1, n_spin
             call  solve_dyson_equation_re_v2 (         &
                  self_energy_freq      (:,:,:,i_spin), &
                  xc_matr               (:,:,i_spin),   &
                  hartree_pot       (:,:),          &
                  hartree_pot0      (:,:),          &
                  exchange_self_energy0 (:,:,i_spin),   &
                  exchange_self_energy  (:,:,i_spin),   &
                  number_reiterations,                  &
                  scgw_converged,                       &
                  i_spin)
           enddo
        endif

        if(allocated(self_energy_freq))then
           deallocate(self_energy_freq)
        endif
        if(allocated(hartree_pot))then
           deallocate(hartree_pot)
        endif
        if(allocated(exchange_self_energy))then
           deallocate(exchange_self_energy)
        endif

        if (myid.eq.0) then
         write(use_unit,*) " ----------------------------------------  "
         write(use_unit,*) "     End of the SCF-GW loop      "
         write(use_unit,*) " ----------------------------------------  "
         write(use_unit,*) "   "
         write(use_unit,*) "   Number of loops = ",number_reiterations
        endif

        !evaluate dipole moment with self-consistent density matrix
        if(.not. allocated (green_fn))then
           allocate(green_fn(n_basis,n_basis))
        endif
        green_fn(:,:)=0.d0
        if (n_spin .eq.1)then
          green_fn(:,:) = green_fn_time(:,:,0,1)
        elseif (n_spin .eq.2)then
          green_fn(:,:) = green_fn_time(:,:,0,1) + green_fn_time(:,:,0,2)
        endif
        call compute_scgw_dipole ()

       if (allocated(ovlp_3KS)) then
         deallocate( ovlp_3KS )
       endif
       if(allocated(n_homo)) then
        deallocate(n_homo)
       endif
       if (allocated(green_fn_freq)) then
         deallocate(green_fn_freq)
       endif
       if (allocated(green_fn_time)) then
         deallocate(green_fn_time)
       endif
       if (allocated(xc_matr)) then
         deallocate( xc_matr )
       endif

       if (myid.eq.0) then
         write(use_unit,*) " --- scGW deallocations: completed --- "
       endif

      end subroutine self_consistent_gw
