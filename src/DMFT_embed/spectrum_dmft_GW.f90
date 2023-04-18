            subroutine spectrum_dmft_GW (new_green_fn_freq,&
                                         hamiltonian_GW,&
                                         new_chemical_potential_cluster,&
                                         new_chemical_potential,&
                                         embed_part_number_LDA)
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
              use scgw_grid
              use poles_fit
              use localorb_io, only: use_unit


            integer number_reiterations
            complex*16 new_green_fn_freq (n_basis, n_basis, nomega, n_k_points)
            complex*16, dimension(:,:,:), allocatable:: diagonal_green_fn
!            complex*16, dimension(:,:,:,:), allocatable:: diagonal_green_fn
            complex*16, dimension(:,:,:), allocatable:: green_fn_par
            complex*16 hamiltonian_GW(n_basis, n_basis, nomega, n_k_points)
            integer aux_nomega
            real*8  aux_omegamax
            real*8  embed_part_number_LDA
            real*8  n_particles
            real*8  new_chemical_potential
            real*8  new_chemical_potential_cluster
            real*8 , dimension (:), allocatable :: aux_omega
            real*8 , dimension (:), allocatable :: aux_womega
            real*8 , dimension (:), allocatable :: spectrum
            real*8 , dimension (:), allocatable :: add_spectrum
            complex*16 w

            integer i_spin,i_k_point, sign_freq,i_freq
            character*4 iter
            character*25 filename
            character*25 filename_k
            character*25 filename_k_1
            character*25 filename_k_2
            character*25 filename_k_3
            character*25 filename_k_4
            real*8 new_chem_pot

              if(.not.allocated(diagonal_green_fn))then
                allocate(diagonal_green_fn(n_states, n_states, nomega))
                !allocate(diagonal_green_fn(n_states, n_states, nomega,n_spin))
                !diagonal_green_fn(:,:,:,:) = (0.d0,0.d0)
              endif
              if(.not.allocated(green_fn_par))then
                allocate(green_fn_par(n_max_par, n_states, n_states))
              endif

!          new_chemical_potential = -0.85/hartree
          aux_nomega = 15000
          aux_omegamax = 6.50
          !write(use_unit,*) new_chemical_potential, new_chemical_potential*hartree 
              allocate(aux_omega(aux_nomega))
              allocate(aux_womega(aux_nomega))
              allocate(spectrum(-aux_nomega:aux_nomega))
              allocate(add_spectrum(-aux_nomega:aux_nomega))
              !allocate(spectrum_partial(-aux_nomega:aux_nomega))

      call gauleg(0d0, aux_omegamax, aux_omega, aux_womega, aux_nomega)




!            do i_spin = 1, n_spin
!              if(scgw_converged) then
                 if(n_spin .eq. 1 )then
                 !  filename = "sp_ImG"//iter//".dat"
                 filename = 'spectrum_sc.dat'
                 else
                     if (i_spin ==1)  filename = "spect_sc_up.dat"
                     if (i_spin ==2)  filename = "spect_sc_do.dat"
                 endif

            !  else
                 if (number_reiterations.lt.10)then
                   write(iter,'(A,I1)') "0", number_reiterations
                 else
                   write(iter,'(I2)') number_reiterations
                 endif
                 
                 if(n_spin .eq. 1 )then
                   filename = "sp_ImG.dat"
                 else
                   if (i_spin ==1)  filename = "sp_ImG"//iter//"_SU.dat"
                   if (i_spin ==2)  filename = "sp_ImG"//iter//"_SD.dat"
                 endif
!              endif
         add_spectrum(:) = 0.d0
              do i_k_point = 1, n_k_points,1
!              do i_k_point = 1,1
              !transformation in the KS basis
              diagonal_green_fn(:,:,:) = (0.d0,0.d0)
              green_fn_par(:,:,:) = (0.d0,0.d0)
              call diagonalize_green_fn_dmft_GW &
                  ( new_green_fn_freq(:,:,:,i_k_point), &
                   hamiltonian_GW(:,:,:,i_k_point),&
                   diagonal_green_fn, i_k_point,&
                   new_chemical_potential)

              !analytic contitnuation (preferably with a 2 poles fit)
              call analy_continue_green_fn &
                  (anacon_type,&
                   nomega, &
                   n_max_par, &
                   green_fn_par,omega, &
                   !diagonal_green_fn(:,:,:,i_spin), n_states)
                   diagonal_green_fn, n_states)

              call get_spectrum_dmft_GW (anacon_type, green_fn_par, n_max_par, &
                     omega, nomega  , filename, n_states,new_chem_pot,spectrum,&
                     aux_nomega, aux_omega , aux_womega) 

        if (myid.eq.0)then
!        if (.false.)then

                 if (i_k_point.lt.10)then
                   write(iter,'(A,I1)') "000", i_k_point
                 elseif(i_k_point.lt.100)then
                   write(iter,'(A,I2)') "00", i_k_point
                 elseif(i_k_point.lt.1000)then
                   write(iter,'(A,I3)') "0", i_k_point
                 else  
                   write(iter,'(I4)') i_k_point
                 endif

          filename_k = "sp_ImG"//iter//".dat"
          open(55, file= filename_k)
          filename_k_1 = "sp_ImG_1_"//iter//".dat"
          open(56, file= filename_k_1)
          filename_k_2 = "sp_ImG_2_"//iter//".dat"
          open(57, file= filename_k_2)
          filename_k_3 = "sp_ImG_3_"//iter//".dat"
          open(58, file= filename_k_3)
          filename_k_4 = "sp_ImG_4_"//iter//".dat"
          open(59, file= filename_k_4)
           do sign_freq = -1, 1, 2
             do i_freq = 1, aux_nomega, 1
               w = sign_freq*aux_omega(i_freq)
           if(add_spectrum (sign_freq*i_freq).lt.0) then

              write(55, *)(real(w)-(new_chemical_potential-chemical_potential))*hartree,&
              - spectrum (sign_freq*i_freq)
              write(56, *)(real(w)-(chemical_potential))*hartree,&
              - spectrum (sign_freq*i_freq)
              write(57, *)(real(w)-(new_chemical_potential_cluster+chemical_potential))*hartree,&
              - spectrum (sign_freq*i_freq)
              write(58, *)(real(w)+(new_chemical_potential_cluster+chemical_potential))*hartree,&
              - spectrum (sign_freq*i_freq)
              write(59, *)(real(w))*hartree,&
!              write(77, *)(real(w)-(new_chemical_potential))*hartree,&
!              write(55, *)(real(w))*hartree,&
!               write(77, *)(real(w)-(chemical_potential+new_chemical_potential_cluster-new_chemical_potential))*hartree,&
              ! write(77, *)(real(w))*hartree-5.3,&
                         - spectrum (sign_freq*i_freq)
           else
               write(55, *)(real(w)-(new_chemical_potential-chemical_potential))*hartree,&
               spectrum (sign_freq*i_freq)
               write(56, *)(real(w)-(chemical_potential))*hartree,&
               spectrum (sign_freq*i_freq)
               write(57, *)(real(w)-(new_chemical_potential_cluster+chemical_potential))*hartree,&
               spectrum (sign_freq*i_freq)
               write(58, *)(real(w)+(new_chemical_potential_cluster+chemical_potential))*hartree,&
               spectrum (sign_freq*i_freq)
               write(59, *)(real(w))*hartree,&
!               write(77, *)(real(w)-(new_chemical_potential))*hartree,&
!              write(55, *)(real(w))*hartree,&
!               write(77, *)(real(w)-(chemical_potential+new_chemical_potential_cluster-new_chemical_potential))*hartree,&
               !write(77, *)real(w)*hartree-5.3,&
                         spectrum (sign_freq*i_freq)
!                        spectrum_partial (sign_freq*i_freq)
           endif
                enddo
                enddo
               close(55)
               close(56)
               close(57)
               close(58)
               close(59)
              endif

               add_spectrum(:) = add_spectrum(:) + spectrum(:)

             enddo ! i_k_points

               add_spectrum(:) = add_spectrum(:)*(1./n_k_points)

!write(use_unit,*) , add_spectrum(:)
     !if(.true.)then
     if(.false.)then
         n_particles = 0.d0
         do i_freq= 1, aux_nomega, 1
          if( aux_omega(i_freq) .lt. -new_chemical_potential )then
          !if( aux_omega(i_freq) .lt. 0 )then
            n_particles = n_particles + add_spectrum(i_freq)*aux_womega(i_freq)
          else
            n_particles = n_particles +&!  spectrum(i_freq)*aux_womega(i_freq)+&
              add_spectrum(-i_freq)*aux_womega(i_freq)
          endif
         enddo

        if (myid.eq.0)then
          write(use_unit,*) " "
          write(use_unit,*) "n_particles as calculated from the spectrum = ", &
            n_particles
          write(use_unit,*) " "
        endif
      endif

      if(.true.)then
!      if(.false.)then
!    chemical_potential  = -10.d0
         n_particles = 0.d0
         do i_freq= 1, aux_nomega, 1
          if( aux_omega(i_freq) .lt. chemical_potential )then
            n_particles = n_particles + add_spectrum(i_freq)*aux_womega(i_freq)
          else
            n_particles = n_particles +&!  spectrum(i_freq)*aux_womega(i_freq)+&
              add_spectrum(-i_freq)*aux_womega(i_freq)
          endif

         enddo

        if (myid.eq.0)then
          write(use_unit,*) " "
          write(use_unit,*) "n_particles as calculated from the spectrum = ", &
            n_particles
          write(use_unit,*) " "
        endif
      endif

!stop
        if (myid.eq.0)then

          open(77, file='spectrum_sc.dat')
           do sign_freq = -1, 1, 2
             do i_freq = 1, aux_nomega, 1
               w = sign_freq*aux_omega(i_freq)
!           if(add_spectrum (sign_freq*i_freq).lt.0) then

!              write(77, *)(real(w)-(new_chemical_potential-new_chemical_potential_cluster-chemical_potential))*hartree,&
              write(77, *)(real(w))*hartree,&
!              write(77, *)(real(w))*hartree,&
!               write(77, *)(real(w)-(chemical_potential+new_chemical_potential_cluster-new_chemical_potential))*hartree,&
              ! write(77, *)(real(w))*hartree-5.3,&
                          add_spectrum (sign_freq*i_freq)
!           else
!               write(77, *)(real(w)-(new_chemical_potential-new_chemical_potential_cluster-chemical_potential))*hartree,&
!               write(77, *)(real(w)+(new_chemical_potential))*hartree,&
!              write(77, *)(real(w))*hartree,&
!               write(77, *)(real(w)-(chemical_potential+new_chemical_potential_cluster-new_chemical_potential))*hartree,&
               !write(77, *)real(w)*hartree-5.3,&
!                         add_spectrum (sign_freq*i_freq)
!             write(45, *)(real(w)+e_fermi)*hartree,&
!                        spectrum_partial (sign_freq*i_freq)
!           endif
                enddo
                enddo
               close(77)
              endif
              if(allocated(diagonal_green_fn))then
                deallocate(diagonal_green_fn)
              endif
              if(allocated(green_fn_par))then
                deallocate(green_fn_par)
              endif
 
          end subroutine spectrum_dmft_GW
