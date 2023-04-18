            subroutine spectrum (new_green_fn_freq, number_reiterations, scgw_converged)
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

            implicit none

            integer number_reiterations
            logical scgw_converged
            real*8 inv_overlap_matrix (n_basis, n_basis)
            complex*16 new_green_fn_freq (n_basis, n_basis, nomega, n_spin)
            complex*16, dimension(:,:,:), allocatable:: diagonal_green_fn
!            complex*16, dimension(:,:,:,:), allocatable:: diagonal_green_fn
            complex*16, dimension(:,:,:), allocatable:: green_fn_par

            integer i_spin
            character*2 iter
            character*15 filename
            real*8 new_chem_pot

              if(.not.allocated(diagonal_green_fn))then
                allocate(diagonal_green_fn(n_states, n_states, nomega))
                !allocate(diagonal_green_fn(n_states, n_states, nomega,n_spin))
                !diagonal_green_fn(:,:,:,:) = (0.d0,0.d0)
              endif
              if(.not.allocated(green_fn_par))then
                allocate(green_fn_par(n_max_par, n_states, n_states))
              endif

            do i_spin = 1, n_spin
              if(scgw_converged) then
                 if(n_spin .eq. 1 )then
                 !  filename = "sp_ImG"//iter//".dat"
                 filename = 'spectrum_sc.dat'
                 else
                     if (i_spin ==1)  filename = "spect_sc_up.dat"
                     if (i_spin ==2)  filename = "spect_sc_do.dat"
                 endif

              else
                 if (number_reiterations.lt.10)then
                   write(iter,'(A,I1)') "0", number_reiterations
                 else
                   write(iter,'(I2)') number_reiterations
                 endif
                 
                 if(n_spin .eq. 1 )then
                   filename = "sp_ImG"//iter//".dat"
                 else
                     if (i_spin ==1)  filename = "sp_ImG"//iter//"_SU.dat"
                     if (i_spin ==2)  filename = "sp_ImG"//iter//"_SD.dat"
                 endif
              endif


              !transformation in the KS basis
              diagonal_green_fn(:,:,:) = (0.d0,0.d0)
              call diagonalize_green_fn &
                  ( new_green_fn_freq(:,:,:,i_spin), &
                   diagonal_green_fn(:,:,:), i_spin)

              !analytic contitnuation (preferably with a 2 poles fit)
              call analy_continue_green_fn &
                  (anacon_type,&
                   nomega, &
                   n_max_par, &
                   green_fn_par,omega, &
                   !diagonal_green_fn(:,:,:,i_spin), n_states)
                   diagonal_green_fn(:,:,:), n_states)

              call get_spectrum (anacon_type, green_fn_par, n_max_par, &
                     omega, nomega  ,inv_overlap_matrix, filename, n_states,new_chem_pot)

             enddo ! I SPIN
              if(allocated(diagonal_green_fn))then
                deallocate(diagonal_green_fn)
              endif
              if(allocated(green_fn_par))then
                deallocate(green_fn_par)
              endif
 
          end subroutine spectrum
