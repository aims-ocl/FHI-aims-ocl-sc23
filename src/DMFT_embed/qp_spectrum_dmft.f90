            subroutine qp_spectrum_dmft (self_energy_freq, &
                                         DMFT_ham_NAO,&
                                         full_ovlp_matr_sqrt,&
                                         full_ovlp_matr)
                                         
                                      



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

      integer i_freq
      integer i_k_point
      integer i_state
      integer i_k


      character*17 filename

      complex*16, dimension(n_basis, n_basis, nomega) :: self_energy_freq
!      complex*16, dimension(n_basis, n_basis, nomega,n_k_points) :: self_energy_freq
      complex*16, dimension(:,:,:),allocatable :: DMFT_ham_KS
      complex*16, dimension(n_basis, n_basis, n_k_points) :: DMFT_ham_NAO
      complex*16, dimension(n_basis, n_basis, n_k_points) :: full_ovlp_matr_sqrt
      complex*16, dimension(n_basis, n_basis, n_k_points) :: full_ovlp_matr
      !complex*16, dimension(n_states, n_states, n_k_points) :: DMFT_ham_NAO
      complex*16, dimension(:,:,:), allocatable :: self_energy_omega
      complex*16, dimension(:,:,:,:), allocatable :: self_energy_omega_KS
      complex*16, dimension(:,:,:,:), allocatable :: sigma_par_loc
!      complex*16, dimension(:,:,:), allocatable :: sigma_par_loc
      complex*16, dimension(:,:,:), allocatable :: aux_matr
      complex*16, dimension(:,:,:), allocatable :: KS_eigenvector_tmp





              if(myid.eq.0)then
                write(use_unit,'(A)') " -------------------------------------"
              endif

        if(.not.allocated(self_energy_omega))then
          allocate(self_energy_omega (n_basis,n_basis, nomega))
        endif
        if(.not.allocated(self_energy_omega_KS))then
          allocate(self_energy_omega_KS (n_states,n_states, nomega,n_k_points))
!          allocate(self_energy_omega_KS (n_basis,n_basis, nomega,n_k_points))
        endif
        if(.not.allocated(DMFT_ham_KS))then
          allocate(DMFT_ham_KS (n_states,n_states, n_k_points))
        endif
        if(.not.allocated(aux_matr))then
          allocate(aux_matr(n_states,n_states, nomega))
!          allocate(aux_matr(n_basis,n_basis, nomega))
        endif
       if (.not.allocated(sigma_par_loc)) then
          allocate(sigma_par_loc(n_max_par, n_states,n_states,n_k_points))
!          allocate(sigma_par_loc(n_max_par, n_states,n_states))
!          allocate(sigma_par_loc(n_max_par, n_basis,n_basis,n_k_points))
!          allocate(sigma_par_loc(n_max_par, n_basis,n_basis))
       end if


          aux_matr(:,:,:) = (0.d0,0.d0)

   i_k = 0
              do i_k_point = 1, n_k_points,1

      if(.not. allocated(KS_eigenvector_tmp))then
        allocate(KS_eigenvector_tmp(n_basis,n_states,n_spin))
      endif

           if (myid .eq. MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then

            i_k = i_k + 1

             if(real_eigenvectors)then

               KS_eigenvector_tmp(:,:,:) = DCMPLX(KS_eigenvector(:,:,:,i_k))

              else

               KS_eigenvector_tmp(:,:,:) = KS_eigenvector_complex(:,:,:,i_k)

             endif

           else
             ! zero temp. KS eigenvector on all other threads, prior to allreduce
             KS_eigenvector_tmp(:,:,:) = (0.d0,0.d0)
          end if
           call sync_eigenvector_complex(KS_eigenvector_tmp)




              do i_freq = 1, nomega,1


                 call multiply_ovlp_matr(self_energy_freq(:,:,i_freq),&
                                         full_ovlp_matr_sqrt(:,:,i_k_point),&
                                         !full_ovlp_matr(:,:,i_k_point),&
                                         self_energy_omega(:,:,i_freq))
                                         !self_energy_omega_KS(:,:,i_freq,i_k_point))


                call transform_to_KS_basis_complex(&
                                       self_energy_omega(:,:,i_freq),i_k_point,&
                                       !self_energy_freq(:,:,i_freq),i_k_point,&
                                       self_energy_omega_KS(:,:,i_freq,i_k_point),&
                                       KS_eigenvector_tmp)

!               do i_state = 1, n_states,1
!                 self_energy_omega(i_state,i_freq,i_k_point)=&
!                 self_energy_omega_KS(i_state,i_state,i_freq,1)
!               enddo
          aux_matr(:,:,i_freq) =  aux_matr(:,:,i_freq) + &
                                  self_energy_omega_KS(:,:,i_freq,i_k_point) !&
              enddo

           
!       if(.false.)&
                call analy_continue_self_energy_dmft &
                      (anacon_type, &
                       1,n_states,nomega, &
                       !1,n_basis,nomega, &
                       n_max_par, &
                       chemical_potential, &
                       sigma_par_loc(:,:,:,i_k_point),omega, &
                       !self_energy_omega(:,:,:))
                       self_energy_omega_KS(:,:,:,i_k_point),&
                       i_k_point)

              
          !aux_matr(:,:,:) =  aux_matr(:,:,:) + &
          !                    self_energy_omega_KS(:,:,:,i_k_point) + &
          !                    DMFT_ham_NAO(:,:,i_k_point) 

      if(allocated(KS_eigenvector_tmp))then
         deallocate(KS_eigenvector_tmp)
      endif


             enddo
          aux_matr(:,:,:) = aux_matr(:,:,:)*(1.d0/n_k_points)

!                call analy_continue_self_energy_dmft &
!                      (anacon_type, &
!                       1,n_states,nomega, &
!                       !1,n_basis,nomega, &
!                       n_max_par, &
!                       chemical_potential, &
!                       sigma_par_loc(:,:,:),omega, &
!                       aux_matr(:,:,:))
!                       !self_energy_freq(:,:,:))


                  filename = "spec_sc_new.dat"
                  call get_qp_spectrum_dmft &
                   ( anacon_type, &
                    n_max_par, n_low_state, n_states, &
                    nomega, omega, &
                    !sigma_par_loc, DMFT_ham_KS, &
                    sigma_par_loc, DMFT_ham_NAO, &
                    full_ovlp_matr, &
                    filename)
!                endif

              end subroutine qp_spectrum_dmft 
