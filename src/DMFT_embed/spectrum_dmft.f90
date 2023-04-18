            subroutine spectrum_dmft (new_k_ham,&
                                      spectrum,&
                                      aux_omegamax,&
                                      aux_omega, aux_womega,& 
                                      aux_nomega,&
                                      green_fn_par,&
                                      i_k_point,&
                                      diag_ham_k,&
                                      DMFT_dos_k,&
                                      ovlp_matr)



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
              use arch_specific, only: arch_erf
              use localorb_io, only: OL_high, output_priority
              use pbc_lists, only: k_weights

              implicit none


            integer number_reiterations
            complex*16 new_k_ham (n_basis, n_basis)
            complex*16 ovlp_matr (n_basis, n_basis)
            complex*16, dimension(:,:,:), allocatable:: diagonal_green_fn
            complex*16, dimension(n_max_par, n_basis, n_basis):: green_fn_par
            complex*16, dimension(n_states, n_states):: diag_ham_k 

            integer i_a, i_freq,n ,l
            integer i_k_point,i_e
            integer i_spin, i_state
            character*2 iter
            real*8 new_chem_pot
            integer   n_nonsingular

            real*8  :: de,  E1, E2
            real*8  :: en

            integer :: output_priority_old

      integer aux_nomega
      real*8  aux_omegamax
      real*8, dimension(-aux_nomega:aux_nomega) :: spectrum
      real*8 , dimension (aux_nomega) :: aux_omega
      real*8 , dimension (aux_nomega) :: aux_womega
      real*8, allocatable :: dmft_eigenvalues(:)
      complex*16, dimension(:,:,:), allocatable:: dmft_eigenvectors
      complex*16, dimension(:), allocatable:: ham_triangle
      complex*16, dimension(:), allocatable:: ovlp_triangle
!  real*8, dimension(:), allocatable :: KS_dos
      real*8, dimension(dos_n_en_points) :: DMFT_dos_k



              if(.not.allocated(diagonal_green_fn))then
                allocate(diagonal_green_fn(n_states, n_states, nomega))
!                allocate(diagonal_green_fn(n_basis, n_basis, nomega))
              endif


              if(.not.allocated(dmft_eigenvectors))then
                !allocate(dmft_eigenvectors(n_basis, n_basis,n_k_points))
                allocate(dmft_eigenvectors(n_states, n_states,n_k_points))
              endif

              if(.not.allocated(dmft_eigenvalues))then
                !allocate(dmft_eigenvalues(n_basis))
                allocate(dmft_eigenvalues(n_states))
              endif

allocate(ham_triangle   (n_basis*(n_basis+1)/2))
allocate(ovlp_triangle   (n_basis*(n_basis+1)/2))


                    n = 0
                        do l = 1, n_basis

                           ham_triangle(n+1:n+l) = new_k_ham(1:l,l)
                           ham_triangle(n+1:n+l-1) = new_k_ham(l,1:l-1)

                           ovlp_triangle(n+1:n+l) = ovlp_matr(1:l,l)
                           ovlp_triangle(n+1:n+l-1) = ovlp_matr(l,1:l-1)

                           n = n+l
                        enddo

                output_priority_old = output_priority
                output_priority = OL_high


                   call improve_complex_eigenfunctions &
                        ( ovlp_triangle, ham_triangle(:),  &
                        dmft_eigenvalues, dmft_eigenvectors (:,:,1),1)

                output_priority = output_priority_old





              diagonal_green_fn(:,:,:) = (0.d0, 0.d0)

!if(i_k_point .eq. 1) then

             diag_ham_k(:,:) = (0.d0, 0.d0)

              do i_a =1, n_states

       diag_ham_k(i_a,i_a) = diag_ham_k(i_a,i_a) + &
                       DCMPLX(dmft_eigenvalues(i_a))

              do i_freq =1, nomega

      diagonal_green_fn(i_a,i_a,i_freq) = diagonal_green_fn(i_a,i_a,i_freq) +&
                                             1./(((0.d0,1.d0)*omega(i_freq)+ &
                                                       chemical_potential) - &
                                        dmft_eigenvalues(i_a)+(0.d0,0.0005d0))


              enddo
              enddo

       de= (dos_high_energy - dos_low_energy)/dble(dos_n_en_points-1)
       DMFT_dos_k =0d0
           do i_e = 1, dos_n_en_points
              en= dos_low_energy + dble(i_e-1)*de
              E1 = en - de/2d0
              E2 = en + de/2d0
              do i_state = 1, n_states, 1
                 DMFT_dos_k(i_e) =  DMFT_dos_k (i_e) + k_weights(i_k_point)&
                       *(arch_erf((E2-(dmft_eigenvalues(i_state))*hartree)*&
                                                  one_over_sqrt2/dos_alpha)&
                        -arch_erf((E1-(dmft_eigenvalues(i_state))*hartree)*&
                                                  one_over_sqrt2/dos_alpha)&
                                                                  )/(2d0*dE)
              enddo
           enddo




              call analy_continue_green_fn_dmft &
                  (anacon_type,&
                   nomega, &
                   n_max_par, &
                   green_fn_par,omega, &
                   diagonal_green_fn, n_states, i_k_point)

              call get_spectrum_k_dmft (anacon_type, &
                                       green_fn_par, &
                                          n_max_par, &
                                              omega, &
                                             nomega, & 
                                            n_basis, &
                                           spectrum, &
                                       aux_omegamax, &
                                          aux_omega, &
                                         aux_womega, &
                                         aux_nomega)

      !            enddo ! I SPIN

              if(allocated(diagonal_green_fn))then
                deallocate(diagonal_green_fn)
              endif

          end subroutine spectrum_dmft
