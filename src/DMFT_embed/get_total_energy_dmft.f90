  subroutine get_total_energy_dmft(on_site_gf,&
                                   loc_self_enrg,&
                                   inner_loop_reiterations,&
                                   free_cluster_ham,&
                                   hartree_pot_LDA,&
                                   on_site_xc_matr,&
                                   hartree_pot_HF,&
                                   exchange_self_energy,&
                                   on_site_PBE_x_matr,& 
                                   on_site_PBE_c_matr,& 
                                   !full_k_xc_matr,&
                                   full_hamiltonian,&
                                   full_ovlp_matrix_sqrt,&
                                   !inv_full_ovlp_matrix_sqrt,&
                                   full_ovlp_matrix,&
                                   full_k_PBE_xc_matr,&
                                   dens_matr_LDA,&
                                   dens_matr,&
                                   !inv_full_ovlp_matrix,&
                                   embed_part_number_LDA,&
                                   self_energy_freq)




  use dimensions
  use runtime_choices
  use species_data
  use pbc_lists 
  use physics
  use prodbas
  use scgw_grid
  use constants
  use mpi_tasks
  use synchronize_mpi 
!  use hartree_fock
  use localized_basbas
  use gw_para
  use dmft_para
  use poles_fit 
  use gt
  use timing
  use localorb_io

   implicit none


! local parameters
   integer l,a,b,n, i_a, j_a,i_states, j_states 
   integer hamiltonian_size
   integer ovlp_matrix_size
   integer i_hamilton_size, i_xc_matrix_size, j_xc_matrix_size
   integer inner_loop_reiterations 
   integer i_matrix_size
   integer j_matrix_size
   integer i_spin
   integer i_state
   integer j_basis
   integer k_basis
   integer i_basis
   integer i_basis_1
   integer i_freq
   integer i_k_point, i_k, i
!   integer nomega
   character*2 iter
   character*17 filename_RE
   character*17 filename_OUT
   character*17 filename_IM
   complex*16, dimension(n_basis,n_basis,nomega), intent(in) :: on_site_gf 
   complex*16, dimension(n_basis,n_basis,n_spin,n_k_points), &
   intent(in) :: full_hamiltonian 
   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 

   real*8, dimension(n_basis,n_basis), intent(in) :: hartree_pot_HF
   real*8, dimension(n_basis,n_basis), intent(in) :: exchange_self_energy
   real*8, dimension(n_basis,n_basis), intent(in) :: loc_self_enrg
   real*8, dimension(n_basis,n_basis), intent(in) :: free_cluster_ham 
   real*8, dimension(n_basis,n_basis), intent(in) :: hartree_pot_LDA 
   real*8, dimension(n_basis,n_basis), intent(in) :: on_site_xc_matr 
   real*8, dimension(n_basis,n_basis), intent(in) :: on_site_PBE_x_matr 
   real*8, dimension(n_basis,n_basis), intent(in) :: on_site_PBE_c_matr 
   real*8 :: total_enrg_dmft
   real*8 :: total_X_enrg_dmft
   real*8 :: total_hartree_enrg_dmft
   real*8 :: total_kinetic_enrg_dmft
   real*8 :: embed_part_number_LDA


   complex*16, dimension(:,:), allocatable :: total_enrg_hartree_dmft
   real*8, dimension(:,:), allocatable :: total_enrg_fock_dmft
   complex*16, dimension(:,:), allocatable :: total_enrg_PBE_x_dmft
   complex*16, dimension(:,:), allocatable :: total_enrg_PBE_x
   complex*16, dimension(:,:), allocatable :: total_enrg_PBE_c_dmft
   complex*16, dimension(:,:), allocatable :: total_enrg_PBE_c_lda
   complex*16, dimension(:,:), allocatable :: E_tot_hartree_lda
   complex*16, dimension(:,:), allocatable :: E_tot_lda
   complex*16, dimension(:,:), allocatable :: E_tot_xc_non_loc


   real*8, dimension(:,:), allocatable:: on_site_dens_matr 
   real*8, dimension(:,:), allocatable:: hartree_pot_hybrid 
   real*8, dimension(n_basis,n_basis):: dens_matr 
   real*8, dimension(n_basis,n_basis):: dens_matr_LDA 
   real :: embed_part_number
   complex :: integral_test 
!   complex*16, dimension(n_basis,n_basis,n_spin,n_k_points) :: full_k_xc_matr
   complex*16, dimension(n_basis,n_basis,n_spin,n_k_points) :: &
   full_k_PBE_xc_matr
   complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) :: &
   full_ovlp_matrix_sqrt
   complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) :: & 
   full_ovlp_matrix
!   complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) :: inv_full_ovlp_matrix
!   complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) :: inv_full_ovlp_matrix_sqrt
   complex*16, dimension(:,:), allocatable:: free_ham
   complex*16, dimension(:,:), allocatable:: free_ham_dmft
   complex*16, dimension(:,:), allocatable:: kin_enrg_surround_matr
   complex*16, dimension(:,:), allocatable:: kin_enrg_dmft_matr
   real*8 :: kin_enrg_surround
   real*8 :: kin_enrg_dmft

      real*8, dimension(:,:,:), allocatable :: xc_matr

   logical  :: output

   integer, allocatable ::  myid_tmp(:)
   character*150 :: info_str

   real*8 :: sum_eigen_val
   real*8 :: total_hartree_enrg_LDA
   real*8 :: total_XC_enrg_LDA
   real*8 :: total_X_enrg_LDA
   complex*16 , dimension (n_basis,n_basis,nomega) :: self_energy_freq
   real*8 :: total_C_enrg_LDA
   real*8 :: trace_XC_pot_LDA
   real*8 :: total_X_enrg_surround
   real*8 :: total_X_enrg_surround_LDA
   real*8 :: total_C_enrg_surround
   real*8 :: total_C_enrg_dmft
   real*8 :: trace
   real*8 :: delta_hartree
   real*8 :: delta_kin
   real*8 :: delta_dens_matr
   real*8 :: c_energy
   real*8 :: gw_correlation_energy
   real*8, dimension(n_spin) :: x_energy
   real*8 :: xc_energy
   real*8 :: c_energy_lda
   real*8, dimension(n_spin) :: x_energy_lda


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------




      if (.not. allocated(E_tot_lda))then
          allocate(E_tot_lda(n_basis,n_basis))
       endif
      if (.not. allocated(E_tot_xc_non_loc))then
          allocate(E_tot_xc_non_loc(n_basis,n_basis))
       endif
      if (.not. allocated(E_tot_hartree_lda))then
          allocate(E_tot_hartree_lda(n_basis,n_basis))
       endif
      if (.not. allocated(total_enrg_fock_dmft))then
          allocate(total_enrg_fock_dmft(n_basis,n_basis))
       endif
      if (.not. allocated(total_enrg_hartree_dmft))then
          allocate(total_enrg_hartree_dmft(n_basis,n_basis))
       endif
      if (.not. allocated(total_enrg_PBE_x))then
          allocate(total_enrg_PBE_x(n_basis,n_basis))
       endif
      if (.not. allocated(total_enrg_PBE_x_dmft))then
          allocate(total_enrg_PBE_x_dmft(n_basis,n_basis))
       endif
      if (.not. allocated(total_enrg_PBE_c_dmft))then
          allocate(total_enrg_PBE_c_dmft(n_basis,n_basis))
       endif
      if (.not. allocated(total_enrg_PBE_c_lda))then
          allocate(total_enrg_PBE_c_lda(n_basis,n_basis))
       endif
      if (.not. allocated(free_ham_dmft))then
          allocate(free_ham_dmft(n_basis,n_basis))
       endif
      if (.not. allocated(free_ham))then
          allocate(free_ham(n_basis,n_basis))
       endif
      if (.not. allocated(kin_enrg_dmft_matr))then
          allocate(kin_enrg_dmft_matr(n_basis,n_basis))
       endif
      if (.not. allocated(kin_enrg_surround_matr))then
          allocate(kin_enrg_surround_matr(n_basis,n_basis))
       endif
      if (.not. allocated(hartree_pot_hybrid))then
          allocate(hartree_pot_hybrid(n_basis,n_basis))
       endif
       
      if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


!       if(myid.eq.0) then
!         on_site_dens_matr_temp(:,:) = 0.d0

       hartree_pot_hybrid(:,:) = (1.d0-hybrid_alpha)*(hartree_pot_LDA(:,:))+&
       hybrid_alpha*(hartree_pot_HF(:,:))

         !dens_matr_k_LDA(:,:) = 0.d0
       E_tot_lda(:,:) = (0.d0,0.d0)
       E_tot_xc_non_loc(:,:) = (0.d0,0.d0)
       E_tot_hartree_lda(:,:) = (0.d0,0.d0)
       total_enrg_hartree_dmft(:,:) = (0.d0,0.d0)
       total_enrg_fock_dmft(:,:) = (0.d0,0.d0)
!       total_enrg_PBE0_exx_dmft(:,:) = (0.d0,0.d0)
       total_enrg_PBE_x_dmft(:,:) = (0.d0,0.d0)
       total_enrg_PBE_x(:,:) = (0.d0,0.d0)
       total_enrg_PBE_c_dmft(:,:) = (0.d0,0.d0)
       total_enrg_PBE_c_lda(:,:) = (0.d0,0.d0)
       kin_enrg_dmft_matr(:,:) = (0.d0,0.d0)
       kin_enrg_surround_matr(:,:) = (0.d0,0.d0)


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------






! write(use_unit,*) 'dens_matr', gf_non_loc
          i_k = 0
     do i_k_point = 1, n_k_points

         free_ham(:,:) = full_hamiltonian(:,:,1,i_k_point) - &
                         1.d0*DCMPLX(hartree_pot_LDA(:,:)) - &
                         DCMPLX(on_site_xc_matr(:,:))
         free_ham_dmft(:,:) = full_hamiltonian(:,:,1,i_k_point) - &
                              2.d0*DCMPLX(hartree_pot_LDA(:,:)) - &
                              DCMPLX(on_site_xc_matr(:,:))




           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            DCMPLX(dens_matr_LDA(:,:)), n_basis, &
            free_ham(:,:), &
            n_basis,(1.d0,0.d0),&
            kin_enrg_surround_matr(:,:),n_basis)

           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            DCMPLX(dens_matr(:,:)), n_basis, &
            free_ham(:,:), &
            n_basis,(1.d0,0.d0),&
            kin_enrg_dmft_matr(:,:),n_basis)



!!
           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            DCMPLX(dens_matr(:,:)), n_basis, &
            full_hamiltonian(:,:,1,i_k_point), &
            n_basis,(1.d0,0.d0),&
            E_tot_lda(:,:),n_basis)

           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            DCMPLX(dens_matr_LDA(:,:)), n_basis, &
            !full_k_xc_matr(:,:,1,i_k_point), &
            DCMPLX(on_site_xc_matr(:,:)), &
            n_basis,(1.d0,0.d0),&
            E_tot_xc_non_loc(:,:),n_basis)

           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            DCMPLX(dens_matr_LDA(:,:)), n_basis,&
            DCMPLX(hartree_pot_LDA(:,:)),&
            n_basis,(0.d0,0.d0),&
            E_tot_hartree_lda(:,:),n_basis)

           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            DCMPLX(dens_matr(:,:)), n_basis, & 
            DCMPLX(hartree_pot_LDA(:,:)), &
            n_basis,(0.d0,0.d0),&
            total_enrg_hartree_dmft(:,:),n_basis)

           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            DCMPLX(dens_matr_LDA(:,:)), n_basis, &
            DCMPLX(on_site_PBE_x_matr(:,:)), &
            n_basis,(1.d0,0.d0),&
            total_enrg_PBE_x_dmft(:,:),n_basis)
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            DCMPLX(dens_matr(:,:)), n_basis, &
            DCMPLX(on_site_PBE_x_matr(:,:)), &
            n_basis,(1.d0,0.d0),&
            total_enrg_PBE_x(:,:),n_basis)

           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            DCMPLX(dens_matr(:,:)), n_basis, &
            DCMPLX(on_site_PBE_c_matr(:,:)), &
            n_basis,(1.d0,0.d0),&
            total_enrg_PBE_c_dmft(:,:),n_basis)

           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            DCMPLX(dens_matr_LDA(:,:)), n_basis, &
            DCMPLX(on_site_PBE_c_matr(:,:)), &
            n_basis,(1.d0,0.d0),&
            total_enrg_PBE_c_lda(:,:),n_basis)



!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

           call dgemm ('N','N', n_basis,&
            n_basis,n_basis,1.d0,&
            dens_matr(:,:), n_basis, &
            exchange_self_energy, &
            n_basis,0.d0,&
            total_enrg_fock_dmft(:,:),n_basis)



     enddo

       E_tot_xc_non_loc(:,:)= (1./n_k_points)*E_tot_xc_non_loc(:,:)
       E_tot_lda(:,:) = E_tot_lda(:,:)*(1./n_k_points)
       total_enrg_PBE_x_dmft(:,:) = &
       total_enrg_PBE_x_dmft(:,:)*(1./n_k_points)
       total_enrg_PBE_x(:,:) = total_enrg_PBE_x(:,:)*(1./n_k_points)
       total_enrg_PBE_c_dmft(:,:) =&
       total_enrg_PBE_c_dmft(:,:)*(1./n_k_points)
       total_enrg_PBE_c_lda(:,:) = total_enrg_PBE_c_lda(:,:)*(1./n_k_points)
       kin_enrg_dmft_matr(:,:) = kin_enrg_dmft_matr(:,:)*(1./n_k_points)
       kin_enrg_surround_matr(:,:) = &
       kin_enrg_surround_matr(:,:)*(1./n_k_points)

       gw_correlation_energy = 0.d0

         do i_basis = 1, n_basis, 1
           do j_basis  = 1, n_basis , 1
             do i_freq = 2, nomega, 1

               gw_correlation_energy = gw_correlation_energy +&
               self_energy_freq(i_basis, j_basis, i_freq) *&
               on_site_gf (i_basis, j_basis, i_freq)* womega (i_freq) &
               / (2.d0 * pi)

             enddo

             gw_correlation_energy = gw_correlation_energy +&
               self_energy_freq(i_basis, j_basis,1) *&
               on_site_gf (i_basis, j_basis,1 )* womega (1) &
               / (2.d0 * pi)/2.d0 

           enddo
         enddo

       gw_correlation_energy = gw_correlation_energy*(2.d0/n_spin)

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  total_enrg_dmft =0.d0
  total_X_enrg_dmft =0.d0
  total_hartree_enrg_dmft =0.d0
  total_kinetic_enrg_dmft =0.d0
  sum_eigen_val =0.d0
  total_hartree_enrg_LDA =0.d0
  total_XC_enrg_LDA =0.d0
  total_X_enrg_LDA=0.d0
  total_C_enrg_LDA=0.d0
  trace_XC_pot_LDA =0.d0
  total_X_enrg_surround =0.d0
  total_X_enrg_surround_LDA =0.d0
  total_C_enrg_surround =0.d0 
  total_C_enrg_dmft =0.d0 
  delta_dens_matr =0.d0
  kin_enrg_dmft =0.d0
  kin_enrg_surround =0.d0

  do i_a = 1, n_basis,1

    total_X_enrg_dmft = &
    total_X_enrg_dmft - total_enrg_fock_dmft(i_a,i_a)
    total_hartree_enrg_dmft = &
    total_hartree_enrg_dmft + total_enrg_hartree_dmft(i_a,i_a)
    total_hartree_enrg_LDA =&
    total_hartree_enrg_LDA + E_tot_hartree_lda(i_a,i_a)
    total_XC_enrg_LDA = &
    total_XC_enrg_LDA + E_tot_xc_non_loc(i_a,i_a)
    total_X_enrg_surround =&
    total_X_enrg_surround +  total_enrg_PBE_x_dmft(i_a,i_a)
    total_X_enrg_surround_LDA =&
    total_X_enrg_surround_LDA +  total_enrg_PBE_x(i_a,i_a)
    total_C_enrg_surround = &
    total_C_enrg_surround +  total_enrg_PBE_c_lda(i_a,i_a)
    total_C_enrg_dmft =&
    total_C_enrg_dmft +  total_enrg_PBE_c_dmft(i_a,i_a)
    kin_enrg_dmft = &
    kin_enrg_dmft + kin_enrg_dmft_matr(i_a,i_a) 
    kin_enrg_surround = &
    kin_enrg_surround + kin_enrg_surround_matr(i_a,i_a) 

    delta_dens_matr = delta_dens_matr + dens_matr(i_a,i_a)&
     - dens_matr_LDA(i_a,i_a)

  enddo

  do i_a = 1, n_basis,1
!   do j_a = 1, n_basis,1
    sum_eigen_val = sum_eigen_val + E_tot_lda(i_a,i_a)
!   enddo
  enddo

  delta_hartree = total_hartree_enrg_dmft - total_hartree_enrg_LDA
  delta_kin = kin_enrg_dmft - kin_enrg_surround

!---------------------------------------------------------------------------
!-------------------------------------------------------------------------

          call integrate_xc_energy ( &
                                     partition_tab,&
                                     rho,rho_gradient,&
                                     kinetic_density,&
                                     xc_energy,&
                                     x_energy,&
                                     c_energy,&
                                     x_energy_lda,&
                                     c_energy_lda&
                                     )


total_enrg_dmft =  total_energy     - 0.25*en_xc&
                                    + 0.25*(total_X_enrg_dmft)  &
                                    + 0.5*(en_pot_xc-2*total_XC_enrg_LDA) &
                                    + 0.25*(gw_correlation_energy)&
                                    + 1*(delta_hartree + delta_kin)


if(.false.) then
total_enrg_dmft =  total_energy     - 0.25*en_xc&
                                    + 0.25*(hybrid_alpha*total_X_enrg_dmft& 
                                    + (1.d0-hybrid_alpha)*1*x_energy(1))&
                                    + 0.25*(1*c_energy)&
                                    + 1*(delta_hartree + delta_kin)
endif


if(.false.) then
total_enrg_dmft =  2*sum_eigen_val  - 2*total_XC_enrg_LDA&
                                    - 0.5*(hartree_energy_free &
                                    + 2*hartree_delta_energy)  &
                                    + 1*delta_hartree&
                                    + total_X_enrg_dmft&
                                    + gw_correlation_energy
 endif


         write(info_str, '(A)') ''
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X)') &
              "DMFT-embedding Total energy components:"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Sum of eigenvalues                       :", 2*sum_eigen_val, " Ha", &
              2*sum_eigen_val * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Kin energy surround                      :", kin_enrg_surround, " Ha", &
              kin_enrg_surround * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Kin energy DMFT                          :", kin_enrg_dmft, " Ha", &
              kin_enrg_dmft * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| DELTA-Kin energy                         :", (delta_kin), " Ha", &
              (delta_kin) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| EXX energy correction from DMFT          :", total_X_enrg_dmft, " Ha", &
              total_X_enrg_dmft * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Hybrid-EXX energy correction from DMFT   :"&
              , hybrid_alpha*total_X_enrg_dmft+&
               (1.d0-hybrid_alpha)*total_X_enrg_surround, " Ha", &
              (hybrid_alpha*total_X_enrg_dmft+&
              (1.d0-hybrid_alpha)*total_X_enrg_surround) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| X energy from environement               :"&
              , total_X_enrg_surround_LDA, " Ha", &
              total_X_enrg_surround_LDA * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| DELTA-X energy from environement         :"&
              , 2*( hybrid_alpha*en_xc -  (1.d0-hybrid_alpha)*en_pot_xc), " Ha", &
              2*(hybrid_alpha*en_xc - (1.d0-hybrid_alpha)*en_pot_xc) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| embedding DELTA-X energy                 :"&
              , ((hybrid_alpha*total_X_enrg_dmft+&
                (1.d0-hybrid_alpha)*total_X_enrg_surround_LDA)-&
                total_X_enrg_surround_LDA), " Ha", &
              ((hybrid_alpha*total_X_enrg_dmft+&
              (1.d0-hybrid_alpha)*total_X_enrg_surround_LDA)-&
              total_X_enrg_surround_LDA) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| C energy from environement               :", 0.5*total_C_enrg_surround, " Ha", &
              total_C_enrg_surround*0.5 * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| C energy from DMFT                       :", 0.5*total_C_enrg_DMFT, " Ha", &
              total_C_enrg_DMFT *0.5* hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Delta C energy from environement         :"&
            , 0.5*(total_C_enrg_dmft-total_C_enrg_surround), " Ha", &
              0.5*(total_C_enrg_dmft-total_C_enrg_surround) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| XC potential from environement           :",  total_XC_enrg_LDA, &
              " Ha",  total_XC_enrg_LDA * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Hartree potential DMFT                   :",  total_hartree_enrg_dmft, &
              " Ha",  total_hartree_enrg_dmft * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Hartree potential KS-scf                 :",  total_hartree_enrg_LDA, &
              " Ha",  total_hartree_enrg_LDA * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| DELTA DENSITY from scDMFT                :", delta_dens_matr, &
               " Ha", delta_dens_matr * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| DMFT-embedding total energy              :", 2*total_enrg_dmft, &
              " Ha", 2*total_enrg_dmft * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| DELTA Hartree                            :", ( 1*delta_hartree ), &
              " Ha", ( 1*delta_hartree ) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "|  sum of eingenvalues                     :", (sum_eigen_val), &
              " Ha", (sum_eigen_val) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )


        if (myid.eq.0) then
         write(use_unit,*) " ----------------------------------------  "
         write(use_unit,*) "  The embedded total Energy Contributions  "
         write(use_unit,*) " ----------------------------------------  "
         write(use_unit,*) "   "
        endif



         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| embedding-XX energy Contribution         :"&
            , 0.25*(total_X_enrg_dmft+gw_correlation_energy-en_xc), " Ha", &
              0.25*(total_X_enrg_dmft+gw_correlation_energy-en_xc) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Hartree  Contribution                    :"&
              ,1*((delta_hartree )+(delta_kin)), &
              " Ha",1*(delta_hartree+1*(delta_kin)) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| XC  Correction                           :"&
              ,0.5*(en_pot_xc-2*total_XC_enrg_LDA), &
              " Ha",0.5*(en_pot_xc-2*total_XC_enrg_LDA) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )


        if (myid.eq.0) then
         write(use_unit,*) " ----------------------------------------  "
         write(use_unit,*) "  The embedded total Energy for Testing  "
         write(use_unit,*) " ----------------------------------------  "
         write(use_unit,*) "   "
        endif



         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Single Particle Energy                   :"&
              , 2*sum_eigen_val* hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Exchange Energy  DFT                     :"& 
              ,  2*total_XC_enrg_LDA* hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Exchange Energy                          :"&
              ,  total_X_enrg_dmft* hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| GW correlation Energy                    :"&
              , gw_correlation_energy* hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Hartree  Contribution                    :"&
              , 0.5*((delta_hartree )+(delta_kin)), &
              " Ha",0.5*(delta_hartree+1*(delta_kin)) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )


        if (myid.eq.0) then
         write(use_unit,*) " ----------------------------------------  "
         write(use_unit,*) "  The embedded total Energy for Testing  "
         write(use_unit,*) " ----------------------------------------  "
         write(use_unit,*) "   "
        endif



         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Single Particle Energy                   :"&
              , 2*sum_eigen_val* hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Exchange Energy  DFT                     :"&
              , 2*total_XC_enrg_LDA* hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Exchange Energy                          :"&
              ,  total_X_enrg_dmft* hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| GW correlation Energy                    :"&
              , gw_correlation_energy* hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Hartree energy from DFT density          :"&
              , total_hartree_enrg_LDA * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Hartree energy from GW density           :"&
              , total_hartree_enrg_dmft * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Hartree free Energy  DFT                 :"&
              , (0.5*hartree_energy_free ) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Delta Hartree Energy  DFT                :"&
              , (hartree_delta_energy) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Delta Hartree Energy  GW                 :"&
              , (delta_hartree) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Delta XC Energy  GW                      :"&
             , ( 2*total_XC_enrg_LDA-&
               (total_X_enrg_dmft+gw_correlation_energy)) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Galitskii-Migdal Total Energy            :"&
             , total_enrg_dmft * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )



       if(myid.eq.0) then

     output = .false.
      if(output) then
          open(57, file='dens_matr.dat')
            do i_a = 1, n_basis, 1
              
              write(57,*) i_a, &
                       (on_site_dens_matr(i_a,i_a)) 
            enddo
            do i_a = 1, n_basis, 1

              write(57,*) i_a, &
                       (on_site_dens_matr(1,i_a)) !,  embed_gf_FT(i_a,i_a,0)
            enddo
             write(57,*) 'particle numbers, embed and on_site = ', embed_part_number

  !        do i_a = 1, n_basis, 1
             !write(57,*) 'k_summed_overlap_matr = ', k_summed_overlap_matr!(i_a,i_a)
  !        enddo

!             write(57,*) 'womega = ', womega(1:10)!, lda_part_number
          close(57)
      endif


      endif


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!-----------------------------------------------------------------------
!--------------------------------------------------------------------

      if ( allocated(E_tot_lda))then
          deallocate(E_tot_lda)
       endif
      if ( allocated(E_tot_xc_non_loc))then
          deallocate(E_tot_xc_non_loc)
       endif
      if ( allocated(E_tot_hartree_lda))then
          deallocate(E_tot_hartree_lda)
       endif
      if (allocated(total_enrg_fock_dmft))then
          deallocate(total_enrg_fock_dmft)
       endif
      if ( allocated(total_enrg_hartree_dmft))then
          deallocate(total_enrg_hartree_dmft)
       endif
!      if (.not. allocated(total_enrg_PBE0_exx_dmft))then
!          allocate(total_enrg_PBE0_exx_dmft(n_basis,n_basis))
!       endif
!      if ( allocated(total_enrg_PBE_xc_dmft))then
!          deallocate(total_enrg_PBE_xc_dmft)
!       endif

      if ( allocated(ipiv))then
          deallocate(ipiv)
       endif

       if ( allocated(work))then
          deallocate(work)
       endif






!deallocate(on_site_dens_matr_temp)


  end subroutine get_total_energy_dmft 
