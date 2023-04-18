  subroutine update_chemical_pot_dmft_cluster(embed_gf,& !k_summed_overlap_matr,&
                                      inner_loop_reiterations, &
                                      number_reiterations,&
                                      !inv_k_summed_overlap_matr,&
                                      full_ovlp_matrix_sqrt,&
                                      inv_full_ovlp_matrix_sqrt,&
                                      full_ovlp_matrix,free_cluster_ham,&
                                      hybrid_func,loc_self_enrg,&
                                      hartree_pot_LDA,&
                                      embed_part_number_LDA,&
                                      chemical_potential_new,&
                                      k_summed_overlap_matr,&
                                      self_energy_freq) ! (KS_eigenvalue)

  use dimensions
  use runtime_choices
  use species_data
  use pbc_lists 
  use physics
  use prodbas
  use scgw_grid
  use constants
  use mpi_tasks
!  use synchronize_mpi 
!  use hartree_fock
!  use localized_basbas
  use gw_para
  use poles_fit 
  use gt
  use timing
  use localorb_io, only: use_unit

   implicit none


! local parameters
   integer l,a,b,n, i_a, j_a 
   integer hamiltonian_size
   integer ovlp_matrix_size
   integer i_hamilton_size, i_xc_matrix_size, j_xc_matrix_size
   integer inner_loop_reiterations , number_reiterations
   integer i_matrix_size
   integer j_matrix_size
   integer i_spin
   integer i_state
   integer j_basis
   integer k_basis
   integer i_basis
   integer i_basis_1
   integer i_freq
   integer i_k_point
!   integer nomega
   character*2 iter
   character*17 filename_RE
   character*17 filename_OUT
   character*17 filename_IM
!      real*8  :: womega(nomega) 
! subroutine actual parameters
   complex*16, dimension(n_basis,n_basis,nomega), intent(in) :: embed_gf 
   complex*16, dimension(n_basis,n_basis,nomega)  :: self_energy_freq
   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 

   real*8, dimension(n_basis,n_basis) :: k_summed_overlap_matr
   real*8, dimension(:,:), allocatable:: embed_dens_matr 
   real*8, dimension(:,:), allocatable:: sin_anteil_temp 
   real*8, dimension(:,:), allocatable:: cos_anteil_temp 
   real*8, dimension(:,:), allocatable:: embed_dens_matr_temp_2
   real :: cos_anteil 
   real :: sin_anteil 
   real :: embed_part_number
   real :: embed_part_number_test
   real :: embed_part_number_test_2
   complex :: integral_test 
       complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) ::&
       full_ovlp_matrix_sqrt
       complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) ::&
       inv_full_ovlp_matrix_sqrt
       complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) ::&
       full_ovlp_matrix
       complex*16, dimension(n_basis,n_basis,nomega), intent(in) :: &
       hybrid_func

!       real*8, dimension(n_basis,n_basis) ::  inv_k_summed_overlap_matr
       real*8, dimension(n_basis,n_basis) ::  free_cluster_ham
       real*8, dimension(n_basis,n_basis) ::  loc_self_enrg
       real*8, dimension(n_basis,n_basis) ::  hartree_pot_LDA
  !     complex*16, dimension(n_basis,n_basis, n_k_points), intent(in) :: full_ovlp_matrix

   real*8  :: embed_part_number_LDA
   real*8 :: chemical_potential_new
   real*8 :: embed_part_number_new

   logical  :: output


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------



      if (.not. allocated(embed_dens_matr_temp_2))then
          allocate(embed_dens_matr_temp_2(n_basis,n_basis))
       endif
      if (.not. allocated(sin_anteil_temp))then
          allocate(sin_anteil_temp(n_basis,n_basis))
       endif
      if (.not. allocated(cos_anteil_temp))then
          allocate(cos_anteil_temp(n_basis,n_basis))
       endif
      if (.not. allocated(embed_dens_matr))then
          allocate(embed_dens_matr(n_basis,n_basis))
       endif
       
      if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif




!       if(myid.eq.0) then
         embed_dens_matr_temp_2(:,:) = 0.d0
         cos_anteil_temp(:,:) =0.d0
         sin_anteil_temp(:,:) =0.d0
        do i_freq = 1, nomega
            do j_matrix_size = 1, n_basis, 1
              do i_matrix_size = 1, n_basis, 1
                  embed_dens_matr_temp_2(i_matrix_size,j_matrix_size) = &
                  embed_dens_matr_temp_2(i_matrix_size,j_matrix_size) + &
                  womega(i_freq)*&
                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*&
                  cos((omega(i_freq))*(4.1023E-5))&! for grid with 200 points
                 -aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*&
                  sin((omega(i_freq))*(4.1023E-5)))!*exp(-chemical_potential*(4.1023E-5)) ! for grid with 200 points


                  cos_anteil_temp(i_matrix_size,j_matrix_size)= &
                  cos_anteil_temp(i_matrix_size,j_matrix_size) + &
                  womega(i_freq)*&
                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*&
                  cos((omega(i_freq))*(4.1023E-5)))!*(2.545E-4)))! for grid with 80 points

                  sin_anteil_temp(i_matrix_size,j_matrix_size)= &
                  sin_anteil_temp(i_matrix_size,j_matrix_size) + &
                  womega(i_freq)*&
                  aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*&
                  sin((omega(i_freq))*(4.1023E-5)) ! for grid with 200 points


                 enddo
             enddo
          enddo




!embed_dens_matr_temp_2(:,:)= (1./(pi))*embed_dens_matr_temp_2(:,:)

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------





embed_dens_matr(:,:) = (1./(pi))*embed_dens_matr_temp_2(:,:) 
!embed_dens_matr =embed_dens_matr_temp_2(:,:) 

cos_anteil_temp(:,:) = (1./(pi))*cos_anteil_temp(:,:)
sin_anteil_temp(:,:) = (1./(pi))*sin_anteil_temp(:,:)



 embed_part_number =0.d0
 cos_anteil = 0.d0
 sin_anteil = 0.d0
 do i_a = 1, n_basis,1
! do j_a = 1, n_basis,1

    !embed_part_number= embed_part_number + embed_dens_matr_ovlp(i_a,i_a)
    embed_part_number= embed_part_number + (embed_dens_matr(i_a,i_a))

    cos_anteil = cos_anteil+ cos_anteil_temp(i_a,i_a) 
    sin_anteil = sin_anteil+ sin_anteil_temp(i_a,i_a) 
     
enddo

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

call  update_chemical_pot_cluster(inv_full_ovlp_matrix_sqrt, &
                                  full_ovlp_matrix_sqrt, &
                                  full_ovlp_matrix,&
                                  free_cluster_ham,&
                                  hybrid_func,&
                                  loc_self_enrg,&
                                  hartree_pot_LDA,&
                                  KS_eigenvalue(:,:,:),&
                                  embed_part_number_LDA,&
                                  .false.,&
                                  chemical_potential_new,&
                                  embed_part_number_new,&
                                  k_summed_overlap_matr,&
                                  self_energy_freq)


       if(myid.eq.0) then
              write(use_unit,*) 'chemical_pot NEW  = '&
                          , chemical_potential_new*hartree, 'eV'

       endif



!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------





      if ( allocated(embed_dens_matr_temp_2))then
          deallocate(embed_dens_matr_temp_2)
       endif
      if ( allocated(sin_anteil_temp))then
          deallocate(sin_anteil_temp)
       endif
      if ( allocated(cos_anteil_temp))then
          deallocate(cos_anteil_temp)
       endif
      if ( allocated(embed_dens_matr))then
          deallocate(embed_dens_matr)
       endif

      if ( allocated(ipiv))then
          deallocate(ipiv)
       endif

       if ( allocated(work))then
          deallocate(work)
       endif




  end subroutine update_chemical_pot_dmft_cluster 
