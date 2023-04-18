  subroutine update_chemical_pot_dmft_spectrum(embed_gf,& !k_summed_overlap_matr,&
                                      inner_loop_reiterations, &
                                      number_reiterations,&
                                      !inv_k_summed_overlap_matr,&
                                      omega_term,full_ovlp_matrix_sqrt,&
                                      inv_full_ovlp_matrix_sqrt,&
                                      full_ovlp_matrix,free_cluster_ham,&
                                      hybrid_func,loc_self_enrg,&
                                      hartree_pot_LDA,&
                                      embed_part_number_LDA,&
                                      chemical_potential_new) ! (KS_eigenvalue)

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
   complex*16, dimension(n_basis,n_basis,nomega), intent(out) :: omega_term 
   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 

!   real*8, dimension(n_basis,n_basis), intent(in) :: k_summed_overlap_matr
   real*8, dimension(:,:), allocatable:: embed_dens_matr 
   real*8, dimension(:,:), allocatable:: sin_anteil_temp 
   real*8, dimension(:,:), allocatable:: cos_anteil_temp 
   !complex*16, dimension(:,:), allocatable:: embed_dens_matr_ovlp 
!   complex*16, dimension(:,:,:,:), allocatable:: aux_embed_gf
!   complex*16, dimension(:,:,:,:), allocatable:: embed_gf_temp_k
   real*8, dimension(:,:), allocatable:: embed_dens_matr_temp_2
!   complex*16, dimension(:,:), allocatable:: inv_k_summed_overlap_matr_sqrt 
!   complex*16, dimension(:,:), allocatable:: embed_dens_matr_temp 
   real :: cos_anteil 
   real :: sin_anteil 
   real :: embed_part_number
   real :: embed_part_number_test
   real :: embed_part_number_test_2
   complex :: integral_test 
       complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) :: full_ovlp_matrix_sqrt
       complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) :: inv_full_ovlp_matrix_sqrt
       complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) :: full_ovlp_matrix
       complex*16, dimension(n_basis,n_basis,nomega), intent(in) :: hybrid_func

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


!      if (.not. allocated(aux_embed_gf))then
!          allocate(aux_embed_gf(n_basis,n_basis,nomega,n_k_points))
!       endif
!      if (.not. allocated(inv_k_summed_overlap_matr_sqrt))then
!          allocate(inv_k_summed_overlap_matr_sqrt(n_basis,n_basis))
!       endif
!       if (.not. allocated(embed_part_number_LDA))then
!          allocate(embed_part_number_LDA)
!       endif








!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!         inv_k_summed_overlap_matr_sqrt(:,:) = k_summed_overlap_matr(:,:)

!     call power_genmat_lapack(n_basis,inv_k_summed_overlap_matr_sqrt(:,:), -1, &
!                             safe_minimum, 1.d-5 , '')

    !     k_summed_overlap_matr_sqrt(:,:) = k_summed_overlap_matr(:,:)

    ! call power_genmat_lapack(n_basis,k_summed_overlap_matr_sqrt(:,:), 0.5, &
    !                         safe_minimum, 1.d-5 , '')




!------------------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------------------

!aux_embed_gf (:,:,:,:)= (0.d0,0.d0)
!embed_gf_temp_k(:,:,:,:)= (0.d0,0.d0)
!embed_gf_temp(:,:,:) =(0.d0,0.d0)

!!do i_k_point =1 , n_k_points
!do i_freq = 1,nomega,1

!           call zgemm ('N','N', n_basis,&
!             n_basis,n_basis,(1.d0,0.d0),&
!             DCMPLX(k_summed_overlap_matr(:,:)), n_basis, &
!             embed_gf(:,:,i_freq), &
!             n_basis,(0.d0,0.d0), &
!             embed_gf_temp (:,:, i_freq) , n_basis )

!          call zgemm ('N','N', n_basis,&
!             n_basis,n_basis,(1.d0,0.d0),&
!             aux_embed_gf(:,:,i_freq,i_k_point), n_basis, &
!             inv_full_ovlp_matrix_sqrt(:,:,i_k_point), &
!             n_basis,(0.d0,0.d0), &
!             embed_gf_temp_k (:,:,i_freq,i_k_point) , n_basis )

!enddo

!     embed_gf_temp(:,:,i_freq) = embed_gf_temp(:,:,i_freq) + embed_gf_temp_k (:,:,i_freq,i_k_point)


!enddo
!enddo
!embed_gf_temp(:,:,:) = embed_gf_temp(:,:,:)*(1./n_k_points)
!------------------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


!       if(myid.eq.0) then

!embed_dens_matr_temp(:,:)= (1./(pi))*embed_dens_matr_temp(:,:)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------



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
!                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(0))&!*(2.545E-4))&! for grid with 80 points
!                 -aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin(omega(i_freq)*(2.545E-4))) ! for grid with 80 points
!                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega_term(i_matrix_size,j_matrix_size,i_freq))*(4.1023E-5))&! for grid with 200 points
!                 -aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega_term(i_matrix_size,j_matrix_size,i_freq))*(4.1023E-5)))*exp(-chemical_potential*(4.1023E-5)) ! for grid with 200 points
                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(4.1023E-5))&! for grid with 200 points
                 -aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega(i_freq))*(4.1023E-5)))!*exp(-chemical_potential*(4.1023E-5)) ! for grid with 200 points
!                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(1.00552E-3))&! for grid with 40 points
!                 -aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin(omega(i_freq)*(1.00552E-3))) ! for grid with 40 points


                  cos_anteil_temp(i_matrix_size,j_matrix_size)= &
                  cos_anteil_temp(i_matrix_size,j_matrix_size) + &
                  womega(i_freq)*&
                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(4.1023E-5)))!*(2.545E-4)))! for grid with 80 points

                  sin_anteil_temp(i_matrix_size,j_matrix_size)= &
                  sin_anteil_temp(i_matrix_size,j_matrix_size) + &
                  womega(i_freq)*&
!                  (aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega(i_freq))*(2.545E-4)))! for grid with 80 points
                  aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega(i_freq))*(4.1023E-5)) ! for grid with 200 points


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
!if(inner_loop_reiterations .eq. 0) then
!if(number_reiterations .eq. 0) then

!embed_part_number_LDA= embed_part_number
!endif
!endif

!if(inner_loop_reiterations .ge. 1) then

!write(use_unit,*) inner_loop_reiterations, number_reiterations
!write(use_unit,*) 'embed_part_number_LDA= ', embed_part_number_LDA
!write(use_unit,*) 'embed_part_number from update= ', embed_part_number

!endif
!if(inner_loop_reiterations .eq. 2) stop

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!if(inner_loop_reiterations.eq.0) chemical_potential_new = chemical_potential
!chemical_potential_new = chemical_potential
!if(inner_loop_reiterations.ge.1) then
call  update_chemical_pot_cluster(inv_full_ovlp_matrix_sqrt, full_ovlp_matrix_sqrt, full_ovlp_matrix,&
                          free_cluster_ham, hybrid_func,omega_term,loc_self_enrg,hartree_pot_LDA,&
                          KS_eigenvalue(:,:,:), embed_part_number_LDA, .false., chemical_potential_new,embed_part_number_new)
!                          KS_eigenvalue(:,:,1), 1.d0, .false., chemical_potential_new,embed_part_number_new)


!write(use_unit,*) 'embed_part_number_LDA= ', embed_part_number_LDA
       if(myid.eq.0) then

              write(use_unit,*) 'chemical_pot NEW  = ', chemical_potential_new*hartree, 'eV'

       endif


!endif
!       if(myid.eq.0) then

!do i_basis =1 , n_basis,1
!!write(use_unit,*) "omega_term", inner_loop_reiterations, omega_term(i_basis,i_basis,100)
!write(use_unit,*) "hartree_pot_LDA",i_basis, hartree_pot_LDA(i_basis,i_basis) 
!enddo
!endif
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------









! enddo
! enddo

       if(myid.eq.0) then
!
      !        write(use_unit,*) 'particle number_OLD (up-date chem_pot)',  embed_part_number
      !        write(use_unit,*) 'particle number_NEW (up-date chem_pot)',  embed_part_number_new
      !        write(use_unit,*) 'cos_anteil ',  cos_anteil 
      !        write(use_unit,*) 'sin_anteil ',  sin_anteil 
!stop

     output = .false.
      if(output) then
          open(57, file='dens_matr.dat')
            do i_a = 1, n_basis, 1
              
              write(57,*) i_a, &
                       (embed_dens_matr(i_a,i_a)) 

            enddo
            do i_a = 1, n_basis, 1

              write(57,*) i_a, &
                       (embed_dens_matr(1,i_a)) !,  embed_gf_FT(i_a,i_a,0)
            enddo
             write(57,*) 'particle numbers, embed and on_site = ', embed_part_number

  !        do i_a = 1, n_basis, 1
!             write(57,*) 'k_summed_overlap_matr = ', k_summed_overlap_matr!(i_a,i_a)
  !        enddo

!             write(57,*) 'womega = ', womega(1:10)!, lda_part_number
          close(57)
      endif


      endif


!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------








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




  end subroutine update_chemical_pot_dmft_spectrum 
