  subroutine get_particle_number(embed_gf,&
                              embed_part_number_LDA)                              

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
  use hartree_fock
  use localized_basbas
  use gw_para
  use poles_fit 
  use gt
  use timing

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
   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 

   real*8, dimension(:,:), allocatable:: embed_dens_matr 
   real*8, dimension(:,:), allocatable:: sin_anteil_temp 
   real*8, dimension(:,:), allocatable:: cos_anteil_temp 
   !complex*16, dimension(:,:), allocatable:: embed_dens_matr_ovlp 
   real*8, dimension(:,:), allocatable:: embed_dens_matr_temp_2
   real*8, dimension(:,:,:), allocatable:: embed_gf_time
!   complex*16, dimension(:,:), allocatable:: inv_k_summed_overlap_matr_sqrt 
!   complex*16, dimension(:,:), allocatable:: embed_dens_matr_temp 
   complex*16, dimension(:,:,:), allocatable:: ovlp_embed_gf 
   complex*16, dimension(:,:,:), allocatable:: ovlp_embed_gf_ovlp
   real :: cos_anteil 
   real :: sin_anteil 
   real :: embed_part_number
   real :: embed_part_number_test
   real :: embed_part_number_test_2
   complex :: integral_test 

  !     complex*16, dimension(n_basis,n_basis, n_k_points), intent(in) :: full_ovlp_matrix

   real*8 :: embed_part_number_LDA
   real*8 :: chemical_potential_new
   real*8 :: embed_part_number_new
   real :: sin_anteil_calibrate

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
      if (.not. allocated(embed_gf_time))then
          allocate(embed_gf_time(n_basis,n_basis,-ntau:ntau))
       endif
       
      if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif


!      if (.not. allocated(inv_k_summed_overlap_matr_sqrt))then
!          allocate(inv_k_summed_overlap_matr_sqrt(n_basis,n_basis))
!       endif
!       if (.not. allocated(embed_part_number_LDA))then
!          allocate(embed_part_number_LDA)
!       endif



!      if (.not. allocated(ovlp_embed_gf))then
!          allocate(ovlp_embed_gf(n_basis,n_basis,nomega!))
!       endif
!      if (.not. allocated(ovlp_embed_gf_ovlp))then
!          allocate(ovlp_embed_gf_ovlp(n_basis,n_basis,nomega))
!       endif





!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------



!embed_dens_matr_temp(:,:)= (1./(pi))*embed_dens_matr_temp(:,:)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!      do i_basis_1 = 1, n_basis, 1
!            write(use_unit,*)'embed_gf from gent_density',  i_basis_1, &
!             embed_gf(i_basis_1, i_basis_1,100)
!      enddo


!       if(myid.eq.0) then
         embed_dens_matr_temp_2(:,:) = 0.d0
         cos_anteil_temp(:,:) =0.d0
         sin_anteil_temp(:,:) =0.d0
         sin_anteil_calibrate =0.d0
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

!                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(4.1023E-5))&! for grid with 200 points
!                 -aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega(i_freq))*(4.1023E-5)))!*exp(-chemical_potential*(4.1023E-5)) ! for grid with 200 points
                  !(real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(5.1023E-4))&! for grid with 200 points
                 !-aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega(i_freq))*(5.1023E-4)))!*exp(-chemical_potential*(4.1023E-5)) ! for grid with 200 points
                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(3.3523E-4))&! for grid with 200 points
                 -aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega(i_freq))*(3.3523E-4)))!*exp(-chemical_potential*(4.1023E-5)) ! for grid with 200 points

                 ! (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(4.0923E-3))&! for grid with 200 points
                 !-aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega(i_freq))*(4.0923E-3)))!*exp(-chemical_potential*(4.1023E-5)) ! for grid with 200 points
!                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(1.00552E-3))&! for grid with 40 points
!                 -aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin(omega(i_freq)*(1.00552E-3))) ! for grid with 40 points


                  cos_anteil_temp(i_matrix_size,j_matrix_size)= &
                  cos_anteil_temp(i_matrix_size,j_matrix_size) + &
                  womega(i_freq)*&
                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(4.1023E-5)))!*(2.545E-4)))! for grid with 80 points

                  sin_anteil_temp(i_matrix_size,j_matrix_size)= &
                  sin_anteil_temp(i_matrix_size,j_matrix_size) + &
                  womega(i_freq)*&
                 ! (aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega(i_freq))*(2.545E-4)))! for grid with 80 points
                  aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega(i_freq))*(4.1023E-5)) ! for grid with 200 points


                 enddo
             enddo
                  sin_anteil_calibrate = &
                  sin_anteil_calibrate + &
                  womega(i_freq)*&
!                  (aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega(i_freq))*(2.545E-4)))! for grid with 80 points
                  !aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega(i_freq))*(4.1023E-5)) ! for grid with 200 points
                  !((sin(omega(i_freq)))*(1.d0/omega(i_freq))) ! for grid with 200 points
                  ((sin(omega(i_freq)*(5.1023E-4)))*(omega(i_freq)**(-1))) ! for grid with 200 points

!write(use_unit,*) 'sin_anteil_calibrate',(omega(i_freq)**(-1)),sin(omega(i_freq)), womega(i_freq),omega(i_freq),sin_anteil_calibrate
          enddo
!stop



!embed_dens_matr_temp_2(:,:)= (1./(pi))*embed_dens_matr_temp_2(:,:)

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

             call transform_G (embed_gf , n_basis, &
               n_basis, embed_gf_time(:,:,:))




embed_dens_matr(:,:) = (1./(pi))*embed_dens_matr_temp_2(:,:) 
!embed_dens_matr(:,:) = embed_gf_FT(:,:,0) 
!embed_dens_matr =embed_dens_matr_temp_2(:,:) 

cos_anteil_temp(:,:) = (1./(pi))*cos_anteil_temp(:,:)
sin_anteil_temp(:,:) = (1./(pi))*sin_anteil_temp(:,:)




 embed_part_number =0.d0
 cos_anteil = 0.d0
 sin_anteil = 0.d0
 do i_a = 1, n_basis,1
! do j_a = 1, n_basis,1

    !embed_part_number= embed_part_number + embed_dens_matr_ovlp(i_a,i_a)
    !embed_part_number= embed_part_number + (embed_dens_matr(i_a,i_a))
    embed_part_number= embed_part_number + (embed_gf_time(i_a,i_a,0))

    cos_anteil = cos_anteil+ cos_anteil_temp(i_a,i_a) 
    sin_anteil = sin_anteil+ sin_anteil_temp(i_a,i_a) 
     
enddo

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!if(inner_loop_reiterations .eq. 0) then
!if(number_reiterations .eq. 0) then
embed_part_number_LDA= embed_part_number
!endif
!endif

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------








! enddo
! enddo

!       do i_basis_1 = 1, n_basis, 1
!         do i_basis_2 = 1, i_basis_1, 1
!             write(use_unit,*)'embed_GF from dens_matr',  i_basis_1, &
!               embed_gf(i_basis_1, i_basis_1,1)
!         enddo
!       enddo


       if(myid.eq.0) then
!
              write(use_unit,*) 'particle number from Density-Matrix',  embed_part_number
              write(use_unit,*) 'LDA-particle number from Density-Matrix',  embed_part_number_LDA
!              write(use_unit,*) 'sin() fct. Calibration',  sin_anteil_calibrate
              write(use_unit,*) 'cos_anteil ',  cos_anteil 
              write(use_unit,*) 'sin_anteil ',  sin_anteil
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
  !        enddo

!             write(57,*) 'womega = ', womega(1:10)!, lda_part_number
          close(57)
      endif


      endif


!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------

  !    do i_basis_1 = 1, n_basis, 1
  !          write(use_unit,*)'omega_term from get_dens_matr',  i_basis_1, &
  !           omega_term(i_basis_1, i_basis_1,100)
  !    enddo

      do i_basis_1 = 1, n_basis, 1
  !          write(use_unit,*)'embed_gf from get_dens_matr',  i_basis_1, &
  !           embed_gf(i_basis_1, i_basis_1,100)
      enddo








      if ( allocated(embed_dens_matr_temp_2))then
          deallocate(embed_dens_matr_temp_2)
       endif
      if ( allocated(sin_anteil_temp))then
          deallocate(sin_anteil_temp)
       endif
      if ( allocated(cos_anteil_temp))then
          deallocate(cos_anteil_temp)
       endif

      if ( allocated(ipiv))then
          deallocate(ipiv)
       endif

       if ( allocated(work))then
          deallocate(work)
       endif




  end subroutine get_particle_number 
