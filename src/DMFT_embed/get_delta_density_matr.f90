  subroutine get_delta_density_matr(embed_gf,on_site_gf, &
        k_summed_overlap_matr, embed_dens_matr, inner_loop_reiterations, &
        inv_k_summed_overlap_matr, delta_dens_matr_temp) ! (KS_eigenvalue)

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
   integer i_k_point
!   integer nomega
   character*2 iter
   character*17 filename_RE
   character*17 filename_OUT
   character*17 filename_IM
!      real*8  :: womega(nomega) 
! subroutine actual parameters
   complex*16, dimension(:,:,:), allocatable :: delta_gf
   complex*16, dimension(n_basis,n_basis,nomega), intent(in) :: embed_gf 
   complex*16, dimension(n_basis,n_basis,nomega) :: on_site_gf 
   !complex*16, dimension(n_basis,n_basis,nomega), intent(in) :: on_site_gf 
   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 

   complex*16, dimension(n_basis,n_basis), intent(in) :: k_summed_overlap_matr
   real*8, dimension(n_basis,n_basis), intent(out):: embed_dens_matr 
   real*8, dimension(n_basis,n_basis), intent(out):: delta_dens_matr_temp 
   real*8, dimension(:,:), allocatable:: sin_anteil_temp 
   real*8, dimension(:,:), allocatable:: cos_anteil_temp 
   !complex*16, dimension(:,:), allocatable:: embed_dens_matr_ovlp 
   complex*16, dimension(:,:,:), allocatable:: ovlp_embed_gf
   complex*16, dimension(:,:,:), allocatable:: embed_gf_temp
   real*8, dimension(:,:), allocatable:: embed_dens_matr_temp_2
   complex*16, dimension(:,:), allocatable:: inv_k_summed_overlap_matr_sqrt 
!   complex*16, dimension(:,:), allocatable:: embed_dens_matr_temp 
   real :: cos_anteil 
   real :: sin_anteil 
   real :: embed_part_number
   real :: embed_part_number_test
   real :: embed_part_number_test_2
   complex :: integral_test 

       complex*16, dimension(n_basis,n_basis) ::  inv_k_summed_overlap_matr
  !     complex*16, dimension(n_basis,n_basis, n_k_points), intent(in) :: full_ovlp_matrix


   logical  :: output


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------



       if (.not. allocated(delta_gf))then
          allocate(delta_gf(n_basis,n_basis,nomega))
       endif
      if (.not. allocated(embed_dens_matr_temp_2))then
          allocate(embed_dens_matr_temp_2(n_basis,n_basis))
       endif
      if (.not. allocated(sin_anteil_temp))then
          allocate(sin_anteil_temp(n_basis,n_basis))
       endif
      if (.not. allocated(cos_anteil_temp))then
          allocate(cos_anteil_temp(n_basis,n_basis))
       endif
!      if (.not. allocated(delta_dens_matr_temp))then
!          allocate(delta_dens_matr_temp(n_basis,n_basis))
!       endif
       
      if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif


      if (.not. allocated(ovlp_embed_gf))then
          allocate(ovlp_embed_gf(n_basis,n_basis,nomega))
       endif
      if (.not. allocated(embed_gf_temp))then
          allocate(embed_gf_temp(n_basis,n_basis,nomega))
       endif
      if (.not. allocated(inv_k_summed_overlap_matr_sqrt))then
          allocate(inv_k_summed_overlap_matr_sqrt(n_basis,n_basis))
       endif








!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

if (inner_loop_reiterations.eq.0)then
on_site_gf(:,:,:) = (0.d0,0.d0)
endif




        delta_gf(:,:,:) = (0.d0,0.d0)

        do i_freq = 1, nomega
            do j_matrix_size = 1, n_basis, 1
              do i_matrix_size = 1, n_basis, 1
                 delta_gf(i_matrix_size,j_matrix_size,i_freq) = delta_gf(i_matrix_size,j_matrix_size,i_freq)+ &
                 (on_site_gf(i_matrix_size,j_matrix_size,i_freq) - embed_gf(i_matrix_size,j_matrix_size,i_freq))
                 enddo
             enddo
          enddo
do j_matrix_size = 1, n_basis, 1
write(use_unit,*)"delta_gf", on_site_gf(j_matrix_size,j_matrix_size,1), embed_gf(j_matrix_size,j_matrix_size,1)
enddo
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!   call  transform_G(embed_gf, n_basis, n_basis, embed_gf_FT)

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!if (inner_loop_reiterations.eq.0)then
!on_site_gf(:,:,:) = (0.d0,0.d0)
!endif

!       if(myid.eq.0) then
         delta_dens_matr_temp(:,:) = 0.d0

        do i_freq = 1, nomega
            do j_matrix_size = 1, n_basis, 1
              do i_matrix_size = 1, n_basis, 1
                  delta_dens_matr_temp(i_matrix_size,j_matrix_size) = &
                  delta_dens_matr_temp(i_matrix_size,j_matrix_size) + &
                  womega(i_freq)*&
                  (real(delta_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(4.1023E-5))&
                 -aimag(delta_gf(i_matrix_size,j_matrix_size, i_freq))*sin(omega(i_freq)*(4.1023E-5)))
                 enddo
             enddo
          enddo

do j_matrix_size = 1, n_basis, 1
write(use_unit,*)"delta_dens_matr", delta_dens_matr_temp(j_matrix_size,j_matrix_size)
enddo

!stop
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
                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(4.1023E-5))&! for grid with 200 points
                 -aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin(omega(i_freq)*(4.1023E-5))) ! for grid with 200 points
!                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(1.00552E-3))&! for grid with 40 points
!                 -aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin(omega(i_freq)*(1.00552E-3))) ! for grid with 40 points


                  cos_anteil_temp(i_matrix_size,j_matrix_size)= &
                  cos_anteil_temp(i_matrix_size,j_matrix_size) + &
                  womega(i_freq)*&
                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*cos((omega(i_freq))*(0)))!*(2.545E-4)))! for grid with 80 points

                  sin_anteil_temp(i_matrix_size,j_matrix_size)= &
                  sin_anteil_temp(i_matrix_size,j_matrix_size) + &
                  womega(i_freq)*&
!                  (aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin((omega(i_freq))*(2.545E-4)))! for grid with 80 points
                  aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*sin(omega(i_freq)*(4.1023E-5)) ! for grid with 200 points


                 enddo
             enddo
          enddo




embed_dens_matr_temp_2(:,:)= (1./(pi))*embed_dens_matr_temp_2(:,:)
delta_dens_matr_temp(:,:)= (1./(pi))*delta_dens_matr_temp(:,:)

cos_anteil_temp(:,:) = (1./(pi))*cos_anteil_temp(:,:)
sin_anteil_temp(:,:) = (1./(pi))*sin_anteil_temp(:,:)


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------





!embed_dens_matr = (1./(pi))*embed_dens_matr_temp(:,:) 
embed_dens_matr =embed_dens_matr_temp_2(:,:) 




 embed_part_number =0.d0
 cos_anteil = 0.d0
 sin_anteil = 0.d0
 do i_a = 1, n_basis,1

    !embed_part_number= embed_part_number + embed_dens_matr_ovlp(i_a,i_a)
    embed_part_number= embed_part_number + (embed_dens_matr(i_a,i_a))

    cos_anteil = cos_anteil+ cos_anteil_temp(i_a,i_a) 
    sin_anteil = sin_anteil+ sin_anteil_temp(i_a,i_a) 
     

 enddo



       if(myid.eq.0) then
!
!          do i_a = 1, n_basis, 1
!              write(use_unit,*) i_a,(embed_dens_matr(i_a,i_a)), embed_gf_FT(i_a,i_a,0) 
!          enddo



              write(use_unit,*) 'particle number',  embed_part_number
              write(use_unit,*) 'cos_anteil ',  cos_anteil 
              write(use_unit,*) 'sin_anteil ',  sin_anteil 


     output = .false.
      if(output) then
          open(57, file='dens_matr.dat')
            do i_a = 1, n_basis, 1
              
              write(57,*) i_a, &
                       (embed_dens_matr(i_a,i_a)) 

                      ! k_summed_overlap_matr
                      ! real(embed_dens_matr(i_a,i_a)), &
                      ! aimag(embed_dens_matr(i_a,i_a))
                      ! real( prod(:,:,1) ),
                      ! real( prod(i_a,i_a,2)) , real( prod(i_a,i_a,1)) 
            enddo
            do i_a = 1, n_basis, 1

              write(57,*) i_a, &
                       (embed_dens_matr(1,i_a)) !,  embed_gf_FT(i_a,i_a,0)
            enddo
             write(57,*) 'particle numbers, embed and on_site = ', embed_part_number

  !        do i_a = 1, n_basis, 1
             write(57,*) 'k_summed_overlap_matr = ', k_summed_overlap_matr!(i_a,i_a)
  !        enddo

!             write(57,*) 'womega = ', womega(1:10)!, lda_part_number
          close(57)
      endif

      !do i_a =1 , nomega
      !write(use_unit,*) i_a, omega(i_a), womega (i_a)!, lda_part_number
      !enddo


      endif


!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------







!deallocate(embed_dens_matr_temp)
!deallocate(embed_gf_FT)





  end subroutine get_delta_density_matr 
