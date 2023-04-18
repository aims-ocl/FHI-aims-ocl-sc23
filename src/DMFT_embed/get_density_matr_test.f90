  subroutine get_density_matr_test(embed_gf, embed_dens_matr) ! (KS_eigenvalue)

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
!   integer nomega
   character*2 iter
   character*17 filename_RE
   character*17 filename_OUT
   character*17 filename_IM
!      real*8  :: womega(nomega) 
! subroutine actual parameters
   real*8, dimension(:,:,:), allocatable :: embed_gf_FT 
   complex*16, dimension(n_basis,n_basis,nomega) :: embed_gf 
   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 

   real*8, dimension(n_basis,n_basis):: embed_dens_matr 
   real*8, dimension(:,:), allocatable:: embed_dens_matr_temp 
   real :: embed_part_number

       real*8, dimension(:,:), allocatable ::  density_matrix
       !real*8, dimension(n_spin) ::  n_homo


   logical  :: output

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
nomega=n_freq

write(use_unit,*) 'n_homo', n_homo
write(use_unit,*) 'nomega from get_dens_matr', nomega      

       if (.not. allocated(embed_gf_FT))then
          allocate(embed_gf_FT(n_basis,n_basis,-n_freq:n_freq))
       endif

       if (.not. allocated(embed_dens_matr_temp))then
          allocate(embed_dens_matr_temp(n_basis,n_basis))
       endif
       
      if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif




!!!!!!!!!!!!! Frequency - grid definition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (.not. allocated(wtau))then
          allocate(wtau(n_freq))
       endif

       if (.not. allocated(tau))then
          allocate(tau(n_freq))
       endif

       if (.not. allocated(womega))then
          allocate(womega(n_freq))
       endif

       if (.not. allocated(omega_grid))then
          allocate(omega_grid(n_freq))
       endif

     call tf_ini_trans(n_freq,n_freq, taumax,omegamax, &
                                      tau,omega_grid, &
                                      wtau,womega,.true.)

!!!!!!!!!!!!! Frequency - grid definition END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !  call  transform_G(embed_gf, n_basis, n_basis, embed_gf_FT)



!write(use_unit,*) 'embed_gf_FT', embed_gf_FT(1,1,0)



       if(myid.eq.0) then
         embed_dens_matr_temp(:,:) = 0.d0
 
        do i_freq = 1, n_freq 
            do j_matrix_size = 1, n_basis, 1 
              do i_matrix_size = 1, n_basis, 1 
!                do i_freq = 1, nomega 
                  embed_dens_matr_temp(i_matrix_size,j_matrix_size) = &
                  embed_dens_matr_temp(i_matrix_size,j_matrix_size) + &
                  womega(i_freq)*&!embed_gf(i_matrix_size,j_matrix_size, i_freq) !&
!                  ((real(embed_gf(i_matrix_size,j_matrix_size, i_freq)))*cos((omega_grid(i_freq))*(-1E-6)) &
!                  - (aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq)))*sin(omega_grid(i_freq)*(-1E-6)))
                  (real(embed_gf(i_matrix_size,j_matrix_size, i_freq))&
                  - (aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq)))*sin(omega_grid(i_freq)*(-1E-2)))

                 enddo
             enddo
          enddo
        endif



        !   do i_a = 1, n_basis, 1

        !      write(use_unit,*) i_a,(embed_dens_matr_temp(i_a,i_a)) 
        !    enddo



!embed_dens_matr(:,:)=0.d0


write(use_unit,*) 'chemical_potential',chemical_potential

!            do j_matrix_size = 1, n_basis, 1
!              do i_matrix_size = 1, n_basis, 1

!               if(i_matrix_size.eq.j_matrix_size) then
!                  embed_dens_matr_temp(i_matrix_size,j_matrix_size) = embed_dens_matr_temp(i_matrix_size,j_matrix_size)+(pi)/2.d0
!               endif

!              enddo
!            enddo



!   End : Hier wird die GF zusammengestellt

!        call dgemm ('N','N',n_basis, n_basis, n_basis, &
!           1.d0, k_summed_overlap_matr(:,:), n_basis, &
!          embed_dens_matr_temp(:,:), n_basis, 1.d0,&
!         embed_dens_matr(:,:), n_basis  )



!write(use_unit,*)'inner_loop_reiterations',inner_loop_reiterations

embed_dens_matr = (1./(pi))*(embed_dens_matr_temp) !- chemical_potential
!embed_dens_matr = embed_gf_FT(:,:,0)

          do i_a = 1, n_basis, 1
              write(use_unit,*) i_a,(embed_dens_matr(i_a,i_a)) 
          enddo



 embed_part_number =0.d0

 do i_a = 1, n_basis,1

    embed_part_number= embed_part_number + embed_dens_matr(i_a,i_a)

 enddo


     output = .true.
      if(output) then
          open(57, file='dens_matr_HF_test.dat')
            do i_a = 1, n_basis, 1
              
              write(57,*) i_a, &
                       (embed_dens_matr(i_a,i_a))!,  embed_gf_FT(i_a,i_a,0)
                      ! k_summed_overlap_matr
                      ! real(embed_dens_matr(i_a,i_a)), &
                      ! aimag(embed_dens_matr(i_a,i_a))
                      ! real( prod(:,:,1) ),
                      ! real( prod(i_a,i_a,2)) , real( prod(i_a,i_a,1)) 
            enddo
            do i_a = 1, n_basis, 1

              write(57,*) i_a, &
                       (embed_dens_matr(1,i_a)) 
            enddo
             write(57,*) 'particle numbers, embed and on_site = ', embed_part_number

  !        do i_a = 1, n_basis, 1
  !        enddo

             write(57,*) 'womega = ', womega(1:10)!, lda_part_number
          close(57)
      endif


  end subroutine get_density_matr_test 
