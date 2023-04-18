  !subroutine test_loc_gf(overlap_matrix, hamiltonian, inv_gf_loc) ! (KS_eigenvalue)
  subroutine get_spectral_func(gf_loc,full_ovlp_matrix, on_site_gf, k_summed_overlap_matr, number_reiterations) ! (KS_eigenvalue)

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
   integer i_basis
   integer j_basis
   integer k_basis
   integer i_count
   integer i_state
   integer i_freq
   integer i_k_point
   integer ham_size
   character*2 iter
   character*17 filename
   character*17 filename_IM
   integer i_matrix_size
   integer j_matrix_size
 
! subroutine actual parameters
   complex*16, dimension(n_basis,n_basis,n_k_points) :: full_ovlp_matrix
   complex*16, dimension(:,:), allocatable :: full_inv_ovlp_matrix
   complex*16, dimension(:,:,:,:), allocatable :: full_hamiltonian
   complex*16, dimension(:,:,:,:), allocatable :: aux_gf 
   complex*16, dimension(:,:,:,:), allocatable :: new_loc_gf
   complex*16, dimension(:,:,:), allocatable :: aux_ovlp 
   complex*16, dimension(:,:,:), allocatable :: new_ovlp
   real*8, dimension(n_basis,n_basis) :: k_summed_overlap_matr 
   real*8, dimension(n_basis,n_basis) :: loc_self_enrg
   complex*16, dimension(n_basis,n_basis,nomega) :: inv_on_site_gf
   complex*16, dimension(n_basis,n_basis,nomega) :: inv_gf_test
!   complex*16, dimension(n_basis,n_basis,n_freq,n_k_points) :: inv_gf_loc
   complex*16, dimension(n_basis,n_basis,nomega) :: on_site_gf
   complex*16, dimension(:,:,:,:), allocatable :: inv_gf_loc 
!   complex*16, dimension(:,:,:), allocatable :: inv_on_site_gf 
   complex*16, dimension(n_basis,n_basis,n_k_points,nomega) :: gf_loc

!   complex*16, dimension(:,:,:), allocatable :: on_site_gf 
   complex*16, dimension(:,:,:), allocatable :: product_of_the_two 
   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 

!  real*8, dimension(:), allocatable:: overlap_matrix_loc
  real*8, dimension(:), allocatable:: spectral_func 
  real*8, dimension(:,:,:), allocatable:: on_site_dens_matr_temp
  complex*16, dimension(:,:), allocatable:: on_site_dens_matr_ovlp_k_sum
  complex*16, dimension(:,:), allocatable:: on_site_dens_matr_ovlp_sum
  complex*16, dimension(:,:,:), allocatable:: on_site_dens_matr_ovlp
  real*8, dimension(:,:), allocatable:: on_site_dens_matr
  complex*16, dimension(:,:), allocatable:: trace_KS_GF 
  complex*16, dimension(:,:), allocatable:: trace_loc_GF
!  complex*16, dimension(:), allocatable:: trace_KS_GF 
!  complex*16, dimension(:), allocatable:: trace_loc_GF
  real*8, dimension(:), allocatable:: on_site_part_number_k
  complex*16, dimension(:,:,:), allocatable:: inv_on_site_dens_matr_ovlp_k
  complex*16, dimension(:,:), allocatable:: inv_on_site_dens_matr_ovlp_k_sum
  complex*16, dimension(:,:), allocatable:: on_site_dens_matr_k_sum
  complex*16, dimension(:,:,:), allocatable:: KS_gf 

   real :: on_site_part_number
   real :: on_site_part_number_test



   complex*16, dimension(:,:,:,:), allocatable :: spectral_func_k 


   logical  :: output
   integer number_reiterations

  real*8, dimension(:,:,:), allocatable ::    xc_matr

!e_freq= nomega

!      if(myid.eq.0)then

!write(use_unit,*) shape(xc_pot), shape(hamiltonian)
!stop

!i_k_points=1

!------------------------------------------------------   
      if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif
 
  
       if (.not. allocated(spectral_func_k))then
          allocate(spectral_func_k(n_basis,n_basis,n_k_points,nomega))
       endif
       if (.not. allocated(aux_ovlp))then
          allocate(aux_ovlp(n_basis,n_states,n_k_points))
       endif
       if (.not. allocated(new_ovlp))then
          allocate(new_ovlp(n_states,n_states,n_k_points))
       endif
       if (.not. allocated(aux_gf))then
          allocate(aux_gf(n_basis,n_states,n_k_points,nomega))
       endif
       if (.not. allocated(new_loc_gf))then
          allocate(new_loc_gf(n_states,n_states,n_k_points,nomega))
       endif

       if (.not. allocated(spectral_func))then
          allocate(spectral_func(nomega))
       endif

      if (.not. allocated(on_site_dens_matr_temp))then
          allocate(on_site_dens_matr_temp(n_basis,n_basis, n_k_points))
       endif
      if (.not. allocated(on_site_part_number_k))then
          allocate(on_site_part_number_k(n_k_points))
       endif
      if (.not. allocated(on_site_dens_matr))then
          allocate(on_site_dens_matr(n_basis,n_basis))
       endif

      if (.not. allocated(on_site_dens_matr_ovlp))then
          allocate(on_site_dens_matr_ovlp(n_basis,n_basis, n_k_points))
       endif
      if (.not. allocated(on_site_dens_matr_ovlp_sum))then
          allocate(on_site_dens_matr_ovlp_sum(n_basis,n_basis))
       endif
      if (.not. allocated(on_site_dens_matr_ovlp_k_sum))then
          allocate(on_site_dens_matr_ovlp_k_sum(n_basis,n_basis))
       endif

      if (.not. allocated(inv_on_site_dens_matr_ovlp_k))then
          allocate(inv_on_site_dens_matr_ovlp_k(n_basis,n_basis, n_k_points))
       endif
      if (.not. allocated(inv_on_site_dens_matr_ovlp_k_sum))then
          allocate(inv_on_site_dens_matr_ovlp_k_sum(n_basis,n_basis))
       endif

      if (.not. allocated(on_site_dens_matr_k_sum))then
          allocate(on_site_dens_matr_k_sum(n_basis,n_basis))
       endif
      if (.not. allocated(KS_gf))then
          allocate(KS_gf(n_states,n_k_points,nomega))
       endif
      if (.not. allocated(trace_KS_GF))then
          allocate(trace_KS_GF(n_k_points,nomega))
       endif
      if (.not. allocated(trace_loc_GF))then
          allocate(trace_loc_GF(n_k_points,nomega))
       endif
!      if (.not. allocated(trace_KS_GF))then
!          allocate(trace_KS_GF(nomega))
!       endif
!      if (.not. allocated(trace_loc_GF))then
!          allocate(trace_loc_GF(nomega))
!       endif



!----------------------------------------------------------------------------

!----------------------------------------------------------------------------

 aux_gf(:,:,:,:) = (0.d0,0.d0)
 new_loc_gf(:,:,:,:) = (0.d0,0.d0)
 aux_ovlp(:,:,:) = (0.d0,0.d0)
 new_ovlp(:,:,:) = (0.d0,0.d0)
do i_k_point = 1,n_k_points
!do i_freq = 1,nomega


!        call zgemm ('N','N',n_basis,n_states, n_basis,&
!           (1.d0,0.d0), gf_loc(:,:,i_k_point,i_freq), n_basis, &
!           !(1.d0), aux_hamiltonian (:,:), n_basis, &
!           KS_eigenvector_complex(:,:,1,i_k_point), n_basis, (0.d0,0d0),&
!            aux_gf(:,:,i_k_point,i_freq),n_basis)
!
!        call zgemm ('T','N',n_states,n_states, n_basis,&
!            (1.d0, 0.d0),KS_eigenvector_complex(:,:,1,i_k_point) , n_basis, &
!            aux_gf(:,:,i_k_point,i_freq), n_basis, (0.d0,0.d0), &
!            new_loc_gf (:,:,i_k_point,i_freq) , n_states)



        !call zgemm ('N','N',n_basis,n_basis, n_basis,&
        !   (1.d0,0.d0), gf_loc(:,:,i_k_point,i_freq), n_basis, &
        !   !(1.d0), aux_hamiltonian (:,:), n_basis, &
        !   full_ovlp_matrix(:,:,i_k_point), n_basis, (0.d0,0d0),&
        !    aux_gf(:,:,i_k_point,i_freq),n_basis)

        !call zgemm ('T','N',n_states,n_basis, n_basis,&
        !    (1.d0, 0.d0),full_ovlp_matrix(:,:,i_k_point) , n_basis, &
        !    aux_gf(:,:,i_k_point,i_freq), n_basis, (0.d0,0.d0), &
        !    new_loc_gf (:,:,i_k_point,i_freq) , n_basis)



        call zgemm ('N','N',n_basis,n_states, n_basis,&
           (1.d0,0.d0), KS_eigenvector_complex(:,:,1,i_k_point), n_basis, &
           full_ovlp_matrix(:,:,i_k_point), n_basis, (0.d0,0d0),&
            aux_ovlp(:,:,i_k_point),n_basis)

        call zgemm ('T','N',n_states,n_states, n_basis,&
            (1.d0, 0.d0),KS_eigenvector_complex(:,:,1,i_k_point) , n_basis, &
            aux_ovlp(:,:,i_k_point), n_basis, (0.d0,0.d0), &
            new_ovlp (:,:,i_k_point) , n_states)





!enddo
enddo

!----------------------------------------------------------------------------

!----------------------------------------------------------------------------


 spectral_func(:) = 0.d0

   do i_k_point= 1, n_k_points
   do i_freq= 1, nomega
    do i_basis= 1, n_basis
     spectral_func(i_freq) = spectral_func(i_freq)+aimag(spectral_func_k(i_basis,i_basis,i_k_point,i_freq))
    enddo 
   enddo
   enddo
 
     spectral_func(:) = -(1./pi)*spectral_func(:)
     spectral_func(:) = (1./n_k_points)*(spectral_func(:))
 

!----------------------------------------------------------------------------

!----------------------------------------------------------------------------


         on_site_dens_matr_temp(:,:,:) = 0.d0

   do i_k_point= 1, n_k_points
        do i_freq = 1, nomega
            do j_matrix_size = 1, n_basis, 1
              do i_matrix_size = 1, n_basis, 1
                  on_site_dens_matr_temp(i_matrix_size,j_matrix_size,i_k_point) = &
                  on_site_dens_matr_temp(i_matrix_size,j_matrix_size,i_k_point) + &
                  womega(i_freq)*real(gf_loc(i_matrix_size,j_matrix_size,i_k_point, i_freq)) 
!                  (real(gf_loc(i_matrix_size,j_matrix_size,i_k_point, i_freq))*cos((omega(i_freq))*(0))&! for grid with 80 points
!                  -aimag(gf_loc(i_matrix_size,j_matrix_size,i_k_point, i_freq))*sin(omega(i_freq)*(2.545E-4))) ! for grid with 80 points
!                  (real(gf_loc(i_matrix_size,j_matrix_size,i_k_point, i_freq))*cos((omega(i_freq))*(4.1023E-5))&! for grid with 200 points
!                 -aimag(gf_loc(i_matrix_size,j_matrix_size,i_k_point,i_freq))*sin(omega(i_freq)*(4.1023E-5))) ! for grid with 200 points
!                  (real(gf_loc(i_matrix_size,j_matrix_size,i_k_point, i_freq))*cos((omega(i_freq))*(1.00552E-3))&! for grid with 40 points
!                 -aimag(gf_loc(i_matrix_size,j_matrix_size,i_k_point, i_freq))*sin(omega(i_freq)*(1.00552E-3))) ! for grid with 40 points
                 enddo
             enddo
          enddo
    enddo


         KS_gf(:,:,:) = (0.d0,0.d0)

        do i_freq = 1, nomega
            do i_k_point = 1, n_k_points, 1
              do i_matrix_size = 1, n_states, 1
                  KS_gf(i_matrix_size,i_k_point,i_freq) = &
                  KS_gf(i_matrix_size,i_k_point, i_freq) + &
                  1.d0/(((0.d0,1.d0)*omega(i_freq)+chemical_potential) - KS_eigenvalue(i_matrix_size,1,i_k_point))
                 enddo
             enddo
          enddo





   do i_k_point= 1, n_k_points
!        do i_freq = 1, nomega
              do i_matrix_size = 1, n_states, 1

               !write(use_unit,*) 'The local GF and KS-GF', i_k_point, i_matrix_size, new_loc_gf(i_matrix_size,i_matrix_size,i_k_point, 1), KS_gf(i_matrix_size,i_k_point, 1)

               enddo
!              do i_matrix_size = 1, n_states, 1

!               write(use_unit,*) 'The KS GF', i_k_point, i_matrix_size, KS_gf(i_matrix_size,i_k_point, 1)

!               enddo


        enddo
!    enddo
!stop

trace_KS_GF(:,:) = (0.d0,0.d0)
trace_loc_GF(:,:) = (0.d0,0.d0)

   do i_k_point= 1, n_k_points
        do i_freq = 1, nomega
              do i_matrix_size = 1, n_states, 1
!if (occ_numbers(i_matrix_size,1, i_k_point).ne.0)then

!                trace_KS_GF(i_k_point,i_freq) = trace_KS_GF(i_k_point,i_freq) + KS_gf(i_matrix_size,i_k_point, i_freq)             
                trace_KS_GF(i_k_point,i_freq) = trace_KS_GF(i_k_point,i_freq) + KS_gf(i_matrix_size,i_k_point, i_freq)             
!endif
               enddo
               do j_matrix_size = 1, n_basis, 1
               !do j_matrix_size = 1, n_basis, 1

                trace_loc_GF(i_k_point,i_freq) = &
                     trace_loc_GF(i_k_point,i_freq) &
                   + gf_loc(j_matrix_size,j_matrix_size,i_k_point, i_freq)

               !enddo
               enddo

        enddo
   enddo



   do i_k_point= 1, n_k_points

               write(use_unit,*) 'The trace of the local-GF and the KS-GF', i_k_point, trace_loc_GF(i_k_point, 1), trace_KS_GF(i_k_point,1)
!               write(use_unit,*) 'The trace of the local-GF and the KS-GF', i_k_point, trace_loc_GF( 1), trace_KS_GF(1)
   !            write(use_unit,*) 'The new_ovlp', i_k_point, new_ovlp(:,:,i_k_point) 


!               write(use_unit,*) 'The trace of the KS-GF', i_k_point, trace_KS_GF(i_k_point,1) 
   enddo

!stop



on_site_dens_matr_temp = (1./(pi))*(on_site_dens_matr_temp) + full_ovlp_matrix/2.d0
if(myid.eq.0) then
on_site_part_number =0.d0
do i_k_point = 1, n_k_points, 1
!   write(use_unit,*) "i_k_point:", i_k_point
  on_site_part_number_test =0.d0
  do i_matrix_size = 1, n_basis, 1
   on_site_part_number_test = on_site_part_number_test + &
             on_site_dens_matr_temp(i_matrix_size,i_matrix_size,i_k_point)
!   write(use_unit,'(I4,2f16.8)')  i_matrix_size,  on_site_dens_matr_temp(i_matrix_size,i_matrix_size,i_k_point), on_site_part_number_test/(n_k_points)
    
  enddo 
  on_site_part_number = on_site_part_number + on_site_part_number_test
enddo
  write(use_unit,'(A,f16.8)') "Final electron number:", on_site_part_number/n_k_points
endif

!stop


       on_site_dens_matr_ovlp(:,:,:)=0.d0


do i_k_point = 1,n_k_points

        call zgemm ('N','N',n_basis, n_basis, n_basis, &
           1.d0, on_site_dens_matr_temp(:,:,i_k_point), n_basis, &
          full_ovlp_matrix(:,:,i_k_point), n_basis, 0.d0,&
         on_site_dens_matr_ovlp(:,:,i_k_point), n_basis  )

enddo





!----------------------------------------------------------------------------


!----------------------------------------------------------------------------

on_site_dens_matr_ovlp_sum(:,:) =0.d0

do i_k_point = 1,n_k_points

on_site_dens_matr_ovlp_sum(:,:)=on_site_dens_matr_ovlp_sum(:,:)+ on_site_dens_matr_temp(:,:,i_k_point)

enddo


on_site_dens_matr_ovlp_sum(:,:) = (1./n_k_points)*on_site_dens_matr_ovlp_sum(:,:)


 on_site_part_number =0.d0

 do i_a = 1, n_basis,1

    on_site_part_number= on_site_part_number + on_site_dens_matr_ovlp_sum(i_a,i_a) 

 enddo

 on_site_part_number_k(:) =0.d0
do i_k_point = 1,n_k_points

 do i_a = 1, n_basis,1

    on_site_part_number_k(i_k_point)= on_site_part_number_k(i_k_point) + on_site_dens_matr_temp(i_a,i_a,i_k_point)

 enddo
 enddo





!  open(65, file='spectral_function.dat')
!    write(65,*) 'particle_number_test', on_site_part_number
!  do i_freq= 1, nomega
!    write(65,*) omega(i_freq), spectral_func(i_freq)
!  enddo
!  close(65)



if(myid.eq.0)then
     output = .true.
      if(output) then
           ! do a = 1, n_freq, 1
          if( number_reiterations.lt.10 ) then
             write(iter,'(A,I1)') "0", number_reiterations
          else
             write(iter,'(I2)') number_reiterations
          endif

          filename = "dos_iter_"//iter//".dat"
        !    do a = 1, n_freq, 1
          open(65, file=filename)
              write(65,*) 'particle_number_test', on_site_part_number
!           do a = 1, nomega, 1
           do a = 1, n_k_points, 1
              write(65,*) a,  on_site_part_number_k(a)
            enddo
!              write(65,*) omega(a),  spectral_func(a)
!            enddo
          close(65)
      endif
endif




  !  endif

  end subroutine get_spectral_func 
