  subroutine check_norm_dmft_p2_spectrum(aux_omegamax, &
                                        aux_omega, &
                                        aux_womega, &
                                        aux_nomega, &
                                        add_spectrum, &
                                        elec_N, diff_in_elec_N,&
                                        chemical_pot_i,&
                                        i_counter)


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
   integer i_hamilton_size 
   integer inner_loop_reiterations, number_reiterations
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
   character*2 iter
   character*17 filename_RE
   character*17 filename_OUT
   character*17 filename_IM
   integer aux_nomega
   real*8  aux_omegamax
! subroutine actual parameters

   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 


!   complex*16, dimension(:,:,:), allocatable:: omega_term
   integer, intent(inout) :: i_counter


   real*8, intent(out) :: diff_in_elec_N 
!   real*8, intent(out) :: embed_part_number_new 
   real*8, intent(in) :: elec_N
   real*8 :: embed_part_number
   real*8 ::chemical_pot_i 
  real*8 , dimension (aux_nomega) :: aux_omega
  real*8 , dimension (aux_nomega) :: aux_womega
  real*8 , dimension (-aux_nomega:aux_nomega) :: add_spectrum

      logical  :: output

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------


!          allocate(omega_term(n_basis,n_basis,nomega))
!       endif



!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  i_counter = i_counter + 1

!omega_ovlp(:,:,:,:) =(0.d0,0.d0)
!
!   do i_k_point=1,n_k_points
!       do i_freq=1,nomega
!
!          omega_ovlp(:,:, i_freq,i_k_point) = omega_ovlp(:,:, i_freq, i_k_point) +&
!!                                             ((((0.d0,1.d0)*omega(i_freq)))*&!+chemical_potential)*&
!                                             ((((0.d0,1.d0)*omega(i_freq))+chemical_pot_i)*&
!                                             (full_ovlp_matrix(:,:,i_k_point)))
!
!       enddo
!   enddo
!write(use_unit,*) "nomega fro check_norm", nomega
!write(use_unit,*) "nomega fro check_norm", womega
!stop
!full_omega_term_temp(:,:,:,:) =(0.d0,0.d0)
!full_omega_term(:,:,:,:) =(0.d0,0.d0)

!       do i_k_point=1,n_k_points
!         do i_freq = 1, nomega
!
!           call zgemm ('N','N', n_basis,&
!            n_basis,n_basis,(1.d0,0.d0),&
!            inv_full_ovlp_matrix_sqrt(:,:,i_k_point), n_basis, &
!            omega_ovlp(:,:,i_freq,i_k_point), &
!            n_basis,(0.d0,0.d0), &
!            full_omega_term_temp(:,:,i_freq,i_k_point) , n_basis )
!
!           call zgemm ('N','N', n_basis,&
!            n_basis,n_basis,(1.d0,0.d0),&
!            full_omega_term_temp(:,:,i_freq,i_k_point), n_basis, &
!            inv_full_ovlp_matrix_sqrt(:,:,i_k_point), &
!            n_basis,(0.d0,0.d0), &
!            full_omega_term(:,:,i_freq,i_k_point) , n_basis )
!

!!     if(.true.)then
!         embed_part_number = 0.d0
!         do i_freq= 1, aux_nomega, 1
!          if( aux_omega(i_freq) .lt. chemical_pot_i)then
!          !if( aux_omega(i_freq) .lt. 0 )then
!            embed_part_number = embed_part_number + add_spectrum(i_freq)*aux_womega(i_freq)
!          else
!            embed_part_number = embed_part_number +&!  spectrum(i_freq)*aux_womega(i_freq)+&
!              add_spectrum(-i_freq)*aux_womega(i_freq)
!          endif
!         enddo
!if(chemical_pot_i.gt.0) chemical_pot_i = -chemical_pot_i

         embed_part_number = 0.d0
         do i_freq= 1, aux_nomega, 1
          if( aux_omega(i_freq) .le. chemical_pot_i )then
            embed_part_number = embed_part_number + 0.5*add_spectrum(i_freq)*aux_womega(i_freq)
          else
            embed_part_number = embed_part_number + 0.5*add_spectrum(-i_freq)*aux_womega(i_freq)
          endif
         enddo



!        if (myid.eq.0)then
!          write(use_unit,*) " "
!          write(use_unit,*) "n_particles as calculated from the spectrum = ", &
!            embed_part_number
!          write(use_unit,*) " "
!        endif
!      endif


!
!         enddo
!       enddo


!       if(myid.eq.0) then

 diff_in_elec_N = embed_part_number - elec_N
!              write(use_unit,*) 'diff in particle number',  diff_in_elec_N
!              write(use_unit,*) 'max_zeroin', max_zeroin 
    write(use_unit,*) 'chemical_pot_i, particle number',i_counter,chemical_pot_i,embed_part_number, diff_in_elec_N

!embed_part_number_new = embed_part_number
!stop
!write(use_unit,*) "diff_electrons",i_counter, diff_in_elec_N

!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------


  end subroutine check_norm_dmft_p2_spectrum 
