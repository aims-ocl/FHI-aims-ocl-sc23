!****s* FHI-aims/extrapolar
!  NAME
!  SYNOPSIS

      subroutine  E_kin ( green_fn_time , e_kinetic )

! PURPOSE

      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use physics

! INPUT
! - green_fn_time      is green function in the time domain and NAO basis
! - e_kinetic          the sc-GW kinetic energy 
                    
      implicit none                     
!ARGUMENTS
      real*8  :: green_fn_time        (n_basis,n_basis, n_spin)
      real*8  :: e_kinetic 
!EXTRA
      real*8 sum_eigenvalues
      real*8 energy_from_G
      real*8 e_fermi
      real*8 n_part 
      character*2 iter
      character*17 filename
      logical output 

      real*8, dimension(n_hamiltonian_matrix_size, n_spin) :: kinetic
      real*8  :: kinetic_energy_density (n_basis,n_basis )

!COUNTERS 
      integer i_tau
      integer i_spin
      integer i_freq
      integer i_state
      integer j_basis
      integer i_basis
      integer k_basis
      integer i_index

!      if(myid.eq.0)then
!        write(use_unit,*)"  --- Calculation of single particle energy ---"
!      endif  

      e_kinetic       = 0.d0

      call integrate_real_kinetic_matrix &
               (rho,  &
                partition_tab, l_shell_max, &
                kinetic &
               )

      kinetic_energy_density (:,:) = 0.d0

      do i_spin = 1, n_spin, 1 
        i_index = 0
        do i_basis =1 , n_basis, 1
          do j_basis =1, i_basis, 1
            i_index = i_index+1
            if (i_basis .eq. j_basis)then

                e_kinetic = e_kinetic + &
                    green_fn_time (i_basis,j_basis,i_spin)* &
                    kinetic (i_index,i_spin)

!                kinetic_energy_density (i_basis,j_basis) = &
!                    kinetic_energy_density (i_basis,j_basis) +&
!                    green_fn_time (i_basis,j_basis,i_spin)* &
!                    kinetic (i_index,i_spin)

            else

                e_kinetic = e_kinetic + &
                    2* green_fn_time (i_basis,j_basis,i_spin)* &
                    kinetic (i_index,i_spin)

!               kinetic_energy_density (i_basis,j_basis) = &
!                   kinetic_energy_density (i_basis,j_basis) +&
!                   green_fn_time (i_basis,j_basis,i_spin)* &
!                   kinetic (i_index,i_spin)

!               kinetic_energy_density  (j_basis,i_basis) = &
!                   kinetic_energy_density (i_basis,j_basis)

            endif
          enddo
        enddo
      enddo

      e_kinetic = e_kinetic * (2.0/n_spin)
!      kinetic_energy_density (:,:) = 2.d0 * kinetic_energy_density (:,:)

      end subroutine E_kin 
