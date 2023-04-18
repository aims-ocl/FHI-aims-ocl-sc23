!****s* FHI-aims

      subroutine  single_particle_energy_v2 ( &
         green_fn_time,           &
         hartree_pot,             &
         hartree_pot0,            &
         xc_matr ,                &
         exchange_self_energy,    &
         exchange_self_energy0,   &
         exchange_energy,         &
         exchange_energy0,        &
         hartree_en,              &
         hartree_en0,             &          
         dft_xc_en,               &
         KS_hamilt_en             &
         )

! PURPOSE
     !calculate the total energy contribution due to 
     !the static operator in the Hamiltonian:
     ! the static terms can be computed as:
     ! Tr(G(it=0) h ) 
     !where h refers to a generic static operator      

      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use physics

      implicit none
!ARGUMENTS
      real*8  :: green_fn_time          (n_basis,n_basis,n_spin)
      real*8  :: hartree_pot            (n_basis,n_basis)
      real*8  :: hartree_pot0           (n_basis,n_basis)
      real*8  :: xc_matr                (n_basis,n_basis,n_spin)
      real*8  :: exchange_self_energy   (n_basis,n_basis,n_spin)
      real*8  :: exchange_self_energy0  (n_basis,n_basis,n_spin)
      real*8  :: hartree_en
      real*8  :: hartree_en0 
      real*8  :: exchange_energy
      real*8  :: exchange_energy0
      real*8  :: KS_hamilt_en 
      real*8  :: dft_xc_en
!COUNTERS 
      integer j_basis, i_basis
      integer i_index
      integer i_spin

!start work

      i_index           = 0
      hartree_en        = 0.d0
      hartree_en0       = 0.d0
      exchange_energy   = 0.d0
      exchange_energy0  = 0.d0
      KS_hamilt_en      = 0.d0
      dft_xc_en         = 0.d0

      do i_spin = 1, n_spin
        i_index = 0 
        do i_basis =1 , n_basis, 1
          do j_basis =1, i_basis, 1
            i_index = i_index+1
            if (i_basis .eq. j_basis)then

                KS_hamilt_en = KS_hamilt_en +&
                    green_fn_time (i_basis,j_basis,i_spin)* &
                    hamiltonian   (i_index,i_spin)
       
                dft_xc_en = dft_xc_en +&
                    green_fn_time (i_basis,j_basis,i_spin)*  &
                    xc_matr       (i_basis,j_basis,i_spin)
       
                exchange_energy0 = exchange_energy0 +&
                    green_fn_time         (i_basis,j_basis,i_spin)*  &
                    exchange_self_energy0 (i_basis,j_basis,i_spin)
       
                exchange_energy = exchange_energy +&
                    green_fn_time        (i_basis,j_basis,i_spin)*  &
                    exchange_self_energy (i_basis,j_basis,i_spin)
       
                hartree_en0 = hartree_en0 +&
                    green_fn_time (i_basis,j_basis,i_spin)*  &
                    hartree_pot0  (i_basis,j_basis)
                
                hartree_en = hartree_en +&
                    green_fn_time (i_basis,j_basis,i_spin)*  &
                    hartree_pot   (i_basis,j_basis)
            else
                KS_hamilt_en = KS_hamilt_en +&
                    2* green_fn_time (i_basis,j_basis,i_spin)* &
                    hamiltonian      (i_index,i_spin)
       
                dft_xc_en = dft_xc_en +&
                    2* green_fn_time (i_basis,j_basis,i_spin)*  &
                    xc_matr          (i_basis,j_basis,i_spin)
       
                exchange_energy0 = exchange_energy0 +&
                    2* green_fn_time      (i_basis,j_basis,i_spin)*  &
                    exchange_self_energy0 (i_basis,j_basis,i_spin)
       
                exchange_energy = exchange_energy +&
                    2* green_fn_time     (i_basis,j_basis,i_spin)*  &
                    exchange_self_energy (i_basis,j_basis,i_spin)
       
                hartree_en0 = hartree_en0 +&
                    2* green_fn_time (i_basis,j_basis,i_spin)*  &
                    hartree_pot0     (i_basis,j_basis)
       
                hartree_en = hartree_en +&
                    2* green_fn_time (i_basis,j_basis,i_spin)*  &
                    hartree_pot      (i_basis,j_basis)
            endif
          enddo
        enddo
      enddo

      !add the spin degeneracy
      KS_hamilt_en    = KS_hamilt_en         * (2.d0/n_spin)
      dft_xc_en        = dft_xc_en           * (2.d0/n_spin)
      exchange_energy  = exchange_energy     * (2.d0/n_spin) / 2.d0
      exchange_energy0 = exchange_energy0    * (2.d0/n_spin) / 2.d0
      hartree_en       = hartree_en          * (2.d0/n_spin) / 2.d0
      hartree_en0      = hartree_en0         * (2.d0/n_spin) / 2.d0

      end subroutine single_particle_energy_v2

!--------------------------------


!-----------------------------------------------------
!****s* FHI-aims

      subroutine  single_particle_energy ( &
         n_homo2,&
         green_fn_time, &
         !inv_overlap_matrix,&
         hartree_pot, &
         KS_hamilt_en, &
         xc_matr , exchange_self_energy,&
         exchange_self_energy0 ,&
         delta_hartree, exchange_energy,&
          exchange_energy0,&
         dft_xc_en , full_hartree_pot,&
         hartree_en, hartree_en_GW &!, hartree_dft, hartree_en_ensemble &
         )

! PURPOSE
      !calculate the total energy contribution due to 
      !the static operator in the Hamiltonian:
     ! the static terms can be computed as:

     ! Tr(G(it=0) h ) 

      !where h refers to a generic static operator      

      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use physics

      implicit none
!ARGUMENTS
      integer n_homo2
      real*8  :: green_fn_time        (n_basis,n_basis,n_spin)
!      real*8  :: inv_overlap_matrix   (n_basis,n_basis)
      real*8  :: hartree_pot          (n_basis,n_basis,n_spin)
!      real*8  :: hartree_dft          (n_basis,n_basis,n_spin)
      real*8  :: full_hartree_pot     (n_basis,n_basis,n_spin)
      real*8  :: xc_matr              (n_basis, n_basis,n_spin)
      real*8  :: exchange_self_energy   (n_basis, n_basis,n_spin)
      real*8  :: exchange_self_energy0 (n_basis, n_basis,n_spin)
      real*8  :: KS_hamilt_en 
      real*8  :: exchange_energy
      real*8  :: exchange_energy0
      real*8  :: delta_hartree
      real*8  :: hartree_en
      real*8  :: hartree_en_GW 
!      real*8  :: hartree_en_ensemble 
      real*8  :: dft_xc_en
!EXTRA
!      real*8 sum_eigenvalues
!      real*8 energy_from_G
!      real*8 e_fermi
!      real*8 n_part 
!      character*2 iter
!      character*17 filename
!      logical output 
!COUNTERS 
!      integer i_tau
!      integer i_freq
!      integer i_state
      integer j_basis, i_basis, i_spin
      integer i_index
      integer n, i
!AUX(n_basis, n_basis)
!      real*8 aux_KS_eigenvector (n_basis, n_states)
!      real*8 diagonal_green_fn  (n_states, n_states)
!      real*8 green_tmp (n_basis, n_basis)
!      real*8 aux_matr (n_basis, n_states)
!      real*8 aux_green_fn (n_basis, n_basis)
!      real*8 sum_GH
!      real*8 aux_en_matr (n_basis, n_basis) 
!      real*8 aux_en_matr_KS (n_states, n_states) 
!      real*8 H_times_G (n_states, n_states)
       
!      real*8  :: aux_hartree_pot          (n_basis,n_basis)
!      real*8  :: aux_hartree_en
!start work

      i_index           = 0
      KS_hamilt_en   = 0.d0
      dft_xc_en         = 0.d0
      delta_hartree     = 0.d0
      hartree_en        = 0.d0
     ! aux_hartree_en    = 0.d0
      exchange_energy   = 0.d0
      exchange_energy0 = 0.d0
      hartree_en_GW     = 0.d0
!n      hartree_en_ensemble     = 0.d0

!      if(n_spin.gt.1)then
!        aux_hartree_pot (:,:) = (full_hartree_pot (:,:,1) + full_hartree_pot (:,:,2))/2.d0 
!      endif
  
      do i_spin = 1, n_spin
        i_index = 0 
        do i_basis =1 , n_basis, 1
          do j_basis =1, i_basis, 1
            i_index = i_index+1
       
            if (i_basis .eq. j_basis)then
                KS_hamilt_en = KS_hamilt_en +&
                    green_fn_time (i_basis,j_basis,i_spin)* &
                    hamiltonian (i_index,i_spin)
       
                dft_xc_en = dft_xc_en +&
                    green_fn_time (i_basis,j_basis,i_spin)*  &
                     xc_matr(i_basis,j_basis,i_spin)
       
                exchange_energy0 = exchange_energy0 +&
                    green_fn_time (i_basis,j_basis,i_spin)*  &
                     exchange_self_energy0 (i_basis,j_basis,i_spin)
       
                exchange_energy = exchange_energy +&
                    green_fn_time (i_basis,j_basis,i_spin)*  &
                     exchange_self_energy (i_basis,j_basis,i_spin)
       
                delta_hartree = delta_hartree +&
                    green_fn_time (i_basis,j_basis,i_spin)*  &
                    hartree_pot (i_basis,j_basis,i_spin)
                
!                hartree_en_ensemble = hartree_en_ensemble + & 
!                    green_fn_time (i_basis,j_basis,i_spin)*  &
!                    hartree_dft (i_basis,j_basis,i_spin)
      
                hartree_en_GW = hartree_en_GW +&
                    green_fn_time (i_basis,j_basis,i_spin)*  &
                    ( hartree_pot (i_basis,j_basis,i_spin)+&
                 full_hartree_pot (i_basis,j_basis,i_spin))
       
!                aux_hartree_en = aux_hartree_en +&
!                    green_fn_time (i_basis,j_basis,i_spin)*  &
!                    aux_hartree_pot (i_basis,j_basis)
       
                hartree_en = hartree_en +&
                    green_fn_time (i_basis,j_basis,i_spin)*  &
                    full_hartree_pot (i_basis,j_basis,i_spin)
            else
                KS_hamilt_en = KS_hamilt_en +&
                   2* green_fn_time (i_basis,j_basis,i_spin)* &
                    hamiltonian (i_index,i_spin)
       
                dft_xc_en = dft_xc_en +&
                   2* green_fn_time (i_basis,j_basis,i_spin)*  &
                     xc_matr(i_basis,j_basis,i_spin)
       
                exchange_energy0 = exchange_energy0 +&
                   2* green_fn_time (i_basis,j_basis,i_spin)*  &
                     exchange_self_energy0 (i_basis,j_basis,i_spin)
       
                exchange_energy = exchange_energy +&
                   2* green_fn_time (i_basis,j_basis,i_spin)*  &
                     exchange_self_energy (i_basis,j_basis,i_spin)
       
                delta_hartree = delta_hartree +&
                   2* green_fn_time (i_basis,j_basis,i_spin)*  &
                     hartree_pot (i_basis,j_basis,i_spin)
       
                hartree_en_GW = hartree_en_GW +&
                   2* green_fn_time (i_basis,j_basis,i_spin)*  &
                    ( hartree_pot (i_basis,j_basis,i_spin)+&
                 full_hartree_pot (i_basis,j_basis,i_spin))
       
!                hartree_en_ensemble = hartree_en_ensemble +&
!                   2* green_fn_time (i_basis,j_basis,i_spin)*  &
!                    ( hartree_pot (i_basis,j_basis,i_spin)+&
!                 full_hartree_pot (i_basis,j_basis,i_spin))
       
!                aux_hartree_en = aux_hartree_en +&
!                   2* green_fn_time (i_basis,j_basis,i_spin)*  &
!                    aux_hartree_pot (i_basis,j_basis)
       
                hartree_en = hartree_en +&
                   2* green_fn_time (i_basis,j_basis,i_spin)*  &
                    full_hartree_pot (i_basis,j_basis,i_spin)
       
                
            endif
          enddo
        enddo
      enddo

      !add the spin degeneracy
      KS_hamilt_en = KS_hamilt_en     * (2.d0/n_spin)
      dft_xc_en       = dft_xc_en           * (2.d0/n_spin)
      exchange_energy = exchange_energy     * (2.d0/n_spin)
      exchange_energy0 = exchange_energy0 * (2.d0/n_spin)
      delta_hartree   = delta_hartree       * (2.d0/n_spin)
      hartree_en      = hartree_en          * (2.d0/n_spin)
      hartree_en_GW   = hartree_en_GW       * (2.d0/n_spin)
!      hartree_en_ensemble   = hartree_en_ensemble       * (2.d0/n_spin)

      end subroutine single_particle_energy

!--------------------------------
!****s* FHI-aims

      subroutine  single_particle_energy_v1 ( &
         green_fn_time,  &
         matrix, &
         exp_value )
!         hartree_pot, &
!         KS_hamilt_en, &
!         xc_matr , exchange_self_energy,&
!         exchange_self_energy0 ,&
!         delta_hartree, exchange_energy,&
!          exchange_energy0,&
!         dft_xc_en , full_hartree_pot,&
!         hartree_en, hartree_en_GW &
!         )

! PURPOSE
      !calculate the total energy contribution due to 
      !the static operator in the Hamiltonian:
     ! the static terms can be computed as:

     ! Tr(G(it=0) h ) 

      !where h refers to a generic static operator      

      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use physics

      implicit none
!ARGUMENTS
      real*8  :: green_fn_time        (n_basis,n_basis)
      real*8  :: matrix               (n_basis,n_basis)
      real*8  :: exp_value 

!COUNTERS 
      integer j_basis, i_basis
      integer i_spin

!start work

!      do i_spin = 1, n_spin
        do i_basis =1 , n_basis, 1
          do j_basis =1, i_basis, 1
            if (i_basis .eq. j_basis)then
                exp_value = exp_value +&
                    green_fn_time (i_basis,j_basis)*  &
                     matrix(i_basis,j_basis)
            else
                exp_value = exp_value +&
                   2* green_fn_time (i_basis,j_basis)*  &
                     matrix(i_basis,j_basis)
            endif
          enddo
        enddo
!      enddo

      !add the spin degeneracy
!      KS_hamilt_en = KS_hamilt_en     * (2.d0/n_spin)
!      dft_xc_en       = dft_xc_en           * (2.d0/n_spin)
!      exchange_energy = exchange_energy     * (2.d0/n_spin)
!      exchange_energy0 = exchange_energy0 * (2.d0/n_spin)
!      delta_hartree   = delta_hartree       * (2.d0/n_spin)
!      hartree_en      = hartree_en          * (2.d0/n_spin)
!      hartree_en_GW   = hartree_en_GW       * (2.d0/n_spin)

      end subroutine single_particle_energy_v1



