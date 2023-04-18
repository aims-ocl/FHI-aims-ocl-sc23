      subroutine get_GW_self_energy (embed_gf, embed_gf_time, &
                                     self_energy_freq,&
                                     k_summed_overlap_matr)


      use dimensions
      use runtime_choices
      use species_data
      use physics
      use prodbas
      use hartree_fock
      use gw_para
      use dmft_para
      use constants
      use mpi_tasks
      use synchronize_mpi
      use scgw_grid
      use poles_fit
      use gt
      use localorb_io, only: use_unit


!  begin variables

      implicit none
!  Output

!for self- consistency
      integer i_matrix_size, i_basis, j_basis, i_freq, i_basbas, j_basbas
      integer number_reiterations, i_spin, i_index
!      real*8 chemical_pot_c 
      logical output
      logical print_spectrum
      logical compute_spectral_func 
      complex*16, dimension(n_basis,n_basis,nomega) :: embed_gf



      real*8, dimension(n_basis,n_basis,-ntau:ntau,n_spin) :: embed_gf_time
      complex*16, dimension(n_basis,n_basis) :: k_summed_overlap_matr
      real*8, dimension(:,:,:), allocatable :: polarizability_freq
      real*8, dimension(:,:,:), allocatable :: polarizability_time
      real*8, dimension(:,:,:), allocatable :: screened_coul_int_freq
      real*8, dimension(:,:,:), allocatable :: screened_coul_int_time
      real*8, dimension(:,:,:), allocatable :: self_energy_time
      complex*16, dimension(n_basis,n_basis,nomega) :: self_energy_freq
      complex*16, dimension(:,:,:), allocatable :: sumedup_self_energy_gw
!      complex*16, dimension(n_basis,n_basis,n_k_points) :: full_ovlp_matrix_sqrt




      if (.not.allocated(sumedup_self_energy_gw)) then
         allocate(sumedup_self_energy_gw(n_basis,n_basis,nomega))
      endif





             call transform_G (embed_gf , n_basis, &
               n_basis, embed_gf_time(:,:,:,1))

!------------POLARIZABILITY

        if (.not.allocated(polarizability_time)) then
           allocate(polarizability_time(n_basbas,n_loc_prodbas,0:ntau))
        endif

        call get_polar (embed_gf_time, tau, &
           ntau, wtau, polarizability_time)


!---------------FOURIER TRANSFORM THE POLARIZABILITY

      if (.not.allocated(polarizability_freq)) then
         allocate(polarizability_freq(n_basbas,n_loc_prodbas,nomega))
      endif

      call transform_polar(polarizability_time,  &
       n_basbas, n_loc_prodbas, &
       polarizability_freq )




      if (allocated(polarizability_time)) then
          deallocate(polarizability_time)
      end if


!----------EVALUATE THE SCREENED COULOMB INTERACTION

      if (.not.allocated(screened_coul_int_freq)) then
         allocate(screened_coul_int_freq(n_basbas,n_loc_prodbas,nomega))
      endif
      screened_coul_int_freq(:,:,:) = 0.d0

      if (myid.eq.0) then
        write(use_unit,*) "  --- Evaluating the Screened Coulomb interaction"
      endif

      do i_freq = 1, nomega,1

        call screened_coulomb_interaction     &
             (polarizability_freq(:,:,i_freq),  &
              screened_coul_int_freq(:,:,i_freq))
!    subtract the bare Coulomb potential from W, Wc = W-V 

        do j_basbas = 1, n_loc_prodbas, 1
          i_index = map_prodbas(j_basbas,myid+1)
         do i_basbas = 1, n_basbas, 1
            if(i_index.eq.i_basbas) then
               screened_coul_int_freq(i_basbas, j_basbas,i_freq) = &
                screened_coul_int_freq(i_basbas, j_basbas,i_freq) - 1.d0
            endif
          enddo
       enddo
      enddo
      if (allocated(polarizability_freq)) then
         deallocate(polarizability_freq)
      end if
!stop
!!----------------TRANSFORM THE SCI TO TIME  


!stop
        if (.not.allocated(screened_coul_int_time)) then
           allocate(screened_coul_int_time(n_basbas,n_loc_prodbas,0:ntau))
        endif

        call transform_W (screened_coul_int_freq,&
           n_basbas, n_loc_prodbas,&
           screened_coul_int_time )


        if (allocated(screened_coul_int_freq)) then
           deallocate(screened_coul_int_freq)
        endif


!------------- EVALUATE THE SELF ENERGY          

        if(myid.eq.0)then
          write(use_unit,*)"  --- Evaluating the Self Energy "
        endif
        if (.not.allocated(self_energy_time)) then
           allocate(self_energy_time(n_basis,n_basis,-ntau:ntau))
        endif

!        do i_spin = 1, n_spin
          call get_self_energy (  embed_gf_time(:,:,:,1), &
             screened_coul_int_time, &
             self_energy_time  )
!        enddo

          if (allocated(screened_coul_int_time)) then
             deallocate(screened_coul_int_time)
          endif

!----------------FOURIER TRANSFORM THE SELF ENERGY 
      if(myid.eq.0) then
        write(use_unit,*)&
     "         | Fourier transform of the Self Energy from time to frequency "
      endif

          self_energy_freq(:,:,:) =(0.d0,0.d0)
           call transform_sigma (&
                 self_energy_time, n_basis, n_basis,  &
                 self_energy_freq)


!           call multiply_ovlp_self_enrg(self_energy_freq,&
!                                     full_ovlp_matrix_sqrt,&
!                                     k_summed_overlap_matr,&
!                                     self_energy_gw_k)




!     endif

      end subroutine get_GW_self_energy
