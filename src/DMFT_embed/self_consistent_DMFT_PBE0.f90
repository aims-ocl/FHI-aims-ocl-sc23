      subroutine self_consistent_DMFT_PBE0 &
       ( )


      use dimensions
      use timing
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
      use localorb_io, only: localorb_info, OL_norm, use_unit


!  begin variables

      implicit none
!  Output

!for self- consistency
      logical scdmft_converged
      logical inner_loop_converged
      integer i_matrix_size, i_basis, j_basis, i_freq, i_basbas, j_basbas
      integer max_reiteration
      integer max_inner_loop_reiteration
      integer number_reiterations, i_spin, i_index
      integer :: inner_loop_reiterations
!      real*8 chemical_pot_c 
      real*8  threshold_green
      real*8  max_error
      logical output
      logical print_spectrum
      logical compute_spectral_func 

      character*200 :: info_str
      character(*), parameter :: func = 'self_consistent_dmft'

      character*8  :: cdate
      character*10 :: ctime


      real*8, dimension(:,:), allocatable :: free_cluster_ham
      real*8, dimension(:,:), allocatable :: on_site_xc_matr
      real*8, dimension(:,:), allocatable :: on_site_PBE_xc_matr
      real*8, dimension(:,:), allocatable :: on_site_PBE_x_matr
      real*8, dimension(:,:), allocatable :: on_site_PBE_c_matr
      complex*16, dimension(:,:,:), allocatable :: inv_on_site_gf
      complex*16, dimension(:,:,:), allocatable :: on_site_gf
      complex*16, dimension(:,:,:), allocatable :: new_on_site_gf
      real*8, dimension(:,:), allocatable :: embed_dens_matr 
      complex*16, dimension(:,:,:), allocatable :: hybrid_func
      real*8, dimension(:,:), allocatable :: loc_self_enrg 
      real*8, dimension(:,:), allocatable :: loc_self_enrg_ovlp 
      real*8, save, dimension(:,:), allocatable :: hartree_pot_LDA 
      real*8, dimension(:,:), allocatable :: loc_self_enrg_new 
      complex*16, dimension(:,:,:), allocatable :: inv_embed_gf
      complex*16, dimension(:,:,:), allocatable :: embed_gf
      complex*16, dimension(:,:,:), allocatable :: new_embed_gf
      complex*16, dimension(:,:), allocatable :: k_summed_overlap_matr
      complex*16, dimension(:,:,:), allocatable :: full_ovlp_matrix
      complex*16, dimension(:,:,:), allocatable :: inv_full_ovlp_matrix
      complex*16, dimension(:,:,:), allocatable :: full_ovlp_matrix_sqrt
      complex*16,dimension(:,:,:),allocatable::inv_full_ovlp_matrix_sqrt
      complex*16, dimension(:,:,:,:),allocatable :: full_hamiltonian
      complex*16, dimension(:,:,:,:), allocatable :: full_k_xc_matr
      complex*16, dimension(:,:,:,:), allocatable :: full_k_PBE_xc_matr 
      real*8, dimension(:,:), allocatable :: exchange_self_energy 
      real*8, dimension(:,:), allocatable :: hartree_pot_HF
      real*8,dimension(:,:),allocatable:: dens_matr
      real*8,dimension(:,:),allocatable:: dens_matr_LDA

      character*30 name_of_quantity_embed
      character*30 name_of_quantity_on_site
      real*8 :: embed_part_number_LDA

      real*8  error_GF_embed
      real*8  error_GF_on_site
      real*8  new_chemical_potential 
      real*8  new_chemical_potential_cluster 


      real*8, dimension(:,:,:,:), allocatable :: embed_gf_time
      complex*16, dimension(:,:,:), allocatable :: self_energy_freq
      complex*16, dimension(:,:,:), allocatable :: self_energy_freq_old

          nomega = n_full_freq
          ntau = n_full_time

!       if (.not. allocated(inv_on_site_gf))then
!          allocate(inv_on_site_gf(n_basis,n_basis,nomega))
!       endif

!       if (.not. allocated(on_site_gf))then
!          allocate(on_site_gf(n_basis,n_basis,nomega))
!       endif

!       if (.not. allocated(hybrid_func))then
!          allocate(hybrid_func(n_basis,n_basis, nomega))
!       endif

!       if (.not. allocated(embed_gf))then
!          allocate(embed_gf(n_basis,n_basis,nomega))
!       endif

!       if (.not. allocated(new_embed_gf))then
!          allocate(new_embed_gf(n_basis,n_basis,nomega))
!       endif
       
!       if (.not. allocated(inv_embed_gf))then
!          allocate(inv_embed_gf(n_basis,n_basis,nomega))
!       endif

!       if (.not. allocated(loc_self_enrg_new))then
!          allocate(loc_self_enrg_new(n_basis,n_basis))
!       endif

!       if (.not. allocated(hartree_pot_LDA))then
!          allocate(hartree_pot_LDA(n_basis,n_basis))
!       endif

!       if (.not. allocated(loc_self_enrg))then
!          allocate(loc_self_enrg(n_basis,n_basis))
!       endif

!       if (.not. allocated(embed_dens_matr))then
!          allocate(embed_dens_matr(n_basis,n_basis))
!       endif

!       if (.not. allocated(exchange_self_energy))then
!          allocate(exchange_self_energy(n_basis,n_basis))
!       endif
!       if (.not. allocated(hartree_pot_HF))then
!          allocate(hartree_pot_HF(n_basis,n_basis))
!       endif

!       if (.not. allocated(embed_gf_time))then
!          allocate(embed_gf_time(n_basis,n_basis,-ntau:ntau,n_spin))
!       endif


        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "----------------------------&
        &--------------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "Self-Consistent DMFT &
         &embedding starts ..."
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate,&
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "-------------------------------&
        &-----------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )


        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "Initialization of &
         &the frequency grid ..."
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

       if (.not. allocated(womega))then
          allocate(womega(nomega))
       endif

         call init_grid()
               womega(:) = womega1(:)



       if (.not. allocated(full_hamiltonian))then
          allocate(full_hamiltonian(n_basis,n_basis,n_spin,n_k_points))
       endif

       if (.not. allocated(full_ovlp_matrix))then
          allocate(full_ovlp_matrix(n_basis,n_basis,n_k_points))
       endif

!       if (.not. allocated(inv_full_ovlp_matrix_sqrt))then
!         allocate(inv_full_ovlp_matrix_sqrt(n_basis,n_basis,n_k_points))
!       endif
!!       if (.not. allocated(inv_full_ovlp_matrix))then
!         allocate(inv_full_ovlp_matrix(n_basis,n_basis,n_k_points))
!       endif

       if (.not. allocated(free_cluster_ham))then
          allocate(free_cluster_ham(n_basis,n_basis))
       endif

       if (.not. allocated(k_summed_overlap_matr))then
          allocate(k_summed_overlap_matr(n_basis,n_basis))
       endif

       if (.not. allocated(on_site_xc_matr))then
          allocate(on_site_xc_matr(n_basis,n_basis))
       endif

       if (.not. allocated(on_site_PBE_xc_matr))then
          allocate(on_site_PBE_xc_matr(n_basis,n_basis))
       endif

       if (.not. allocated(on_site_PBE_x_matr))then
          allocate(on_site_PBE_x_matr(n_basis,n_basis))
       endif

       if (.not. allocated(on_site_PBE_c_matr))then
          allocate(on_site_PBE_c_matr(n_basis,n_basis))
       endif

       if (.not. allocated(full_ovlp_matrix_sqrt))then
          allocate(full_ovlp_matrix_sqrt(n_basis,n_basis,n_k_points))
       endif

       if (.not. allocated(full_k_PBE_xc_matr))then
          allocate(full_k_PBE_xc_matr(n_basis,n_basis,n_spin,n_k_points))
       endif

!       if (.not. allocated(full_k_xc_matr))then
!          allocate(full_k_xc_matr(n_basis,n_basis,n_spin,n_k_points))
!       endif
 


        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "evaluating the PBC and Cluster &
         &quantities..."
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )


         call get_pbc_and_cluster_quantities(full_hamiltonian, &
                                             full_ovlp_matrix, &
                                             !inv_full_ovlp_matrix, &
                                             free_cluster_ham, &
                                             on_site_xc_matr, &
                                             full_ovlp_matrix_sqrt,&
                                             !inv_full_ovlp_matrix_sqrt,&
                                             on_site_PBE_xc_matr,&
                                             on_site_PBE_x_matr,&
                                             on_site_PBE_c_matr,&
                                             !full_k_xc_matr,&
                                             full_k_PBE_xc_matr,&
                                             k_summed_overlap_matr)


       number_reiterations = 0
       scdmft_converged =  .false.
       max_reiteration =100
       max_inner_loop_reiteration =200
       threshold_green = 1.d-4
 
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "---------------------------------------&
        &---------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  Self-Consistent main loop &
         &initialization ..."
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, &
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "---------------------------------------&
        &---------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )




      do while(.not. scdmft_converged)


       if (.not. allocated(loc_self_enrg))then
          allocate(loc_self_enrg(n_basis,n_basis))
       endif

       if (.not. allocated(loc_self_enrg_new))then
          allocate(loc_self_enrg_new(n_basis,n_basis))
       endif

       if (.not. allocated(hartree_pot_LDA))then
          allocate(hartree_pot_LDA(n_basis,n_basis))
       endif

      if (.not.allocated(self_energy_freq)) then
         allocate(self_energy_freq(n_basis,n_basis,nomega))
      endif

           if (number_reiterations .eq. 0) loc_self_enrg = on_site_xc_matr
           if (number_reiterations .eq. 0) hartree_pot_LDA = (0.d0,0.d0) 
           if (number_reiterations >= 1) loc_self_enrg=loc_self_enrg_new
           if (number_reiterations .eq. 0) self_energy_freq(:,:,:) =&
           (0.d0,0.d0)



        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "Calculating the on-Site Green &
        & Function iteration", number_reiterations
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate,&
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )


       if (.not. allocated(inv_on_site_gf))then
          allocate(inv_on_site_gf(n_basis,n_basis,nomega))
       endif

       if (.not. allocated(on_site_gf))then
          allocate(on_site_gf(n_basis,n_basis,nomega))
       endif

       if (.not. allocated(dens_matr_LDA))then
          allocate(dens_matr_LDA(n_basis,n_basis))
       endif


      call get_on_site_gf_freq_PBE0(full_hamiltonian,&
                                 full_ovlp_matrix,&
                                 inv_on_site_gf,&
                                 on_site_gf,&
                                 on_site_xc_matr,&
                                 loc_self_enrg,&
                                 number_reiterations,&
                                 max_reiteration,&
                                 hartree_pot_LDA,&
                                 full_ovlp_matrix_sqrt,&
                                 !inv_full_ovlp_matrix_sqrt,&
                                 dens_matr_LDA,&
                                 on_site_PBE_xc_matr,& 
                                 !inv_full_ovlp_matrix,& 
                                 new_chemical_potential)!,&
                                 !full_k_xc_matr) 



        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "Calculating the Hybridization &
         &Self-Energy iteration"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, &
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )



       if (.not. allocated(hybrid_func))then
          allocate(hybrid_func(n_basis,n_basis, nomega))
       endif

       if (.not. allocated(free_cluster_ham))then
          allocate(free_cluster_ham(n_basis,n_basis))
       endif
         
         if(number_reiterations.eq.0)  new_chemical_potential=&
           chemical_potential    

          call get_embedding_self_enrg(free_cluster_ham, &
                                       inv_on_site_gf,&
                                       on_site_xc_matr,&
                                       k_summed_overlap_matr,&
                                       hybrid_func,&
                                       loc_self_enrg, &
                                       self_energy_freq, &
                                       hartree_pot_LDA, &
                                       number_reiterations,&
                                       new_chemical_potential)

       if (allocated(inv_on_site_gf))then
          deallocate(inv_on_site_gf)
       endif




       if(number_reiterations.eq.0) then
         call initialize_hartree_fock()
       endif


       inner_loop_converged =  .false.
       inner_loop_reiterations = 0


          if(number_reiterations.eq.0) then

        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  Calculating the embedded &
         &Green function for Particle number Determination..."
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate,&
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )


           loc_self_enrg = on_site_xc_matr
           new_chemical_potential_cluster=  chemical_potential

       if (.not. allocated(embed_gf))then
          allocate(embed_gf(n_basis,n_basis, nomega))
       endif

          call  get_embed_gf_freq(free_cluster_ham, &
                                  hybrid_func,&!, on_site_xc_matr,&
                                  loc_self_enrg,&
                                  embed_gf,&
                                  hartree_pot_LDA,&
                                  inner_loop_reiterations,&
                                  number_reiterations,&
                                  self_energy_freq,&
                                  k_summed_overlap_matr,&
                                  new_chemical_potential_cluster)

        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  Particle Number Calculation ..."
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate,&
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )

        call get_particle_number(embed_gf, &
                               embed_part_number_LDA)
        endif

        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "--------------------------------------------&
        &----------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  Self-Consistent Self-Consistent &
         &DMFT inner-loop starts ..."
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, &
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "--------------------------------------------&
        &----------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )

          

      do while(.not. inner_loop_converged)


          loc_self_enrg = loc_self_enrg_new
          if(number_reiterations.eq.0)then
           if(inner_loop_reiterations.eq.0)then
            loc_self_enrg = on_site_xc_matr    
            hartree_pot_LDA = 0.d0 
           endif
          endif
         

          if(number_reiterations.eq.0) then
          if(inner_loop_reiterations.eq.0) new_chemical_potential_cluster =&
                                            chemical_potential
          endif
         !                                   new_chemical_potential_cluster =&
         !                                   chemical_potential
!         if(.false.)then
         if(.true.)then
         if(number_reiterations.gt.0)then
           if(inner_loop_reiterations.eq.0)then
             new_chemical_potential_cluster = new_chemical_potential
          if(.false.)then
            call  update_chemical_pot_dmft_cluster(embed_gf,& 
                                      inner_loop_reiterations, &
                                      number_reiterations,&
                                      full_ovlp_matrix_sqrt,&
                                      !inv_full_ovlp_matrix_sqrt,&
                                      full_ovlp_matrix,free_cluster_ham,&
                                      hybrid_func,loc_self_enrg,&
                                      hartree_pot_LDA,&
                                      embed_part_number_LDA,&
                                      new_chemical_potential_cluster,&
                                      k_summed_overlap_matr,&
                                      self_energy_freq)
            endif
           endif
         endif
        endif


        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "---------------------------------------&
        &---------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "Calculating the embedded Green &
         &function for iteration", inner_loop_reiterations
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, &
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

       if (.not. allocated(embed_gf))then
          allocate(embed_gf(n_basis,n_basis, nomega))
       endif


          call  get_embed_gf_freq(free_cluster_ham, &
                                  hybrid_func,&
                                  loc_self_enrg,&
                                  embed_gf,&
                                  hartree_pot_LDA,&
                                  inner_loop_reiterations,&
                                  number_reiterations,&
                                  self_energy_freq,&
                                  k_summed_overlap_matr,&
                                  new_chemical_potential_cluster)


        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "Local Exact Exchange &
         &Self-Energy Caluclation starts ..."
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, &
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )


        if(.true.) then
        !if(.false.) then

       if (.not. allocated(embed_gf_time))then
          allocate(embed_gf_time(n_basis,n_basis,-ntau:ntau,n_spin))
       endif

       if (.not. allocated(hartree_pot_HF))then
          allocate(hartree_pot_HF(n_basis,n_basis))
       endif

       if (.not. allocated(exchange_self_energy))then
          allocate(exchange_self_energy(n_basis,n_basis))
       endif

        call transform_G (embed_gf , n_basis, &
                          n_basis, embed_gf_time(:,:,:,1))

        call evaluate_hartree_and_exchange_embed_PBE0(&
                                          embed_gf_time(:,:,0,1), &
                                          loc_self_enrg_new, &
                                          number_reiterations, &
                                          inner_loop_reiterations, &
                                          hartree_pot_LDA,&
                                          hartree_pot_HF,&
                                          exchange_self_energy,&
                                          on_site_xc_matr, &
                                          full_ovlp_matrix_sqrt,&
                                          !inv_full_ovlp_matrix_sqrt,&
                                          on_site_PBE_xc_matr,&
                                          k_summed_overlap_matr)
        endif
 

!     if(.false.)then
     if(.true.)then

        if (number_reiterations.eq.0)then
         if (inner_loop_reiterations.ge.1)then
!        self_energy_freq(:,:,:) = dmft_alpha*self_energy_freq(:,:,:) + &
!                            (1.d0-dmft_alpha)*self_energy_freq_old(:,:,:)
        loc_self_enrg_new(:,:) = dmft_alpha*loc_self_enrg_new(:,:) + &
                            (1.d0-dmft_alpha)*loc_self_enrg(:,:)
         endif
       else
!        self_energy_freq(:,:,:) = dmft_alpha*self_energy_freq(:,:,:) + &
!                            (1.d0-dmft_alpha)*self_energy_freq_old(:,:,:)
        loc_self_enrg_new(:,:) = dmft_alpha*loc_self_enrg_new(:,:) + &
                            (1.d0-dmft_alpha)*loc_self_enrg(:,:)
       endif

     endif


        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "----------------------------------&
        &--------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  Updating the Chemical & 
         &potential for the inner-loop iteration", inner_loop_reiterations
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, &
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "-----------------------------------&
        &-------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )


 
         if(.false.)&
          call  update_chemical_pot_dmft_cluster(embed_gf,& 
                                      inner_loop_reiterations, &
                                      number_reiterations,&
                                      full_ovlp_matrix_sqrt,&
                                      !inv_full_ovlp_matrix_sqrt,&
                                      full_ovlp_matrix,free_cluster_ham,&
                                      hybrid_func,loc_self_enrg_new,&
                                      hartree_pot_LDA,&
                                      embed_part_number_LDA,&
                                      new_chemical_potential_cluster,&
                                      k_summed_overlap_matr,&
                                      self_energy_freq)

       if (.not. allocated(new_embed_gf))then
          allocate(new_embed_gf(n_basis,n_basis,nomega))
       endif

      
          call  get_embed_gf_freq(free_cluster_ham, &
                                  hybrid_func,&
                                  loc_self_enrg_new,&
                                  new_embed_gf,&
                                  hartree_pot_LDA,&
                                  inner_loop_reiterations,&
                                  number_reiterations,&
                                  self_energy_freq,&
                                  k_summed_overlap_matr,&
                                  new_chemical_potential_cluster)



         inner_loop_reiterations =  inner_loop_reiterations +1


       do i_freq = 1, nomega, 1
       do i_basis = 1, n_basis, 1
       do j_basis = 1, n_basis, 1

        new_embed_gf(i_basis,j_basis,i_freq) = DCMPLX(dmft_alpha)*&
                              new_embed_gf(i_basis,j_basis,i_freq) +&
                              ((1.d0,0.d0)-DCMPLX(dmft_alpha))*&
                              embed_gf(i_basis,j_basis,i_freq)

       enddo
       enddo
       enddo

       

      name_of_quantity_embed = "new_embed_GF"

        call check_the_error_dmft(embed_gf, &
                                new_embed_gf, &
                                name_of_quantity_embed,&
                                error_GF_embed, max_error, .false.)

        if ( error_GF_embed.lt. threshold_green)then
         
         inner_loop_converged = .true.
!       endif

        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "---------------------------------------&
        &---------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  Self-Consistent DMFT Embedding &
         &inner-loop converged!"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, &
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "----------------------------------------&
        &--------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )


        elseif(inner_loop_reiterations.gt.max_inner_loop_reiteration) then
         
         inner_loop_converged = .true.


        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "----------------------------------------&
        &--------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  Self-Consistent DMFT &
         &Embedding inner-loop didn't converge!"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  at iteration # ",&
         inner_loop_reiterations
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, &
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "----------------------------------------&
        &--------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )


        endif 
!stop


        enddo ! while(.not. inner_loop_converged), end of the self consistent inner-loop


!--------------------------------------------------------------------
!--------------------inner LOOP deallocations!!!---------------------
!--------------------------------------------------------------------

       if (allocated(hybrid_func))then
          deallocate(hybrid_func)
       endif

       if (allocated(embed_gf))then
          deallocate(embed_gf)
       endif

       if ( allocated(new_embed_gf))then
          deallocate(new_embed_gf)
       endif

!       if ( allocated(inv_embed_gf))then
!          deallocate(inv_embed_gf)
!       endif

       if (allocated(embed_gf_time))then
          deallocate(embed_gf_time)
       endif


!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------


         number_reiterations =  number_reiterations +1

       if (.not. allocated(new_on_site_gf))then
          allocate(new_on_site_gf(n_basis,n_basis,nomega))
       endif

       if (.not. allocated(dens_matr))then
          allocate(dens_matr(n_basis,n_basis))
       endif


     if(.false.) then
!       if(.true.) then
!       if(number_reiterations .eq. 1) then
         !embed_part_number_LDA = 14.d0
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "-------------------------------------&
        &--------------------------------------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  Updating the Chemical &
         &potential for the main-loop iteration", number_reiterations
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, &
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "--------------------------------------&
        &-------------------------------------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )


         call update_chemical_pot_dmft(embed_gf, &
                               inner_loop_reiterations, &
                               number_reiterations, &
                               full_ovlp_matrix_sqrt,&
                               !inv_full_ovlp_matrix_sqrt,&
                               full_ovlp_matrix,free_cluster_ham,&
                               loc_self_enrg_new,&
                               hartree_pot_LDA,&
                               embed_part_number_LDA,&
                               !14.d0,&
                               new_chemical_potential,&
                               !on_site_xc_matr,&
                               on_site_xc_matr,&
                               full_hamiltonian,&
                               on_site_PBE_xc_matr,&
                               self_energy_freq)
!       endif
       endif

        new_chemical_potential = chemical_potential

        call get_NEW_on_site_gf_freq_PBE0(full_hamiltonian, &
                                     full_ovlp_matrix, &
                                     !inv_full_ovlp_matrix, &
                                     new_on_site_gf,&
                                     on_site_xc_matr, &
                                     loc_self_enrg_new, &
                                     hartree_pot_LDA,&
                                     number_reiterations,&
                                     scdmft_converged, &
                                     full_ovlp_matrix_sqrt,&
                                     !inv_full_ovlp_matrix_sqrt,&
                                     on_site_PBE_xc_matr,&
                                     dens_matr,&
                                     new_chemical_potential,&
                                     new_chemical_potential_cluster,&
                                     embed_part_number_LDA)!,&
                                     !full_k_xc_matr)

        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "--------------------------------&
        &----------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  Galitskii-Migdal Total &
         &Energy Calculation iteration", number_reiterations
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  "&
        , cdate, ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "-----------------------------------&
        &-------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )


           call  get_total_energy_dmft_PBE0(new_on_site_gf,&
                                  loc_self_enrg_new,&
                                  inner_loop_reiterations,&
                                  free_cluster_ham,&
                                  hartree_pot_LDA,&
                                  on_site_xc_matr,&
                                  hartree_pot_HF,&
                                  exchange_self_energy,&
                                  on_site_PBE_x_matr,&
                                  on_site_PBE_c_matr,&
                                  !full_k_xc_matr,&
                                  full_hamiltonian,&
                                  full_ovlp_matrix_sqrt,&
                                  !inv_full_ovlp_matrix_sqrt,&
                                  full_ovlp_matrix,&
                                  full_k_PBE_xc_matr,&
                                  dens_matr_LDA,&
                                  dens_matr,&
                                  !inv_full_ovlp_matrix,&
                                  embed_part_number_LDA)




      new_on_site_gf(:,:,:) = DCMPLX(dmft_alpha)*on_site_gf(:,:,:) +&
                       DCMPLX(1-dmft_alpha)*new_on_site_gf(:,:,:)

          name_of_quantity_on_site = "new on-site Green's function "
      if(.true.) &
          call check_the_error_dmft(on_site_gf,&
                                   new_on_site_gf, &
                                   name_of_quantity_on_site,&
                                   error_GF_on_site, max_error, .false.)

 
         if ( error_GF_on_site .lt. threshold_green)then
       
            scdmft_converged = .true.

        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "-----------------------------------&
        &-------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  Self-Consistent DMFT &
         &Embedding converged!"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  at iteration #"&
        , number_reiterations
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, &
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "-----------------------------------&
        &-------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )


       if(.false.) then
!       if(scdmft_converged) then

         call update_chemical_pot_dmft(embed_gf, &
                               inner_loop_reiterations, &
                               number_reiterations, &
                               full_ovlp_matrix_sqrt,&
                               !inv_full_ovlp_matrix_sqrt,&
                               full_ovlp_matrix,free_cluster_ham,&
                               inv_on_site_gf,loc_self_enrg_new,&
                               hartree_pot_LDA,&
                               embed_part_number_LDA,&
                               !14.d0,&
                               new_chemical_potential,&
                               !on_site_xc_matr,&
                               on_site_xc_matr,&
                               full_hamiltonian,&
                               on_site_PBE_xc_matr,&
                               self_energy_freq)

!       endif
       endif

        call get_NEW_on_site_gf_freq_PBE0(full_hamiltonian, &
                                     full_ovlp_matrix, &
                                     !inv_full_ovlp_matrix, &
                                     new_on_site_gf,&
                                     on_site_xc_matr, &
                                     loc_self_enrg_new, &
                                     hartree_pot_LDA,&
                                     number_reiterations,&
                                     scdmft_converged,&
                                     full_ovlp_matrix_sqrt,&
                                     !inv_full_ovlp_matrix_sqrt,&
                                     on_site_PBE_xc_matr,&
                                     dens_matr,&
                                     new_chemical_potential,&
                                     new_chemical_potential_cluster,&
                                     embed_part_number_LDA,&
                                     full_k_xc_matr)

!        endif
          elseif(number_reiterations.gt.max_reiteration) then
          
            scdmft_converged = .true.

        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "------------------------------&
        &------------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  Self-Consistent &
         &DMFT Embedding didn't converge!"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  at iteration #", number_reiterations
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, &
        ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "------------------------------&
        &------------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )


       if(.false.) then
!       if(scdmft_converged) then
         call update_chemical_pot_dmft(embed_gf, &
                               inner_loop_reiterations, &
                               number_reiterations, &
                               full_ovlp_matrix_sqrt,&
                               !inv_full_ovlp_matrix_sqrt,&
                               full_ovlp_matrix,free_cluster_ham,&
                               inv_on_site_gf,loc_self_enrg_new,&
                               hartree_pot_LDA,&
                               embed_part_number_LDA,&
                               new_chemical_potential,&
                               on_site_xc_matr,&
                               full_hamiltonian,&
                               on_site_PBE_xc_matr,&
                               self_energy_freq)
!       endif
       endif


!       if (.not. allocated(gf_non_loc))then
!          allocate(gf_non_loc(n_basis,n_basis,n_k_points,nomega))
!       endif

        call get_NEW_on_site_gf_freq_PBE0(full_hamiltonian, &
                                     full_ovlp_matrix, &
                                     !inv_full_ovlp_matrix, &
                                     new_on_site_gf,&
                                     on_site_xc_matr, &
                                     loc_self_enrg_new, &
                                     hartree_pot_LDA,&
                                     number_reiterations,&
                                     scdmft_converged,&
                                     full_ovlp_matrix_sqrt,&
                                     !inv_full_ovlp_matrix_sqrt,&
                                     on_site_PBE_xc_matr,&
                                     dens_matr,&
                                     new_chemical_potential,&
                                     new_chemical_potential_cluster,&
                                     embed_part_number_LDA)!!,&
                                     !full_k_xc_matr)


         endif


        enddo ! while(.not. scgw_convergedd), end of the self consistent loop
        
        write(info_str,'(A)') ''
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "---------------------------------&
        &---------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(10X,A,1X,I4)') "  Galitskii-Migdal Total &
         &Energy Calculation iteration", number_reiterations
        call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  "&
        , cdate, ", Time     :  ", ctime
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "---------------------------------&
        &---------------------------"
        call localorb_info ( info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_sc_loop, clock_time_sc_loop )


           call  get_total_energy_dmft_PBE0(new_on_site_gf,&
                                  loc_self_enrg_new,&
                                  inner_loop_reiterations,&
                                  free_cluster_ham,&
                                  hartree_pot_LDA,&
                                  on_site_xc_matr,&
                                  hartree_pot_HF,&
                                  exchange_self_energy,&
                                  on_site_PBE_x_matr,&
                                  on_site_PBE_c_matr,&
                                  !full_k_xc_matr,&
                                  full_hamiltonian,&
                                  full_ovlp_matrix_sqrt,&
                                  !inv_full_ovlp_matrix_sqrt,&
                                  full_ovlp_matrix,&
                                  full_k_PBE_xc_matr,&
                                  dens_matr_LDA,&
                                  dens_matr,&
                                  !inv_full_ovlp_matrix,&
                                  embed_part_number_LDA)


      if (allocated(full_k_PBE_xc_matr))then
          deallocate(full_k_PBE_xc_matr)
       endif

           call cleanup_hartree_fock()

!      endif
 

   





       if (myid.eq.0) then
         write(use_unit,*) " ---scDMFT deallocations: completed  "
      ! Still to do for DMFT!!!!! 
       endif



!     endif

      end subroutine self_consistent_DMFT_PBE0
