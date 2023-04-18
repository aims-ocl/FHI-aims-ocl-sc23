!****s* FHI-aims/evaluate_gw_selfenergy_band_kpoint
!  NAME
!   evaluate_gw_selfenergy_band_kpoint
!  SYNOPSIS

      subroutine evaluate_gw_selfenergy_band_kpoint &
           (n_low_state, n_high_state, &
            i_q_point, n_freq, n_full_freq, &
            omega, omega_full, womega_full, &
            chemical_potential_spin, mode, &
            screened_coulomb, gw_selfe_band, &
            delta_gw_selfe_band &
           )

!  PURPOSE
!  Subroutine "evaluate_gw_selfenergy_band_kpoint" evaluates the 
!  the GW self-energy for a periodic system.
!

! USES
      use dimensions
      use prodbas
      use hartree_fock
      use physics, only : overlap_matrix, hamiltonian, n_electrons
      use lvl_triples
      use tight_binding_auxmat
      use constants
      use mpi_tasks
      use synchronize_mpi_basic
      use lapack_wrapper
      use timing
      use g_times_w_p0
      use crpa_blacs
      use exchange_trico
      use exchange_ev
      use lvl_tricoeff
      use pbc_lists

      implicit none

! ARGUMENTS 

      integer :: n_low_state
      integer :: n_high_state
      integer :: n_freq
      integer :: n_full_freq
      integer :: i_q_point
      real*8  :: omega(n_freq)
      real*8  :: omega_full(n_full_freq)
      real*8  :: womega_full(n_full_freq)
      real*8  :: chemical_potential_spin(n_spin)
      integer :: mode
      complex*16  :: screened_coulomb(lbb_row:ubb_row, lbb_col:ubb_col, n_freq)

!     output
      complex*16  :: gw_selfe_band(n_freq,n_low_state:n_high_state,n_spin,n_band_kpoints)
      complex*16  :: delta_gw_selfe_band(n_freq,n_low_state:n_high_state,n_spin,n_band_kpoints)
      integer tasks

! INPUTS
! o  n_freq -- integer number, the number of frequency points for the self energy
! o  n_full_freq -- integer number,
!            the number of frequency points for the screened Coulomb interaction W
! o  n_low_state  -- integer number,
!            the lowest KS/HF eigenstate taken into account in the GW self-energy calculations
! o  n_high_state -- integer number,
!            the highest KS/HF eigenstate taken into account
! o  n_electrons -- real number
!            the total number of electrons in the system
! o  omega(n_freq) -- real array
!            the frequency grid for the self energy
! o  omega_full(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  womega_full(n_freq) -- real array
!            the weigth of the Gauss-Legendre frequency grid for the screened Coulomb 
!            in teraction
! o  chemical_potential_spin -- real number, the chemical potential of the system
! o  screened_coulomb -- the screened coulomb matrix at a set of irreducible (regular) k grid points

!
! OUTPUT
! o  gw_selfe_band -- complex array, the calculated GW self energy at a special set of k points for
!             band structure plotting
! o  delta_gw_selfe_band -- a correction term of the calculated GW self energy at a special set of k points for
!             band structure plotting
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
 
!  local variables

      integer :: output_priority_old
      integer :: info
      integer :: n_k_points_band
      integer :: n_k_points_band_task

      real*8 :: k_band(3)
      real*8, dimension(3):: k_minus_q_vector

      integer win_ev

      integer dcsize, count, mpierr, k_task
      integer(kind=MPI_ADDRESS_KIND):: nbytes, offset

      integer, dimension(n_k_points):: n_states_k_saved
      complex*16, dimension(:,:), allocatable :: kphase_saved, k_phase_band
      complex*16, dimension(:,:,:,:), allocatable :: lvl_tricoeff_recip_r_k
      complex*16, dimension(:,:,:,:), allocatable :: lvl_r_kq_current

      real*8, dimension(:), allocatable :: k_weights_band
      integer, dimension(n_spin):: n_homo_kq_aux

      real*8, dimension(:,:,:), allocatable :: KS_eigenvalue_kq_aux
      complex*16, dimension(:,:,:,:), allocatable :: KS_eigenvector_complex_kq_aux
      real*8, dimension(:,:,:), allocatable :: occ_numbers_kq_aux

      real*8, dimension(:,:,:), allocatable :: KS_eigenvalue_k
      complex*16, dimension(:,:,:,:), allocatable :: KS_eigenvector_complex_k
      real*8, dimension(:,:,:), allocatable :: occ_numbers_k

      complex*16, dimension(:,:,:), allocatable :: KS_eigenvector_k_current, KS_eigenvector_kq_current

      character(*), parameter :: &
              func='evaluate_gw_selfenergy_band_kpoint'
!     timing
      
      real*8  time_lvl_coeff, clock_time_lvl_coeff
      real*8  tot_time_lvl_coeff, tot_clock_time_lvl_coeff
      real*8  time_gw_selfe, clock_time_gw_selfe
      real*8  tot_time_gw_selfe, tot_clock_time_gw_selfe

!     parameters of the fitting tails
!       real*8  s1, s2, omega_1, omega_2
!       real*8  alpha, beta, a, b

      character*100 :: info_str
      integer, dimension(4):: mem
!     counters

      integer i_band
      integer i_k_point_band
      integer i_k_point
      integer i_index
      integer i_index_local
      integer i_k_point_local
      integer i_cell_n
      integer i_state
      integer i_spin
      integer i_k_ev
      integer n_qs_points_task, i_qt_point

      integer:: mirror
      integer:: status(MPI_STATUS_SIZE)
      integer, dimension(2):: k_point_band_loc

!     begin work
      call perfon('gwse')

      output_priority_old = output_priority
      output_priority = OL_high

      allocate(kphase_saved(n_cells,n_k_points),stat=info) 
      call check_allocation(info, 'kphase_saved',func)

      allocate(KS_eigenvector_k_current(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_k_current',func)

      allocate(KS_eigenvalue_kq_aux(n_states,n_spin,n_k_points),stat=info) 
      call check_allocation(info, 'KS_eigenvalue_kq_aux',func)

      n_qs_points_task = 0
      do i_k_point = 1, n_k_points
         if(myid_bl.eq.mod(i_k_point,n_tasks_bl) .and.(myid_bl .le. n_k_points)) then
            n_qs_points_task = n_qs_points_task + 1
         endif
      enddo

      allocate(KS_eigenvector_complex_kq_aux(n_basis,n_states,n_spin,n_qs_points_task),stat=info) 
      call check_allocation(info, 'KS_eigenvector_complex_kq_aux',func)

      allocate(occ_numbers_kq_aux(n_states,n_spin,n_k_points),stat=info) 
      call check_allocation(info, 'occ_numbers_kq_aux',func)

      kphase_saved(:,:) = k_phase(:,:)
      n_states_k_saved = n_states_k

!      tot_time_lvl_coeff = 0.d0
!      tot_clock_time_lvl_coeff = 0.d0
!      tot_time_gw_selfe = 0.d0
!      tot_clock_time_gw_selfe = 0.d0

!      if(myid.eq.0) then
!         write(use_unit,*) 
!         write(use_unit,'(2X,2A)') "Evaluating the GW self-energy along", &
!              " the high-symmetry directions in the Brillouin zone "
!      endif

      allocate(lvl_tricoeff_recip_r_k(lbb_row:ubb_row,max_n_basis_sp,n_states,n_spin),stat=info)
      call check_allocation(info, 'lvl_tricoeff_recip_r_k', func)

      allocate(lvl_r_kq_current(lbb_row:ubb_row,max_n_basis_sp,n_states,n_spin),stat=info)
      call check_allocation(info, 'lvl_r_kq_current', func)

      if(n_tasks_bl.gt.1) then
         !need arrays to collect data from other ks
         allocate(KS_eigenvector_kq_current(n_basis,n_states,n_spin),stat=info) 
         call check_allocation(info, 'KS_eigenvector_kq_current',func)
      end if

      !      call initialize_lvl_triples(OVLP_TYPE_COULOMB)
      i_index = 0

      do i_band = 1, n_plot_band
         
         n_k_points_band =  n_points_in_band(i_band)

         n_k_points_band_task = 0
         do i_k_point = 1, n_k_points_band, 1
            if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points_band )then
               n_k_points_band_task = n_k_points_band_task + 1
            end if
         end do
         
         !for the moment, some global variables have to be redefined because they are used in some subroutines
         n_k_points=n_k_points_band
         n_k_points_task=n_k_points_band_task
         
!         if(myid.eq.0) then
!            write(use_unit,'(2X,A,I4)') "i_band : ", i_band
!         endif
         
         ! First compute the LVL coefficients and KS eigenvectors at these special k points for band plotting
         deallocate(k_weights) 
         allocate(k_weights(n_k_points_band))
         k_weights(:)=1.d0/dble(n_k_points_band)

         
         deallocate(k_phase) 
         allocate(k_phase(n_cells,n_k_points_band))
         do i_k_point = 1, n_k_points_band
            k_band(:) = band_begin(i_band,:) +  &
                 real(i_k_point-1)/real(n_k_points_band-1) &
                 *( band_end(i_band,:) -  band_begin(i_band,:))
            k_band(:) = modulo(k_band(:),1.d0)
            
            !           if(i_k_point .eq. 1) then
            !            write(use_unit,'(A,I4,3f16.4)') "i_band_k_point first", i_band, k_band(:)
            !           endif
            do i_cell_n = 1, n_cells
               k_phase( i_cell_n, i_k_point) = exp((0,2)*pi*sum(k_band(:)*cell_index(i_cell_n,:)))
            enddo
         enddo

         allocate(k_phase_band(n_cells,n_k_points_band))
         k_phase_band=k_phase

         allocate(KS_eigenvalue_k(n_states,n_spin,n_k_points_band),stat=info) 
         call check_allocation(info, 'KS_eigenvalue_k',func)
         
         allocate(KS_eigenvector_complex_k(n_basis,n_states,n_spin,n_k_points_band_task),stat=info) 
         call check_allocation(info, 'KS_eigenvector_complex_k',func)
         
         allocate(occ_numbers_k(n_states,n_spin,n_k_points_band),stat=info) 
         call check_allocation(info, 'occ_numbers_k',func)

         deallocate(flag_KS_k_points) 
         allocate(flag_KS_k_points(n_k_points_band))
         deallocate(n_states_k)
         allocate(n_states_k(n_k_points_band))
         call check_allocation(info, 'flag_KS_k_points',func)
         flag_KS_k_points(:)=BASIS_SINGULAR_NOT_TESTED

         call get_KS_orbitals_bandplot(overlap_matrix, hamiltonian, n_electrons, &
              KS_eigenvalue_k, KS_eigenvector_complex_k, occ_numbers_k, .false.) 

         ! Now compute the LVL coefficients and KS eigenvectors at the k-dependent k-q mesh

         !create windows to access distributed arrays         
         call init_access_ev_complex(n_k_points_band,n_k_points_band_task,KS_eigenvector_complex_k,win_ev)


         ! change everything back to normal setting
         n_k_points=n_q_points
         n_k_points_task=n_q_points_task

         deallocate(k_phase) 
         allocate(k_phase(n_cells,n_k_points))
         deallocate(flag_KS_k_points) 
         allocate(flag_KS_k_points(n_k_points))
         call check_allocation(info, 'flag_KS_k_points',func)
         flag_KS_k_points(:)=BASIS_SINGULAR_NOT_TESTED
         deallocate(n_states_k)
         allocate(n_states_k(n_k_points))

         deallocate(k_weights) 
         allocate(k_weights(n_k_points))
         k_weights(:)=1.d0/dble(n_k_points)

         allocate(k_weights_band(n_k_points_band))
         k_weights_band(:)=1.d0/dble(n_k_points)

call perfon('seloop')
         do i_k_point_band = 1, n_k_points_band

!            call get_timestamps(time_lvl_coeff, clock_time_lvl_coeff)
            
            if (mode.eq.-1) then
               call sync_ev(win_ev)
!               call sync_trico(win_tri_k)
               cycle
            end if

            i_index = i_index+1
            i_k_point_local = (i_k_point_band-1)/n_tasks + 1
            i_index_local = (i_index-1)/n_tasks + 1
            k_point_band_loc(1)=mod(i_k_point_band,n_tasks)
            k_point_band_loc(2)=(i_k_point_band-1)/n_tasks + 1
!            call access_ev(win_ev,k_point_loc(:,i_k_point_band),KS_eigenvector_k_current)
            call access_ev(win_ev,k_point_band_loc,KS_eigenvector_k_current)

!            call get_timestamps(rtime, clock_rtime)
!            time_lvl_coeff = rtime - time_lvl_coeff
!            clock_time_lvl_coeff = clock_rtime - clock_time_lvl_coeff
!            tot_time_lvl_coeff = tot_time_lvl_coeff + time_lvl_coeff
!            tot_clock_time_lvl_coeff = tot_clock_time_lvl_coeff + clock_time_lvl_coeff
            
            call get_timestamps(time_gw_selfe, clock_time_gw_selfe)
              
            call gw_get_single_lvl_tricoeff(n_cells, k_phase_band(:,i_k_point_band), &
                 KS_eigenvector_k_current, lvl_tricoeff_recip_r_k)
            
            n_k_points_task=n_qs_points_task
            
            k_band(:) = band_begin(i_band,:) +  &
                 real(i_k_point_band-1)/real(n_k_points_band-1) &
                 *( band_end(i_band,:) -  band_begin(i_band,:))
            k_band(:) = modulo(k_band(:),1.d0)
            do i_qt_point = 1, n_k_points
               k_minus_q_vector(:) = &
                    modulo(k_band(:)-k_point_list(i_qt_point,:), 1.d0)
               do i_cell_n = 1, n_cells
                  k_phase(i_cell_n, i_qt_point) = exp((0,2)*pi* &
                       sum(k_minus_q_vector(:)*cell_index(i_cell_n,:)))
               enddo
             enddo
call perfon('gkso')
             call get_KS_orbitals_bandplot(overlap_matrix, hamiltonian, n_electrons, &
                  KS_eigenvalue_kq_aux, KS_eigenvector_complex_kq_aux, occ_numbers_kq_aux, .true.) 
call perfoff
             do i_spin = 1, n_spin, 1
                do i_state = 1, n_states, 1
                   if(occ_numbers_kq_aux(i_state,i_spin,i_q_point).gt.1.e-6) then
                      n_homo_kq_aux(i_spin) = i_state
                   endif
                enddo
             enddo
            n_k_points_task=n_q_points_task

            if (n_tasks_bl.gt.1) then

               k_task=mod(i_q_point,n_tasks_bl)
               i_k_ev=(i_q_point-1)/n_tasks_bl + 1
               if (myid_bl.eq.k_task) KS_eigenvector_kq_current=KS_eigenvector_complex_kq_aux(:,:,:,i_k_ev)

               call mpi_bcast(KS_eigenvector_kq_current,n_basis*n_states*n_spin, MPI_DOUBLE_COMPLEX, &
                    k_task,comm_blacs,mpierr)
               
               call gw_get_single_lvl_tricoeff(n_cells, k_phase(:,i_q_point), &
                    KS_eigenvector_kq_current, lvl_r_kq_current)

               call kq_pair(KS_eigenvector_k_current, KS_eigenvector_kq_current, &
                    KS_eigenvalue_kq_aux(:,:,i_q_point), k_weights_band(i_k_point_band),&
                    lvl_tricoeff_recip_r_k, lvl_r_kq_current, &
                    mode, n_homo_kq_aux, chemical_potential_spin, screened_coulomb, &
                    gw_selfe_band(:,:,:,i_index), delta_gw_selfe_band(:,:,:,i_index))               

            else
               !no n_basbas parallelization
 
               call gw_get_single_lvl_tricoeff(n_cells, k_phase(:,i_q_point), &
                    KS_eigenvector_complex_kq_aux(:,:,:,i_q_point), lvl_r_kq_current)               
               
               call kq_pair(KS_eigenvector_k_current, KS_eigenvector_complex_kq_aux(:,:,:,i_q_point), &
                    KS_eigenvalue_kq_aux(:,:,i_q_point), k_weights_band(i_k_point_band),&
                    lvl_tricoeff_recip_r_k, lvl_r_kq_current, &
                    mode, n_homo_kq_aux, chemical_potential_spin, screened_coulomb, &
                    gw_selfe_band(:,:,:,i_index), delta_gw_selfe_band(:,:,:,i_index))

            end if
         !  end of loop over i_k_point_band
         end do
call perfoff
         deallocate(k_phase_band)

        deallocate(KS_eigenvalue_k) 
        deallocate(KS_eigenvector_complex_k) 
 
        deallocate(occ_numbers_k) 

        deallocate(k_weights_band) 

       call finalize_access_ev(win_ev)
!  end of loop over i_band
      enddo

      deallocate(lvl_tricoeff_recip_r_k)
      deallocate(lvl_r_kq_current) 
      if (n_tasks_bl.gt.1) then
         deallocate(KS_eigenvector_kq_current)
      end if

      deallocate(KS_eigenvalue_kq_aux) 
      deallocate(KS_eigenvector_complex_kq_aux) 
      deallocate(occ_numbers_kq_aux)

!  change the content of "k_phase" back to the regular k grid
 
      k_phase(:,:) = kphase_saved(:,:)
      n_states_k=n_states_k_saved


      if (allocated(kphase_saved)) then
        deallocate(kphase_saved)
      endif

      if (allocated(KS_eigenvector_k_current)) then
        deallocate(KS_eigenvector_k_current)
      endif

      output_priority = output_priority_old
      call perfoff
      return

    end subroutine evaluate_gw_selfenergy_band_kpoint
