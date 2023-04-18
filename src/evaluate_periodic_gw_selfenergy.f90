!****s* FHI-aims/evaluate_periodic_gw_selfenergy
!  NAME
!   evaluate_periodic_gw_selfenergy
!  SYNOPSIS

      subroutine evaluate_periodic_gw_selfenergy &
           ( n_low_state_polar, n_low_state, n_high_state, &
            occ_numbers, n_freq, n_full_freq, &
            omega, omega_full, womega_full, &
            chemical_potential_spin, dielec_func_imagfreq, &
            KS_eigenvalue, KS_eigenvector, &
            KS_eigenvector_complex, &
            KS_eigenvector_irk,KS_eigenvector_complex_irk, &
            gw_selfenergy, gw_selfe_band, &
            out_self_energy &
           )

!  PURPOSE
!  Subroutine "evaluate_periodic_gw_selfenergy" evaluates the 
!  the GW self-energy for a periodic system.
!

! USES
      use dimensions
      use prodbas
      use hartree_fock
      use constants
      use mpi_tasks
      use synchronize_mpi
      use timing
      use geometry, only: recip_lattice_vector
      use runtime_choices
      use tight_binding_auxmat, only: singularity_lifting_chi
      use evaluate_polarisability_kspace_mod
      use g_times_w_p0
      use crpa_blacs
      use pbc_lists

      implicit none

! ARGUMENTS 

      integer :: n_low_state_polar
      integer :: n_low_state
      integer :: n_high_state
      integer :: n_freq
      integer :: n_full_freq
      real*8  :: occ_numbers(n_states,n_spin,n_k_points)
      real*8  :: omega(n_freq)
      real*8  :: omega_full(n_full_freq)
      real*8  :: womega_full(n_full_freq)
      real*8  :: chemical_potential_spin(n_spin)
      real*8  :: KS_eigenvalue(n_states,n_spin,n_k_points)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin,n_k_points_task)
      complex*16  :: dielec_func_imagfreq(n_full_freq)
      complex*16  :: KS_eigenvector_complex(n_basis,n_states,n_spin,n_k_points_task)
      real*8  :: KS_eigenvector_irk(n_basis,n_states,n_spin,n_irk_points_task)
      complex*16  :: KS_eigenvector_complex_irk(n_basis,n_states,n_spin,n_irk_points_task)

      logical :: out_self_energy
!     output                      
      complex*16  :: gw_selfenergy(n_freq,n_low_state:n_high_state,n_spin,n_irk_points_task)
      complex*16  :: gw_selfe_band(n_freq,n_low_state:n_high_state,n_spin,n_band_kpoints_task)
      integer  :: mpierr
      real*8:: checkloc, checkglob


! INPUTS
! o  n_freq -- integer number,
!            the number of frequency points for the self energy
! o  n_full_freq -- integer number,
!            the number of frequency points for the screened Coulomb interaction W
! o  n_low_state_polar  -- integer number,
!            the lowest KS/HF eigenstate taken into account in the polarizability calculation -- frozen-core
!            approximation
! o  n_low_state  -- integer number,
!            the lowest KS/HF eigenstate taken into account in the GW self-energy calculations
! o  n_high_state -- integer number,
!            the highest KS/HF eigenstate taken into account
! o  n_electrons -- real number
!            the total number of electrons in the system
! o  occ_numbers -- real 2-dimentianal array of length (n_states, n_spin)
!            the occupation number of the electrons for each eigenstate and each spin
! o  omega(n_freq) -- real array
!            the frequency grid for the self energy
! o  omega_full(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  womega_full(n_freq) -- real array
!            the weigth of the Gauss-Legendre frequency grid for the screened Coulomb 
!            in teraction
! o  chemical_potential_spin -- real number, the chemical potential of the system
! o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
! o  KS_eigenvector -- real array,
!            the eigenvector of the single-particle calculation
! o  KS_eigenvector_complex -- complex array,
!            the complex eigenvector of the single-particle calculation,
!            used when "real_eigenvectors == .false."
! o  KS_eigenvector_irk -- the real eigenvector of the single-particle (KS/HF) self-consistent calculation
!            on an irreducible set of k vectors
! o  KS_eigenvector_complex_irk -- the complex eigenvector of the single-particle (KS/HF) self-consistent 
!            calculation on an irreducible set of k vectors

! o  out_self_energy -- if true, print out the self-energy
!           
!
! OUTPUT
! o  GW_selfenergy -- complex array, the calculated GW self energy on the regualr k grid
! o  GW_selfe_band -- complex array, the calculated GW self energy on specific k points for band plotting
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

      real*8  ::  qsq
      real*8   k_lattvec(3)

      complex*16  det_v_times_polar
      complex*16  trace_v_times_polar

!    local timing
      real*8  temp_time_polar
      real*8  temp_clock_time_polar
      real*8  temp_time_self_energy
      real*8  temp_clock_time_self_energy

!     auxiliary matrices for Level 3 Blas matrix multiplications

      integer :: ipiv(n_basbas)
      integer :: info
      integer :: n_nonsingular

      integer, dimension(:,:), allocatable :: n_homo_k
      complex*16, dimension(:,:,:), allocatable :: polar_kspace_complex
      real*8,dimension(:,:,:), allocatable:: polar_kspace_real
      complex*16, dimension(:,:,:), allocatable :: screened_coulomb_full
      complex*16, dimension(:,:), allocatable :: v_times_polar
      complex*16, dimension(:,:), allocatable :: dielec_gamma
      real*8, dimension(:), allocatable :: coulomb_eigenvalues
      complex*16, dimension(:,:), allocatable :: coulomb_eigenvectors
      complex*16, dimension(:,:), allocatable :: sqrtv_eigenvectors

      complex*16, dimension(:,:,:,:), allocatable :: delta_gw_selfenergy, gw_selfenergy1, delta_gw_selfenergy1, gw_selfenergy2
      complex*16, dimension(:,:,:,:), allocatable :: gw_selfe_band1, gw_selfe_band2, delta_gw_selfe_band1

!      real*8, dimension(:), allocatable :: dielec_func_imagfreq
      real*8 :: omega_tmp

      character(*), parameter :: func='evaluate_periodic_gw_selfenergy'
      character*100  filename
      character*80  tmp_str

!     counters

      integer :: i_state
      integer :: i_freq
      integer :: i_spin
      integer :: i_index
      integer :: i_index_local
      integer :: i_band
      integer :: i_k_point
      integer :: i_k_point_band
      integer :: i_k_point_local
      integer :: i_irk_point
      integer :: i_irk_point_local
      integer :: i_prodbas_1,i_prodbas_2
      integer :: i_basis_fn_1,i_basis_fn_2
      integer :: i_task
      integer :: id_recv
      integer :: id_send
      integer :: lc
      integer(kind=MPI_ADDRESS_KIND) :: nbytes, offset
      integer :: max_irkq_points_task, k_task
      integer :: i_irkq_point_local, i_irkq_point, count
      integer, dimension(n_irk_points):: pcount
      integer, dimension(:,:), allocatable:: mode_k
      integer :: irkp, ind, looprange, mode
      integer :: i_row, i_col

!     begin work

      call perfon('ev_gw')
      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"---------------------------------------------------", &
                  "---------------------------------"
        write(use_unit,'(2X,A)') &
              "Start to calculate the periodic GW self energy  ... "
      endif

!     determine the highest occupied orbital level
!     such complication occurs when some of the orbitals are
!     either not fully occupied or not fully empty
!      n_first(:) = 1
!      do i_spin = 1, n_spin
!       do i_state = 1, n_states
!        if (abs(occ_numbers(i_state,i_spin)-dble(2/n_spin)) &
!                         .lt.1.d-8) then
!         n_first(i_spin)= i_state + 1
!        endif
!       enddo
!       if(n_first(i_spin) .gt. n_states) then
!         n_first(i_spin) = n_states
!       endif
!      enddo

      if(flag_KS_eigenfunc_conjg) then
         KS_eigenvector_complex = conjg(KS_eigenvector_complex)
      endif

!     if(myid.eq.0) then
!      write(use_unit,'(2X, A,A,4I5)') "HOMO and first non-fully-occupied", &
!                   " orbitals:", n_homo(:), n_first(:)
!      write(use_unit,*)
!     endif

      !because of the global mpi_fence operations, all mpi tasks have to run all iteration so of the irkq loop
      if (mod(n_irk_points,n_tasks_irkq).eq.0) then    
         max_irkq_points_task=n_irk_points/n_tasks_irkq
      else
         max_irkq_points_task=n_irk_points/n_tasks_irkq+1
      end if

      allocate(n_homo_k(n_spin,n_k_points),stat=info) 
      call check_allocation(info, 'n_homo_k')

      if (real_eigenvectors) then
         allocate(polar_kspace_real(lbb_row:ubb_row, lbb_col:ubb_col, n_full_freq),stat=i_index)
      else
         allocate(polar_kspace_real(0, lbb_col:ubb_col, n_full_freq),stat=i_index)
      end if
      call check_allocation(i_index, 'polar_kspace_real              ')

! full polar_kspace_real path (for n_k_points=1) not yet implemented, need to convert to complex for some computations 
!      if(n_k_points.gt.1) then
         allocate(polar_kspace_complex(lbb_row:ubb_row, lbb_col:ubb_col, n_full_freq),stat=i_index)
!      else
!         allocate(polar_kspace_complex(0, lbb_col:ubb_col, n_full_freq),stat=i_index)
!      end if
      call check_allocation(i_index, 'polar_kspace_complex       ')

      allocate(v_times_polar(lbb_row:ubb_row, lbb_col:ubb_col),stat=info)
      call check_allocation(info, 'v_times_polar                   ')

      allocate(delta_gw_selfenergy(n_freq,n_low_state:n_high_state,n_spin, n_irk_points_task),stat=info)
      call check_allocation(info, 'delta_gw_selfenergy             ')

      allocate(gw_selfenergy1(n_freq,n_low_state:n_high_state,n_spin, n_irk_points),stat=info)
      call check_allocation(info, 'gw_selfenergy1             ')

      allocate(gw_selfenergy2(n_freq,n_low_state:n_high_state,n_spin, n_irk_points),stat=info)
      call check_allocation(info, 'gw_selfenergy2             ')

      allocate(delta_gw_selfenergy1(n_freq,n_low_state:n_high_state,n_spin, n_irk_points),stat=info)
      call check_allocation(info, 'delta_gw_selfenergy1             ')

      if(use_gw_gamma_corr) then
         allocate(dielec_gamma(lbb_row:ubb_row, lbb_col:ubb_col),stat=info)
         call check_allocation(info, 'dielec_gamma                   ')
      endif
      allocate(coulomb_eigenvalues(n_basbas),stat=info)
      call check_allocation(info, 'coulomb_eigenvalues             ')
      allocate(coulomb_eigenvectors(n_bb_row,n_bb_col),stat=info)
      call check_allocation(info, 'coulomb_eigenvectors             ')
      allocate(sqrtv_eigenvectors(n_bb_row,n_bb_col),stat=info)
      call check_allocation(info, 'sqrtv_eigenvectors             ')

      gw_selfenergy1=0.
      delta_gw_selfenergy1=0.

      if(out_band) then
         allocate(gw_selfe_band1(n_freq,n_low_state:n_high_state,n_spin, n_band_kpoints),stat=info)
         call check_allocation(info, 'gw_selfe_band1',func)
         gw_selfe_band1=0.
         allocate(gw_selfe_band2(n_freq,n_low_state:n_high_state,n_spin, n_band_kpoints),stat=info)
         call check_allocation(info, 'gw_selfe_band2',func)
         allocate(delta_gw_selfe_band1(n_freq,n_low_state:n_high_state,n_spin, n_band_kpoints),stat=info)
         call check_allocation(info, 'delta_gw_selfe_band1',func)
         delta_gw_selfe_band1=0.
      endif

      time_polar = 0.d0
      clock_time_polar = 0.d0
      time_self_energy = 0.d0
      clock_time_self_energy = 0.d0
      time_band_self_energy = 0.d0
      clock_time_band_self_energy = 0.d0

!      open(unit=8, file="dielec_imagfreq.dat", action='WRITE')
!!      read(8,*) tmp_str
!      do i_freq = 1, n_full_freq
!!       read(8,*) omega_tmp, dielec_func_imagfreq(i_freq)
!       if(myid.eq.0) then
!         write(8,'(f18.8,2f18.10)') omega_full(i_freq), dielec_func_imagfreq(i_freq)
!       endif
!      enddo
!      close(8)

      n_homo_k(:,:) = 0
      do i_k_point = 1, n_k_points, 1
       do i_spin = 1, n_spin, 1
         do i_state = 1, n_states
          if(occ_numbers(i_state,i_spin,i_k_point) .gt. 1.e-12) then
            n_homo_k(i_spin,i_k_point) = i_state
          endif
         enddo
     
       enddo
      enddo

      ! first compute the square_root of the Coulomb matrix   
      do i_irkq_point_local = 1, n_irkq_points_task
         i_irk_point=n_tasks_irkq*(i_irkq_point_local-1) + myid_irkq + 1
         i_k_point = inv_irk_point_mapping(i_irk_point)

         if (all(abs(k_point_list(i_k_point,:)) .lt. 1.e-10)) then
            call diagonalize_auxmat_scalapack_complex &
                 (n_basbas, coulomb_matr_blacs(:,:,i_irkq_point_local), prodbas_threshold, n_nonsingular, &
                  coulomb_eigenvalues, coulomb_eigenvectors, sqrtv_eigenvectors)

         else
            call power_auxmat_scalapack_complex(n_basbas,0.5d0,coulomb_matr_blacs(:,:,i_irkq_point_local),'')
         endif
      enddo

      do i_irkq_point_local = 1, n_irkq_points_task
           call power_auxmat_scalapack_complex(n_basbas,0.5d0,coulomb_cut_blacs(:,:,i_irkq_point_local),'')
      enddo
      
      checkloc=0.

      call init_compute_g_times_w_p0( n_low_state, n_high_state, &
           n_full_freq, n_freq, n_full_freq, omega_full, womega_full, omega)

!      omega_full(:) = omega_full(1)
      call perfon('loop')
      if (myid.eq.0) then
         write(use_unit,*) 
         write(use_unit,'(2X,A)') &
           'Starting the main loop of GW computation'
         write(use_unit,*) 
      end if
      do i_irkq_point_local=1, max_irkq_points_task
         call get_timestamps(temp_time_polar, temp_clock_time_polar )

         ! evaluate the polarisability
         call  evaluate_polarisability_kspace_list &
              (n_low_state_polar, i_irkq_point_local, n_full_freq, omega_full, KS_eigenvalue, KS_eigenvector, &
              KS_eigenvector_complex, occ_numbers, polar_kspace_real, polar_kspace_complex)      

         !the rest of this subroutine currently only works for complex polar_kspace, need to convert here
         if(n_k_points.eq.1) then
            polar_kspace_complex=cmplx(polar_kspace_real)
         endif

         if (irkblacs_member.and.(myid_irkq.lt.n_irk_points)) checkloc = checkloc + sum(abs(polar_kspace_complex))
         call perfon('epgwla')
         if(i_irkq_point_local.le.n_irkq_points_task) then
            !!$OMP PARALLEL DO private(v_times_polar,i_prodbas_1)

            i_irk_point=n_tasks_irkq*(i_irkq_point_local-1) + myid_irkq + 1
            i_k_point = inv_irk_point_mapping(i_irk_point)

            do i_freq = 1, n_full_freq
               
               !  Multiply sqrt(v) from both sides of \chi_0 
               !  Note that now "coulomb_matrix_recip" contains the square root of the Coulomb matrix
               if (all(abs(k_point_list(i_k_point,:)) .lt. 1.e-10)) then

                  call pzgemm('C', 'N', n_nonsingular, n_basbas, n_basbas, (1.d0,0.d0), &
                          sqrtv_eigenvectors, 1, 1, bb2desc, &
                          polar_kspace_complex(lbb_row,lbb_col,i_freq), 1, 1, bb2desc, (0.d0, 0.d0), &
                          v_times_polar, 1, 1, bb2desc)
               
                  call pzgemm('N', 'N', n_nonsingular, n_nonsingular, n_basbas, (1.d0,0.d0), &
                          v_times_polar, 1, 1, bb2desc, &
                          sqrtv_eigenvectors, 1, 1, bb2desc, (0.d0,0.d0), &
                          polar_kspace_complex(lbb_row,lbb_col,i_freq), 1, 1, bb2desc) 

! a crucial step to modify the head matrix element of the dielectric function (here v^1/2 chi_0 v^1/2)	
                  if(all(n_k_points_xyz(:) .gt. 1)) then
! We switch this correction off for 1D system
                    do i_col = lbb_col, ubb_col, 1
                      do i_row = lbb_row, ubb_row, 1
                        if(i_row .eq. 1 .and. i_col .eq. 1) then
                           polar_kspace_complex(i_row, i_col, i_freq)=  - dielec_func_imagfreq(i_freq)                                       
                         endif
!                        if(i_col .eq. i_row) then
!                           write(use_unit,'(3I4,3f16.8)') i_k_point, i_row, i_freq, polar_kspace_complex(i_row, i_col, i_freq), dielec_func_imagfreq(i_freq)
!                        endif
                       enddo
                     enddo 
                   endif

                  call pzgemm('N', 'N', n_basbas, n_nonsingular, n_nonsingular, (1.d0,0.d0), &
                          coulomb_eigenvectors, 1, 1, bb2desc, &
                          polar_kspace_complex(lbb_row,lbb_col,i_freq), 1, 1, bb2desc, (0.d0, 0.d0), &
                          v_times_polar, 1, 1, bb2desc)
               
                  call pzgemm('N', 'C', n_basbas, n_basbas, n_nonsingular, (1.d0,0.d0), &
                          v_times_polar, 1, 1, bb2desc, &
                          coulomb_eigenvectors, 1, 1, bb2desc, (0.d0,0.d0), &
                          polar_kspace_complex(lbb_row,lbb_col,i_freq), 1, 1, bb2desc) 

!                  write(use_unit,*) "i_k_point", i_k_point
!                  do i_prodbas_2 = 1, n_basbas, 1
!                     do i_prodbas_1 = 1, n_basbas, 1
!                         write(use_unit,'(3I4,2f16.8)') i_k_point, i_prodbas_1, i_prodbas_2, polar_kspace_complex(i_prodbas_1, i_prodbas_2, i_freq)
!                     enddo
!                  enddo
               else

                  call pzgemm('N', 'N', n_basbas, n_basbas, n_basbas, (1.d0,0.d0), &
                          coulomb_matr_blacs(lbb_row,lbb_col,i_irkq_point_local), 1, 1, bb2desc, &
                          polar_kspace_complex(lbb_row,lbb_col,i_freq), 1, 1, bb2desc, (0.d0, 0.d0), &
                          v_times_polar, 1, 1, bb2desc)
               
                  call pzgemm('N', 'N', n_basbas, n_basbas, n_basbas, (1.d0,0.d0), &
                          v_times_polar, 1, 1, bb2desc, &
                          coulomb_matr_blacs(lbb_row,lbb_col,i_irkq_point_local), 1, 1, bb2desc, (0.d0,0.d0), &
                          polar_kspace_complex(lbb_row,lbb_col,i_freq), 1, 1, bb2desc) 

!               if (all(abs(k_point_list(i_k_point,:)) .lt. 1.e-10)) then
!                  write(use_unit,*) "i_k_point", i_k_point
!                  do i_prodbas_2 = 1, n_basbas, 1
!                     do i_prodbas_1 = 1, n_basbas, 1
!                        write(use_unit,'(3I4,2f16.8)') i_k_point, i_prodbas_1, i_prodbas_2, polar_kspace_complex(i_prodbas_1, i_prodbas_2, i_freq)
!                     enddo
!                  enddo
!               endif

               endif
               
                ! computing the inverse of 1-sqrt(v)*\chi_0*sqrt(v)
                v_times_polar(:,:) = - polar_kspace_complex(:,:,i_freq)
                !            write(use_unit,*) "1-v\chi-0"

               do i_prodbas_1 = lbb_row, ubb_row
                  if((i_prodbas_1.ge.lbb_col).and.(i_prodbas_1.le.ubb_col)) then
                     v_times_polar(i_prodbas_1,i_prodbas_1) = &
                          v_times_polar(i_prodbas_1,i_prodbas_1) + (1.d0,0.d0)
                  end if
               enddo
               call perfon('epgwsolv')
               call power_auxmat_scalapack_complex(n_basbas, -1.d0, v_times_polar, '')
               call perfoff
               polar_kspace_complex(:,:,i_freq) = v_times_polar(:,:)
               do i_prodbas_1 = lbb_row, ubb_row
                   if((i_prodbas_1.ge.lbb_col).and.(i_prodbas_1.le.ubb_col)) then
                      polar_kspace_complex(i_prodbas_1,i_prodbas_1,i_freq) = &
                           v_times_polar(i_prodbas_1,i_prodbas_1) - (1.d0,0.d0)
                   end if
               enddo
!               if(use_gw_gamma_corr) then
!                  dielec_gamma(:,:) = polar_kspace_complex(:,:,i_freq)
!
!                  call pzgemm('N', 'N', n_basbas, n_basbas, n_basbas, (1.d0,0.d0), &
!                         coulomb_coeff_blacs(lbb_row,lbb_col,i_irkq_point_local), 1, 1, bb2desc, &
!                         dielec_gamma(lbb_row,lbb_col), 1, 1, bb2desc, (0.d0,0.d0), &
!                         v_times_polar, 1, 1, bb2desc)
!               
!                  call pzgemm('N', 'N', n_basbas, n_basbas, n_basbas, (1.d0,0.d0), &
!                         v_times_polar, 1, 1, bb2desc, &
!                         coulomb_coeff_blacs(lbb_row,lbb_col,i_irkq_point_local), 1, 1, bb2desc, (0.d0,0.d0), &
!                         dielec_gamma(lbb_row,lbb_col), 1, 1, bb2desc) 
!               endif

                ! multiplying 1/(1-v*\chi_0) with sqrt(v) from left-hand and right-hand sides, 
                ! and still puting back into polar_kspace_complex
               call pzgemm('N', 'N', n_basbas, n_basbas, n_basbas, (1.d0,0.d0), &
                     coulomb_cut_blacs(lbb_row,lbb_col,i_irkq_point_local), 1, 1, bb2desc, &
                     polar_kspace_complex(lbb_row,lbb_col,i_freq), 1, 1, bb2desc, (0.d0,0.d0), &
                     v_times_polar, 1, 1, bb2desc)
               
               call pzgemm('N', 'N', n_basbas, n_basbas, n_basbas, (1.d0,0.d0), &
                     v_times_polar, 1, 1, bb2desc, &
                     coulomb_cut_blacs(lbb_row,lbb_col,i_irkq_point_local), 1, 1, bb2desc, (0.d0,0.d0), &
                     polar_kspace_complex(lbb_row,lbb_col,i_freq), 1, 1, bb2desc) 
               
!                  write(use_unit,*) "screened coulomb", i_k_point
!                  do i_prodbas_2 = 1, n_basbas, 1
!                     do i_prodbas_1 = 1, n_basbas, 1
!                        write(use_unit,'(3I4,2f16.8)') i_k_point, i_prodbas_1, i_prodbas_2, polar_kspace_complex(i_prodbas_1, i_prodbas_2, i_freq)
!                     enddo
!                  enddo
!                if(use_gw_gamma_corr) then
!                   i_irk_point=n_tasks_irkq*(i_irkq_point_local-1) + myid_irkq + 1
!                   i_k_point = inv_irk_point_mapping(i_irk_point)
!                   k_lattvec(1:3) = matmul(recip_lattice_vector,k_minus_q_point_list(i_k_point,1:3))
!                   qsq=k_lattvec(1)**2 + k_lattvec(2)**2 + k_lattvec(3)**2
!
!                   if( (abs(k_minus_q_point_list(i_k_point,1)).le.1.e-10 .and. abs(k_minus_q_point_list(i_k_point,2)).le.1.e-10 .and. qsq .gt. 1.e-10 )) then
!                       do i_prodbas_1 = lbb_col, ubb_col, 1
!                         i_basis_fn_1=basbas_fn(i_prodbas_1)
!                         do i_prodbas_2 = lbb_row, ubb_row, 1
!                           i_basis_fn_2=basbas_fn(i_prodbas_2)
!                           if(basbas_l(i_prodbas_1).eq.0 .and. basbas_l(i_prodbas_2).eq.0 .and. multipole_basbas_fn(i_basis_fn_1) .gt. 1.e-10 &
!                         .and. multipole_basbas_fn(i_basis_fn_2).gt.1.e-10) then
!                              write(use_unit,'(A, 2I5,2f13.6,4f18.8)')"w", i_prodbas_2, i_prodbas_1, abs(k_minus_q_point_list(i_k_point,3)), sqrt(qsq), &
!                                       polar_kspace_complex(i_prodbas_2,i_prodbas_1,i_freq), dielec_gamma(i_prodbas_2,i_prodbas_1)/qsq
!                            endif
!                         enddo
!                       enddo
!                    endif
!                    if(qsq.gt.1.e-10) then
!                      polar_kspace_complex(:,:,i_freq)=polar_kspace_complex(:,:,i_freq)-dielec_gamma(:,:)/qsq
!                    endif
!
!                 endif

               !       end of loop over i_freq
             enddo
            !!$OMP END PARALLEL DO
         end if
         call perfoff
         call get_timestamps(rtime, clock_rtime)
         time_polar = time_polar + rtime - temp_time_polar
         clock_time_polar = clock_time_polar + clock_rtime - temp_clock_time_polar
         
         !  Now what is actually contained in "polar_kspace" is the screened coulomb interaction 
         call get_timestamps(temp_time_self_energy, temp_clock_time_self_energy )         
         call compute_g_times_w_p0 &
              ( n_homo_k, i_irkq_point_local, KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
              KS_eigenvector_irk, KS_eigenvector_complex_irk, &
              occ_numbers, chemical_potential_spin, &
              polar_kspace_complex, gw_selfenergy1, &
              delta_gw_selfenergy1)

         call get_timestamps(rtime, clock_rtime)
         time_self_energy = time_self_energy + rtime - temp_time_self_energy
         clock_time_self_energy = clock_time_self_energy + clock_rtime - temp_clock_time_self_energy

         if(out_band) then
            i_irkq_point=(i_irkq_point_local-1)*n_tasks_irkq+myid_irkq+1

            !get the maximum number of k points mapping to one i_irk_point            
            pcount=0
            do i_k_point=1,n_k_points
               irkp=irk_point_mapping(i_k_point)
               pcount(irkp)=pcount(irkp)+1
            end do
            looprange=maxval(pcount)

            allocate(mode_k(2,looprange))
            mode_k=-1
            ind=0
            if (irkblacs_member) then
               do i_k_point=1,n_k_points
                  if (irk_point_mapping(i_k_point).eq.i_irkq_point) then
                     ind=ind+1
                     mode_k(1,ind)=i_k_point
                     if (irk_point_included(i_k_point)) then
                        mode_k(2,ind)=0
                     else
                        mode_k(2,ind)=1
                     end if
                  end if
               end do
            end if

            do ind=1,looprange
               i_k_point=mode_k(1,ind)
               mode=mode_k(2,ind)
               call  evaluate_gw_selfenergy_band_kpoint &
                    ( n_low_state, n_high_state, i_k_point, &
                    n_freq, n_full_freq, &
                    omega, omega_full, womega_full, &
                    chemical_potential_spin, mode, &
                    polar_kspace_complex, gw_selfe_band1, &
                    delta_gw_selfe_band1 & 
                    )
            end do
            deallocate(mode_k)          
            call get_timestamps(temp_time_self_energy, temp_clock_time_self_energy )         

            time_band_self_energy = time_band_self_energy + temp_time_self_energy - rtime 
            clock_time_band_self_energy = clock_time_band_self_energy + temp_clock_time_self_energy - clock_rtime
!         else
!            gw_selfe_band1=0.
!            delta_gw_selfe_band1=0.
         end if

         if(myid.eq.0) then
            write(use_unit,'(2X,A,F5.1,A)') "finished ",100.*i_irkq_point_local/max_irkq_points_task, &
                 '% of GW computation'
         endif
         !     end of loop over i_irkq_point_local
      enddo
      if(myid.eq.0) write(use_unit,*) 

      call mpi_allreduce(checkloc, checkglob,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr)
      if (myid.eq.0) print*,'CPOUT: ', checkglob

      call finalize_compute_g_times_w_p0      
      count=n_freq*(n_high_state-n_low_state+1)*n_spin*n_irk_points

      if (irkblacs_member.and.(myid_bl.eq.0)) then
         call mpi_reduce(gw_selfenergy1,gw_selfenergy2,count,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm_irkq,mpierr)
         if(myid.eq.0) gw_selfenergy1=gw_selfenergy2
         call mpi_reduce(delta_gw_selfenergy1,gw_selfenergy2,count,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm_irkq,mpierr)
         if(myid.eq.0) delta_gw_selfenergy1=gw_selfenergy2
      end if
      call mpi_bcast(gw_selfenergy1,count,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,mpierr)
      call mpi_bcast(delta_gw_selfenergy1,count,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,mpierr)
      call perfoff

      do i_irk_point_local=1,n_irk_points_task
         i_irk_point=n_tasks*(i_irk_point_local-1)+mod(myid+n_tasks-1,n_tasks)+1      
         gw_selfenergy(:,:,:,i_irk_point_local)=gw_selfenergy1(:,:,:,i_irk_point)/2.d0/pi + delta_gw_selfenergy1(:,:,:,i_irk_point)
      end do

      if (myid.eq.1) then
         checkloc=sum(abs(gw_selfenergy(:,:,:,1)))
         print*,'GWOUT: ', checkloc
      end if

      if(out_band) then
         count=n_freq*(n_high_state-n_low_state+1)*n_spin*n_band_kpoints
         if (irkblacs_member.and.(myid_bl.eq.0)) then
            deallocate(gw_selfenergy2)
            allocate(gw_selfenergy2(n_freq,n_low_state:n_high_state,n_spin, n_band_kpoints))
            call mpi_reduce(gw_selfe_band1,gw_selfe_band2,count,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm_irkq,mpierr)
            if(myid.eq.0) gw_selfe_band1=gw_selfe_band2

            call mpi_reduce(delta_gw_selfe_band1,gw_selfe_band2,count,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm_irkq,mpierr)
            if(myid.eq.0) delta_gw_selfe_band1=gw_selfe_band2

         end if
         call mpi_bcast(gw_selfe_band1,count,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,mpierr)
         call mpi_bcast(delta_gw_selfe_band1,count,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,mpierr)

         do i_irk_point_local=1,n_band_kpoints_task
            i_irk_point=n_tasks*(i_irk_point_local-1)+mod(myid+n_tasks-1,n_tasks)+1      
            gw_selfe_band(:,:,:,i_irk_point_local)=gw_selfe_band1(:,:,:,i_irk_point)/2.d0/pi &
                 + delta_gw_selfe_band1(:,:,:,i_irk_point)
         end do

         if (myid.eq.1) then
            checkloc=sum(abs(gw_selfe_band(:,:,:,1)))
            print*,'GWBOUT: ', checkloc
         end if
      endif

!    printing out  
      if (out_self_energy) then

        do i_spin = 1, n_spin, 1
          do i_irk_point = 1, n_irk_points, 1

            if(myid.eq.mod(i_irk_point,n_tasks)) then
              i_irk_point_local = (i_irk_point-1)/n_tasks + 1
              i_k_point = inv_irk_point_mapping(i_irk_point)
              do i_state = n_low_state, n_high_state, 1
  
               if(n_spin.eq.2) then
                  write(filename,'(A,I0,A,I0,A,I0,A)') &
                     "self_energy/Sigma.omega.n_",  &
                        i_state,".s_",i_spin, ".k_", i_k_point, &
                        ".dat"
               else
                  write(filename,'(A,I0,A,I0,A)') &
                    "self_energy/Sigma.omega.n_", &
                        i_state, ".k_", i_k_point, ".dat"
               endif

               open(102, file=filename)
               do i_freq = 1, n_freq, 1
                  write(102,'(3f16.8)')omega(i_freq), gw_selfenergy(i_freq,i_state,i_spin,i_irk_point_local)
               enddo
               close(102)

              enddo
! end of if myid
           endif
          enddo
       
          i_index = 0
          do i_band = 1, n_plot_band, 1
             do i_k_point_band = 1, n_points_in_band(i_band), 1
               i_index = i_index + 1
               i_index_local = (i_index-1)/n_tasks + 1
               if(myid .eq. mod(i_index, n_tasks)) then
                 do i_state = n_low_state, n_high_state, 1
     
                  if(n_spin.eq.2) then
                     write(filename,'(A,I0,A,I0,A,I0,A,I0,A)') &
                        "self_energy/Sigma.omega.n_",  &
                           i_state,".s_",i_spin, ".band_", i_band, ".k_", i_k_point_band, &
                           ".dat"
                  else
                     write(filename,'(A,I0,A,I0,A,I0,A)') &
                       "self_energy/Sigma.omega.n_", &
                           i_state, ".band_", i_band, ".k_", i_k_point_band, ".dat"
                  endif
   
                  open(102, file=filename)
                  do i_freq = 1, n_freq, 1
                     write(102,'(3f16.8)')omega(i_freq), gw_selfe_band(i_freq,i_state,i_spin,i_index_local)
                  enddo
                  close(102)
   
                 enddo
               endif
             enddo
          enddo
! end of loop over i_spin
        enddo
!   end if out_self_energy
     endif

      if (allocated (n_homo_k)) then
        deallocate (n_homo_k)
      endif

      deallocate (polar_kspace_real)
      deallocate (polar_kspace_complex)
      
      if (allocated (coulomb_eigenvalues)) then
        deallocate (coulomb_eigenvalues)
      endif
      if (allocated (coulomb_eigenvectors)) then
        deallocate (coulomb_eigenvectors)
      endif
      if (allocated (sqrtv_eigenvectors)) then
        deallocate (sqrtv_eigenvectors)
      endif
      if (allocated (v_times_polar)) then
        deallocate (v_times_polar)
      endif
      if (allocated (gw_selfenergy1)) then
        deallocate (gw_selfenergy1)
      endif
      if (allocated (gw_selfenergy2)) then
         deallocate (gw_selfenergy2)
      endif
      if (allocated (delta_gw_selfenergy)) then
        deallocate (delta_gw_selfenergy)
      endif
      if (allocated (gw_selfe_band1)) then
        deallocate (gw_selfe_band1)
      endif
      if (allocated (gw_selfe_band2)) then
        deallocate (gw_selfe_band2)
      endif
      if (allocated (delta_gw_selfe_band1)) then
        deallocate (delta_gw_selfe_band1)
      endif
      call perfoff

      return

      end subroutine evaluate_periodic_gw_selfenergy
