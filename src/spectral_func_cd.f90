!****h* FHI-aims/spectral_func_cd
!  NAME
!   Routines for calculating the spectral function 
!  SYNOPSIS

module spectral_func_cd

   use constants,                    only: hartree
   use dimensions
   use physics,                      only: chemical_potential
   use runtime_choices
   use mpi_tasks,                    only: aims_stop,&
                                           aims_warn,&
                                           n_tasks,&
                                           myid
   use localorb_io,                  only: use_unit
   use prodbas                   
   use timing,                       only: time_spec_func_cd, &
                                           clock_time_spec_func_cd,&
                                           get_timestamps,&
                                           get_times
   use evaluate_self_energy_freq,    only: evaluate_self_energy_cd
   use contour_def_gw_types,         only: cd_environment_type,&
                                           spectral_environment_type,&
                                           deallocate_cd_environment
   use contour_def_gw_environment,   only: basic_init_contour_def_env,& 
                                           init_contour_def_env 
   use synchronize_mpi,              only: sync_vector,&
                                           sync_vector_complex
   use output_handling,              only: open_file, close_file
   use gw_para,                      only: use_hedin_shift

   implicit none
   private

   public :: calc_spectral_func
 
contains

! **************************************************************************************************
!> brief calculation of spectral function for GW using the contour deformation
!  o gw_cd  --  contour deformation environment
!  o Wmn_freq_cd -- Wmn matrix elements
!                dim on entry: (ndim2_o3KS,ndim1_o3KS,n_full_freq,n_spin)
!                dim on exit:  (n_high_states, n_states, n_full_freq, n_spin)
!  o KS_eigenvalue -- KS/HF eigenvalues of the single-particle calculation, dim: (n_states,n_spin)
!  o exchange_self_energy -- the exact exchange part of the self-energy (KS/HF)
!                            dim: (n_high_state,n_spin) 
!  o xc_KS_matr -- the matrix elements of the exchange correlation (KS/HF) 
!                  dim: (n_states,n_states,n_spin) 
!  o x_KS_array -- the matrix elements of the exchange (KS/HF), dim: (n_states,n_spin)
!  o c_KS_array -- the matrix elements of the correlation (KS/HF), dim: (n_states,n_spin)
!  o n_high_state -- highest KS/HF eigenstate for self-energy correction
!  o n_full_freq -- the number of frequency points for numerical integration
!  o n_homo -- the HOMO level, i.e., the number of occupied state, dim: (n_spin)
!  o omega_full(n_freq) -- the Gauss-Legendre frequency grid for numerical integration,
!                          dim:n_full_freq
!  o womega_full(n_freq) -- the weigth of the Gauss-Legendre frequency grid for the numerical
!                           interaction, dim: (n_full_freq)
!  o chemical_potential -- the chemical potential of the system dim: (n_spin)
!  o ovlp_3KS -- transformed 3-center integrals, 
!                dim on entry: (n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
!                dim on exit:  (n_basbas, n_states, n_high_states, n_spin)
! **************************************************************************************************
  subroutine calc_spectral_func(gw_cd, Wmn_freq_cd, KS_eigenvalue, KS_eigenvalue_last, &
                                exchange_self_energy, xc_KS_matr, x_KS_array, c_KS_array, occ_numbers, &
                                n_high_state, n_full_freq, n_homo, omega_full, womega_full, &
                                chemical_potential_spin, ovlp_3KS)

     use constants, only: pi
     implicit none
     type(cd_environment_type)                                :: gw_cd
     complex(kind=8), dimension(:,:,:,:), allocatable, &
       intent(inout)                                          :: Wmn_freq_cd
     real(kind=8), dimension(:,:), intent(in)                 :: KS_eigenvalue
     real(kind=8), dimension(:,:), intent(in)                 :: KS_eigenvalue_last
     real(kind=8), dimension(:,:), intent(in)                 :: exchange_self_energy
     real(kind=8), dimension(:,:,:), intent(in)               :: xc_KS_matr
     real(kind=8), dimension(:,:), intent(in)                 :: x_KS_array      
     real(kind=8), dimension(:,:), intent(in)                 :: c_KS_array      
     real(kind=8), dimension(:,:), intent(in)                 :: occ_numbers
     integer, intent(in)                                      :: n_high_state,&
                                                                 n_full_freq
     integer, dimension(:), intent(in)                        :: n_homo
     real(kind=8), dimension(:), intent(in)                   :: omega_full
     real(kind=8), dimension(:), intent(in)                   :: womega_full
     real(kind=8), dimension(:), intent(in)                   :: chemical_potential_spin   
     real(kind=8), dimension(:,:,:,:), allocatable,&
       intent(inout)                                          :: ovlp_3KS

     character(*), parameter :: func = 'calc_spectral_func_cd'

     character(len=40)                                        :: filename, &
                                                                 spin_tag, method,&
                                                                 state_tag
     integer                                                  :: i_spin, i_freq,  &
                                                                 i_freq_global,&
                                                                 index_contour_def, &
                                                                 i_state, unit_number,&
                                                                 npoints_global
     integer, dimension(n_spin)                               :: n_first
     real(kind=8)                                             :: max_freq, min_freq,&
                                                                 freq_real, sgn, resolution,&
                                                                 my_freq_real
     real(kind=8)                                             :: spectrum
     real(kind=8), dimension(:,:), allocatable                :: my_exchange_self_energy, &
                                                                 xc_energy, my_hedin_shift
     complex(kind=8)                                          :: im_unit
     complex(kind=8), dimension (:), allocatable              :: aux_G 
  
     call get_timestamps(time_spec_func_cd, clock_time_spec_func_cd)

     !if(myid.eq.0)then
     !   write(use_unit,"(T3,A)") "Start calculating the spectral function" 
     !   write(use_unit,*)
     !   write(use_unit,'(2A)')"------------------------------------------", &
     !     "------------------------------------------------------------"
     !   write(use_unit,*)
     !endif
     !*** print some info
     call write_header(gw_cd,n_spin)

     im_unit = (0.0d0, 1.0d0)   

     allocate(my_exchange_self_energy(n_high_state,n_spin)) 
     if(use_hartree_fock .and. (.not.use_screx) &
           .and. .not. use_gw_and_hse) then
       my_exchange_self_energy(:,:) = exchange_self_energy * & 
                                (1.0d0-hybrid_coeff)
     else
       my_exchange_self_energy(:,:) = exchange_self_energy
     endif
     !
     allocate(my_hedin_shift(n_states,n_spin))
     call preserve_hedin_shift(gw_cd,my_hedin_shift,n_spin)
     call deallocate_cd_environment(gw_cd)
     if(.not.allocated(gw_cd%contour_def_start)) then
       allocate(gw_cd%contour_def_start(n_spin))
     endif
     if(.not.allocated(gw_cd%contour_def_end)) then
       allocate(gw_cd%contour_def_end(n_spin))
     endif
     gw_cd%contour_def_start(:) = 1
     gw_cd%contour_def_end(:) = n_states
     call basic_init_contour_def_env(gw_cd, n_low_state=1, n_high_state=n_states,&
          n_states=n_states,setup_qp_iteration=.false.,n_spin=n_spin)

     ! read in DFT exchange-correlation energy for all states
     allocate(xc_energy(n_states,n_spin))
     if(use_split_xc_gw) then
       do i_spin = 1, n_spin, 1
        do i_state = 1, n_states, 1
           xc_energy(i_state,i_spin) = &
           x_KS_array(i_state,i_spin)+c_KS_array(i_state,i_spin)
        enddo
       enddo
     else
       do i_spin = 1, n_spin, 1
        do i_state = 1, n_states, 1
           xc_energy(i_state,i_spin) = &
           xc_KS_matr(i_state,i_state,i_spin)
         enddo
       enddo
     endif

     allocate(aux_G(n_states))
     !*** determine the highest occupied orbital level
     n_first(:) = 1
     do i_spin = 1, n_spin
        do i_state = 1, n_states
           if (abs(occ_numbers(i_state,i_spin)-dble(2/n_spin)) &
                           .lt.1.d-8) then
            n_first(i_spin)= i_state + 1
           endif 
        enddo 
        if(n_first(i_spin) .gt. n_states) then
          n_first(i_spin) = n_states
        endif
     enddo

     resolution = gw_cd%spec_resolution/hartree
     max_freq = gw_cd%max_freq_spec/hartree
     min_freq = gw_cd%min_freq_spec/hartree

     call setup_spectral_env(gw_cd%spectral_env,max_freq, min_freq, resolution, KS_eigenvalue_last)
     call switch_parallelization(ovlp_3KS, Wmn_freq_cd, n_full_freq, n_high_state)

     do i_spin = 1, n_spin
       
        do i_freq = 1, gw_cd%spectral_env(i_spin)%npoints_local
           i_freq_global = gw_cd%spectral_env(i_spin)%freq_index(i_freq)
           freq_real = gw_cd%spectral_env(i_spin)%freqs(i_freq_global)
           sgn = sign(1.0d0,freq_real-chemical_potential_spin(i_spin))
 
           spectrum = 0.0d0   

           do i_state = 1, n_states
              if(gw_cd%calc_single_state_spec) then
                if(i_state /= gw_cd%state_spec) cycle
              endif
              index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
              my_freq_real = freq_real - my_hedin_shift(index_contour_def,i_spin)
              call init_contour_def_env(gw_cd,my_freq_real, KS_eigenvalue_last, &
                                        chemical_potential_spin, i_state, i_spin,&
                                        n_states, n_homo)
              call evaluate_self_energy_cd(gw_cd, Wmn_freq_cd, i_state, i_spin,&
                                           n_homo, n_first, occ_numbers, n_full_freq,&
                                           omega_full, womega_full, chemical_potential_spin, &
                                           my_freq_real, KS_eigenvalue, KS_eigenvalue_last, &
                                           ovlp_3KS, do_scalapack=.FALSE.)
              aux_G(i_state) = & 
                    -(my_exchange_self_energy(i_state, i_spin) &
                    - xc_energy(i_state, i_spin) &
                    + gw_cd%self_energy%complx(index_contour_def,i_spin)) &
                    + freq_real &
                    - KS_eigenvalue(i_state, i_spin) !&
                   ! - im_unit*gw_cd%eta &
                   !   *sign(1.0d0,chemical_potential_spin(i_spin)&
                   !         -KS_eigenvalue(i_state,i_spin))
              aux_G(i_state) =  1.0d0/aux_G(i_state)
              spectrum =  spectrum - 2.0d0*sgn*aimag(aux_G(i_state))/pi/real(n_spin,kind=8)
              deallocate(gw_cd%residue_from_freq,gw_cd%real_freq)
           enddo
           gw_cd%spectral_env(i_spin)%spectrum(i_freq_global) = spectrum
        enddo
        npoints_global = gw_cd%spectral_env(i_spin)%npoints_global
        call sync_vector(gw_cd%spectral_env(i_spin)%spectrum, npoints_global) 
     enddo

     !*** write to file
     spin_tag = ""
     if(gw_cd%calc_single_state_spec) then
       write(state_tag,*) gw_cd%state_spec
       method = "cd_state_" // trim(adjustl(state_tag))
     else
       method =  "cd"
     endif

     do i_spin = 1, n_spin
        if(n_spin > 1) then
          if(i_spin == 1) spin_tag = "up"
          if(i_spin == 2) spin_tag = "down"
        endif
        if (myid == 0)  then
           call open_file(filename,&
                          file_action='write',&
                          front_name=spin_tag,&
                          middle_name='spectrum',&
                          end_name=method,&
                          extension='dat',&
                          unit_number=unit_number)

           npoints_global = gw_cd%spectral_env(i_spin)%npoints_global
           do i_freq = 1, npoints_global
              freq_real = gw_cd%spectral_env(i_spin)%freqs(i_freq)        
              spectrum = gw_cd%spectral_env(i_spin)%spectrum(i_freq) 
              write(unit_number,'(F20.10, T25, F20.10)') (freq_real)*hartree, spectrum
           enddo
           call close_file(unit_number)
        endif
     enddo
     deallocate(aux_G)
     deallocate(my_exchange_self_energy, xc_energy, my_hedin_shift)

     call get_times(time_spec_func_cd, clock_time_spec_func_cd)

  end subroutine calc_spectral_func

! **************************************************************************************************
!> brief write header and write warnings oder stop the code
!  o gw_cd  --  contour deformation environment
!  o n_spin --  number of spins
! **************************************************************************************************
  subroutine write_header(gw_cd,n_spin)

     type(cd_environment_type)                                :: gw_cd
     integer                                                  :: n_spin

     character(len=40)                                        :: state_string
     integer                                                  :: i_spin,&
                                                                 i_state,&
                                                                 index_contour_def

     if(.not.gw_cd%calc_single_state_spec.and.use_hedin_shift.and.&
        (gw_cd%num_levels_tot /= n_states)) then
        call aims_warn("* Warning! Total spectral function is calculated in combination with" // &
                       " Hedin shift. Hedin shift not computed for all levels, which might be" // &
                       " inconsistent")
     endif

     if(gw_cd%calc_single_state_spec.and.use_hedin_shift) then
       do i_spin = 1, n_spin
          i_state = gw_cd%state_spec
          index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
          if(index_contour_def == 0) then
            write(state_string,*) i_state
            call aims_stop("No Hedin shift calculated for state " // trim(adjustl(state_string)))
          endif
       enddo
     endif
  
     if(myid.eq.0) then
        if(gw_cd%calc_single_state_spec) then
          i_state = gw_cd%state_spec
          write(use_unit,"(T3,A50,I3)") "Start calculating the spectral function for state:", i_state 
          if(use_hedin_shift) then
            write(use_unit,*)
            do i_spin = 1, n_spin
               index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
               write(use_unit,"(T3,A37,F12.5,T53,A2)") "Spectral function printed with shift:", &
                     gw_cd%shift(index_contour_def,i_spin)*hartree, "eV"
            enddo
          endif
        else
          write(use_unit,"(T3,A)") "Start calculating the spectral function" 
        endif
        write(use_unit,*)
        write(use_unit,'(2A)')"------------------------------------------", &
          "------------------------------------------------------------"
        write(use_unit,*)
     endif

  end subroutine write_header

! **************************************************************************************************
!> brief preserve Hedin-Shift
!  o gw_cd  --  contour deformation environment
!  o my_hedin_shift -- contains the Hedin shift obtained from previous QP calculations
!  o n_spin --  number of spins
! **************************************************************************************************
  subroutine preserve_hedin_shift(gw_cd,my_hedin_shift,n_spin)

     type(cd_environment_type)                                :: gw_cd
     real(kind=8), dimension(:,:), intent(inout)              :: my_hedin_shift
     integer                                                  :: n_spin

     integer                                                  :: i_spin,&
                                                                 i_state,&
                                                                 index_contour_def
     my_hedin_shift = 0.0d0

     if(use_hedin_shift) then
       do i_spin = 1, n_spin
          do i_state =1, n_states
             index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
             if(index_contour_def /= 0) then
               my_hedin_shift(i_state,i_spin) = gw_cd%shift(index_contour_def,i_spin)
             endif
          enddo 
       enddo 
     endif

  end subroutine preserve_hedin_shift

! **************************************************************************************************
!> brief switch parallelization scheme. we want to parallelize over frequenies points. However,
!>       at this point the ovlp_3KS and Wmn_freq are distributed (2D scheme). Sync these matrices
!>       and reset the 2D parallelization variables
!  o ovlp_3KS -- transformed 3-center integrals, 
!                dim on entry: (n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
!                dim on exit:  (n_basbas, n_states, n_high_states, n_spin)
!  o Wmn_matrix_freq -- Wmn matrix elements
!                dim on entry: (ndim2_o3KS,ndim1_o3KS,n_full_freq,n_spin)
!                dim on exit:  (n_high_states, n_states, n_full_freq, n_spin)
!  o n_full_freq -- the number of frequency points for numerical integration
!  o n_high_state -- highest KS/HF eigenstate for self-energy correction
! **************************************************************************************************
  subroutine switch_parallelization(ovlp_3KS, Wmn_freq, n_full_freq, n_high_states)

     real(kind=8), dimension(:,:,:,:), allocatable,&
        intent(inout)                                         :: ovlp_3KS
     complex(kind=8), dimension(:,:,:,:), allocatable,&
        intent(inout)                                         :: Wmn_freq
     integer, intent(in)                                      :: n_full_freq
     integer, intent(in)                                      :: n_high_states    

     character(*), parameter :: func = 'switch_parallelization'

     integer                                                  :: i, i_spin, &
                                                                 m1, m2, m3, m4
     integer                                                  :: j_state, k_state, &
                                                                 j_state_loc, k_state_loc 
     real(kind=8), dimension(:,:,:,:), allocatable            :: tmp_ovlp
     complex(kind=8), dimension(:,:,:,:), allocatable         :: tmp_Wmn

     !*** reallocate ovlp_3KS
     m1 = SIZE(ovlp_3KS,1)
     m2 = SIZE(ovlp_3KS,2)
     m3 = SIZE(ovlp_3KS,3)
     m4 = SIZE(ovlp_3KS,4)
   
     allocate(tmp_ovlp(m1,m2,m3,m4))
     tmp_ovlp(:,:,:,:) = ovlp_3KS
     deallocate(ovlp_3KS)
     allocate(ovlp_3KS(n_basbas, n_states, n_high_states, n_spin))
     ovlp_3KS(:,:,:,:) = 0.0d0 

     do i_spin = 1, n_spin
        do j_state = 1, n_high_states
           if(own_dim2_o3KS(j_state) /= myp2_o3KS) cycle
           j_state_loc = loc_dim2_o3KS(j_state)
           do k_state = 1, n_states
              if(own_dim1_o3KS(k_state) /= myp1_o3KS) cycle 
              k_state_loc = loc_dim1_o3KS(k_state)
              ovlp_3KS(:,k_state,j_state,i_spin) = &
                tmp_ovlp(:,k_state_loc,j_state_loc,i_spin)
           enddo
        enddo
     enddo

     call sync_vector(ovlp_3KS,n_basbas*n_states*n_high_states*n_spin)  
     deallocate(tmp_ovlp)

     !*** reallocate Wmn_freq
     m1 = SIZE(Wmn_freq,1)
     m2 = SIZE(Wmn_freq,2)
     m3 = SIZE(Wmn_freq,3)
     m4 = SIZE(Wmn_freq,4)
   
     allocate(tmp_Wmn(m1,m2,m3,m4))
     tmp_Wmn(:,:,:,:) = Wmn_freq
     deallocate(Wmn_freq)
     allocate(Wmn_freq(n_high_states, n_states, n_full_freq, n_spin))
     Wmn_freq(:,:,:,:) = (0.0d0, 0.0d0) 

     do i_spin = 1, n_spin
        do k_state = 1, n_states
           if(own_dim1_o3KS(k_state) /= myp1_o3KS) cycle 
           k_state_loc = loc_dim1_o3KS(k_state)
           do j_state = 1, n_high_states
              if(own_dim2_o3KS(j_state) /= myp2_o3KS) cycle
              j_state_loc = loc_dim2_o3KS(j_state)
              Wmn_freq(j_state,k_state,1:n_full_freq,i_spin) = &
                tmp_Wmn(j_state_loc,k_state_loc,1:n_full_freq,i_spin)
           enddo
        enddo
     enddo
   
     call sync_vector_complex(Wmn_freq,n_high_states*n_states*n_full_freq*n_spin)  

     deallocate(tmp_Wmn)

     !*** reset variables; this is a hack to switch from the 2D distribution

     ndim1_o3KS = n_states
     ndim2_o3KS = n_high_states
     
     max_row_2d = n_basbas 
     max_col_2d = n_basbas

     !*** myp1_o3KS and myp2_o3KS should be always identical to own_dim1_o3KS 
     !*** and own_dim2_o3KS, respectively (so that the cycle cond. are never
     !*** fullfilled)
     myp1_o3KS = 0
     myp2_o3KS = 0

     own_dim1_o3KS(:) = 0
     own_dim2_o3KS(:) = 0

     do i = 1, n_states
        loc_dim1_o3KS(i) = i
        loc_dim2_o3KS(i) = i
     enddo 

  end subroutine switch_parallelization

! **************************************************************************************************
!> brief allocate spectral env, get frequencies points and set up parallelization. Each proc 
!        receives a set of freq points of the spectral function. The self-energy for the given
!        freq point is then calculated locally by the respective proc.
!  o spectral_env spectral environment  
!  o max_freq -- maximal frequency calculated for spectral function
!  o min_freq -- minimal frequency calculated for spectral function
!  o resolution -- resolution of spectral function
!  o KS_eigenvalue -- KS/HF eigenvalues of the single-particle calculation, dim: n_states,n_spin
! **************************************************************************************************
  subroutine setup_spectral_env(spectral_env,max_freq, min_freq, resolution, KS_eigenvalue)

     type(spectral_environment_type), dimension(:), &
      allocatable                                             :: spectral_env 
     real(kind=8), intent(in)                                 :: max_freq,&
                                                                 min_freq,&
                                                                 resolution
     real(kind=8), dimension(:,:), intent(in)                 :: KS_eigenvalue

     character(*), parameter :: func = 'setup_spectral_env'

     logical                                                  :: exclude_omega
     integer                                                  :: i_spin, i_state, &
                                                                 npoints, i_freq, i_index
     real(kind=8)                                             :: freq_real, KS_value
     
     allocate(spectral_env(n_spin)) 

     !*** find the frequency points, where we want to evaluate the function; exclude freqs
     !    that are too close to the KS/HF values; there the contour deformation is not valid
     do i_spin = 1, n_spin
        spectral_env(i_spin)%npoints_global = 0
        freq_real = min_freq - resolution
        do while (freq_real < max_freq - resolution)
           freq_real = freq_real + resolution
           exclude_omega = .false. 
           do i_state = 1, n_states 
              KS_value = KS_eigenvalue(i_state,i_spin) 
              if(ABS(freq_real - KS_value) < 0.005d0/hartree) then
                 exclude_omega = .true.
                 exit
              endif
           enddo
           if(exclude_omega) cycle
           spectral_env(i_spin)%npoints_global = spectral_env(i_spin)%npoints_global + 1
        enddo 
     enddo

     !*** allocate arrays
     do i_spin = 1, n_spin
        npoints = spectral_env(i_spin)%npoints_global
        allocate(spectral_env(i_spin)%freqs(npoints))
        allocate(spectral_env(i_spin)%spectrum(npoints))
        spectral_env(i_spin)%freqs(:) = 0.0d0
        spectral_env(i_spin)%spectrum(:) = 0.0d0
     enddo

     !*** get frequencies
     do i_spin = 1, n_spin
        freq_real = min_freq - resolution
        i_freq = 0
        do while (freq_real < max_freq - resolution)
           freq_real = freq_real + resolution
           exclude_omega = .false. 
           do i_state = 1, n_states 
              KS_value = KS_eigenvalue(i_state,i_spin) 
              if(ABS(freq_real - KS_value) < 0.005d0/hartree) then
                 exclude_omega = .true.
                 exit
              endif
           enddo
           if(exclude_omega) cycle
           i_freq = i_freq + 1
           spectral_env(i_spin)%freqs(i_freq) = freq_real
        enddo 
     enddo

     !*** set-up parallelization
     do i_spin = 1, n_spin
        npoints = spectral_env(i_spin)%npoints_global
        allocate(spectral_env(i_spin)%freq_index(INT(npoints/n_tasks)+1))
        i_index = 0 
        do i_freq = 1, npoints
           if(mod(i_freq, n_tasks) == myid) then
             i_index = i_index + 1
             spectral_env(i_spin)%freq_index(i_index) = i_freq
           endif
        enddo
        spectral_env(i_spin)%npoints_local = i_index
     enddo

  end subroutine setup_spectral_env

end module spectral_func_cd
