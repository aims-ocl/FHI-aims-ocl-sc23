!****h* FHI-aims/evaluate_self_energies
!  NAME
!   Routines for printing self energies
!  SYNOPSIS

module print_self_energies

   use dimensions
   use physics,                      only: chemical_potential
   use gw_para,                      only: use_hedin_shift,&
                                           hedin_shift
                                             
   use runtime_choices
   use mpi_tasks
   use localorb_io,                  only: use_unit
   use prodbas                   
   use timing,                       only: time_print_sigma_cd,&
                                           clock_time_print_sigma_cd,&
                                           get_timestamps,&
                                           get_times
   use evaluate_self_energy_freq,    only: evaluate_self_energy_cd
   use contour_def_gw_types,         only: cd_environment_type
   use contour_def_gw_environment,   only: init_contour_def_env 
   use output_handling,              only: open_file, close_file


   implicit none
   private

   public :: print_self_energy_analytic, print_self_energy_cd
 
contains

! **************************************************************************************************

  subroutine print_self_energy_analytic            &
           (anacon_type,n_max_par, &
            n_freq, omega, sigma_par,&
            state_to_print, print_range )

!  ARGUMENTS

      integer :: anacon_type
      integer :: n_max_par
      integer :: n_freq
      integer :: state_to_print

      real*8                                 :: omega(n_freq)
      real(kind=8), dimension(2), intent(in) :: print_range
      complex*16                             :: sigma_par(n_max_par,n_states,n_spin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      character*40 :: filename
      character*40 :: spin_tag, state_string

      real*8 ims, res

      complex*16 selfe

      integer :: i_spin
      integer :: i_freq
      integer :: unit_number
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
      real*8 interval
      real*8 freq
      real*8 my_freq
      real*8 min_freq
      real*8 max_freq

      if(myid.eq.0)then
         write(use_unit,*) " Printing out the self-energy matrix element:",state_to_print 
         if(use_hedin_shift) then
            write(use_unit,*)
            do i_spin = 1, n_spin
               write(use_unit,"(T3,A31,F12.5,T47,A2)") "Self-energy printed with shift:", &
                     hedin_shift(state_to_print,i_spin)*hartree, "eV" 
            enddo
         endif
        write(use_unit,*)
        write(use_unit,'(2A)')"------------------------------------------", &
          "------------------------------------------------------------"
      endif
      spin_tag = ""
      write(state_string,*) state_to_print

      max_freq = print_range(2)/hartree 
      min_freq = print_range(1)/hartree
      if(max_freq < 6.5d0)  max_freq = 6.5d0 
      if(min_freq > -6.5d0) min_freq = -6.5d0 

      do i_spin = 1, n_spin
        if(n_spin > 1) then
          if(i_spin == 1) spin_tag = "up"
          if(i_spin == 2) spin_tag = "down"
        endif
        if (myid == 0) then
           call open_file(filename,&
                          file_action='write',&
                          front_name=spin_tag,&
                          middle_name='self_energy_analytic_state',&
                          end_name=state_string,&
                          extension='dat',&
                          unit_number=unit_number)
        endif

        interval = 6.5d0 /10000d0
        freq = min_freq - interval - chemical_potential  
        do while (freq < max_freq)
           freq = freq + interval
           my_freq = freq - hedin_shift(state_to_print,i_spin) 
           call get_real_selfenergy ( anacon_type, n_freq , omega, &
                 dcmplx(my_freq, 0.d0), n_max_par, &
                 sigma_par(1:n_max_par,state_to_print , i_spin),selfe)
     
           res = real(selfe) 
           ims = aimag(selfe)

         if(myid.eq.0)then
             write(unit_number,'(3F16.10)') (freq+chemical_potential)*hartree, res, ims
         endif
        enddo

        if (myid == 0) call close_file(unit_number)
 
      enddo!i_spin
!      print*, chemical_potential 

  end subroutine print_self_energy_analytic 

! **************************************************************************************************
!> brief print self energy
!  o gw_cd  --  contour deformation environment
!  o Wmn_freq_cd --  screened Coulomb matrix element in the MO basis, i.e. (mn|W(iomega)|mn)
!  o KS_eigenvalue -- real array,
!           the eigenvalues of the single-particle calculation. For DFT calculation,
!           this is the KS eigenvalue, but for HF calculation, this is then the HF
!  o occ_numbers -- occupation numbers of single-particle energy levels
!  o n_full_freq -- the number of frequency points for the screened Coulomb interaction W
!  o n_homo -- the HOMO level for each spin channel
!  o omega_full -- the Gauss-Legendre frequency grid for the screened Coulomb interaction
!  o womega_full -- the weight of the Gauss-Legendre frequency grid for the self-energy
!  o chemical_potential_spin -- the chemical potential for the i-th spin
!  o ovlp_3KS -- transformed 3-center integrals
!  o state_to_print -- state that should be printed
!  o print_range -- print in range of KS_value +- print_rang
! **************************************************************************************************
  subroutine print_self_energy_cd(gw_cd, Wmn_freq_cd, KS_eigenvalue, KS_eigenvalue_last, occ_numbers,&
                                  n_full_freq, n_homo, omega_full, womega_full, chemical_potential_spin,&
                                  ovlp_3KS, state_to_print, print_range)

     type(cd_environment_type)                                :: gw_cd
     complex(kind=8), dimension(ndim2_o3KS,ndim1_o3KS,&
       n_full_freq, n_spin), intent(in)                       :: Wmn_freq_cd
     real(kind=8), dimension(n_states,n_spin),&
       intent(in)                                             :: KS_eigenvalue
     real(kind=8), dimension(n_states,n_spin),&
       intent(in)                                             :: KS_eigenvalue_last
     real(kind=8), dimension(n_states,n_spin)                 :: occ_numbers
     integer, intent(in)                                      :: n_full_freq
     integer, dimension(n_spin), intent(in)                   :: n_homo
     real(kind=8), dimension(n_full_freq), intent(in)         :: omega_full
     real(kind=8), dimension(n_full_freq), intent(in)         :: womega_full
     real(kind=8), dimension(n_spin), intent(in)              :: chemical_potential_spin   
     real(kind=8), &
       dimension(n_basbas, ndim1_o3KS, ndim2_o3KS, n_spin),&
       intent(in)                                             :: ovlp_3KS
     integer, intent(in)                                      :: state_to_print
     real(kind=8), dimension(2), intent(in)                   :: print_range

     character(*), parameter :: func = 'print_self_energies'

     character(len=40)                                        :: filename, &
                                                                 spin_tag, state_string
     logical, dimension(n_spin)                               :: no_print
     integer                                                  :: i_spin, i_level, &
                                                                 index_contour_def, &
                                                                 i_state, unit_number
     integer, dimension(n_spin)                               :: n_first
     real(kind=8)                                             :: max_freq, min_freq, KS_value, &
                                                                 my_freq_real, freq_real, delta,&
                                                                 dist_KS_val

     call get_timestamps(time_print_sigma_cd, clock_time_print_sigma_cd)

     no_print = .false.
     do i_spin = 1, n_spin
        if(state_to_print < gw_cd%contour_def_start(i_spin) &
           .or. state_to_print > gw_cd%contour_def_end(i_spin)) then
          if(myid.eq.0)then
             write(use_unit,"(T3,2(A,I1))") &
                        "State is not included in contour deformation;"// &
                        " won't print ",state_to_print, " for spin ", i_spin
          endif 
          no_print(i_spin) = .true.
        endif 
     enddo
     if(all(no_print)) return
     if(myid.eq.0)then
        do i_spin = 1, n_spin
           if(no_print(i_spin)) cycle
           write(use_unit,"(T3,2(A,I3))") "Printing out the self-energy matrix element:",&
                                        state_to_print, " for spin", i_spin
           if(use_hedin_shift) then
             write(use_unit,*)
             index_contour_def = gw_cd%index_cd_levels(state_to_print,i_spin)
             write(use_unit,"(T3,A31,F12.5,T47,A2)") "Self-energy printed with shift:", &
                   gw_cd%shift(index_contour_def,i_spin)*hartree, "eV"
           endif
           write(use_unit,*)
           write(use_unit,'(2A)')"------------------------------------------", &
             "------------------------------------------------------------"
           write(use_unit,*)
        enddo
     endif

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

     spin_tag = ""
     write(state_string,*) state_to_print
     delta = 0.01d0/hartree
     dist_KS_val = 0.05d0/hartree
     if(use_ev_scgw.or.use_ev_scgw0) dist_KS_val = dist_KS_val/10.0d0 
     do i_spin = 1, n_spin
        if(no_print(i_spin)) cycle
        if(gw_cd%spin_channel(i_spin) /= i_spin)  cycle
        if(n_spin > 1) then
          if(i_spin == 1) spin_tag = "up"
          if(i_spin == 2) spin_tag = "down"
        endif
        if (myid == 0) then
           call open_file(filename,&
                          file_action='write',&
                          front_name=spin_tag,&
                          middle_name='self_energy_cd_state',&
                          end_name=state_string,&
                          extension='dat',&
                          unit_number=unit_number)
        endif

        KS_value = KS_eigenvalue_last(state_to_print,i_spin)
        index_contour_def = gw_cd%index_cd_levels(state_to_print,i_spin)

        !*** define print range 
        max_freq = print_range(2)/hartree 
        min_freq = print_range(1)/hartree
        freq_real = min_freq - delta 
        do while (freq_real < max_freq)
           freq_real = freq_real + delta
           !*** if the freq is close to KS value, contour deformation is not valid
           if(ABS(freq_real - KS_value) < dist_KS_val) cycle
           my_freq_real = freq_real - gw_cd%shift(index_contour_def,i_spin)
           call init_contour_def_env(gw_cd,my_freq_real,KS_eigenvalue_last,chemical_potential_spin,&
                                     state_to_print, i_spin, n_states, n_homo)
           call evaluate_self_energy_cd(gw_cd, Wmn_freq_cd, state_to_print,i_spin,&
                                        n_homo, n_first, occ_numbers, n_full_freq,&
                                        omega_full, womega_full, chemical_potential_spin, &
                                        my_freq_real, KS_eigenvalue, KS_eigenvalue_last, ovlp_3KS)
           if (myid == 0) then
              write(unit_number,'(3F16.10)') (freq_real)*hartree, &
                    gw_cd%self_energy%re(index_contour_def,i_spin), gw_cd%self_energy%im(index_contour_def,i_spin)
           endif 
           deallocate(gw_cd%residue_from_freq,gw_cd%real_freq)
        enddo

        if (myid == 0) call close_file(unit_number)
     enddo

     call get_times(time_print_sigma_cd, clock_time_print_sigma_cd)

  end subroutine print_self_energy_cd 

end module print_self_energies

