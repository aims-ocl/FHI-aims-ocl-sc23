! ****h* FHI-aims/quasi_particle_energies
!  NAME
!   Routines for evaluating the quasi_particle_energies
!  SYNOPSIS

module quasi_particle_energies

   use dimensions
   use runtime_choices
   use mpi_tasks
   use constants
   use prodbas     
   use timing,                       only: time_qp_energy_cd,&
                                           clock_time_qp_energy_cd,&
                                           get_timestamps,&
                                           get_times,&
                                           warn_qpe_nonconvergence
   use localorb_io,                  only: use_unit, localorb_info
   use contour_def_gw_types,         only: cd_environment_type
   use contour_def_gw_environment,   only: init_contour_def_env,&
                                           print_restart_information,&
                                           read_restart_information
   use evaluate_self_energy_freq,    only: evaluate_self_energy_cd
   use gw_para,                      only: use_hedin_shift,&
                                           hedin_shift, &
                                           state_to_shift,&
                                           z_value

   implicit none

   private
  
   public :: quasi_particle_energy, quasi_particle_energy_zshot,&
             quasi_particle_energy_contour_def,&
             quasi_particle_energy_contour_def_zshot,&
             get_exchange_terms,&
             output_quasi_particle_energies


contains

! **************************************************************************************************
!> brief iterative calcuation of the quasi particle energies for the GW calculation.  Note, now the
!        self-energy has been calculated and  ontinuted to the real frequency axis 
! **************************************************************************************************
  subroutine quasi_particle_energy &
       (anacon_type,n_max_par, &
        n_low_state,n_high_state, &
        n_homo, n_freq, omega, &
        sigma_par,occ_numbers,KS_eigenvalue, &
        chemical_potential_spin, &
        exchange_self_energy, &
        xc_KS_matr,x_KS_array,c_KS_array, &
        qp_energy, correl_energy, qp_non_convergence )

!
!  ARGUMENTS

      integer :: anacon_type 
      integer :: n_max_par
      integer :: n_low_state, n_high_state
      integer :: n_freq
      integer :: n_homo(n_spin)

      real*8 :: omega(n_freq)
      real*8 :: occ_numbers(n_states,n_spin)
      real*8 :: KS_eigenvalue(n_states,n_spin)
      real*8 :: exchange_self_energy(n_high_state,n_spin)
      real*8 :: xc_KS_matr(n_states,n_states,n_spin)
      real*8 :: x_KS_array(n_states,n_spin)
      real*8 :: c_KS_array(n_states,n_spin)
      real*8 :: chemical_potential_spin(n_spin)

      complex*16 :: sigma_par(n_max_par,n_states,n_spin)

      real*8 :: qp_energy(n_low_state:n_high_state,n_spin)
      real(kind=8),  intent(inout), optional  :: correl_energy(n_states,n_spin)
      logical, dimension(:,:), intent(inout), &
        optional :: qp_non_convergence(n_low_state:n_high_state,n_spin)

!  INPUTS
!  o  anacon_type -- integer number, if 0, the two-pole fitting for analytical
!          continuation; if 1, using Pade approximation for ana. cont.        
!  o  n_max_par -- the number of parameters used for analytical continuation 
!          For anacon_type = 0, recommended n_max_par is  4 or 5. If 4, this will 
!          be the normal two-pole fitting, else if 5, it will be two-pole plus a 
!          (small) constant number
!          For anacon_type = 1, recommended n_max_par is the half of n_freq  
!  o  n_low_state  -- integer number,
!          the lowest KS/HF eigenstate for self-energy correction
!  o  n_high_state -- integer number,
!          the highest KS/HF eigenstate for self-energy correction
!  o  n_freq -- integer number, the number of frequency points for the GW self-energy
!  o  n_homo -- integer array, the HOMO level for each spin channel
!  o  omega(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the self-energy
!  o  occ_numbers -- occupation numbers of single-particle energy levels
!  o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
!  o  chemical_potential -- real number, the chemical potential of the system
!  o  exchange_self_energy -- real array the exact exchange part of the self-energy 
!            of each KS/HF state for each spin channel
!  o  xc_KS_matr -- real array, the matrix elements of the exchange correlation 
!            potential witin KS/HF orbitals
!  o  x_KS_array -- real array, the matrix elements of the exchange
!            potential witin KS/HF orbitals
!  o  c_KS_array -- real array, the matrix elements of the correlation
!            potential witin KS/HF orbitals
!  o  sigma_par -- complex array, the fitting parameters from analytical continuation
! OUTPUTS
!  o  qp_energy -- real array, the final calculated quasiparticle energy for each
!  o  correl_energy -- real array, correlation part of self-energy
!            concerned state (from n_low_state to n_high_state).

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

      real*8, allocatable :: xc_energy(:,:)
      real*8, allocatable :: qp_energy_old(:,:)
      real(kind=8), dimension(:,:), allocatable :: my_exchange_self_energy
      real(kind=8), dimension(:,:), allocatable :: my_correl_energy
      logical, dimension(:,:), allocatable :: my_qp_non_convergence

      real*8  e_diff
      real*8  my_shift
      real*8  en
      real*8  mu

!      real*8 :: qpe_sum

      complex*16 selfe

      real*8, parameter ::  qp_energy_thr = 1.d-5

!     counters


      integer :: i_state,i_state_1
      integer :: my_state, my_state_1
      integer :: i_count
      integer :: i_spin

!   external function


!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A)') "Quasi particle energy calculation using "// &
                                 "analytic continuation starts..."
      endif

      allocate(my_exchange_self_energy(n_high_state,n_spin)) 
      if(use_hartree_fock .and. (.not.use_screx) &
            .and. .not. use_gw_and_hse) then
        my_exchange_self_energy(:,:) = exchange_self_energy * & 
                                 (1.0d0-hybrid_coeff)
      else
        my_exchange_self_energy(:,:) = exchange_self_energy
      endif

!  read in DFT exchange-correlation energy for all states
      allocate(xc_energy(n_states,n_spin))
!      open (106, file='xc_energy.dat')
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
!         read(106, *) i, xc_energy(i_state,i_spin)
            xc_energy(i_state,i_spin) = &
            xc_KS_matr(i_state,i_state,i_spin)
            
!         if (i.ne.i_state) then
!          write(use_unit,'(2X,A)') "Wrong in reading the xc energy data file"
!          stop
!         endif
          enddo
        enddo
      endif
!      close(106)

!  quasi particle energy calculation

      allocate(qp_energy_old(n_states,n_spin))
      allocate(my_correl_energy(n_states,n_spin))
      allocate(my_qp_non_convergence(n_low_state:n_high_state,n_spin))
      my_qp_non_convergence(n_low_state:n_high_state,:) = .false.

!   determine homo level
      qp_energy(n_low_state:n_high_state,1:n_spin)= &
              KS_eigenvalue(n_low_state:n_high_state,1:n_spin)
      qp_energy_old(:,:)=KS_eigenvalue(:,:)
      do i_spin = 1, n_spin
        do i_state = n_low_state, n_high_state

          i_state_1 = i_state - n_low_state + 1

          e_diff = 1.d-3
          i_count =0

          !*** calculate Hedin shift
          my_shift = 0.0d0
          if(use_hedin_shift) then
            my_state = i_state
            if(state_to_shift > 0) my_state = state_to_shift
            my_state_1 = my_state - n_low_state + 1
            en =  KS_eigenvalue(my_state,i_spin)-chemical_potential_spin(i_spin)
            call get_real_selfenergy(anacon_type,n_freq,omega, &
                      dcmplx(en,0.d0), n_max_par, &
                      sigma_par(1:n_max_par,my_state_1,i_spin), selfe)
            my_shift = real(selfe) &
                      + my_exchange_self_energy(my_state,i_spin) &
                      - xc_energy(my_state,i_spin)
            hedin_shift(i_state,i_spin) = my_shift
            if(myid == 0) write(use_unit, '(T3,A,T10,F10.4,A13,I3)') "Shift:",&
                          my_shift*hartree, " eV for state", i_state
          endif

          do while (abs(e_diff).gt.qp_energy_thr)
            i_count = i_count +1
            qp_energy(i_state,i_spin) = &
                  qp_energy_old(i_state,i_spin) + 0.5d0* e_diff
            qp_energy_old(i_state,i_spin) = qp_energy(i_state,i_spin)

            mu =  chemical_potential_spin(i_spin)

            en = qp_energy(i_state,i_spin)-mu
            en = en - my_shift

            call get_real_selfenergy(anacon_type,n_freq,omega, &
                      dcmplx(en,0.d0), n_max_par, &
                      sigma_par(1:n_max_par,i_state_1,i_spin), selfe)

            qp_energy(i_state,i_spin)= KS_eigenvalue(i_state,i_spin) &
                      + real(selfe) &
                      + my_exchange_self_energy(i_state,i_spin) &
                      - xc_energy(i_state,i_spin)

            e_diff =  qp_energy(i_state,i_spin) &
                     - qp_energy_old(i_state,i_spin)

            if(i_count .gt. 100) then
               if(myid == 0) then
                  write(use_unit,'(2X,2A,I4,A,I4 )') &
                 "QUASI_PARTILCE_ENERGY: self-consistent quasiparticle ", &
                 "solution can not  be found for state: i_state = ", &
                 i_state, "  i_spin = ", i_spin
               endif
               my_qp_non_convergence(i_state,i_spin)=.true.
               warn_qpe_nonconvergence = .true.

              exit
            endif

! end of do while
          enddo
          my_correl_energy(i_state,i_spin) = real(selfe)
          if(present(qp_non_convergence)) then
             qp_non_convergence(i_state,i_spin) = my_qp_non_convergence(i_state,i_spin)
          endif
! end of do i_state
        enddo
! end of do i_spin
      enddo


      if(.not.(use_contour_def_gw)) then 
        call output_quasi_particle_energies(qp_energy,KS_eigenvalue,my_exchange_self_energy,&
                                            x_KS_array,c_KS_array,my_correl_energy,xc_energy,&
                                            occ_numbers,n_low_state,n_high_state,n_spin,n_states,&
                                            n_homo,qp_non_convergence=my_qp_non_convergence)
      endif

      if(present(correl_energy)) then
        correl_energy = my_correl_energy
      endif

      deallocate(my_exchange_self_energy)
      deallocate(xc_energy)
      deallocate(qp_energy_old)
      deallocate(my_correl_energy)
      deallocate(my_qp_non_convergence)

   end subroutine quasi_particle_energy

! **************************************************************************************************
!> brief evalutes the quasi particle energies based on the tailor expansion of the quasiparticle
!>       energies for the GW calculation. Note, now the self-energy has been calculated and
!>       continued to the real frequency axisi
!  o  anacon_type -- integer number, if 0, the two-pole fitting for analytical
!          continuation; if 1, using Pade approximation for ana. cont.        
!  o  n_max_par -- the number of parameters used for analytical continuation 
!          For anacon_type = 0, recommended n_max_par is  4 or 5. If 4, this will 
!          be the normal two-pole fitting, else if 5, it will be two-pole plus a 
!          (small) constant number
!          For anacon_type = 1, recommended n_max_par is the half of n_freq  
!  o  n_low_state  -- integer number,
!          the lowest KS/HF eigenstate for self-energy correction
!  o  n_high_state -- integer number,
!          the highest KS/HF eigenstate for self-energy correction
!  o  n_freq -- integer number, the number of frequency points for the GW self-energy
!  o  n_homo -- integer array, the HOMO level for each spin channel
!  o  omega(n_freq) -- the Gauss-Legendre frequency grid for the self-energy
!  o  occ_numbers -- occupation numbers of single-particle energy levels
!  o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
!  o  chemical_potential -- real number, the chemical potential of the system
!  o  exchange_self_energy -- real array the exact exchange part of the self-energy 
!            of each KS/HF state for each spin channel
!  o  xc_KS_matr -- real array, the matrix elements of the exchange correlation 
!            potential witin KS/HF orbitals
!  o  x_KS_array -- real array, the matrix elements of the exchange
!            potential witin KS/HF orbitals
!  o  c_KS_array -- real array, the matrix elements of the correlation
!            potential witin KS/HF orbitals
!  o  sigma_par -- complex array, the fitting parameters from analytical continuation
!  o  qp_energy -- real array, the final calculated quasiparticle energy for each
!  o  correl_energy -- real array, correlation part of self-energy
! **************************************************************************************************
  subroutine quasi_particle_energy_zshot &
              (anacon_type,n_max_par, &
               n_low_state,n_high_state, &
               n_homo, n_freq, omega, &
               sigma_par,occ_numbers,KS_eigenvalue, &
               chemical_potential_spin, &
               exchange_self_energy, &
               xc_KS_matr,x_KS_array,c_KS_array, &
               qp_energy, correl_energy )

     integer, intent(in)                                :: anacon_type 
     integer, intent(in)                                :: n_max_par
     integer, intent(in)                                :: n_low_state, n_high_state
     integer, intent(in)                                :: n_freq
     integer, intent(in)                                :: n_homo(n_spin)
     real(kind=8), intent(in)                           :: omega(n_freq)
     real(kind=8), intent(in)                           :: occ_numbers(n_states,n_spin)
     real(kind=8), intent(in)                           :: KS_eigenvalue(n_states,n_spin)
     real(kind=8), intent(inout)                        :: exchange_self_energy(n_high_state,n_spin)
     real(kind=8), intent(in)                           :: xc_KS_matr(n_states,n_states,n_spin)
     real(kind=8), intent(in)                           :: x_KS_array(n_states,n_spin)
     real(kind=8), intent(in)                           :: c_KS_array(n_states,n_spin)
     real(kind=8), intent(in)                           :: chemical_potential_spin(n_spin)
     complex(kind=8), intent(in)                        :: sigma_par(n_max_par,n_states,n_spin)
     real(kind=8), intent(inout)                        :: qp_energy(n_low_state:n_high_state,n_spin)
     real(kind=8),  intent(inout)                       :: correl_energy(n_states,n_spin)

     character(*), parameter :: func = 'quasi_particle_energy_zshot'

     real(kind=8)                                       :: en
     real(kind=8)                                       :: mu
     real(kind=8), dimension(:,:), allocatable          :: my_exchange_self_energy
     real(kind=8), dimension(:,:), allocatable          :: xc_energy
     complex(kind=8)                                    :: selfe
     integer                                            :: i_state,i_state_1
     integer                                            :: i_spin

!    begin work

     if(myid.eq.0) then
       write(use_unit,*)
       write(use_unit,'(2X,A)') "Z-shot quasi particle energy calculation using "// &
                                 "analytic continuation starts..."
     endif
    
     if(.not.allocated(z_value)) then
       allocate(z_value(n_states,n_spin))
       z_value = 0.0d0
     endif
     allocate(my_exchange_self_energy(n_high_state,n_spin)) 
     if(use_hartree_fock .and. (.not.use_screx) &
           .and. .not. use_gw_and_hse) then
       my_exchange_self_energy(:,:) = exchange_self_energy * & 
                                (1.0d0-hybrid_coeff)
     else
       my_exchange_self_energy(:,:) = exchange_self_energy
     endif

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

    ! quasi particle energy calculation

     qp_energy(n_low_state:n_high_state,1:n_spin)= &
             KS_eigenvalue(n_low_state:n_high_state,1:n_spin)
     do i_spin = 1, n_spin
        do i_state = n_low_state, n_high_state

           i_state_1 = i_state - n_low_state + 1


           mu =  chemical_potential_spin(i_spin)

           en = qp_energy(i_state,i_spin)-mu

           call get_real_selfenergy(anacon_type, n_freq, omega, &
                     dcmplx(en,0.d0), n_max_par, &
                     sigma_par(1:n_max_par,i_state_1,i_spin), selfe)

           call get_z_value(anacon_type,n_freq, omega, dcmplx(en,0.d0), n_max_par, &
                     sigma_par(1:n_max_par,i_state_1,i_spin), z_value(i_state,i_spin))

           qp_energy(i_state,i_spin)= KS_eigenvalue(i_state,i_spin) &
                      + z_value(i_state,i_spin)*(real(selfe) &
                      + my_exchange_self_energy(i_state,i_spin) &
                      - xc_energy(i_state,i_spin))

           correl_energy(i_state,i_spin) = z_value(i_state,i_spin)*real(selfe)
           my_exchange_self_energy(i_state,i_spin) = z_value(i_state,i_spin)&
                                                     *my_exchange_self_energy(i_state,i_spin) 
           xc_energy(i_state,i_spin) = z_value(i_state,i_spin)*xc_energy(i_state,i_spin)
        enddo ! end of do i_state
     enddo  ! end of do i_spin

     if(.not.(use_contour_def_gw)) then 
       call output_quasi_particle_energies(qp_energy,KS_eigenvalue,my_exchange_self_energy,&
                                           x_KS_array,c_KS_array,correl_energy,xc_energy,&
                                           occ_numbers,n_low_state,n_high_state,n_spin,n_states,&
                                           n_homo,z_value=z_value)
     endif

     deallocate (xc_energy)

  end subroutine quasi_particle_energy_zshot

! **************************************************************************************************
!> brief get re-normalization factor Z_n for tailor expansion of the the quasiparticle energies
!  o flag_ana -- flag_ana = 0, then P_N(x) = (a1 + a2*x + a3*x^2 + ...) /(1 + a4*x + a5*x^2 + ...)
!                flag_ana = 1, then P_N(x) = a1/(1 + a2(x-x[1])/(1 + a3(x-x[2])/(1 + ...)) )
!  o n_freq -- the number of frequency points for the GW self-energy
!  o omega(n_freq) -- the Gauss-Legendre frequency grid for the self-energy
!  o x -- is an complex input number. 
!  o npar -- is a positive integer input variable set to the number of parameters
!            in the fitting function
!  o par --  is a complex input array of length npar. On input par contains the 
!            values of the fitting parameter
!  o z_value -- renormalization factor
! **************************************************************************************************
  subroutine  get_z_value(flag_ana, n_freq, omega, x, npar, par, z_value)

     integer, intent(in)                           :: flag_ana
     integer, intent(in)                           :: n_freq
     real(kind=8), dimension(:), intent(in)        :: omega
     integer, intent(in)                           :: npar 
     complex(kind=8), intent(in)                   :: x
     complex(kind=8), intent(in)                   :: par(npar)
     real(kind=8), intent(out)                     :: z_value

     character(*), parameter :: func = 'get_z_value'

     integer                                       :: fac_count
     integer                                       :: i_dat, i_par, n_step 
     integer, dimension(:), allocatable            :: dev_factor
     complex(kind=8)                               :: xtmp, re_unit
     complex(kind=8)                               :: pade_numerator
     complex(kind=8)                               :: pade_denominator
     complex(kind=8)                               :: dev_pade_numerator
     complex(kind=8)                               :: dev_pade_denominator
     complex(kind=8)                               :: temp_z_value
     complex(kind=8)                               :: gtmp, dev_gtmp
     complex(kind=8), dimension(:), allocatable    :: xdata

     xtmp = x

     if(flag_ana.eq.0) then

        if(mod(npar,2) .eq.0)  then
          call aims_stop('anacon_type=0: expected uneven number of parameters')
        endif
        allocate(dev_factor(npar))
        pade_numerator = dcmplx(0.d0,0.d0)
        pade_denominator = dcmplx(0.d0,0.d0)
        dev_pade_numerator = dcmplx(0.d0,0.d0)
        dev_pade_denominator = dcmplx(0.d0,0.d0)
        
        fac_count = 0
        do i_par =2, npar
           fac_count = fac_count + 1
           dev_factor(i_par) = fac_count
           if(fac_count == int(npar/2)) fac_count = 0
        enddo
 
        do i_par = int(npar/2)+1, 2, -1 
          pade_numerator = (pade_numerator + par(i_par))*xtmp
        enddo

        do i_par = npar, int(npar/2)+2, -1 
          pade_denominator = (pade_denominator + par(i_par))*xtmp
        enddo

        do i_par = int(npar/2)+1, 3, -1 
          dev_pade_numerator = (dev_pade_numerator + dev_factor(i_par)*par(i_par))*xtmp
        enddo

        do i_par = npar, int(npar/2)+3, -1 
          dev_pade_denominator = (dev_pade_denominator + dev_factor(i_par)*par(i_par))*xtmp
        enddo

        pade_numerator = pade_numerator + par(1)
        pade_denominator = pade_denominator + dcmplx(1.d0,0.d0)
        dev_pade_numerator = dev_pade_numerator + par(2)
        dev_pade_denominator = dev_pade_denominator + par(int(npar/2)+2)

        temp_z_value = dev_pade_numerator/pade_denominator &
                  - (pade_numerator*dev_pade_denominator)/(pade_denominator**2)

        z_value = 1.0d0 - real(temp_z_value)
        z_value = 1.0d0/z_value
      
        deallocate(dev_factor)
 
     else

        allocate(xdata(npar))

        re_unit = (1.0d0, 0.0d0)

        n_step = n_freq/(npar-1)
        i_dat = 1
        do i_par = 1, npar-1, 1
          xdata(i_par) = dcmplx(0.d0,omega(i_dat))
          i_dat = i_dat + n_step
        enddo
        xdata(npar) = dcmplx(0.d0,omega(n_freq))

        gtmp = re_unit
        dev_gtmp = (0.0d0, 0.0d0)
        do i_par = npar, 2, -1
           pade_numerator = par(i_par)*(re_unit*x-xdata(i_par-1))
           dev_pade_numerator = par(i_par)*re_unit
           pade_denominator = gtmp
           dev_pade_denominator = dev_gtmp
           dev_gtmp = dev_pade_numerator/pade_denominator&
                     -(pade_numerator*dev_pade_denominator)/(pade_denominator**2)
           gtmp = re_unit+par(i_par)*(re_unit*x-xdata(i_par-1))/gtmp
        enddo

        dev_gtmp = -1.0d0*par(1)/(gtmp**2)*dev_gtmp

        z_value = 1.0d0 - real(dev_gtmp)
        z_value = 1.0d0/z_value 
       
        deallocate(xdata) 
     endif

     return

  end subroutine get_z_value

! **************************************************************************************************
!> brief evaluates quasi-particle energies using the contour deformation
!  o gw_cd  --  contour deformation environment
!  o Wmn_freq_cd --  screened Coulomb matrix element in the MO basis, i.e. (mn|W(iomega)|mn)
!  o qp_energy -- real array, the final calculated quasiparticle energy for each
!  o correl_energy -- real array, correlation part of self-energy
!  o KS_eigenvalue -- real array,
!           the eigenvalues of the single-particle calculation. For DFT calculation,
!           this is the KS eigenvalue, but for HF calculation, this is then the HF
!           eigenvalue
!  o KS_eigenvalue_last -- real array,
!           initially the eigenvalues of the single-particle calculation. For GW0, these
!           are then overwritten by the qp energies at the respective step
!  o exchange_self_energy -- real array the exact exchange part of the self-energy 
!           of each KS/HF state for each spin channel
!  o xc_KS_matr -- real array, the matrix elements of the exchange correlation 
!           potential witin KS/HF orbitals
!  o x_KS_array -- real array, the matrix elements of the exchange
!           potential witin KS/HF orbitals
!  o c_KS_array -- real array, the matrix elements of the correlation
!           potential witin KS/HF orbitals
!  o occ_numbers -- occupation numbers of single-particle energy levels
!  o n_low_state  -- the lowest KS/HF eigenstate
!  o n_high_state  -- the highest KS/HF eigenstate for self-energy correction
!  o n_states -- number of KS/HF eigenstates
!  o n_homo -- the HOMO level for each spin channel
!  o n_full_freq -- the number of frequency points for the screened Coulomb interaction W
!  o omega_full -- the Gauss-Legendre frequency grid for the screened Coulomb interaction
!  o womega_full -- the weight of the Gauss-Legendre frequency grid for the self-energy
!  o chemical_potential -- the chemical potential of the system
!  o ovlp_3KS -- transformed 3-center integrals
!  o qp_non_convergence -- flag whether the qp calculation from the previous calculation with
!                          analytic continuation has converged 
! **************************************************************************************************
  subroutine quasi_particle_energy_contour_def(gw_cd, Wmn_freq_cd, qp_energy, correl_energy, &
                                               KS_eigenvalue, KS_eigenvalue_last, exchange_self_energy,&
                                               xc_KS_matr,  x_KS_array, c_KS_array, occ_numbers,&
                                               n_low_state, n_high_state, n_homo, n_full_freq, &
                                               omega_full, womega_full, chemical_potential_spin, &
                                               ovlp_3KS, qp_non_convergence, step_sc_loop)
    
     type(cd_environment_type)                                :: gw_cd
     complex(kind=8), dimension(ndim2_o3KS,ndim1_o3KS,&
       n_full_freq, n_spin), intent(in)                       :: Wmn_freq_cd
     real(kind=8), dimension(n_low_state:n_high_state,&
       n_spin), intent(inout)                                 :: qp_energy
     real(kind=8), dimension(n_states,n_spin), intent(inout)  :: correl_energy
     real(kind=8), dimension(n_states,n_spin),&
       intent(in)                                             :: KS_eigenvalue
     real(kind=8), dimension(n_states,n_spin),&
       intent(in)                                             :: KS_eigenvalue_last
     real(kind=8), dimension(n_high_state,n_spin),&
       intent(in)                                             :: exchange_self_energy
     real(kind=8), dimension(n_states,n_states,n_spin), &
       intent(in)                                             :: xc_KS_matr
     real(kind=8), dimension(n_states,n_spin), & 
       intent(in)                                             :: x_KS_array      
     real(kind=8), dimension(n_states,n_spin), & 
       intent(in)                                             :: c_KS_array      
     real(kind=8), dimension(n_states,n_spin)                 :: occ_numbers
     integer, intent(in)                                      :: n_low_state,n_high_state,&
                                                                 n_full_freq
     integer, dimension(n_spin), intent(in)                   :: n_homo
     real(kind=8), dimension(n_full_freq), intent(in)         :: omega_full
     real(kind=8), dimension(n_full_freq), intent(in)         :: womega_full
     real(kind=8), dimension(n_spin), intent(in)              :: chemical_potential_spin   
     real(kind=8), &
       dimension(n_basbas, ndim1_o3KS, ndim2_o3KS, n_spin),&
       intent(in)                                             :: ovlp_3KS
     logical, dimension(n_low_state:n_high_state,n_spin),&
       intent(inout)                                          :: qp_non_convergence
     integer, intent(in)                                      :: step_sc_loop

     character(*), parameter :: func = 'quasi_particle_energy_contour_def'

     character(len=10)                                        :: spin_tag
     logical                                                  :: already_converged
     integer                                                  :: i_spin, i_state, &
                                                                 index_contour_def, &
                                                                 i_count, i_level,&
                                                                 my_state
     integer, dimension(n_spin)                               :: n_first
     real(kind=8)                                             :: e_diff, my_shift, my_freq
     real(kind=8)                                             :: qp_energy_thr = 1.d-5  
     real(kind=8), dimension(:,:), allocatable                :: my_exchange_self_energy
     real(kind=8), dimension(:,:), allocatable                :: xc_energy
     real(kind=8), dimension(:,:), allocatable                :: qp_energy_old
    
     call get_timestamps(time_qp_energy_cd, clock_time_qp_energy_cd)

     if(myid.eq.0) then
       write(use_unit,*)
       write(use_unit,'(2X,A)') "Self-consistent quasi particle energy calculation with contour "//&
                                "deformation starts..."
     endif

     !*** get correct XC contributions
     allocate(my_exchange_self_energy(n_high_state,n_spin))
     allocate(xc_energy(n_states,n_spin))
     call get_exchange_terms(exchange_self_energy,xc_KS_matr,x_KS_array,c_KS_array,&
                             my_exchange_self_energy,xc_energy)

     !*** determine the highest occupied orbital level
     !    such complication occurs when some of the orbitals are 
     !    either not fully occupied or not fully empty
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

     !*** set threshold, lighter for sc since
     !    CD cannot be before if residue become zero
     if((gw_cd%self_consistent.and.step_sc_loop > 1)&
         .or.use_hedin_shift) then
       qp_energy_thr = 1.5d-4
     else
       qp_energy_thr = 1.d-5
     endif

     !*** start iteration
     allocate(qp_energy_old(n_states,n_spin))
     do i_spin = 1, n_spin
        if(gw_cd%spin_channel(i_spin) /= i_spin)  cycle
        spin_tag = ""
        if(n_spin > 1 .and. i_spin == 1) spin_tag ='Spin up' 
        if(n_spin > 1 .and. i_spin == 2) spin_tag ='Spin down' 
        do i_level = 1, gw_cd%num_levels(i_spin)
           if(gw_cd%sc_env(i_spin)%qp_converged(i_level)) cycle
           i_state = gw_cd%corrected_levels(i_level,i_spin)
           qp_non_convergence(i_state,i_spin) = .false.   
          
           qp_energy(i_state, i_spin) = KS_eigenvalue_last(i_state,i_spin)
           qp_energy_old(i_state, i_spin) = KS_eigenvalue_last(i_state,i_spin)
           e_diff = 1.d-3
           i_count =0
           index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
           call read_restart_information(gw_cd, qp_energy, qp_energy_old, qp_non_convergence,&
                e_diff, i_count, i_spin, i_state, already_converged)
           if(already_converged) cycle
 
           if(myid == 0) then
              write(use_unit,*)
              write(use_unit,'(T2,A)')"------------------------------------------------------------"
              write(use_unit, '(T3,A,I6,T50, A)') "Start contour deformation for level: ", i_state,&
                                           trim(adjustl(spin_tag)) 
              write(use_unit,'(T2,A)')"------------------------------------------------------------"
           endif
 
           !*** calculate Hedin shift
           my_shift = 0.0d0
           if(use_hedin_shift) then
              my_state = i_state
              if(state_to_shift > 0) my_state = state_to_shift
              my_freq = KS_eigenvalue_last(my_state,i_spin)+0.05d0/hartree
              call calculate_delta_shift(gw_cd, Wmn_freq_cd, KS_eigenvalue, occ_numbers,&
                                         n_full_freq, n_homo, n_first, omega_full,&
                                         womega_full, chemical_potential_spin,&
                                         ovlp_3KS, my_exchange_self_energy, xc_energy,& 
                                         i_spin, my_state, my_freq, my_shift)
              gw_cd%shift(index_contour_def,i_spin) = my_shift
              if(myid == 0) write(use_unit, '(T3,A,T10,F10.4,A3)') "Shift:", my_shift*hartree, " eV" 
           endif

           !*** Start QP iteration
           do while (abs(e_diff).gt.qp_energy_thr)
              i_count = i_count +1  
              if(myid == 0) then
                 write(use_unit, '(T3,A,T10,I3)') "Step:", i_count 
              endif
              if(i_count ==1) then
                 qp_energy(i_state,i_spin) = qp_energy(i_state, i_spin)+0.002d0
              else
                 qp_energy(i_state,i_spin) = qp_energy_old(i_state,i_spin) + 0.5d0* e_diff
              endif
              qp_energy_old(i_state,i_spin) = qp_energy(i_state,i_spin)
              qp_energy(i_state,i_spin) = qp_energy(i_state,i_spin) - my_shift
              call init_contour_def_env(gw_cd,qp_energy(i_state,i_spin),KS_eigenvalue_last,&
                                        chemical_potential_spin,i_state,i_spin, n_states,&
                                        n_homo)
              call evaluate_self_energy_cd(gw_cd, Wmn_freq_cd, i_state, i_spin,&
                                           n_homo, n_first, occ_numbers, n_full_freq,&
                                           omega_full, womega_full, chemical_potential_spin, &
                                           qp_energy(i_state, i_spin), KS_eigenvalue, &
                                           KS_eigenvalue_last, ovlp_3KS)
              qp_energy(i_state,i_spin)= KS_eigenvalue(i_state,i_spin) &
                                         + gw_cd%self_energy%re(index_contour_def,i_spin) &
                                         + my_exchange_self_energy(i_state,i_spin) &
                                         - xc_energy(i_state,i_spin)
              e_diff =  qp_energy(i_state,i_spin) - qp_energy_old(i_state,i_spin)
              deallocate(gw_cd%residue_from_freq,gw_cd%real_freq)

              call print_restart_information(gw_cd, qp_energy, qp_non_convergence, e_diff,&
                                             i_count, i_state, i_spin)
              if(i_count .gt. 200) then
                if(myid == 0) then
                  write(use_unit,'(2X,2A,I4,A,I4 )') &
                  "QUASI_PARTILCE_ENERGY: self-consistent quasiparticle ", &
                  "solution with contour deformation can not  be found for state: i_state = ", &
                  i_state, "  i_spin = ", i_spin
                endif
                qp_non_convergence(i_state,i_spin)=.true.
                exit
              endif
           enddo 
           correl_energy(i_state,i_spin) = gw_cd%self_energy%re(index_contour_def,i_spin) 
        enddo ! i_state
     enddo    ! i_spin

     deallocate(my_exchange_self_energy)
     deallocate(xc_energy)
     deallocate(qp_energy_old)

     call get_times(time_qp_energy_cd, clock_time_qp_energy_cd)

  end subroutine quasi_particle_energy_contour_def

! **************************************************************************************************
!> brief calculate Hedin shift \Delta E to achieve level alignment. Can be also considered
!        as poor-man's self-consistency similar to ev-scGW0
!  o gw_cd  --  contour deformation environment
!  o Wmn_freq_cd --  screened Coulomb matrix element in the MO basis, i.e. (mn|W(iomega)|mn)
!  o KS_eigenvalue -- real array,
!           the eigenvalues of the single-particle calculation. For DFT calculation,
!           this is the KS eigenvalue, but for HF calculation, this is then the HF
!           eigenvalue
!  o occ_numbers -- occupation numbers of single-particle energy levels
!  o n_full_freq -- the number of frequency points for the screened Coulomb interaction W
!  o n_homo -- the HOMO level for each spin channel
!  o n_first -- the first orbital which is NOT fully occupied 
!  o omega_full -- the Gauss-Legendre frequency grid for the screened Coulomb interaction
!  o womega_full -- the weight of the Gauss-Legendre frequency grid for the self-energy
!  o chemical_potential -- the chemical potential of the system
!  o ovlp_3KS -- transformed 3-center integrals
!  o my_exchange_self_energy -- Sigma_x 
!  o xc_energy -- DFT-XC energy, v_xc
!  o state_to_shift -- state for which the shift is calculated
!  o freq_real -- frequency at which self-energy is calculated, i.e., current epsilon_QP in iterati.
!  o delta_shift -- calculated Hedin shift 
! **************************************************************************************************
  subroutine calculate_delta_shift(gw_cd, Wmn_freq_cd, KS_eigenvalue, occ_numbers, n_full_freq,&
                                   n_homo, n_first, omega_full, womega_full, chemical_potential_spin,&
                                   ovlp_3KS, my_exchange_self_energy, xc_energy, i_spin, &
                                   state_to_shift, freq_real, delta_shift)

     type(cd_environment_type)                                :: gw_cd
     complex(kind=8), dimension(:,:,:,:), intent(in)          :: Wmn_freq_cd
     real(kind=8), dimension(:,:), intent(in)                 :: KS_eigenvalue
     real(kind=8), dimension(:,:)                             :: occ_numbers
     integer, intent(in)                                      :: n_full_freq
     integer, dimension(:), intent(in)                        :: n_homo
     integer, dimension(:), intent(in)                        :: n_first
     real(kind=8), dimension(:), intent(in)                   :: omega_full
     real(kind=8), dimension(:), intent(in)                   :: womega_full
     real(kind=8), dimension(:), intent(in)                   :: chemical_potential_spin   
     real(kind=8), dimension(:,:,:,:),intent(in)              :: ovlp_3KS
     real(kind=8), dimension(:,:), intent(in)                 :: my_exchange_self_energy
     real(kind=8), dimension(:,:), intent(in)                 :: xc_energy
     integer, intent(in)                                      :: i_spin
     integer, intent(in)                                      :: state_to_shift
     real(kind=8), intent(in)                                 :: freq_real
     real(kind=8), intent(out)                                :: delta_shift

     character(*), parameter :: func = 'calculate_delta_shift'

     integer                                                  :: index_contour_def
     real(kind=8)                                             :: KS_value, dist_KS_val
 
     dist_KS_val = 0.005d0/hartree
     KS_value = KS_eigenvalue(state_to_shift,i_spin)
     index_contour_def = gw_cd%index_cd_levels(state_to_shift,i_spin) 

     if(ABS(freq_real - KS_value) < dist_KS_val) then
       if(myid.eq.0) write(use_unit,"(T3,A)") "freq close to KS"
     endif

     call init_contour_def_env(gw_cd,freq_real,KS_eigenvalue,chemical_potential_spin,&
                               state_to_shift, i_spin, n_states, n_homo)
     call evaluate_self_energy_cd(gw_cd, Wmn_freq_cd, state_to_shift,i_spin,&
                                  n_homo, n_first, occ_numbers, n_full_freq,&
                                  omega_full, womega_full, chemical_potential_spin, &
                                  freq_real, KS_eigenvalue, KS_eigenvalue, ovlp_3KS)
     delta_shift = gw_cd%self_energy%re(index_contour_def,i_spin) &
                    + my_exchange_self_energy(state_to_shift,i_spin) &
                    - xc_energy(state_to_shift,i_spin)

     deallocate(gw_cd%residue_from_freq,gw_cd%real_freq)
 end subroutine calculate_delta_shift

! **************************************************************************************************
!> brief evalutes the quasi particle energies based on the tailor expansion of the quasiparticle
! **************************************************************************************************
  subroutine quasi_particle_energy_contour_def_zshot(gw_cd, Wmn_freq_cd, qp_energy, correl_energy, &
                                               KS_eigenvalue, KS_eigenvalue_last, exchange_self_energy,&
                                               xc_KS_matr,  x_KS_array, c_KS_array, occ_numbers,&
                                               n_low_state, n_high_state, n_homo, n_full_freq, &
                                               omega_full, womega_full, chemical_potential_spin, &
                                               ovlp_3KS,qp_non_convergence)

     type(cd_environment_type)                                :: gw_cd
     complex(kind=8), dimension(:,:,:,:), intent(in)          :: Wmn_freq_cd
     real(kind=8), dimension(:,:), intent(inout)              :: qp_energy
     real(kind=8), dimension(:,:), intent(inout)              :: correl_energy
     real(kind=8), dimension(:,:), intent(in)                 :: KS_eigenvalue
     real(kind=8), dimension(:,:), intent(in)                 :: KS_eigenvalue_last
     real(kind=8), dimension(:,:), intent(in)                 :: exchange_self_energy
     real(kind=8), dimension(:,:,:), intent(in)               :: xc_KS_matr
     real(kind=8), dimension(:,:),  intent(in)                :: x_KS_array      
     real(kind=8), dimension(:,:), intent(in)                 :: c_KS_array      
     real(kind=8), dimension(:,:)                             :: occ_numbers
     integer, intent(in)                                      :: n_low_state,n_high_state,&
                                                                 n_full_freq
     integer, dimension(:), intent(in)                        :: n_homo
     real(kind=8), dimension(:), intent(in)                   :: omega_full
     real(kind=8), dimension(:), intent(in)                   :: womega_full
     real(kind=8), dimension(:), intent(in)                   :: chemical_potential_spin   
     real(kind=8), dimension(:,:,:,:), intent(in)             :: ovlp_3KS
     logical, dimension(:,:), intent(in), optional            :: qp_non_convergence


     character(*), parameter :: func = 'quasi_particle_energy_contour_def_zshot'

     integer                                                  :: i_spin, i_state, &
                                                                 index_contour_def
     integer                                                  :: i_level, i_freq
     integer, dimension(n_spin)                               :: n_first
     real(kind=8), dimension(:,:), allocatable                :: my_exchange_self_energy
     real(kind=8), dimension(:,:), allocatable                :: xc_energy
     real(kind=8)                                             :: freq_real
     real(kind=8), dimension(3)                               :: self_energy
   
     if(myid.eq.0) then
       write(use_unit,*)
       write(use_unit,'(2X,A)') "Z-shot quasi particle energy calculation with contour "//&
                                "deformation starts..."
     endif

     !*** get correct XC contributions
     allocate(my_exchange_self_energy(n_high_state,n_spin))
     allocate(xc_energy(n_states,n_spin))
     call get_exchange_terms(exchange_self_energy,xc_KS_matr,x_KS_array,c_KS_array,&
                             my_exchange_self_energy,xc_energy)
 

     !*** determine the highest occupied orbital level
     !    such complication occurs when some of the orbitals are 
     !    either not fully occupied or not fully empty
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

     if(.not.(allocated(z_value))) then
       allocate(z_value(n_states, n_spin))
       z_value = 0.0d0
     endif
     do i_spin = 1, n_spin
        if(gw_cd%spin_channel(i_spin) /= i_spin)  cycle
        do i_level = 1, gw_cd%num_levels(i_spin)
           i_state = gw_cd%corrected_levels(i_level,i_spin)
           index_contour_def = gw_cd%index_cd_levels(i_state,i_spin)
           if(present(qp_non_convergence)) then
             if(.not.qp_non_convergence(i_state,i_spin)) cycle
           endif 
           do i_freq = 1, 2
              if(i_freq == 1) then
                freq_real = KS_eigenvalue_last(i_state,i_spin) + gw_cd%contour_def_offset
              elseif(i_freq == 2) then
                freq_real = KS_eigenvalue_last(i_state,i_spin) - gw_cd%contour_def_offset 
              endif 
              call init_contour_def_env(gw_cd,freq_real,KS_eigenvalue_last,&
                                        chemical_potential_spin,i_state,i_spin,n_states,&
                                        n_homo)
              call evaluate_self_energy_cd(gw_cd, Wmn_freq_cd, i_state, i_spin,&
                                           n_homo, n_first, occ_numbers, n_full_freq,&
                                           omega_full, womega_full, chemical_potential_spin, &
                                           freq_real, KS_eigenvalue, KS_eigenvalue_last, ovlp_3KS) 
              self_energy(i_freq) = gw_cd%self_energy%re(index_contour_def,i_spin) 
              deallocate(gw_cd%residue_from_freq,gw_cd%real_freq)
           enddo
           call compute_z_and_m_contour_def(self_energy, z_value(i_state,i_spin), &
                                            gw_cd%contour_def_offset )        
           qp_energy(i_state,i_spin)= KS_eigenvalue(i_state,i_spin) &
                      + z_value(i_state,i_spin)*(self_energy(3)&
                      + my_exchange_self_energy(i_state,i_spin) &
                      - xc_energy(i_state,i_spin))

           correl_energy(i_state,i_spin) = z_value(i_state,i_spin)*self_energy(3) 
        enddo ! i_state
     enddo    ! i_spin

     deallocate(my_exchange_self_energy)
     deallocate(xc_energy)
  end subroutine quasi_particle_energy_contour_def_zshot

! **************************************************************************************************
!> brief compute z_value for Z-Shot contour deformation
!  o self_energy_real -- self-energy on real-frequency axis
!  o z_value -- renormalization factor 
!  o m_value -- derivative
!  o self_energy_delta -- self-energies for \epsilon_n +- delta
!  o contour_def_offset  -- delta value
! **************************************************************************************************
  subroutine compute_z_and_m_contour_def(self_energy, z_value, contour_def_offset)
     real(kind=8), dimension(3), intent(inout)          :: self_energy
     real(kind=8), intent(out)                          :: z_value
     real(kind=8), intent(in)                           :: contour_def_offset

     real(kind=8)                                       :: derivative

     ! average: sigma_c(en+delta) and sigma(en-delta)
     self_energy(3) = 0.5d0*(self_energy(1)+self_energy(2))
     derivative = 0.5d0*(self_energy(1)-self_energy(2))/contour_def_offset

     z_value = 1.0d0/(1.0d0-derivative)
 
  end subroutine compute_z_and_m_contour_def

! **************************************************************************************************
!> brief get the right exchange contributions 
!  o exchange_self_energy --  exact exchange part of the self-energy 
!  o xc_KS_matr -- exchange correlation potential witin KS/HF orbitals
!  o x_KS_array -- exchange potential witin KS/HF orbitals
!  o c_KS_array -- correlation potential witin KS/HF orbitals
!  o my_exchange_energy -- exchange self-energy added to correlation term in QP equation
!  o xc_energy -- v_xc potential (KS) subtracted in QP equation 
! **************************************************************************************************
  subroutine get_exchange_terms(exchange_self_energy,xc_KS_matr,x_KS_array,c_KS_array,&
                                my_exchange_self_energy,xc_energy)
     real(kind=8), dimension(:,:), intent(in)                 :: exchange_self_energy
     real(kind=8), dimension(:,:,:), intent(in)               :: xc_KS_matr
     real(kind=8), dimension(:,:), intent(in)                 :: x_KS_array      
     real(kind=8), dimension(:,:), intent(in)                 :: c_KS_array  
     real(kind=8), dimension(:,:), intent(inout)              :: my_exchange_self_energy
     real(kind=8), dimension(:,:), intent(inout)              :: xc_energy

     integer                                                  :: i_spin, i_state

     if(use_hartree_fock .and. (.not.use_screx) &
           .and. .not. use_gw_and_hse) then
       my_exchange_self_energy(:,:) = exchange_self_energy * & 
                                (1.0d0-hybrid_coeff)
     else
       my_exchange_self_energy(:,:) = exchange_self_energy
     endif

     !*** read in DFT exchange-correlation energy for all states
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

  end subroutine get_exchange_terms

! **************************************************************************************************
!> brief output quasi-particle energies
!  o qp_energy --  quasi particle energies
!  o KS_eigenvalue -- KS/HF eigenvalues of the single-particle calculation
!  o exchange_self_energy -- real array the exact exchange part of the self-energy 
!           of each KS/HF state for each spin channel
!  o x_KS_array -- real array, the matrix elements of the exchange
!           potential witin KS/HF orbitals
!  o c_KS_array -- real array, the matrix elements of the correlation
!           potential witin KS/HF orbitals
!  o correl_energy -- real array, correlation part of self-energy
!  o xc_energy -- exchange correlation potential witin KS/HF orbitals
!  o occ_numbers -- occupation numbers of single-particle energy levels
!  o n_low_state  -- the lowest KS/HF eigenstate for self-energy correction
!  o n_high_state  -- the highest KS/HF eigenstate for self-energy correction
!  o n_spin -- number of spins
!  o n_states -- number of KS/HF eigenstates
!  o n_homo -- the HOMO level for each spin channel
!  o contour_def_start -- level where contour def starts
!  o contour_def_end -- level where contour def ends
!  o z_value -- 1/1(-derivative)
! **************************************************************************************************
  subroutine output_quasi_particle_energies(qp_energy, KS_eigenvalue, exchange_self_energy,&
                                            x_KS_array, c_KS_array, correl_energy, xc_energy,&
                                            occ_numbers, n_low_state,n_high_state, n_spin, n_states,&
                                            n_homo, contour_def_start, contour_def_end, &
                                            contour_spin_channel,z_value,qp_non_convergence)

     real(kind=8), dimension(n_low_state:n_high_state,n_spin),&
       intent(in)                                             :: qp_energy
     real(kind=8), dimension(n_states,n_spin),&
       intent(in)                                             :: KS_eigenvalue
     real(kind=8), dimension(n_high_state,n_spin),&
       intent(in)                                             :: exchange_self_energy
     real(kind=8), dimension(n_states,n_spin), & 
       intent(in)                                             :: x_KS_array      
     real(kind=8), dimension(n_states,n_spin), & 
       intent(in)                                             :: c_KS_array      
     real(kind=8), dimension(n_states,n_spin), intent(in)     :: correl_energy
     real(kind=8), dimension(n_states,n_spin), &
       intent(in)                                             :: xc_energy
     real(kind=8), dimension(n_states,n_spin)                 :: occ_numbers
     integer, intent(in)                                      :: n_low_state,n_high_state,&
                                                                 n_spin, n_states
     integer, dimension(n_spin), intent(in)                   :: n_homo
     integer, dimension(n_spin), intent(in), optional         :: contour_def_start,&
                                                                 contour_def_end,&
                                                                 contour_spin_channel
     real(kind=8), dimension(:,:), intent(in), optional       :: z_value
     logical, dimension(n_low_state:n_high_state,n_spin), &
       intent(in), optional                                   :: qp_non_convergence


     character(len=20)                                        :: header_end
     integer                                                  :: i_spin, i_state     
     character*120                                            :: info_str
     logical, dimension(:,:), allocatable                     :: my_qp_non_convergence


     allocate(my_qp_non_convergence(n_states, n_spin))
     my_qp_non_convergence(:,:) = .false.
     if(present(qp_non_convergence)) then
        my_qp_non_convergence(n_low_state:n_high_state,:) &
         = qp_non_convergence(n_low_state:n_high_state,:) 
     endif

     if(present(z_value)) then
       header_end = "e_qp      z-factor"
     else
       header_end = "e_qp"
     endif
      
     if(myid.eq.0) then
       write(use_unit,*)
       write(use_unit,'(2A)')"-----------------------------------------", &
        "-------------------------------------------------------------"
       if (use_gw) then
          write(use_unit,'(15X,A)')"GW quasi-particle energy levels"
       else if (use_mp2sf) then
          write(use_unit,'(15X,A)')"MP2 quasi-particle energy levels"
       else
          write(use_unit,'(15X,A)')"Quasi-particle energy levels"
       end if
       
       if(present(contour_def_start) .and. present(contour_def_end).and.&
          present(contour_spin_channel)) then
          if(contour_spin_channel(1) == 1) then
            write(use_unit,*)
            write(use_unit,'(T16,A35,T52,I6)') &
             'Start contour deformation at level:', contour_def_start(1)
            write(use_unit,'(T16,A33,T52,I6)') &
             'End contour deformation at level:', contour_def_end(1)
          endif
       endif
       write(use_unit,*)
       if(use_split_xc_gw) then

         write(use_unit,'(15X,A)')"e_qp = e_gs + e_x^ex - e_x^gs - e_c^gs + e_c^nloc"
         write(use_unit,*)
        write(use_unit, '(2X, A, 5X, A,8X, A, 8X,A,8X,A,8X,A,8X,A,8X, A)') &
            "state", "occ_num", "e_gs", "e_x^ex", "e_x^gs", "e_c^gs", &
             "e_c^nloc", TRIM(ADJUSTL(header_end))
       else
         write(use_unit,'(15X,A)')"e_qp = e_gs + e_x^ex - e_xc^gs + e_c^nloc"
         write(use_unit,*)
         write(use_unit, '(2X, A, 5X, A,8X, A, 8X,A,8X,A,8X,A, 8X, A)') &
            "state", "occ_num", "e_gs", "e_x^ex", "e_xc^gs", &
             "e_c^nloc", TRIM(ADJUSTL(header_end))
       endif
       write(use_unit,'(2A)')"-----------------------------------------", &
        "-------------------------------------------------------------"
       do i_spin = 1, n_spin
         if(n_spin.eq.2.and.i_spin.eq.1) then
           write(use_unit,'(35X, A)') "Spin Up"
         endif

         if(n_spin.eq.2.and.i_spin.eq.2) then
            if(present(contour_def_start) .and. present(contour_def_end).and.&
               present(contour_spin_channel)) then
               if(contour_spin_channel(i_spin) == i_spin) then
                  write(use_unit,'(T16,A35,T52,I6)') &
                   'Start contour deformation at level:', contour_def_start(2)
                  write(use_unit,'(T16,A33,T52,I6)') &
                   'End contour deformation at level:', contour_def_end(2)
                  write(use_unit,'(2A)')"-----------------------------------------", &
                   "-------------------------------------------------------------"
               endif
            endif
           write(use_unit,'(35X, A)') "Spin Down"
         endif

         write(use_unit,'(2A)')"-----------------------------------------", &
         "-------------------------------------------------------------"
         do i_state = n_low_state, n_high_state, 1
              if(my_qp_non_convergence(i_state,i_spin).eqv..true.) then 
                if(use_split_xc_gw) then
                   write (use_unit, '(2X, I6, 2X,F8.4, 4F14.4, A)') &
                         i_state, occ_numbers(i_state,i_spin), &
                         KS_eigenvalue(i_state,i_spin)*hartree, &
                         exchange_self_energy(i_state,i_spin)*hartree, &
                         x_KS_array(i_state,i_spin)*hartree, &
                         c_KS_array(i_state,i_spin)*hartree, &
                                " *** no convergence (see below) ***"
                else
                  write (use_unit, '(2X, I6, 2X, F8.4, 3F14.4, A)') &
                  i_state, occ_numbers(i_state,i_spin), &
                  KS_eigenvalue(i_state,i_spin)*hartree, &
                  exchange_self_energy(i_state,i_spin)*hartree, &
                  xc_energy(i_state,i_spin)*hartree, &
                              " *** no convergence (see below) ***"
                endif !end of if use_split_xc_gw
              else
                if(use_split_xc_gw) then
                  if(present(z_value)) then
                    write (use_unit, '(2X, I6, 2X,F8.4, 7F14.4)') &
                    i_state, occ_numbers(i_state,i_spin), &
                    KS_eigenvalue(i_state,i_spin)*hartree, &
                    exchange_self_energy(i_state,i_spin)*hartree, &
                    x_KS_array(i_state,i_spin)*hartree, &
                    c_KS_array(i_state,i_spin)*hartree, &
                    correl_energy(i_state,i_spin)*hartree, &
                    qp_energy(i_state,i_spin)*hartree,&
                    z_value(i_state,i_spin)
                 else
                   write (use_unit, '(2X, I6, 2X,F8.4, 6F14.4)') &
                    i_state, occ_numbers(i_state,i_spin), &
                    KS_eigenvalue(i_state,i_spin)*hartree, &
                    exchange_self_energy(i_state,i_spin)*hartree, &
                    x_KS_array(i_state,i_spin)*hartree, &
                    c_KS_array(i_state,i_spin)*hartree, &
                    correl_energy(i_state,i_spin)*hartree, &
                    qp_energy(i_state,i_spin)*hartree
                  endif
                else
                  if(present(z_value)) then
                    write (use_unit, '(2X, I6, 2X, F8.4, 6F14.4)') &
                    i_state, occ_numbers(i_state,i_spin), &
                    KS_eigenvalue(i_state,i_spin)*hartree, &
                    exchange_self_energy(i_state,i_spin)*hartree, &
                    xc_energy(i_state,i_spin)*hartree, &
                    correl_energy(i_state,i_spin)*hartree, &
                    qp_energy(i_state,i_spin)*hartree,&
                    z_value(i_state,i_spin)
                  else
                    write (use_unit, '(2X, I6, 2X, F8.4, 5F14.4)') &
                    i_state, occ_numbers(i_state,i_spin), &
                    KS_eigenvalue(i_state,i_spin)*hartree, &
                    exchange_self_energy(i_state,i_spin)*hartree, &
                    xc_energy(i_state,i_spin)*hartree, &
                    correl_energy(i_state,i_spin)*hartree, &
                    qp_energy(i_state,i_spin)*hartree
                  endif
                endif !end of if use_split_xc_gw
!  end of if(qp_non_convergence)
              endif
         enddo


         if(n_spin.eq.2.and.i_spin.eq.1) then
          write(use_unit,'(2A)')"------------------------------------------", &
           "------------------------------------------------------------"
         endif
       enddo ! end of i_spin

       write(use_unit,'(2A)')"------------------------------------------", &
          "------------------------------------------------------------"
       write(use_unit,*)

       if(n_homo(1).gt.0 .and. n_homo(n_spin).gt.0)then ! for H atom
           write(use_unit,'(2X, A, 2f10.4)') &
           " DFT/Hartree-Fock HOMO level (eV): ", &
          KS_eigenvalue(n_homo(1),1)*hartree, &
          KS_eigenvalue(n_homo(n_spin),n_spin)*hartree
          
         write(use_unit,*)
         write(use_unit,'(2X, A, 2f10.4)') &
         " Quasiparticle HOMO level (eV):    ", &
           qp_energy(n_homo(1),1)*hartree, &
           qp_energy(n_homo(n_spin),n_spin)*hartree
       endif

       if (warn_qpe_nonconvergence.and.&
           .not.(use_contour_def_gw)) then
           write (info_str,'(1X,A)') &
           ' '
           call localorb_info(info_str,use_unit,'(A)')
           write (info_str,'(1X,1A)') &
           "* Warning! One or more quasiparticle eigenvalues listed above could not be"
           call localorb_info(info_str,use_unit,'(A)')
           write (info_str,'(1X,1A)') &
           "* determined unambiguously in the list given above. This is related to the"
           call localorb_info(info_str,use_unit,'(A)')
           write (info_str,'(1X,1A)') &
           "* default numerical fit procedure used to determine parameters for analytical"
           call localorb_info(info_str,use_unit,'(A)')
           write (info_str,'(1X,1A)') &
           "* continuation of the self-energy from the imaginary axis to the real axis."
           call localorb_info(info_str,use_unit,'(A)')
           write (info_str,'(1X,1A)') &
           "* More information can be found in Sec. 3.21 in the manual."
           call localorb_info(info_str,use_unit,'(A)')
           write (info_str,'(1X,1A)') &
           "* In particular, the keyword 'anacon_type 0' can be used to change the default"
           call localorb_info(info_str,use_unit,'(A)')
           write (info_str,'(1X,1A)') &
           "* analytical continuation type to the slightly less accurate 'two-pole' self-"
           call localorb_info(info_str,use_unit,'(A)')
           write (info_str,'(1X,1A)') &
           "* energy model, which is more robust than the 16-parameter Pade approximation"
           call localorb_info(info_str,use_unit,'(A)')
           write (info_str,'(1X,1A)') &
           "* used by default. Again, please see the manual for the mathematical definitions."
           call localorb_info(info_str,use_unit,'(A)')
       end if
       write(use_unit,*)
       write(use_unit,'(2A)')"------------------------------------------", &
         "------------------------------------------------------------"
       write(use_unit,*)


!  end of if myid==0
     endif

     deallocate(my_qp_non_convergence) 
  end subroutine output_quasi_particle_energies
  
end module quasi_particle_energies
