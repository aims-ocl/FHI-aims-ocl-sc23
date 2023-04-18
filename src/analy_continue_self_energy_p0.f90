!****s* FHI-aims/analy_continute_self_energy_p0
!  NAME
!   analy_continute_self_energy_p0
!  SYNOPSIS
!  
      subroutine analy_continue_self_energy_p0 &
          (anacon_type, n_recip_points, n_recip_points_task,&
           n_max_par, n_low_state,n_high_state,n_freq, &
           omega, self_energy_omega, &
           sigma_par)

!  PURPOSE
!  Subroutine analy_continue_self_energy_p0  analytically continue 
!  the self energy from imaginary energy axis to real energy axis.
!  Two schemes can be used here: the two-pole fitting or Pade approximation 
!
! USES

      use dimensions
      use mpi_tasks
      use localorb_io,only :use_unit
      use gw_para, only: anacon_type_defined
      
      implicit none

! ARGUMENTS

      integer ::  anacon_type
      integer ::  n_recip_points
      integer ::  n_recip_points_task
      integer ::  n_low_state
      integer ::  n_high_state
      integer ::  n_freq
      integer ::  n_max_par 

      real*8 ::   omega(n_freq)
      complex*16 ::   self_energy_omega  &
               (n_freq,n_low_state:n_high_state,n_spin,n_recip_points_task)

      complex*16 :: sigma_par(n_max_par,n_low_state:n_high_state,n_spin,n_recip_points_task)

! INPUTS 
!         self-energy comes from the calling subroutine
!  o  anacon_type -- integer number, if 0, the two-pole fitting for analytical
!          continuation; if 1, using Pade approximation for ana. cont.
!  o  n_low_state  -- integer number,
!          the lowest KS/HF eigenstate for self-energy correction
!  o  n_high_state -- integer number,
!          the highest KS/HF eigenstate for self-energy correction
! o  n_recip_points -- number of k points (not necessarily on the regular k mesh)
! o  n_recip_points_task -- number of k points per CPU core 
!  o  n_k_points_task -- number of k_points per processor
!  o  n_freq -- integer number, the number of frequency points for the GW self-energy
!  o  n_max_par -- the number of parameters used for analytical continuation
!          For anacon_type = 0, recommended n_max_par is  4 or 5. If 4, this will
!          be the normal two-pole fitting, else if 5, it will be two-pole plus a
!          (small) constant number
!          For anacon_type = 1, recommended n_max_par is the half of n_freq
!  o  omega(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the self-energy
!  o  self_energy_omega  -- complex array, the calculated self-energy on imaginary
!            frequency grid for each state and each spin channel
! OUTPUT
!  o  sigma_par -- complex array, the fitting parameters from analytical continuation

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
      
!  fitting omega grid range 
      
      complex*16, allocatable :: par(:)

!     counters

      integer :: i_state
      integer :: i_freq
      integer :: i_spin
      integer :: i_k_point
      integer :: i_k_point_local


!     begin work

      

!   write out 
    
      ! Error trap in case the analytical continuation type was not defined in control.in.
      if (.not.anacon_type_defined) then
         if (myid.eq.0) then
            write(use_unit,'(2X,A)') "Analytical continuation of the self-energy starts..."
            write(use_unit,'(1X,A)') "* Oops - error. The analytical continuation of the self-energy"
            write(use_unit,'(1X,A)') "* from the imaginary to the real frequency axis was invoked"
            write(use_unit,'(1X,A)') "* although the keyword 'anacon_type' was not defined in control.in."
            write(use_unit,'(1X,A)') "* This should not happen (this omission should have been caught when"
            write(use_unit,'(1X,A)') "* first reading control.in in read_control.f90). Please notify us"
            write(use_unit,'(1X,A)') "* if you do encounter this error message - normally, it should"
            write(use_unit,'(1X,A)') "* never trigger."
            call aims_stop('anacon_type not defined in control.in.','analy_continue_self_energy_p0.f90')
         end if
      end if
      
      if(myid.eq.0) then
        write(use_unit,*) 
        write(use_unit,'(2X,A)') "Analytical continuation of the self-energy starts..."
      endif 

      allocate (par(n_max_par),stat=i_spin)
      call check_allocation(i_spin, 'par                           ')


! perform analytical continuation, and get the fitting parameters 
      do i_k_point = 1, n_recip_points, 1

        i_k_point_local = (i_k_point-1)/n_tasks + 1
        if(myid.ne.mod(i_k_point,n_tasks)) cycle
        do i_spin = 1, n_spin, 1
!         write(use_unit,*) " | i_spin   ", i_spin
           do i_state = n_low_state, n_high_state, 1
!               write(use_unit,*) " |      i_state   ", i_state
 

              if(anacon_type .eq. 0) then
                call mr_min_lsq(n_freq,dcmplx(0.d0,omega), & 
                          self_energy_omega(:,i_state,i_spin,i_k_point_local), &
                          n_max_par,par) 
              elseif(anacon_type .eq. 1) then
                call get_pade_approx_para(n_freq,omega, &
                         self_energy_omega(:,i_state,i_spin,i_k_point_local), &
                         n_max_par,par)
              endif

              sigma_par(1:n_max_par,i_state,i_spin,i_k_point_local)=par(1:n_max_par)


!            do i_freq = n_freq, 1, -1 
!
!              call get_real_selfenergy(anacon_type,n_freq,omega, &
!                       dcmplx(-omega(i_freq),0.d0), n_max_par, &
!                       sigma_par(1:n_max_par,i_state,i_spin),selfe)
! 
!              self_energy_real_omega(-i_freq+1,i_state,i_spin) = selfe
!            enddo
!
!            do i_freq = 1, n_freq, 1 
!
!             call get_real_selfenergy(anacon_type,n_freq,omega,&
!                       dcmplx(omega(i_freq),0.d0),n_max_par,&
!                       sigma_par(1:n_max_par,i_state,i_spin),selfe)
! 
!              self_energy_real_omega(i_freq,i_state,i_spin) = selfe
!            enddo

          enddo
        enddo
      enddo

      deallocate (par)

      return 
      end subroutine analy_continue_self_energy_p0
!******
!---------------------------------------------------------------------
