       subroutine get_qp_spectrum            &
           (anacon_type,n_max_par, &
            n_low_state,n_high_state, &
            n_homo, n_lumo, n_freq, omega, &
            sigma_par,KS_eigenvalue, &
            exchange_self_energy, &
            xc_KS_matr,filename )

! (xc_KS_matr, KS_eigenvalue,&
!                        anacon_type, n_low_state,n_high_state ,   &
!                        sigma_par,omega_grid , n_max_par, n_homo, n_lumo )

      use constants, only: pi
      use dimensions
      use runtime_choices
      use mpi_tasks

      implicit none

!  ARGUMENTS

      integer :: anacon_type
      integer :: n_max_par
      integer :: n_low_state, n_high_state
      integer :: n_freq
      integer :: n_homo(n_spin)
      integer :: n_lumo(n_spin)

      real*8 :: omega(n_freq)
      real*8 :: KS_eigenvalue(n_states,n_spin)
      real*8 :: exchange_self_energy(n_high_state,n_spin)
      real*8 :: xc_KS_matr(n_states,n_states,n_spin)
      character*14 filename
      complex*16 :: sigma_par(n_max_par,n_states,n_spin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable :: xc_energy(:,:)
      real*8, allocatable :: correl_energy(:,:)

      real*8  e_diff
      real*8  en
      real*8  mu
      real*8  delta_mu
 
      real ims, res

      complex*16 selfe, dselfe

      integer :: i_state,i_state_1
      integer :: i, i_count, i_freq
      integer :: i_spin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
!      real*8, dimension (:), allocatable ::  spectral_function
      real*8 ::  spectral_function
      real*8 ::  spectral_function_homo
      real*8 interval
      real*8 freq

!      if (n_spin.ne.1)then
!        write(use_unit,*) " Quasi particle spectrum ",&
!                      " implemented only for spin unpolarized"
!      endif

      do i_spin = 1, n_spin
        if(myid.eq.0 .and. i_spin .eq.1)then
          if (n_spin.gt.1)then
             open(86,file='up'//filename)
             open(87,file='do'//filename)
          else 
             open(86,file=filename)
!             open(87,file='ho'//filename)
          endif
        endif       
      
        interval = 6.5 /10000
        freq = 0.d0
      
        do i_freq = -10000,10000, 1
         freq = i_freq*interval
      
         spectral_function  = 0.d0
         spectral_function_homo  = 0.d0
      
         do i_state = n_low_state, n_high_state, 1
      
           call get_real_selfenergy ( anacon_type, n_freq , omega, &
                 dcmplx(freq, 0.d0), n_max_par, &
                 sigma_par(1:n_max_par,i_state,i_spin),selfe)
     
           res = real(selfe) 
           ims = aimag(selfe)
!           ims = 0.000001 
!           res = aimag(selfe) 
!           ims = real(selfe) 


           if(.not. use_hartree_fock)then
              spectral_function  = &
              spectral_function    &
              + 1.d0/pi * abs(ims) &
              /( (freq - KS_eigenvalue(i_state, i_spin)-&
                (res+exchange_self_energy (i_state,i_spin) &
                -xc_KS_matr(i_state,i_state,i_spin)))**2 &
                + ims*ims )
           elseif (use_hartree_fock)then
              spectral_function  = &
              spectral_function    &
              + 1.d0/pi * abs(ims) &
              /( (freq - KS_eigenvalue(i_state, i_spin)-&
                res)**2 &
                + ims*ims )
           endif
      
!           if(i_state == n_homo(i_spin))then
!             spectral_function_homo  = &
!             spectral_function_homo    &
!             + 1.d0/pi * abs(aimag(selfe)) &
!             /( (freq - KS_eigenvalue(i_state, i_spin)-&
!             (real(selfe)+exchange_self_energy (i_state,i_spin) &
!             -xc_KS_matr(i_state,i_state,i_spin)))**2 &
!             + aimag(selfe)*aimag(selfe) )
!           endif
      
         enddo
         if(myid.eq.0)then
           if (i_spin.eq.1)then
             write(86,*) freq*hartree, spectral_function!, spectral_function_homo
           elseif(i_spin.eq.2)then
             write(87,*) freq*hartree, spectral_function
           endif 
!           if (n_spin == 1)then 
!             write(87,*) freq*hartree, spectral_function_homo
!           endif 
         endif
        enddo

        if(myid.eq.0)then
          if (n_spin.ne.1)then
             close(86)
             close(87)
          else
             close(86)
          endif
        endif
 
      enddo!i_spin

      end subroutine get_qp_spectrum


