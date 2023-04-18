        subroutine get_gw_tot_en (eex_energy, self_energy_omega)

! calculate the total energy from Galitskij-Migdal (GM)
! using the KS Green's function G_KS.
! Note that Exchange en, Hartree en, Kinetic en and Electron-Ion interaction en
! evaluated with GM using G_KS assume the same values as in RPA!
! What changes here is only the evaluation of the correlation en, 
! that is derived directly from the Self-Energy instead that from the
! ACFD theorem as in RPA.   

        use mpi_tasks
        use gw_para
        use physics
        use dimensions
        use localorb_io,only:use_unit
        use constants, only: hartree, pi

        implicit none

        real*8 eex_energy
        complex*16 self_energy_omega (n_states, n_freq, n_spin)

        complex*16, dimension(:,:,:), allocatable :: green_fn

        real*8 GW_corr_fit
        real*8 tot_en_GW_fit
        real*8 GW_corr
        real*8 tot_en_GW

        integer i_spin
        integer i_state
        integer i_freq
        integer i

        integer n_par, i_par
        complex*16,dimension(:), allocatable :: self_en_par 
        complex*16 a,b,w
        complex*16 denominator, numerator
        complex*16, dimension (:) , allocatable :: fitted_sigma 
 
        real*8, dimension (:) , allocatable :: aux_omega
        real*8, dimension (:) , allocatable :: aux_womega
        integer aux_nomega
        real*8  omega0
        real*8 correction
       aux_nomega = 2000
       if(.not. allocated(aux_omega)) then
          allocate(aux_omega(aux_nomega))
       endif
       if(.not. allocated(aux_womega)) then
          allocate(aux_womega(aux_nomega))
       endif
       if(.not. allocated(fitted_sigma)) then
          allocate(fitted_sigma(aux_nomega))
       endif

       n_par = 4

       !auxiliary grid for the correlation term
       omega0 = 0.5d0
       call gauleg(-1.d0,1.d0, aux_omega(1:aux_nomega),aux_womega(1:aux_nomega),aux_nomega)
       do i = 1, aux_nomega
        aux_womega(i) = aux_womega(i) * 2.0d0*omega0/( (1.0d0-aux_omega(i))**2 )
        aux_omega(i) = omega0*(1.0d0+aux_omega(i))/(1.0d0-aux_omega(i))
       enddo
       

       if(.not. allocated(self_en_par)) then
         allocate(self_en_par (n_par))
       endif
       if(.not. allocated(green_fn)) then
         allocate(green_fn (n_states,n_freq, n_spin))
       endif
       !update the Green's function from dyson's equation
       GW_corr = 0.d0
       GW_corr_fit = 0.d0
       correction = 0.d0
       if(myid.eq.0)then 
         open (64,file='self.dat')
       endif
       do i_spin = 1, n_spin, 1
         do i_state =1 , n_states

           ! interpolation of the self_energy (sum of 2 poles)
           call mr_min_lsq (n_freq,dcmplx(0.d0,omega_grid(1:n_freq)), &
               self_energy_omega(i_state,1:n_freq,i_spin), &
               n_par, self_en_par)

        do i_freq= 1, aux_nomega, 1
          w = (0.0,1.0) *aux_omega(i_freq)

         numerator = dcmplx(0.d0,0.d0)
         denominator = dcmplx(0.d0,0.d0)
!         if(mod(n_par,2) .eq.0) then

            do i_par = n_par/2, 2, -1
              numerator = (numerator + self_en_par(i_par))*w
            enddo

            do i_par = n_par, n_par/2+1, -1
              denominator = (denominator + self_en_par(i_par))*w
            enddo

!         else
!
!            do i_par = n_par/2+1, 2, -1
!              numerator = (numerator + self_en_par(i_par))*w
!            enddo
!
!            do i_par = n_par, n_par/2+2, -1
!              denominator = (denominator + self_en_par(i_par))*w
!            enddo
!
!         endif

         numerator = numerator + self_en_par(1)
         denominator = denominator + dcmplx(1.d0,0.d0)

         fitted_sigma (i_freq) = numerator/denominator

!             a = (1.d0/self_en_par(2))
!             b = (self_en_par(1)/self_en_par(2))     

!             if(i_state == 1)then
!               write(64,*) aux_omega(i_freq), &
!!                 real(self_energy_omega(i_state,i_freq,i_spin))-&
!                 real (fitted_sigma(i_freq))
!             endif

!            print *, "ciao"

           !construct the KS Green's function in the KS basis
           !using the expression :
           ! G_n(iw) = 1/(iw - e_n + mu )

!           green_fn (i_state, i_freq, i_spin) = &
!             1.d0/((0.d0,1.d0)*aux_omega(i_freq)-&
!             (KS_eigenvalue(i_state, i_spin,1)-chemical_potential))
 
            GW_corr_fit = GW_corr_fit + &
!            green_fn (i_state, i_freq, i_spin)*&
             1.d0/((0.d0,1.d0)*aux_omega(i_freq)-&
             (KS_eigenvalue(i_state, i_spin,1)-chemical_potential))*&
             fitted_sigma (i_freq) *&
             aux_womega(i_freq)/pi/2.d0/n_spin

           enddo

           do i_freq = 1, n_freq, 1

           green_fn (i_state, i_freq, i_spin) = &
             1.d0/((0.d0,1.d0)*omega_grid(i_freq)-&
             (KS_eigenvalue(i_state, i_spin,1)-chemical_potential))
    

           GW_corr = GW_corr + &
            green_fn (i_state, i_freq, i_spin)*&
            self_energy_omega(i_state, i_freq, i_spin)*&
            womega(i_freq)/pi/2.d0/n_spin

            w = (0.0,1.0) *omega_grid(i_freq)
            numerator = dcmplx(0.d0,0.d0)
            denominator = dcmplx(0.d0,0.d0)
            do i_par = n_par/2, 2, -1
              numerator = (numerator + self_en_par(i_par))*w
            enddo
            do i_par = n_par, n_par/2+1, -1
              denominator = (denominator + self_en_par(i_par))*w
            enddo
            numerator = numerator + self_en_par(1)
            denominator = denominator + dcmplx(1.d0,0.d0)
            correction  = correction -(numerator/denominator &
               - self_energy_omega(i_state, i_freq, i_spin))*&
            green_fn (i_state, i_freq, i_spin)*&
            womega(i_freq)/pi/2.d0/n_spin



          enddo
         enddo
        enddo

       if(myid.eq.0)then 
         close (64)
       endif

        tot_en_GW = total_energy - en_xc + &
         eex_energy + GW_corr

        tot_en_GW_fit = total_energy - en_xc + &
         eex_energy + GW_corr_fit


!       if(myid.eq.0) print *, "exchange_energy1 :", eex_energy, 'Ha', eex_energy*hartree, 'eV'
!       if(myid.eq.0) print *, "GW_corr          :", GW_corr, 'Ha', GW_corr*hartree, 'eV'
!       if(myid.eq.0) print *, "GW_Tot           :", tot_en_GW, 'Ha', tot_en_GW*hartree, 'eV'

       if (myid.eq.0)then
        write(use_unit,*) " "
        write(use_unit,'(A)') "   --- GW Total Energy Calculation"
        write(use_unit,'(A)') "          |"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | DFT Total Energy                  :"&
                 ,total_energy * hartree, " eV   ", total_energy , " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Galitskij-Migdal Total Energy     :",&
          tot_en_GW* hartree, " eV   ", &
          tot_en_GW, " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Galitskij-Migdal Total Energy Fit :",&
          tot_en_GW_fit* hartree, " eV   ", &
          tot_en_GW_fit, " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | GW correlation Energy             :", &
            GW_corr* hartree, " eV   ", &
            GW_corr , " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | GW correlation Energy Fit         :", &
            GW_corr_fit* hartree, " eV   ", &
            GW_corr_fit , " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Exchange Energy                   :", &
            eex_energy* hartree, " eV   ", &
            eex_energy, " Ha"
        write(use_unit,'(A,F16.6,A,F16.6,A)'  ) "          | Correction                        :", &
            correction* hartree, " eV   ", &
            correction, " Ha"
       endif


       if( allocated(green_fn)) then
         deallocate(green_fn)
       endif

       return
       end subroutine get_gw_tot_en 
