      subroutine get_spectrum_dmft_GW (anacon_type, green_fn_par, n_max_par, &
             omega, nomega,filename, n_matrix , new_chem_pot,&
             spectrum,aux_nomega, aux_omega , aux_womega)

!get the spectrum from the (analytically continued) diagonalized Green's function


      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use physics
      use localorb_io, only: use_unit

!ARGUMENTS
      implicit none
      integer n_matrix
      integer n_max_par
      integer nomega, i_k_points
      integer anacon_type
      integer aux_nomega
      complex*16 green_fn_par (n_max_par, n_matrix, n_matrix)
!      complex*16 aux_green_KS (n_max_par, n_states)
      real*8 omega(nomega)
      real*8  spectrum(-aux_nomega:aux_nomega)
!      real*8  spectrum_partial(-aux_nomega:aux_nomega)
      real*8  aux_omega(aux_nomega)
      real*8  aux_womega(aux_nomega)
      !real*8,  dimension (:), allocatable ::  spectrum
      real*8,  dimension (:), allocatable ::  spectrum_partial
      character*15 filename
      real*8 new_chem_pot

!internal
      complex*16 diagonal_green (n_matrix)
      complex*16 aux_green (n_matrix, n_matrix)
      complex*16 denominator
      complex*16 numerator
!      real*8     w
      complex*16 w 
      complex*16 xdata (n_max_par)
      complex*16 gtmp
      real*8 n_particles
!the local grid
!      integer aux_nomega
!      real*8  aux_omegamax
      !real*8 , dimension (:), allocatable :: aux_omega
      !real*8 , dimension (:), allocatable :: aux_womega

!counter
      integer :: i_basis
      integer :: j_basis
      integer :: i_freq
      integer :: i_par
      integer :: i_state
      integer :: sign_freq
      integer :: n_step
      integer :: i_dat
      integer i_index
      real*8 e_fermi
      integer n_homo
 
!define a grid on which the spectrum is calculated
      e_fermi = chemical_potential
!      e_fermi = 0.d0
!      e_fermi =   new_chem_pot

!      aux_nomega = 15000
!      aux_omegamax = 6.50
!      allocate(aux_omega(aux_nomega))
!      allocate(aux_womega(aux_nomega))
!      allocate(spectrum(-aux_nomega:aux_nomega))
      allocate(spectrum_partial(-aux_nomega:aux_nomega))

!      call gauleg(0d0, aux_omegamax, aux_omega, aux_womega, aux_nomega)

      do i_state = 1, n_states
        if(occ_numbers(i_state,1,1).gt.1.d-6) then
          n_homo=i_state
        endif
      enddo 

!start work
!      if(myid.eq.0)then
!       write(use_unit,*) "Calculating the Spectrum"
!      endif

      spectrum (:) = 0.d0

      if (myid.eq.0)then
        open (44,file=filename)
      !  open (45,file='spectrum_homo.dat')
      endif

      if (anacon_type.eq.0)then
       do sign_freq = -1, 1, 2
        do i_freq = 1, aux_nomega, 1
         do i_basis=1, n_matrix, 1

         j_basis = i_basis         
         numerator = (0.d0,0.d0)
         denominator = (0.d0,0.d0)
         w = sign_freq*aux_omega(i_freq) 

         if(mod(n_max_par,2) .eq.0) then

            do i_par = n_max_par/2, 2, -1 
              numerator = (numerator + &
               green_fn_par (i_par,i_basis,j_basis))&
              *(w)
            enddo

            do i_par = n_max_par, n_max_par/2+1, -1 
              denominator = (denominator + &
               green_fn_par (i_par,i_basis, j_basis))&
              *(w)
            enddo

         else

            do i_par = n_max_par/2+1, 2, -1 
              numerator = (numerator + &
               green_fn_par (i_par,i_basis,j_basis))&
              *(w)
            enddo

            do i_par = n_max_par, n_max_par/2+2, -1 
              denominator = (denominator + &
               green_fn_par (i_par,i_basis,j_basis))&
              *(w)
            enddo

         endif
         
         numerator = numerator +  green_fn_par (1,i_basis,j_basis)
         denominator = denominator + dcmplx(1.d0,0.d0)
         aux_green (i_basis, j_basis ) = numerator/denominator

        enddo !i_basis

        i_index = 0 
        do i_basis = 1, n_matrix , 1
          j_basis = i_basis

          i_index = i_index+1

          spectrum (sign_freq*i_freq) = spectrum (sign_freq*i_freq) -&
             (aimag (aux_green (j_basis,i_basis)))*&
             2./pi
 
          if(i_basis==n_homo)then
             spectrum_partial (sign_freq*i_freq) = &
             spectrum_partial (sign_freq*i_freq) -&
             (aimag (aux_green (j_basis,i_basis)))*&
             2./pi
          endif
        enddo


        if (myid.eq.0)then
             write(44, *)(real(w)+e_fermi)*hartree,&
                        spectrum (sign_freq*i_freq)
!             write(45, *)(real(w)+e_fermi)*hartree,&
!                        spectrum_partial (sign_freq*i_freq)
        endif

       enddo  ! i_freq
      enddo !sign_freq

      elseif(anacon_type.eq.1)then
         do sign_freq = -1, 1, 2
          do i_freq= 1, aux_nomega, 1
           do i_basis=1, n_matrix, 1
             j_basis = i_basis

            w = sign_freq *aux_omega(i_freq)! *(0.d0,1.d0)
 
            n_step = aux_nomega/(n_max_par-1)
            i_dat = 1
            do i_par = 1, n_max_par-1, 1
              !xdata(i_par) = dcmplx(0.d0,aux_omega(i_dat))
              xdata(i_par) = dcmplx(aux_omega(i_dat),0.d0)
              i_dat = i_dat + n_step
            enddo
            xdata(n_max_par) = dcmplx(aux_omega(aux_nomega),0.d0)

            gtmp = dcmplx(1.d0,0.d0)
            do i_par = n_max_par, 2, -1
              gtmp = 1.d0 + green_fn_par(i_par,i_basis,j_basis)*(w-xdata(i_par-1))/gtmp
            enddo

            aux_green (i_basis,j_basis) = green_fn_par(1,i_basis,j_basis)/gtmp

           enddo
         
           i_index = 0
           do i_basis = 1, n_matrix , 1
             j_basis = i_basis

             i_index = i_index+1

              spectrum (sign_freq*i_freq) = spectrum (sign_freq*i_freq) +&
                 abs(aimag (aux_green (i_basis,j_basis)))*&
                 2./pi

             if(i_basis==n_homo)then
               spectrum_partial (sign_freq*i_freq) = &
               spectrum_partial (sign_freq*i_freq) +&
               abs(aimag (aux_green (j_basis,i_basis)))*&
               2./pi
             endif
           enddo

 
           if (myid.eq.0)then
             write(44, *)(real( w)+e_fermi)*hartree,&
                        spectrum (sign_freq*i_freq)
!             write(45, *)(real( w)+e_fermi)*hartree,&
!                        spectrum_partial (sign_freq*i_freq)
           endif

          enddo
!           close(777)
         enddo
      endif !anacon_type

      if (myid.eq.0)then
        close(44)
!        close(45)
      endif

! the integral of the spectrum should give the number of particle!
! let's check

      if(.false.)then
         n_particles = 0.d0
         do i_freq= 1, aux_nomega, 1
          if( aux_omega(i_freq) .lt. chemical_potential )then
            n_particles = n_particles + spectrum(i_freq)*aux_womega(i_freq)
          else
            n_particles = n_particles +&!  spectrum(i_freq)*aux_womega(i_freq)+&
              spectrum(-i_freq)*aux_womega(i_freq)
          endif
            
         enddo     
       
        if (myid.eq.0)then
          write(use_unit,*) " "
          write(use_unit,*) "n_particles as calculated from the spectrum = ", &
            n_particles
          write(use_unit,*) " "
        endif
      endif
!stop
      return

      end subroutine get_spectrum_dmft_GW
