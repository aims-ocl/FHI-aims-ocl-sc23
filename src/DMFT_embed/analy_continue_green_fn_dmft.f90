!****s* FHI-aims/analy_continute_self_energy
!  NAME
!   analy_continute_self_energy
!  SYNOPSIS
!  
      subroutine analy_continue_green_fn_dmft &
          (anacon_type,&
           n_freq, &
           n_max_par, &
           green_fn_par,omega, &
           green_fn_freq, n_matrix, i_k_point)

!  PURPOSE
!  Subroutine analy_continue_self_energy  analytically continue 
!  the self energy from imaginary energy axis to real energy axis.
!  Two schemes can be used here: the two-pole fitting or Pade approximation 
!
! USES

      use dimensions
      use mpi_tasks

      implicit none

! ARGUMENTS
      integer :: n_matrix
      integer ::  anacon_type
      integer ::  n_freq
      integer ::  n_max_par 

      real*8 ::   omega(n_freq)
      complex*16 :: green_fn_par (n_max_par,n_matrix,n_matrix)
      complex*16 :: green_fn_freq (n_matrix,n_matrix,n_freq)
!      real*8 inv_overlap_matrix(n_basis,n_basis)

!extra
      complex*16, dimension(:,:), allocatable ::  green_tmp
      complex*16, dimension(:,:), allocatable ::  aux_ovlp
      character*40 filename 
      character*2  basis_el
      character*2 k_point 
! INPUTS 
!  o  flag_qpe -- logical, if true, read the self-energy from files, if false, the
!         self-energy comes from the calling subroutine
!  o  anacon_type -- integer number, if 0, the two-pole fitting for analytical
!          continuation; if 1, using Pade approximation for ana. cont.
!  o  n_low_state  -- integer number,
!          the lowest KS/HF eigenstate for self-energy correction
!  o  n_high_state -- integer number,
!          the highest KS/HF eigenstate for self-energy correction
!  o  n_freq -- integer number, the number of frequency points for the GW self-energy
!  o  n_max_par -- the number of parameters used for analytical continuation
!          For anacon_type = 0, recommended n_max_par is  4 or 5. If 4, this will
!          be the normal two-pole fitting, else if 5, it will be two-pole plus a
!          (small) constant number
!          For anacon_type = 1, recommended n_max_par is the half of n_freq
!  o  omega(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the self-energy
!  o  chemical_potential -- real number, the chemical potential of the system
!  o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
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

!internal
      complex*16 fitted_green (n_freq)
      complex*16 denominator
      complex*16 numerator
      complex*16 w      
      complex*16 xdata (n_max_par)
      complex*16 gtmp
      logical test
 
! counters
      integer :: i_freq, i_k_point
      integer :: i_basis
      integer :: j_basis
      integer :: i_spin
      integer :: i_par 
      integer :: n_step
      integer :: i_dat

!     begin work
! perform analytical continuation, and get the fitting parameters 
!      do i_spin = 1, n_spin, 1


      do j_basis = 1, n_matrix, 1
        do i_basis = 1, n_matrix, 1

            if(anacon_type .eq. 0) then
             call mr_min_lsq_dmft(n_freq,dcmplx(0.d0,omega), & 
                        green_fn_freq(j_basis, i_basis, :), &
                        n_max_par, green_fn_par(:,j_basis,i_basis), i_k_point) 

            elseif(anacon_type .eq. 1) then

             call get_pade_approx_para(n_freq,omega, &
                        green_fn_freq(j_basis, i_basis, :), &
                        n_max_par, green_fn_par(:,j_basis,i_basis))
            endif

        enddo
      enddo


!test if it works, this section compare the numerical function (on the imaginary axis)
! with its analytical fits, check the output, if everything is fine they should overlap

      test = .false.
      if (test)then

        fitted_green(:) = (0.d0,0.d0)

        if (anacon_type.eq.0)then

        do i_basis = 1, n_matrix, 1

         if (i_basis .lt.10) then
             write(basis_el, '(A,I1)')'0',i_basis
         else
             write(basis_el, '(I2)')i_basis
         endif


         if (i_k_point .lt.10) then
            write(k_point, '(A,I1)')'0',i_k_point
         else
            write(k_point, '(I2)')i_k_point
         endif 

            filename = 'G_an_con_'//basis_el//'k_'//k_point//'.dat'
          
!         filename = 'G_an_con_'//basis_el//'.dat'

         if(myid.eq.0)then
           open (25, file=filename)
         endif
 
        do i_freq= 1, n_freq, 1     
          w = (0.0,1.0) *omega(i_freq)

         numerator = dcmplx(0.d0,0.d0)
         denominator = dcmplx(0.d0,0.d0)
         if(mod(n_max_par,2) .eq.0) then

            do i_par = n_max_par/2, 2, -1
              numerator = (numerator + green_fn_par(i_par,1,i_basis))*w
            enddo

            do i_par = n_max_par, n_max_par/2+1, -1
              denominator = (denominator + green_fn_par(i_par,1,i_basis))*w
            enddo

         else

            do i_par = n_max_par/2+1, 2, -1
              numerator = (numerator + green_fn_par(i_par,1,i_basis))*w
            enddo

            do i_par = n_max_par, n_max_par/2+2, -1
              denominator = (denominator + green_fn_par(i_par,1,i_basis))*w
            enddo

         endif

         numerator = numerator + green_fn_par(1,1,i_basis)
         denominator = denominator + dcmplx(1.d0,0.d0)

         fitted_green (i_freq) = numerator/denominator
 

         if(myid.eq.0)then       
          write(25,*) omega(i_freq), real(green_fn_freq (1,i_basis,i_freq)),&
          real( fitted_green (i_freq))
         endif

        enddo
       if(myid.eq.0)then       
        close(25)
       endif
        enddo

        elseif(anacon_type.eq.1)then
          do i_basis = 1, n_matrix, 1
            if (i_basis .lt.10) then
                write(basis_el, '(A,I1)')'0',i_basis
            else 
                write(basis_el, '(I2)')i_basis
            endif
                write(k_point, '(I6)')i_k_point
                 
            filename = 'G_an_con_'//basis_el//'k_'//k_point//'.dat'
  
            if(myid.eq.0)then       
              open (25, file=filename)
            endif
      
            do i_freq= 1, n_freq, 1

              w = (0.0,1.0) *omega(i_freq)  
              n_step = n_freq/(n_max_par-1)
              i_dat = 1
       
              do i_par = 1, n_max_par-1, 1
                xdata(i_par) = dcmplx(0.d0,omega(i_dat))
                i_dat = i_dat + n_step
              enddo
       
              xdata(n_max_par) = dcmplx(0.d0,omega(n_freq))  
              gtmp = dcmplx(1.d0,0.d0)
       
              do i_par = n_max_par, 2, -1
                gtmp = 1.d0 + green_fn_par(i_par,i_basis,i_basis)&
                             *(w-xdata(i_par-1))/gtmp
              enddo
       
              fitted_green (i_freq) = green_fn_par(1,i_basis,i_basis)/gtmp  

              if(myid.eq.0)then       
                write(25,*) omega(i_freq),&
                 real(green_fn_freq (i_basis,i_basis,i_freq)),&
                 real( fitted_green (i_freq))
              endif
            enddo
          enddo
          if(myid.eq.0)then       
            close(25)
          endif
        endif
      endif
      !endif

      return 
      end subroutine analy_continue_green_fn_dmft
!******
!---------------------------------------------------------------------
