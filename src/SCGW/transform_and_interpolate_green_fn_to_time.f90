subroutine  transform_and_interpolate_green_fn_to_time (green_fn_freq,  &
         tau,ntau,wtau, omega, nomega, womega,  &
         green_fn_time, inv_overlap_matrix , omegamax,  ovlp_NAO_KS) 

!transform the Green's function (expressed in the NAO basis) from the frequency
!to the time domain. To help in treating the slow decaying of the green function 
!(in frequency) the last 2 points of the tail are interpolated with a function
!
! g(w) = a/(b + i w)  
!
!where the coefficients a and b are complex. 
!Hence we calculate F(G(w) - g(w)) + F(g(w))

!This is formally exact and helps it treating the tail.

      use dimensions
!      use runtime_choices
!      use species_data
      use physics
      use prodbas
!      use hartree_fock
!      use gw_para
      use synchronize_mpi
      use constants
      use mpi_tasks
      use localorb_io, only: use_unit

 implicit none

!INPUT 
 integer ntau, nomega
 real*8  tau(-ntau:ntau)
 real*8 wtau(-ntau:ntau)
 real*8  omega(nomega)
 real*8 womega(nomega)
 real*8 green_fn_time (n_basis,n_basis,-ntau:ntau)
 complex*16 green_fn_freq (n_basis,n_basis,nomega)
 complex*16 tmp_green_fn_freq (n_states,n_states,nomega) 
 real*8 inv_overlap_matrix (n_basis, n_basis)
 real*8 ovlp_NAO_KS(n_states,n_basis)
 real*8 omegamax
! EXTRA STUFFi
 real*8, dimension(:,:), allocatable :: green_tmp
 integer :: i, j, k, n_max
 integer :: i_basis, j_basis, k_basis
 integer :: i_state, j_state
 integer :: i_tau
 real*8 norm
 real*8 aux_green_fn_time(n_states,n_states, -ntau:ntau)
 logical output
 character*14 filename
 character*2 iter

!the coefficient for the tail interpolation
 complex*16 a,b, a1,b1, a2,b2 
 real*8     x1,x2
 complex*16 G1,G2
 real*8     threshold

!auxiliary function for the interpolation
! complex*16 aux_green_fn_freq (nomega)
 complex*16 aux_inv_overlap_matrix (n_basis, n_basis)
 integer i_freq
 logical tail_fit_2
 real*8 gre1
 integer window_points 
!the local grid
 integer i_index
 integer i_spin
 integer n_fit
 integer n_max_par
 real*8 point
 complex*16 par(2)

 if (myid.eq.0)then
   write(use_unit,*)"  --- Fourier transform of the Green's function from frequency to time ---"
 endif

 aux_green_fn_time(:,:,:) = 0.d0
 threshold = 1.d-15

 ! this transformt the G(w) in the KS basis, 
 ! where it has a better analytical structure and it's
 ! easier to do the analytical tail fit
 call diagonalize_green_fn &
      (green_fn_freq, &
      tmp_green_fn_freq)

 !first evaluate the interpolation of the tail of the green's
 ! function from the last 2 grid points

  a = 0.d0
  b = 0.d0
  i_index = 0
  n_max_par = 2
 
  if(myid.eq.0)then
    open(55, file = 'test_fit.dat')
  endif

  do i_basis =1, n_states, 1
    do j_basis =1, n_states, 1
 
       G1 = tmp_green_fn_freq (i_basis,j_basis,nomega-1)
       G2 = tmp_green_fn_freq (i_basis,j_basis,nomega)
!       x1 = omega(nomega-1)
!       x2 = omega(nomega)

      if(abs(G2).gt.1.d-15)then  

       n_fit = nomega/4
       call mr_min_lsq(n_fit,dcmplx(0.d0,omega(nomega-n_fit:nomega)), &
            tmp_green_fn_freq(i_basis,j_basis,nomega-n_fit:nomega), &
            n_max_par,par)

       b = real (1.d0/par(2))
       a = real(par(1)/par(2)) 

       if(.true.)then ! this block print out the fit
        if(myid.eq.0 )then
         do i_freq = 1,nomega,1
          point = omega(i_freq)
          write(55, *) point,&! real(par(1)/(1.d0+(0.d0,1.d0)*par(2)*point)), &
                     real( a/(b+(0.d0,1.d0)*point)), real(tmp_green_fn_freq(i_basis,j_basis,i_freq))
         enddo
         do i_freq = 1, 100
         point = omega(nomega)+0.1*i_freq  
          write(55, *) point, real( a/(b+(0.d0,1.d0)*point))  ! real(par(1)/(1.d0+(0.d0,1.d0)*par(2)*point))
         enddo
        endif
       endif

!         b = -real((((0.d0,1.d0)*(G1*x1-G2*x2))/(G1-G2)))
!         a =  -real((((0.d0,1.d0)*(G1*G2*x1 - G1*G2*x2))/(G1-G2)))

!       if(myid.eq.0.and. i_basis .eq.1 .and.j_basis .eq.1)then
!         print *, a, b 
!       endif 

         if(abs(real(b)).lt.1.d-20)then
          if(myid.eq.0)then
           write(use_unit,*) " ***** WARNING: the fitting function could have a pole! *****"
          endif
         endif

! subtract the interpolation to the green's function
         do k =1, nomega, 1
           tmp_green_fn_freq (i_basis,j_basis,k) =&
           tmp_green_fn_freq (i_basis,j_basis,k)&
           - a/(b+(0.d0,1.d0)*omega(k))
         enddo
   
      endif ! abs(G2).gt.1.d-18

      if(.false.)then
      window_points = 20
      if(nomega.gt.window_points+10)then
        i_index = 0
        do k = nomega - window_points , nomega, 1
          i_index=i_index+1
          tmp_green_fn_freq (i_basis,j_basis,k)  =&
          tmp_green_fn_freq (i_basis,j_basis,k) * 0.5d0* &
          (1.+dcos(2*pi *i_index / (2*window_points-1)))
        enddo
      endif
      endif

         do i=-ntau, ntau, 1
           do k=1, nomega, 1
             norm = womega(k)/pi

             aux_green_fn_time(i_basis,j_basis,i) =  aux_green_fn_time(i_basis,j_basis,i) + &
             tmp_green_fn_freq(i_basis,j_basis,k)* exp(cmplx(0,-omega(k)*tau(i)))*norm

           enddo

        !some output 
        output = .false.
        if(myid.eq.0)then
          if(output)then
            open(242,file="FT_0.dat")
            open(243,file="FT_1.dat")
            if(abs(G2).gt.1.d-18)then
              if (i_basis .eq.9 .and. j_basis .eq.2)then

                write(242,*) tau(i), aux_green_fn_time(i_basis,j_basis,i)

                if(real(b).ge.0)then
                  if (i.le.0)then

                    write(243,*) tau(i), &
                      real(a*exp(b*tau(i)))
                  endif
                elseif(real(b).lt.0)then
                  if (i.gt.0)then
                    write(243,*)  tau(i),&
                      real(-a*exp(b*tau(i)))
                  endif
                endif
              endif
            endif
            close(242)
            close(243)
          endif
        endif

! add the FT of the interpolation of the tail
         if(abs(G2).gt.1.d-15)then
           if(real(b).ge.0)then
             if (i.le.0)then
               aux_green_fn_time(i_basis,j_basis,i) = aux_green_fn_time(i_basis,j_basis,i) + &
                a*exp(b*tau(i))
             endif
           elseif(real(b).lt.0)then
             if (i.gt.0)then
               aux_green_fn_time(i_basis,j_basis,i) = aux_green_fn_time(i_basis,j_basis,i) - &
                a*exp(b*tau(i)) 
             endif
           endif
         endif

         enddo! i_tau
      
     
     enddo
   enddo
   if(myid.eq.0)then
     close(55)
   endif

    green_fn_time(:,:,:) = 0.d0

     ! back to the NAO basis
    if(.not.allocated(green_tmp))then
       allocate(green_tmp (n_basis,n_states))
    endif
    do i_tau = -ntau, ntau, 1
      green_tmp = 0.d0
      call dgemm ('N','N', n_basis, n_states, n_states, &
        1.d0, KS_eigenvector, n_basis,& 
        aux_green_fn_time (1,1,i_tau) , n_states, 0.d0,& 
        green_tmp, n_basis )
      call dgemm ('N','T', n_basis, n_states, n_states, &
        1.d0, green_tmp, n_basis , KS_eigenvector, n_basis,&
         0.d0, green_fn_time(1,1,i_tau), n_basis )
    enddo
    if(allocated(green_tmp))then
       deallocate(green_tmp)
    endif


   ! print out the Fourier Trasnsform
   if(myid.eq.0)then
   output = .true.
   if( output.and.myid.eq.0)then
     do i_basis = 1, n_basis, 1
       if( i_basis.lt.10 ) then
          write(iter,'(A,I1)') "0",i_basis
       else
          write(iter,'(I2)') i_basis
       endif
       filename = "green_FT"//iter//".dat"
       open(77, file=filename)
         do i_tau = -ntau, ntau, 1
           write(77,*) tau(i_tau),green_fn_time(1,i_basis,i_tau)
         enddo
       close(77)
     enddo
   endif
   endif

  return
end subroutine transform_and_interpolate_green_fn_to_time
