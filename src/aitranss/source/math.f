c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           march/april 2008
c     f-functions update: sept 2008
c     last revision:  jan 2012
c###########################################################

      module math_functions
c      ************************************************************
c      some math functions needed to calculate overlap integrals
c      original version: ferdinand evers
c      updated for f-functions and adapted for f90 : alexej bagrets
c      ************************************************************
       implicit none
c      http://www.mathwithmrherte.com/pi_digits.htm
       double precision, parameter :: pi=3.14159265358979323846d0

      contains

      double precision function r00(mxi,nxi,r2)
c     ******************************************************
c     the prefactor coming from convolution of two gaussians      
c     ******************************************************
       implicit none
       double precision mxi, nxi, r2
       double precision xisum, gnorm, g3norm
ccc    normalization is approved!  # f.evers had the same norm
       xisum  = mxi + nxi
       gnorm  = sqrt(pi/xisum)
       g3norm = gnorm*gnorm*gnorm
       r00 = g3norm*exp(-mxi*nxi*r2/xisum)

      end function r00

      double precision function p(xia,xib,x_a,x_b,n,m)
c     *************************************************
c     p integrates a polynom weighted with the gaussian
c     *************************************************
       implicit none
       double precision xia, xib,  ! exponents
     &                  x_a, x_b   ! positions
       integer          n,m        ! polynoms pows
ccc                     caution! a <--> n, b <--> m
       double precision r, g1, g2, g3, a, b, a3b, 
     &                  xia2, xib2, tmp, tmpsum  
   
       r  = x_a-x_b  ! caution, it should be so  
       g1 = xia+xib
       g2 = g1*g1 
       g3 = g2*g1
       if (m+n == 0 ) then
        p = 1.0d0
       else if (m+n == 1) then
        p = -(1-m)*xib*r/g1 + m*xia*r/g1
       else if (n+m==2) then
        select case (n)
         case (2) ; tmp = xib*r/g1 
	            p = 0.5d0/g1 + tmp*tmp  
	 case (1) ; tmp = r/g1  
	            p = 0.5d0/g1 - xia*xib*tmp*tmp
	 case (0) ; tmp = xia*r/g1
	            p = 0.5d0/g1 + tmp*tmp 
	 case default;	 p = 1.0d0
	           print *
  	           stop '[FUNCTION p]: n or m are negative???' 
         end select
       else if (n+m==3) then
        select case (n)
         case (3) ; tmp = xib*r 
                    p = -(1.5d0 + tmp*tmp/g1)*tmp/g2
	 case (2) ; tmp = xib*r 
	            p = (-xib + 0.5d0*xia + xia*tmp*tmp/g1)*r/g2 
	 case (1) ; tmp = xia*r
	            p = ( xia - 0.5d0*xib - xib*tmp*tmp/g1)*r/g2 
         case (0) ; tmp = xia*r 
	            p = (1.5d0 + tmp*tmp/g1)*tmp/g2
         case default ;  p = 1.0d0
	           print *
		   stop '[FUNCTION p]: n or m are less than 0 ???'
 	end select
       else if (n+m==4) then
        select case (n) 
         case (1)    
	    tmp = r*r/g1
	    a = xia*(xia-xib)
            b = xia*xia*xia*xib
	    p = (3.0d0/4.0d0 + 1.5d0*a*tmp - b*tmp*tmp)/g2
         case (2) 
	    a=(xia*xia + xib*xib)/2.0d0 - 2.0d0*xia*xib
	    b=xia*xia*xib*xib 
	    tmp = r*r/g1
            p = (3.0d0/4.0d0 + a*r*r/g1 + b*tmp*tmp)/g2
         case (3)
            tmp = r*r/g1
	    a = xib*(xib-xia)
            b = xib*xib*xib*xia
	    p = (3.0d0/4.0d0 + 1.5d0*a*tmp - b*tmp*tmp)/g2
         case default 
	    p = 1.0d0
	    print *
	    stop '[FUNCTION p]: n is neither 1, 2, or 3 ???' 
 	end select
       else if (n+m==5) then
        select case (n)
	 case (3)     
	    tmp = r*r/g1 ;  xia2 = xia*xia ; xib2 = xib*xib
	    tmpsum=0.5d0*(xib2+3.0d0*xia2-6.0d0*xia*xib)+xia2*xib2*tmp
            p = (0.75d0*(2.0d0*xia-3.0d0*xib) - tmpsum*tmp*xib)*r/g3
         case (2)
	    tmp = r*r/g1 ; xia2 = xia*xia ; xib2 = xib*xib
            tmpsum=0.5d0*(xia2+3.0d0*xib2-6.0d0*xia*xib)+xib2*xia2*tmp
            p =-(0.75d0*(2.0d0*xib-3.0d0*xia) - tmpsum*tmp*xia)*r/g3
	 case default
	    p = 1.0d0
	    print *
	    stop '[FUNCTION p]: n+m==5, but n or m is out of range!' 
         end select
       else if (n+m==6) then
        select case (n)
	 case (3)
	    tmp = r*r/g1 ;  xia2 = xia*xia ; xib2 = xib*xib
	    a3b = xia2 + xib2 - 3.0d0*xia*xib
	    tmpsum =  1.5d0*a3b*(1.5d0 - xia*xib*tmp)  
     &               -xia2*xia*xib2*xib*tmp*tmp  
            p = (15.0d0/8.0d0 + tmp*tmpsum)/g3
	 case default
            p = 1.0d0
	    print *
	    stop '[FUNCTION p]: n+m==6, but n is out of range!' 
        end select
       else
            p = 1.0d0
            print *
            stop '[FUNCTION p]: case n+m>6 is not accepted!' 
       end if

      end function p
  
      end module math_functions

