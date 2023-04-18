C
CVB Modified solver for radial equation, derived from D.R. Hamanns dftseq solver.
C   The present version integrates the radial equation for fixed energy only, 
C   i.e. does not solve an actual eigenvalue problem. 
C
C   Ultimately, this is duplication of code. Should rip out the actual integration part, and 
C   use the same in both here (once) and in dftseq (multiple times)
C
c  Header:
c***********************************************************************
c
c integrates radial (pauli-type scalar-relativistic w/ spin-orbit av.) 
c equation on a c logarithmic mesh w/ fixed/variable effective mass 
c
c input
c mode ...... 1 is for full potential bound state (scalar relativistic)
c             2 is for nonrelativistic bound state 
c z ......... nuclear charge
c mmax ...... number of radial grid points
c r() ....... radial grid
c n ......... principal quantum number
c l ......... angular momentum qn
c efm ....... effective mass
c uld ....... logarithmic derivative at mch
c v() ....... potential v*Psi
c v1() ...... potential v1*grad(Psi)
c v2() ...... potential 1/v2 * laplace(Psi)
c
c output
c nin ....... outermost gridpoint used
c mch ....... classical turning point
c e ......... eigenvalue
c u() ....... wave function
c up() ......               1st derivative on log grid
c upp() .....               2nd derivative on log grid
c
c vb's variables
c 
c r_outer ... radius for outermost grid point for integration -- unless zero.
c
c original version by D.R. Hamann, gncpp
c
c revised version from fhiPP code by Martin Fuchs FHI Berlin 
c
c Volker Blum, FHI 2004: 
c adapted dftseq to my needs. Original version
c assumed known limiting behavior of H-like wave functions as r -> infty,
c but we have confining potential.
C FIXME: We use the incorrect outer boundary conditions of the free atom (exponential
C        decay of WF
C FIXME: We should use WKB-like boundary condition when integrating towards the inside.
c
c***********************************************************************
c
      subroutine solve_radial_eq
     + ( mode,z,mmax,r,n,l,efm,v,v1,v2,
     +   nin,mch,uld,e,u,up,upp, r_outer )
c
      implicit real*8 (a-h,o-z)
c
      include 'parameter.h'
      include 'default.h'

      parameter (itmax=200)
      parameter (fss_i=137.0359896d0)
c
      logical trap,trel
      common/err/trap(10)
c
      dimension r(mx),u(mx),up(mx),upp(mx)
      dimension cf(mx),dv(mx),fr(mx),frp(mx),v(mx),v1(mx),v2(mx),
     1          x1(mx),x2(mx),x4(mx),erec(itmax)

cvb
      real*8 r_outer
      integer i_outer

C  begin work

C  -------------------------------------------------------------
C  begin preparations
C  -------------------------------------------------------------

      amesh=r(2)/r(1)
      al=log(amesh)
      sls=l*(l+1)

ctest
c      write(6,*) "Eigenvalue, incoming: ", e
ctest end

C     ics is the actual end point of radial grid (here synonymous with mmax)
C     ico is index of the classical turning point (it's really set below)

      ics=mmax
      ico=mch
        
C     initialize relativistic or non-relativistic

      if(mode .eq. 1) then
        trel=.true.
        fss=(1.0d0/fss_i)**2
        if(l .eq. 0) gamma=sqrt(1.0d0-fss*z**2)
        if(l .gt. 0) gamma=(l*sqrt(l**2-fss*z**2) +
     &   (l+1)*sqrt((l+1)**2-fss*z**2))/(2*l+1)
      else
        trel=.false.
        fss=1.0d-20
        gamma=l+1
      end if

C     alias the mode flag into imod

      imod=mode

C       Verify preset integration energy against specified potential

C       find highest possible bound state energy emax
C       find lowest  possible bound state energy emin

        if (r_outer.eq.0.) then
cvb     screened coulomb potential, always negative

          emin=0.0d0
          emax=v(mmax)+0.5d0*sls/r(mmax)**2

        else
cvb     we have a cutoff potential, so the outermost integration point
cvb     r_outer needs to contain the wave function

          i_outer = 1 + ( log( r_outer/r(1) ) )/al

          emin = v(i_outer)+0.5d0*sls/r(i_outer)**2
          emax = v(i_outer)+0.5d0*sls/r(i_outer)**2

        end if

c       now find minimum effective potential
        do 6 i=1,mmax
    6     emin=dmin1(emin,v(i)+0.5d0*sls/r(i)**2)

ctest
c        write(6,*) "emin, emax, e: ", emin, emax, e
ctest end

        if ( (e .gt. emax).or.(e.lt.emin) ) then
          write(6,'(A,A)') " * The specified eigenenergy lies",
     +      " outside the possible ranges for bound states. "
          write(6,'(A,G15.6,A)') " * Specified energy:", 
     +      e, " Ha"
          write(6,'(A,G15.6,A)') " * Minimum potential: ",
     +      emin, " Ha"
          write(6,'(A,G15.6,A)') " * Maximum potential: ",
     +      emax, " Ha"
          stop
        end if

c     initialize arrays and remove leftover garbage
      do i=1,4
        u(i)=0.d0
        up(i)=0.0d0
        upp(i)=0.0d0
      enddo

      do i=1,mmax
        x1(i)=al*r(i)
        x2(i)=x1(i)*x1(i)
      enddo

      aux=al*al*sls

      if(trel) then
c       calculate dv/dr for darwin correction
        dv(1)=(-50.d0*v(1)+96.d0*v(2)-72.d0*v(3)+32.d0*v(4)
     &         -6.d0*v(5))/(24.d0*x1(1))
        dv(2)=(-6.d0*v(1)-20.d0*v(2)+36.d0*v(3)-12.d0*v(4)
     &         +2.d0*v(5))/(24.d0*x1(2))

        do i=3,ics
          dv(i)=(2.d0*v(i-2)-16.d0*v(i-1)+16.d0*v(i+1)
     &           -2.d0*v(i+2))/(24.d0*x1(i))
        enddo
      else
c       nonrelativistic coefficient arrays for u (fr) and up (frp)
        do i=1,ics
          fr(i)=-2*x2(i)*v1(i)*v2(i)
          frp(i)=2*x1(i)*v1(i)*v2(i)
        enddo
      endif

C  -------------------------------------------------------------
c  original return point for bound state convergence 
C  -------------------------------------------------------------

C     nint is leftover counter for integration cycles
      nint=0

C     FIXME - The following should be comprised into integration subroutines of their own.
C     * prepare
C     * integrate_out
C     * integrate_in
C     * normalize

C  -------------------------------------------------------------
C  prepare all integration prefactors which depend on e
C     * cf(i)
C     * fr(i)  (relativistic)
C     * frp(i) (relativistic)
C     * classical turning point mch
C  -------------------------------------------------------------

c       coefficient array for u in differential eq.
        do i=1,ics
          cf(i)=aux+2*(v(i)-e)*v2(i)*x2(i)
        enddo

c       relativistic coefficient arrays for u (fr) and up (frp)
        if(trel) then
          do i=1,ics
            fr(i)=x2(i)*fss*(-(v(i)-e)**2 + 0.5d0*dv(i)/
     &      (r(i)*(1.0d0+0.5d0*fss*(e-v(i)))))
            frp(i)=-x1(i)*0.5d0*fss*dv(i)/(1.0d0+0.5d0*fss* (e-v(i)))
          enddo
        endif

c         find classical turning point for matching
          do i=ics,2,-1
            if(cf(i-1) .le. 0.d0 .and. cf(i) .gt. 0.d0) then
              mch=i
              ico=mch
              go to 40
            end if
          enddo

c         we reach this point only if no turning point was found

          write(ie,*) '& dftseq - stop: no classical turning point'
          write(ie,*) '&    nint n l e                      ',nint,n,l,e
          write(ie,*) '&    0-order coefficient to file err.dftseq_tp'
          open(10,file='err.dftseq_tp',status='unknown')
          do i=1,ics
            write(10,*) r(i),cf(i)/r(i)**2
          enddo
          close(10)
          trap(5)=.true.
          stop 

   40     continue

C  -------------------------------------------------------------
c  outward integration of the wave function
C  -------------------------------------------------------------

c       using predictor once, corrector twice (cf abramowitz & stegun)
c       integrate outward only up to classical turning point, ico.

c       initialize innermost four points of wave function
C       using limiting behavior of Coulomb-like wave function for r -> 0 (power series)
cmf     l=0: hydrogen-like start (improves eigenvalues and wavefunctions)

        do i=1,4

          if(l .eq. 0 .and. imod .eq. 2 .or. imod .eq. 1) then
            u(i)=r(i)**gamma*(1.d0-z*r(i)+0.5d0*(z*r(i))**2)
            up(i)=al*gamma*r(i)**gamma*((1.d0-z*r(i)+0.5d0*(z*r(i))**2))
     &          +al*r(i)**gamma*(-z*r(i)+(z*r(i))**2)
            upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
          else 
            u(i)=r(i)**gamma
            up(i)=al*gamma*r(i)**gamma
            upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
          endif

        enddo

C       node counts the number of nodes in the wave function
        node=0

        do i=4,ico-1
          u(i+1)=u(i)+aeo(up,i)
          up(i+1)=up(i)+aeo(upp,i)
          do it=1,2
            upp(i+1)=(al+frp(i+1))*up(i+1)+(cf(i+1)+fr(i+1))*u(i+1)
            up(i+1)=up(i)+aio(upp,i)
            u(i+1)=u(i)+aio(up,i)
          enddo
          if(u(i+1)*u(i) .le. 0.0d0) node=node+1
        enddo

c       store wave function value, derivative at the classical turning point
        uout=u(mch)
        upout=up(mch)

        if (node-n+l+1 .lt. 0) then
c         too few nodes - would need to raise energy

          write(6,*) "* After outward integration: "
          write(6,*) "* The wavefunction contains fewer nodes",
     +      " than an eigenfunction would."

        else if ( (node-n+l+1) .gt. 0) then
c         too many nodes - would need to lower energy

          write(6,*) "* After outward integration: "
          write(6,*) "* The wavefunction contains more nodes",
     +      " than an eigenfunction would."

        end if

C  -------------------------------------------------------------
C  inward integration of the wave function
C  -------------------------------------------------------------

c       set outer boundary conditions
cvb     FIXME: The outer boundary conditions are simple exponential decay on four
cvb     FIXME: points! This is correct for Coulombic potential, but NOT AT ALL for
cvb     FIXME: for our cutoff potential. Hence, our cutoff behavior in detail is likely
cvb     FIXME: not what the radial equation really prescribes!

C       determine outermost integration point
c       maximum radius: 15*classical turning point 

        nin=mch+2.7d0/al

        if(nin .gt. mmax-4) nin=mmax-4

        if ( (r_outer.gt.0.).and.(r_outer.le.r(nin)) ) then
          nin = 1 + ( log( r_outer/r(1) ) )/al
        end if

C       initialize wave function values on outermost four points
C       xkap is the decay constant of the wave function
C       FIXME - by shifting xkap into the loop, we would get proper WKB - right?

        xkap=sqrt(sls/r(nin)**2 + 2.0d0*(v(nin)-e)*efm)

        do i=nin,nin+4
              u(i)=exp(-xkap*(r(i)-r(nin)))
              up(i)=-r(i)*al*xkap*u(i)
              upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
        enddo

c       integrate inward, up to classical turning point
c       same principle as outward integration

        do i=nin,mch+1,-1
          u(i-1)=u(i)+aei(up,i)
          up(i-1)=up(i)+aei(upp,i)
          do it=1,2
            upp(i-1)=(al+frp(i-1))*up(i-1)+(cf(i-1)+fr(i-1))*u(i-1)
            up(i-1)=up(i)+aii(upp,i)
            u(i-1)=u(i)+aii(up,i)
          enddo
        enddo

c  match outside wf to inside wf for continuity

        sc=uout/u(mch)
c
        do i=mch,nin
          upp(i)=sc*upp(i)
          up(i)=sc*up(i)
          u(i)=sc*u(i)
        enddo

C       store derivative at turning point, coming from the outside
        upin=up(mch)

C  -------------------------------------------------------------
C  Normalize entire radial function to unity integral
C  -------------------------------------------------------------

        ro=r(1)/sqrt(amesh)
        sn=ro**(2.0d0*gamma+1.0d0)/(2.0d0*gamma+1.0d0)

C       integral, except for the outermost three points    
        do i=1,nin-3
          sn=sn+al*r(i)*u(i)**2
        enddo

C       treat the outermost three points explicitly (trapezoid rule?)
        sn=sn + al*(23.0d0*r(nin-2)*u(nin-2)**2
     &              + 28.0d0*r(nin-1)*u(nin-1)**2
     &              +  9.0d0*r(nin  )*u(nin  )**2)/24.0d0

C       cn becomes the normalisation sum   
        cn=1.0d0/sqrt(sn)

C       normalize the radial function
        do i=1,nin
          upp(i)=cn*upp(i)
          up(i)=cn*up(i)
          u(i)=cn*u(i)
        enddo

C       normalize stored wave function and derivative at classical turning point
        uout=cn*uout
        upout=cn*upout
        upin=cn*upin

C       zero wave function outside the outermost grid point
        do i=nin+1,mmax
          upp(i)=0.d0
          up(i)=0.d0
          u(i)=0.0d0
        enddo

      return
      end
