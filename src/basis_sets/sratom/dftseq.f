c  Header:
c***********************************************************************
c
c integrates radial (pauli-type scalar-relativistic w/ spin-orbit av.) 
c equation on a c logarithmic mesh w/ fixed/variable effective mass 
c
c input
c mode ...... 1 is for full potential bound state (scalar relativistic)
c             2 is for pseudopotential bound state (nonrelativistic)
c             3 is for full potential to find log derivative (sc.rel.)
c             4 is for pseudopotential to find energy which produces
c               specified log derivative (fake bound state)
c             5 is for pseudopotential to produce wavefunction beyond
c               radius used for pseudopotential construction
c             6 like 3 but nonrelativistic
c             7 like 4 
c             8 like 2 but not terminating if no bound state is found
c               returning zero eigenvalue instead
c             9 like 1 with ZORA-scalar relativistic 
c            10 like 1 with eigenvalue corrected ZORA-scalar relativistic 
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
c vb variables
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
      subroutine dftseq(mode,z,mmax,r,n,l,efm,v,v1,v2,
     &           nin,mch,uld,e,u,up,upp, r_outer)
c
      implicit real*8 (a-h,o-z)
      integer mode
      integer mmax
      integer n,l
      integer nin, mch
c
      include './parameter.h'
      include './default.h'

      integer itmax
      parameter (itmax=200)
      parameter (fss_i=137.0359896d0)
c     upper bound for eigenenergy of current radial function: 
c     Any eigenenergy beyond 1d3 Hartrees can not lead to an acceptable 
C     basis function, can it?
      parameter (e_bound=1d3)
c
      logical trap,trel_KH,trel_Z,trel_own
      common/err/trap(10)
c
      dimension r(mmax),u(mmax),up(mmax),upp(mmax)
      dimension cf(mmax),dv(mmax),fr(mmax),frp(mmax),v(mmax),
     +          v1(mmax),v2(mmax),
     1          x1(mmax),x2(mmax),x4(mmax),erec(itmax+1)

cvb
      real*8 r_outer
      integer i_outer
      integer i, i_grid, it
      integer ico, ics

      integer imod, nint, node

      real*8 up_lo, up_hi

c  begin work

ctest
c      write(6,*) "dftseq: n = ", n
c      write(6,*) "dftseq: l = ", l
c      write(6,*) "dftseq: e = ", e
ctest end

c convergence factor for solution of schroedinger eq.  if calculated
c correction to eigenvalue is smaller in magnitude than epsae times
c the magnitude of the current guess, the current guess is not changed.
c
c     epsae=1.0d-9
c
c relativistic - non-relativistic switch
c
! this turns off the zora.
!      mode = 2

      amesh=r(2)/r(1)
      al=log(amesh)
      sls=l*(l+1)

      if(mode .eq. 6 .or. mode .eq. 7) then
        ics=max(1,mch)
        ico=max(1,mch)
      else
        ics=mmax
        ico=mch
      endif

      if((mode.eq.1).or.(mode.eq.3).or.(mode.eq.9).or.(mode.eq.10)) then
        if(mode .eq. 9) then
          trel_Z=.true.
          trel_KH=.false.
          trel_own = .false.
        else if (mode.eq.10) then
          trel_Z=.false.
          trel_KH=.false.
          trel_own = .true.
        else
          trel_Z=.false.
          trel_KH=.true.
          trel_own = .false.
        endif  
        fss=(1.0d0/fss_i)**2
        if(l .eq. 0) gamma=sqrt(1.0d0-fss*z**2)
        if(l .gt. 0) gamma=(l*sqrt(l**2-fss*z**2) +
     &   (l+1)*sqrt((l+1)**2-fss*z**2))/(2*l+1)
      else
        trel_KH=.false.
        trel_Z =.false.
        trel_own =.false.
        fss=1.0d-20
        gamma=l+1
      end if
      imod=mode
      if(imod .eq. 6) imod = 3
      if(imod .eq. 7) imod = 4
      if(imod .eq. 8) imod = 2
c     if(mode .ne. 1 .and. mode .ne. 2) then
c       write(ie,*) '& dftseq- mode,imod,fss,gamma',mode,imod,fss,gamma
c     endif
c   
      if((imod.eq.1).or.(imod.eq.2).or.(imod.eq.9).or.(imod.eq.10)) then
        emax=min( v(mmax)+0.5d0*sls/r(mmax)**2, e_bound )
        if (r_outer.eq.0.) then
cvb       screened coulomb potential, always negative
          emin=0.0d0
        else
cvb       we have a cutoff potential so emin is the minimum radial potential 
cvb      (effective pot'l + angular momentum barrier)
          i_outer = int(1d0 + ( log( r_outer/r(1) ) )/al)
          emin = v(i_outer)+0.5d0*sls/r(i_outer)**2
        end if
        do 6 i=1,mmax
    6     emin=dmin1(emin,v(i)+0.5d0*sls/r(i)**2)
ctest
c        write(6,*) "emin = ", emin
c        write(6,*) "emax = ", emax
ctest end
        if(e .gt. emax) e=1.25d0*emax
        if((e .lt. emin).and.(emin.le.0)) then
          e=0.75d0*emin
        else if((e .lt. emin).and.(emin.gt.0)) then
          e=1.25d0*emin
        end if
        if(e .gt. emax) e=0.5d0*(emax+emin)
      else if(imod .eq. 4) then
        emax=e +  10.0d0
        emin=e - 10.0d0
      end if
c
c initialize arrays and remove leftover garbage
      do i=1,4
        u(i)=0.d0
        up(i)=0.0d0
        upp(i)=0.0d0
      enddo
      do i=1,mmax
        x1(i)=al*r(i)
        x2(i)=x1(i)*x1(i)
      enddo

      if (trel_KH .or. trel_Z .or. trel_own ) then
c calculate dv/dr for darwin correction
        if (ics.lt.5) then
          write(6,*) "* Integration of radial basis functions: "
          write(6,*) "* Need more than four logarithmic grid points."
          stop
        end if
        ! see Eq (17) at http://mathworld.wolfram.com/FiniteDifference.html

        dv(1)=(-50.d0*v(1)+96.d0*v(2)-72.d0*v(3)+32.d0*v(4)
     &         -6.d0*v(5))/(24.d0*x1(1))
        dv(2)=(-6.d0*v(1)-20.d0*v(2)+36.d0*v(3)-12.d0*v(4)
     &         +2.d0*v(5))/(24.d0*x1(2))

        do i=3,ics-2
          dv(i)=(2.d0*v(i-2)-16.d0*v(i-1)+16.d0*v(i+1)
     &           -2.d0*v(i+2))/(24.d0*x1(i))
        enddo

        dv(ics-1) = (+6.d0*v(ics)+20.d0*v(ics-1)-36.d0*v(ics-2)
     &         +12.d0*v(ics-3)
     &         -2.d0*v(ics-4))/(24.d0*x1(2))

        dv(ics) = (50.d0*v(ics)-96.d0*v(ics-1)+72.d0*v(ics-2)
     &         -32.d0*v(ics-3)
     &         +6.d0*v(ics-4))/(24.d0*x1(ics))

      else
c nonrelativistic coefficient arrays for u (fr) and up (frp)
        do i=1,ics
          fr(i)=-2*x2(i)*v1(i)*v2(i)
          frp(i)=2*x1(i)*v1(i)*v2(i)
        enddo
      endif

c return point for bound state convergence
      aux=al*al*sls
      nint=0
   10 nint=nint+1
ctest
c      write(6,'(A,3I5,1X,E16.8)') 
c     +  "dftseq: iteration, n, l, e ", nint, n, l, e
ctest end
      erec(nint)=e
        if(nint .gt. itmax) then
          if(mode .eq. 8) then
            write(ie,*) 
     1        '& dftseq - no bound state found (iter): e=>0 n l',n,l
            e=0.d0
            return
          endif
          write(ie,*) '& subroutine dftseq - stop: convergence error'
          write(ie,*) '&  maximum number of iterations exceeded.'
          write(ie,*) '&  nint n l                ',nint,n,l
          write(ie,*) '&  eigenvalue log written to file err.dftseq_e'
          write(ie,*) '&  potential  log written to file err.dftseq_v'
          write(ie,*) '*  '
          write(ie,*) '*  This error may happen if a free atom state'
          write(ie,*) '*  is unbound in a given functional.'
          write(ie,*) '*  '
          write(ie,*) '*  A possible remedy is to define a faraway'
          write(ie,*) '*  bounding potential using the cut_free_atom'
          write(ie,*) '*  keyword, e.g.,'
          write(ie,*) '*  '
          write(ie,*) '*    cut_free_atom finite 20.'
          write(ie,*) '*  '
          write(ie,*) '*  in the species definition of the offending'
          write(ie,*) '*  element.'
          write(ie,*) '*  '
          open(10,file='err.dftseq_e',status='unknown')
          do i=1,10
            write(10,*) i,erec(i)
          enddo
          do i=nint,1,-1
            write(10,*) i,erec(i)
          enddo
          close(10)
          open(10,file='err.dftseq_v',status='unknown')
          do i=1,ics
            write(10,*) r(i),v(i)
          enddo
          close(10)
          trap(4)=.true.
          stop
        endif

c coefficient array for u in differential eq.
        do i=1,ics
          cf(i)=aux+2*(v(i)-e)*v2(i)*x2(i)
        enddo

c KH relativistic coefficient arrays for u (fr) and up (frp)
        if(trel_KH) then
          do i=1,ics
            fr(i)=x2(i)*fss*(-(v(i)-e)**2 + 0.5d0*dv(i)/
     &      (r(i)*(1.0d0+0.5d0*fss*(e-v(i)))))
            frp(i)=-x1(i)*0.5d0*fss*dv(i)/(1.0d0+0.5d0*fss* (e-v(i)))
          enddo
        endif

c ZORA relativistic coefficient arrays for u (fr) and up (frp)
        if(trel_Z) then
          do i=1,ics
            fr(i)=x2(i)*(-(v(i)-e)*fss*v(i) + dv(i)/
     &      (r(i)*(2.0d0/fss-v(i))))
            frp(i)=-x1(i)/(2.0d0/fss-v(i))*dv(i)
          enddo
        endif

c "own" relativistic coefficient arrays for u (fr) and up (frp)
c
c as programmed here, it's _exactly_ Koelling-Harmon.
c
        if(trel_own) then
          do i=1,ics
!
!  This is the Dirac eqn rewritten, with the eigenvalue left in the denominator,
!  and with the classical kinetic energy term already subtracted.
!  Entertainingly, it's exactly Koelling-Harmon as programmed by Martin Fuchs? above.
!
            fr(i)=x2(i)*(-(v(i)-e)*fss*(v(i)-e) + dv(i)/
     &      (r(i)*(2.0d0/fss+e-v(i))))
            frp(i)=-x1(i)/(2.0d0/fss+e-v(i))*dv(i)
          enddo
        endif

c
c find classical turning point for matching
        if((imod.eq.1).or.(imod.eq.2).or.(imod.eq.9).or.(imod.eq.10))
     +     then
          do 30 i=ics,2,-1
            if(cf(i-1) .le. 0.d0 .and. cf(i) .gt. 0.d0) then
              mch=i
              ico=mch
              go to 40
            end if
   30     continue
c
c error trap
          if(mode .eq. 8) then
            write(ie,*) 
     1        '& dftseq - no bound state found (ctp): e=>0 n l',n,l
            e=0.d0
            return
          endif
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
        else
          nin=mch
        end if


c start wavefunction with series
        do 50 i=1,4
          u(i)=r(i)**gamma
          up(i)=al*gamma*r(i)**gamma
   50     upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)


c
cmf l=0: hydrogen-like start (improves eigenvalues and wavefunctions)
c
        if((l.eq.0).and.(imod.eq.2).or.(imod.eq.1)
     &     .or.(imod.eq.9).or.(imod.eq.10)) then
cmf          z=z*efm
          do i=1,4
            u(i)=r(i)**gamma*(1.d0-z*r(i)+0.5d0*(z*r(i))**2)
            up(i)=al*gamma*r(i)**gamma*((1.d0-z*r(i)+0.5d0*(z*r(i))**2))
     &          +al*r(i)**gamma*(-z*r(i)+(z*r(i))**2)
            ! (Rundong)finite nucleus test:
           !call fn_dftseq_u(r(i),z,fss,l,al,u(i))
           !call fn_dftseq_up(r(i),z,fss,l,al,up(i))
            upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
          enddo
cmf         z=z/efm
        endif
c

ctest
c        write(6,*) "After inner initialization."
ctest end

c outward integration using predictor once, corrector
c twice
        node=0
c
        do 70 i=4,ico-1
          u(i+1)=u(i)+aeo(up,i,mmax)
          up(i+1)=up(i)+aeo(upp,i,mmax)
          do 60 it=1,2
            upp(i+1)=(al+frp(i+1))*up(i+1)+(cf(i+1)+fr(i+1))*u(i+1)
            up(i+1)=up(i)+aio(upp,i,mmax)
   60     u(i+1)=u(i)+aio(up,i,mmax)
          if(u(i+1)*u(i) .le. 0.0d0) node=node+1
ctest
c          write(6,*) i+1, u(i+1)
ctest end          
          if (dabs(u(i+1)).gt.1d+50) then
            go to 71
          end if
   70   continue
c
        uout=u(mch)
        upout=up(mch)

        go to 72

ctest
c        write(6,*) "After outward integration: uout = ", uout
c        write(6,*) "After outward integration: node = ", node
ctest end

 71     continue

        if (node.le.n+l-1) then
          write(6,'(1X,A,A)') 
     +      "* dftseq: Diverging solution u(r), but not too ",
     +      "many nodes. I do no know what to do."
          stop
        end if

 72      continue
c
        if(node-n+l+1 .eq. 0 .or. imod .eq. 3 .or. imod .eq. 5) then
c

ctest
c        write(6,*) "Correct number of nodes."
ctest end

          if((imod.eq.1).or.(imod.eq.2).or.(imod.eq.9).or.(imod.eq.10)) 
     +      then

cvb <----------------- begin setting outer boundary conditions

c             no cutoff:
c             start inward integration at 15*classical turning
c             point with simple exponential

ctest
c        write(6,*) "Initialize outer boundary cond."
ctest end

            nin=int(mch*1d0 + 2.7d0/al)

            if(nin .gt. mmax-4) nin=mmax-4

            if ( (r_outer.gt.0.).and.(r_outer.le.r(nin)) ) then
              nin = int(1d0 + ( log( r_outer/r(1) ) )/al)
            end if

            if (r_outer.eq.0.) then
C             no confining potential: initialize for simple exponential decay only - 
C             H-atom style

              xkap=sqrt(sls/r(nin)**2 + 2.0d0*(v(nin)-e)*efm)

              do 110 i=nin,nin+4
                u(i)=exp(-xkap*(r(i)-r(nin)))
                up(i)=-r(i)*al*xkap*u(i)
  110           upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)

            else
cvb           confining potential: initialize wave function WKB-like
cvb           do my own lumberjack-style initialization and pray

c             u(i) goes first (WKB!)
              u(nin-2) = 1.d0
              do i_grid=nin-1,nin+4
                xkap=sqrt(sls/r(i_grid)**2 + 2.0d0*(v(i_grid)-e)*efm)
                u(i_grid)=u(i_grid-1) *
     +            exp(-xkap*(r(i_grid)-r(i_grid-1)))
              enddo

C             up(i) is next 
              do i_grid=nin-1,nin+3
                up_lo = (u(i_grid)-u(i_grid-1))
     +                / (r(i_grid)-r(i_grid-1))
                up_hi = (u(i_grid+1)-u(i_grid))
     +                / (r(i_grid+1)-r(i_grid))
                up(i_grid) = 0.5d0*(up_lo + up_hi)
              enddo
              up(nin+4) = (u(nin+4)-u(nin+3)) / (r(nin+4)-r(nin+3))

C             upp(i) is next 
              do i_grid=nin,nin+3
                up_lo = (up(i_grid)-up(i_grid-1))
     +                / (r(i_grid)-r(i_grid-1))
                up_hi = (up(i_grid+1)-up(i_grid))
     +                / (r(i_grid+1)-r(i_grid))
                upp(i_grid) = 0.5d0*(up_lo + up_hi)
              enddo
              upp(nin+4) = (up(nin+4)-up(nin+3)) / (r(nin+4)-r(nin+3))

            end if

cvb <---------------- end setting outer boundary conditions

ctest
c        write(6,*) "Begin inward integration."
ctest end

c
c integrate inward
c
            do 130 i=nin,mch+1,-1
              u(i-1)=u(i)+aei(up,i,mmax)
              up(i-1)=up(i)+aei(upp,i,mmax)
              do 130 it=1,2
                upp(i-1)=(al+frp(i-1))*up(i-1)+(cf(i-1)+fr(i-1))*u(i-1)
                up(i-1)=up(i)+aii(upp,i,mmax)
                u(i-1)=u(i)+aii(up,i,mmax)
ctest
c                write(6,*) i-1, u(i-1)
ctest end
 130        continue
c
c scale outside wf for continuity

ctest
c            write(6,*) "Matching radius: u(mch) = ", u(mch)
ctest end

            sc=uout/u(mch)
c
            do 150 i=mch,nin
              upp(i)=sc*upp(i)
              up(i)=sc*up(i)
  150         u(i)=sc*u(i)
c
            upin=up(mch)
c
          else
c
            upin=uld*uout
c
          end if
c
c perform normalization sum
c
          ro=r(1)/sqrt(amesh)
          sn=ro**(2.0d0*gamma+1.0d0)/(2.0d0*gamma+1.0d0)
c
          do 160 i=1,nin-3
  160       sn=sn+al*r(i)*u(i)**2
c
          sn=sn + al*(23.0d0*r(nin-2)*u(nin-2)**2
     &              + 28.0d0*r(nin-1)*u(nin-1)**2
     &              +  9.0d0*r(nin  )*u(nin  )**2)/24.0d0
c

ctest
c          write(6,*) "Normalization sum: ", sn
ctest end
          cn=1.0d0/sqrt(sn)
          uout=cn*uout
          upout=cn*upout
          upin=cn*upin
c
          do 180 i=1,nin
            upp(i)=cn*upp(i)
            up(i)=cn*up(i)
  180       u(i)=cn*u(i)
          do 190 i=nin+1,mmax
            upp(i)=0.d0
            up(i)=0.d0
  190       u(i)=0.0d0
c
c exit for fixed-energy calculation
c
          if(imod .eq. 3 .or. imod .eq. 5) return
c
c perturbation theory for energy shift
          de=0.5d0*uout*(upout-upin)/(al*r(mch))

c
c convergence test and possible exit
c
          if(abs(de) .ge. max(abs(e),0.2d0)*epsae) then
c
ctest
c            write(6,*) "  dftseq, iter. ", nint, ", de = ", de
ctest end
            if(de .gt. 0.0d0) then 
              emin=e
            else
              emax=e
            end if
            e=e+de
            if(e .gt. emax .or. e .lt. emin) e=0.5d0*(emax+emin)
c
c loop back to converge e
c
            go to 10

          else

ctest
c            write(6,*) "  dftseq, iter. ", nint, ", de converged"
ctest end

            return

          endif
c
        else if(node-n+l+1 .lt. 0) then
c too few nodes
ctest
c        write(6,*) "Too few nodes."
ctest end

          emin=e
          e=0.5d0*(emin+emax)
          go to 10
c
        else
c too many nodes
ctest
c        write(6,*) "Too many nodes."
ctest end

          emax=e
          e=0.5d0*(emin+emax)
          go to 10
        end if
      continue
c
c
      end
c
