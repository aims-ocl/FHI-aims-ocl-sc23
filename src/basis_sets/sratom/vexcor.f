c Header:
c***********************************************************************
c output exchange-correlation potential
c 
c input
c    iexc        xc scheme
c    mmax        max. index to evaluate
c    r()         radial grid
c    rho()       density array
c    rhop()      1st derivative of density w.r.t. r
c    rhopp()     2st derivative of density w.r.t. r
c    flag        .true. : evaluate xc-energy 
c
c output 
c    vpxc()      xc potential 
c    evxc        xc potential energy
c    ex          exchange energy
c    ec          correlation energy
c
c Martin Fuchs, FHI der MPG, Berlin, 01-1993
c***********************************************************************
c
      subroutine vexcor(iexc,mmax,r,rho,rhop,rhopp,vpxc,evxc,ex,ec,flag)

ctest
c      use xc
ctest end
c
      implicit none
      include  'parameter.h'

      logical flag
      integer iexc,mmax,i,j,ix
      real*8  al,ecp,evxc,exc,fxc,pi4,rh,rs,cen,xen
      real*8  ec,fx,ex,aln,x,vcuplsd,vcdnlsd,vxup,vcup,dmelm,fex,fvx 
      real*8  r(mmax),d(2),dp(2),dpp(2),xpot(2),cpot(2)
      real*8  rho(mmax),rhop(mmax),rhopp(mmax),vpxc(mmax)
      real*8  dex(mmax),dec(mmax),dork(3)

      real*8  dexc(mmax)

c meta gga
      real*8  dummy_xxx(2),rho_q,rhop_q,tkin_q,ec_q
      real*8  fx_mgga_pk,fc_mgga_pk
      real*8  stat(mmax),stat_core(mmax)
c no longer possible like this      common/use_mgga/stat(mmax),stat_core(mmax)

      external dmelm,fx_mgga_pk,fc_mgga_pk
c
      pi4=16.0d0*atan(1.0d0)
      do i = 1,2
        d(i) = 1.d-10
        dp(i) = 0.d0
        dpp(i) = 0.d0
      enddo
c
c exchange-correlation potential 
c
c x only
      if(iexc.eq.0 .or. iexc.eq.12 .or. iexc.eq.13) then
        do i=1,mmax
          vpxc(i)=0.d0
          dex(i)=0.d0
          dec(i)=0.d0
        enddo
      else if(iexc .eq. 11) then
        do i=1,mmax
          rh=rho(i)/pi4
          call xlda(rh,vpxc(i),dex(i))
          dec(i)=0.d0
        enddo
c          
c C Wigner 
      else if(iexc .eq. 1) then
        do i=1,mmax
          rh=rho(i)/pi4
          call wigner(rh,ex,fx,exc,fxc)
          vpxc(i)=fxc
          dex(i)=ex
          dec(i)=exc-ex
        enddo  
c
c C Hedin-Lundqvist
      else if(iexc .eq. 2) then
        do i=1,mmax
          rh=rho(i)/pi4
          if(rh .ne. 0.0d0) then
            rs=0.62035049d0*rh**(-0.3333333333333333d0)
            x=rs/21.0d0
            aln=dlog(1.0d0 + 1.0d0/x)
            ecp = aln+(x**3*aln-x*x)+x/2-1.0d0/3.0d0
            dex(i)=-0.458175d0/rs - 0.0225d0*ecp
            vpxc(i)=-0.6109d0/rs - 0.0225d0*aln
          else
            dex(i)=0.d0
            vpxc(i)=0.d0
          endif
        enddo  
c
c XC Ceperley - Alder 
c as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)
      else if(iexc .eq. 3) then
        do i=1,mmax
          rh=rho(i)/pi4
          call cepal(rh,ex,fx,exc,fxc)
          vpxc(i)=fxc
          dex(i)=ex
          dec(i)=exc-ex
        enddo
c
c X Perdew, C Perdew, generalized gradient approximation 1992
      else if(iexc .eq. 4) then
        do i=1,mmax
          rh=rho(i)/pi4
          do j=1,2
            d(j)=.5d0*rho(i)/pi4 
            dp(j)=.5d0*rhop(i)/pi4
            dpp(j)=.5d0*rhopp(i)/pi4
          enddo
          call ggaxrad(3,r(i),d,dp,dpp,xpot,dex(i))
          call ggacrad(2,r(i),d,dp,dpp,cpot,dec(i))
C          call ggaxrad(5,r(i),d,dp,dpp,xpot,dex(i))
C          call ggacrad(5,r(i),d,dp,dpp,cpot,dec(i))
          vpxc(i)=xpot(1)+cpot(1)
        enddo
c
c X Becke, C Perdew, generalized gradient approximation
      else if(iexc .eq. 5) then
        do i=1,mmax
          do j=1,2
            d(j)=0.5d0*rho(i)/pi4
            dp(j)=0.5d0*rhop(i)/pi4
            dpp(j)=0.5d0*rhopp(i)/pi4
          enddo
          call ggaxrad(2,r(i),d,dp,dpp,xpot,dex(i))
          call ggacrad(3,r(i),d,dp,dpp,cpot,dec(i))
          vpxc(i)=xpot(1)+cpot(1)
        enddo
c
c XC Ceperley-Alder in Perdew-Wang parametrization of 1991
c + MacDonald-Vosko relativistic correction to exchange
c vosko/wilk/nussair parametrization can be checked against
c http://physics.nist.gov/PhysRefData/DFTdata/
      else if(iexc .eq. 7) then 
        do i=1,mmax
          d(1)=0.5d0*rho(i)/pi4
          d(2)=0.5d0*rho(i)/pi4
          dork(1)=d(1)
          dork(2)=d(2)
          dork(3)=d(1)+d(2)
c      produces LSDA exchange
          call ggaxrad(1,r(i),d,dp,dpp,xpot,dex(i))
c          call ggacrad(1,r(i),d,dp,dpp,cpot,dec(i))
          call cepvwn(dork,ec,cpot,0)
          d(1)=2.d0*d(1)
c      I don't want relativistic corrections here - MR.
c          call relxc(1,d(1),fex,fvx)
          dex(i)=dex(i)
          dec(i)=ec
          vpxc(i)=xpot(1)+cpot(1)
        enddo

c vosko/wilk/nussair parametrization to RPA data
c like it is used in GAUSSIAN
      else if(iexc .eq. 16) then 
        do i=1,mmax
          d(1)=0.5d0*rho(i)/pi4
          d(2)=0.5d0*rho(i)/pi4
          dork(1)=d(1)
          dork(2)=d(2)
          dork(3)=d(1)+d(2)
c      produces LSDA exchange
          call ggaxrad(1,r(i),d,dp,dpp,xpot,dex(i))
c          call ggacrad(1,r(i),d,dp,dpp,cpot,dec(i))
          call cepvwn_gauss(dork,ec,cpot,0)
          d(1)=2.d0*d(1)
c      I don't want relativistic corrections here - MR.
c          call relxc(1,d(1),fex,fvx)
          dex(i)=dex(i)
          dec(i)=ec
          vpxc(i)=xpot(1)+cpot(1)
        enddo
c
c XC Ceperley-Alder in Perdew-Wang parametrization of 1991
      else if(iexc .eq. 8) then

        do i=1,mmax
          d(1)=0.5d0*rho(i)/pi4
          d(2)=0.5d0*rho(i)/pi4
          call ggaxrad(1,r(i),d,dp,dpp,xpot,dex(i))
          call ggacrad(1,r(i),d,dp,dpp,cpot,dec(i))
          vpxc(i)=xpot(1)+cpot(1)
        enddo

c C Lee-Yang-Parr
      else if(iexc .eq. 9 .or. iexc .eq. 10 .or. iexc .eq. 29) then
c X Becke gga, C Lee-Yang-Parr gga
        ix=2
c X Perdew-Wang gga, C Lee-Yang-Parr gga
        if(iexc .eq. 10) ix=3 
        do i=1,mmax
          do j=1,2
            d(j)=0.5d0*rho(i)/pi4
            dp(j)=0.5d0*rhop(i)/pi4
            dpp(j)=0.5d0*rhopp(i)/pi4
          enddo
          call ggaxrad(ix,r(i),d,dp,dpp,xpot,dex(i))
          call ggacrad(4,r(i),d,dp,dpp,cpot,dec(i))
          vpxc(i)=xpot(1)+cpot(1)
        enddo
c
c XC burke-perdew-ernzerhof gga 1996
      else if(iexc .eq. 6) then
        do i=1,mmax
          do j=1,2
            d(j)=0.5d0*rho(i)/pi4
            dp(j)=0.5d0*rhop(i)/pi4
            dpp(j)=0.5d0*rhopp(i)/pi4
          enddo
          call ggaxrad(5,r(i),d,dp,dpp,xpot,dex(i))
          call ggacrad(5,r(i),d,dp,dpp,cpot,dec(i))
C          call ggaxrad(3,r(i),d,dp,dpp,xpot,dex(i))
C          call ggacrad(2,r(i),d,dp,dpp,cpot,dec(i))
          vpxc(i)=xpot(1)+cpot(1)

ctest
c          if ((i.eq.8).or.(i.eq.9)) then
c            write(6,*) "PBE: vpxc(",i,") = ", vpxc(i)
c            write(6,*) "PBE: xpot    = ", xpot(1)
c            write(6,*) "PBE: cpot    = ", cpot(1)
c            write(6,*) "d(1)   = ", d(1)
c            write(6,*) "d(2)   = ", d(2)
c            write(6,*) "dp(1)  = ", dp(1)
c            write(6,*) "dp(2)  = ", dp(2)
c            write(6,*) "dpp(1) = ", dpp(1)
c            write(6,*) "dpp(2) = ", dpp(2)
c            write(6,*) "rhopp(",i," = ", rhopp(i)
c          end if
ctest end
        enddo
c
c XC hammer/norskov pbe
      else if(iexc .eq. 14) then
        do i=1,mmax
          do j=1,2
            d(j)=0.5d0*rho(i)/pi4
            dp(j)=0.5d0*rhop(i)/pi4
            dpp(j)=0.5d0*rhopp(i)/pi4
          enddo
          call ggaxrad(6,r(i),d,dp,dpp,xpot,dex(i))
          call ggacrad(5,r(i),d,dp,dpp,cpot,dec(i))
          vpxc(i)=xpot(1)+cpot(1)
        enddo
c
c XC zhang/wang pbe
      else if(iexc .eq. 15) then
        do i=1,mmax
          do j=1,2
            d(j)=0.5d0*rho(i)/pi4
            dp(j)=0.5d0*rhop(i)/pi4
            dpp(j)=0.5d0*rhopp(i)/pi4
          enddo
          call ggaxrad(7,r(i),d,dp,dpp,xpot,dex(i))
          call ggacrad(5,r(i),d,dp,dpp,cpot,dec(i))
          vpxc(i)=xpot(1)+cpot(1)
        enddo
c
c@@edu>
      else if(iexc .eq. 21) then
        do i=1,mmax
          do j=1,2
            d(j)=0.5d0*rho(i)/pi4
            dp(j)=0.5d0*rhop(i)/pi4
            dpp(j)=0.5d0*rhopp(i)/pi4
          enddo
          call ggaxrad(8,r(i),d,dp,dpp,xpot,dex(i))
          call ggacrad(6,r(i),d,dp,dpp,cpot,dec(i))
          vpxc(i)=xpot(1)+cpot(1)
        enddo
c@@edu<

c meta gga perdew kurth PRL 92 (12) 1999
      else if(iexc .eq. 16) then
        do i=1,mmax
          do j=1,2
            d(j)=0.5d0*rho(i)/pi4
            dp(j)=0.5d0*rhop(i)/pi4
            dpp(j)=0.5d0*rhopp(i)/pi4
          enddo
          call ggaxrad(5,r(i),d,dp,dpp,xpot,dummy_xxx(i))   ! pbe gga for x potential

          call ggaxrad(1,r(i),d,dp,dpp,dummy_xxx,dex(i))    ! mgga for x energy
          rho_q  = rho(i)/pi4
          rhop_q = rhop(i)/pi4
          tkin_q = stat(i)/pi4
          dex(i) = dex(i)
     1      *fx_mgga_pk(rho_q,rhop_q,tkin_q)

          call ggacrad(5,r(i),d,dp,dpp,cpot,ec)       ! pbe gga for c potential
          d(2)   = 1.d-12                             ! mgga c enhancement factor
          dp(2)  = 1.d-12
          dpp(2) = 1.d-12
          call ggacrad(5,r(i),d,dp,dpp,dummy_xxx,ec_q)
          tkin_q = 0.5d0*tkin_q
          dec(i) = ec
     1       *fc_mgga_pk(d(1),dp(1),tkin_q,d(1),dp(1),tkin_q
     1                              ,ec,ec_q,ec_q)

          vpxc(i)=xpot(1)+cpot(1)
        enddo
c
c     PBE GGA exchange + LDA correlation
      else if(iexc .eq. 17) then

        do i=1,mmax
          do j=1,2
            d(j)=0.5d0*rho(i)/pi4
            dp(j)=0.5d0*rhop(i)/pi4
            dpp(j)=0.5d0*rhopp(i)/pi4
          enddo
          call ggaxrad(5,r(i),d,dp,dpp,xpot,dex(i))
          call ggacrad(5,r(i),d,dp,dpp,cpot,dec(i))
          vpxc(i)=xpot(1)+cpot(1)
        enddo
c
c     LDA correlation
      else if(iexc .eq. 18) then

        do i=1,mmax
          d(1)=0.5d0*rho(i)/pi4
          d(2)=0.5d0*rho(i)/pi4
          call ggacrad(1,r(i),d,dp,dpp,cpot,dec(i))
          vpxc(i)=cpot(1)
        enddo
c
c     PBE GGA correlation
      else if(iexc .eq. 19) then

        do i=1,mmax
          do j=1,2
            d(j)=0.5d0*rho(i)/pi4
            dp(j)=0.5d0*rhop(i)/pi4
            dpp(j)=0.5d0*rhopp(i)/pi4
          enddo
          call ggacrad(5,r(i),d,dp,dpp,cpot,dec(i))
          vpxc(i)=cpot(1)
        enddo
c    AM05 xc
c      else if(iexc .eq. 20) then
c
c        do i=1,mmax
c          
c          call am05(rho(i), s(i), u(i),t(i),ex(i)
c     +                ,ec(i),vx(i),vc(i),1)
c          vpxc(i)= vx(i)+vc(i)
c        enddo
c      else
c 
        write(ie,*) '& vexcor - stop: xc option not implemented',iexc
        stop

      endif
c    
c total energy components
      if(flag) then
        do i=1,mmax
          dexc(i)=dex(i)+dec(i)
        enddo
        al=log(r(2)/r(1))
        ex=dmelm(mmax,al,r,rho,dex)
        ec=dmelm(mmax,al,r,rho,dec)
        evxc=dmelm(mmax,al,r,rho,vpxc)
      endif

ctest
c        write(6,*)
c        write(6,*) "Total atomic ex: ", ex
c        write(6,*) "Total atomic ec: ", ec
c        write(6,*) "Total atomic evxc: ", evxc
ctest end 

      return
      end
c
