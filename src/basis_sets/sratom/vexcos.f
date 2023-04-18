c Header:
c***********************************************************************
c output exchange-correlation energy & potential
c spherical density
c spin-unpolarized
c 
c c input
c    iexc        xc scheme
c    nmx         max. index to evaluate
c    mmax        grid dimension
c    r()         radial grid
c    rho()       spin density array, rho(i,1) == spin up
c    rhop()      1st derivative of spin density w.r.t. r
c    rhopp()     2st derivative of spin density w.r.t. r
c    flag        .true. : evaluate xc-energy
c
c output
c    vs()        xc potential, vs(i,1) == spin up
c    evxc        xc potential energy
c    ex          exchange energy
c    ec          correlation energy
c
c ixc=  3      lsda          xc: perdew-wang 1991
c ixc = 4      lsda          xc: perdew-zunger (from Karsten Freihube) <disabled>
c ixc = 5      gga           x: becke c: perdew 1986
c ixc=  6      gga           xc: perdew-wang 1991
c ixc=  9      gga           xc: perdew/burke/ernzerhof
c ixc=  15                   x: lsda c: vwn5
c ixc = 73                   x: lsda c: colle-salvetti                 <disabled>
c ixc = 74                   x: lsda c: colle-salvetti                 <disabled>
c ixc = 75                   x: becke c: colle-salvetti                <disabled>
c ixc = 76                   x: perdew-wang c: colle-salvetti          <disabled>
c ixc = 83                   x: lsda c: lee-yang-parr
c ixc = 84                   x: lsda c: lee-yang-parr
c ixc = 85                   x: becke c: lee-yang-parr
c ixc = 86                   x: perdew-wang c: lee-yang-parr
c ixc = 30     lsda-x        like 3, exchange only
c ixc = 40     lsda-x        like 4,    "      "  
c ixc = 50     gga-x         like 5,    "      "
c ixc = 60     gga-x         like 6,    "      "
c ixc = 90     gga-x         like 9,    "      "
c ixc = 12     RPBE
c ixc = 11     revPBE
c
c tene = .true.  evaluate ex,ec,eexc
c***********************************************************************
      subroutine vexcos(ixc,nmx,mmax,mx,r,rho,rhop,rhopp,
     1   vs,eexc,ex,ec,tene)
c
      implicit real*8 (a-h,o-z)
      integer ixc, nmx, mmax
      logical tene,tloc
c
      include 'parameter.h'
c
      parameter(c0=0.d0,c1=1.d0)
      integer mx
      dimension rho(mx,2),difxc(mx),d(2),dec(mx),dex(mx),dp(2),dpp(2),
     1          pxc(mx),rhop(mx,2),rhopp(mx,2),vs(mx,2),xpot(2),
     1          zet(mx),cpot(2),wksp1(mx),wksp2(mx),r(mmax),excden(mx)
c xc density
      real*8  dexc(mx), dork(3)
c kinetic energy density etc.
      real*8 stat(mx,2), vs1(mx,2), vs2(mx,2), vc_cs(mx)
c no longer possible      common/tkin/stat(mx,2),vs1(mx,2),vs2(mx,2),vc_cs(mx)
      external  fmom,dmelm
      integer ixca
      integer i, j
c
      pi4=16.0d0*atan(1.0d0)
      al=log(r(2)/r(1))
      tloc=.true.
      ixca=ixc
      if(ixc .eq. 73) ixca=30
      if(ixc .eq. 74) ixca=40
      if(ixc .eq. 75) ixca=50
      if(ixc .eq. 76) ixca=60
      if(ixc .eq. 83) ixca=30
      if(ixc .eq. 84) ixca=40
      if(ixc .eq. 85) ixca=50
      if(ixc .eq. 86) ixca=60
c gga's are used
      if(ixca .eq. 6 .or. ixca .eq. 60) tloc=.false.
      if(ixca .eq. 5 .or. ixca .eq. 50) tloc=.false.
      if(ixca .eq. 9 .or. ixca .eq. 90) tloc=.false.
      if(ixca .eq. 11) tloc=.false.
      if(ixca .eq. 12) tloc=.false.
c
c exchange-correlation potential
c
c Perdew's gga (PW91+). J.P. Perdew in "Electronic Structure of Solids 
c '91", eds. P. Zische and H. Eschrig (Akademie Verlag, Berlin, 1991).
c w/out gradients it amounts to Ceperley-Alder
c
      if(ixca .eq. 3 .or. ixca .eq. 30) then
        do i=1,nmx
          do j=1,2
            d(j)=rho(i,j)/pi4
          enddo
          call ggaxrad(1,r(i),d,dp,dpp,xpot,xen)
          vs(i,1)=xpot(1)
          vs(i,2)=xpot(2)
          dex(i)=xen
          dec(i)=c0
          excden(i)=xen
          if(ixca .eq. 3) then
            call ggacrad(1,r(i),d,dp,dpp,cpot,cen)
            vs(i,1)=vs(i,1)+cpot(1)
            vs(i,2)=vs(i,2)+cpot(2)
            vc_cs(i)=cpot(1)
            dec(i)=cen
            excden(i)=excden(i)+cen
          endif
        enddo
c vwn5-lda
      else if(ixca .eq. 15) then
        do i=1,nmx
          do j=1,2
            d(j)=rho(i,j)/pi4
          enddo
          dork(1)=d(1)
          dork(2)=d(2)
          dork(3)=d(1)+d(2)
          call ggaxrad(1,r(i),d,dp,dpp,xpot,xen)
          vs(i,1)=xpot(1)
          vs(i,2)=xpot(2)
          dex(i)=xen
          dec(i)=c0
          excden(i)=xen
          call cepvwn(dork,ec,cpot,1)
          vs(i,1)=vs(i,1)+cpot(1)
          vs(i,2)=vs(i,2)+cpot(2)
          vc_cs(i)=cpot(1)
          dec(i)=ec
          excden(i)=excden(i)+ec
        enddo
c vwn-lda RPA parametrization (as in GAUSSIAN)
      else if(ixca .eq. 16) then
        do i=1,nmx
          do j=1,2
            d(j)=rho(i,j)/pi4
          enddo
          dork(1)=d(1)
          dork(2)=d(2)
          dork(3)=d(1)+d(2)
          call ggaxrad(1,r(i),d,dp,dpp,xpot,xen)
          vs(i,1)=xpot(1)
          vs(i,2)=xpot(2)
          dex(i)=xen
          dec(i)=c0
          excden(i)=xen
          call cepvwn_gauss(dork,ec,cpot,1)
          vs(i,1)=vs(i,1)+cpot(1)
          vs(i,2)=vs(i,2)+cpot(2)
          vc_cs(i)=cpot(1)
          dec(i)=ec
          excden(i)=excden(i)+ec
        enddo
c becke-perdew
      else if(ixca.eq.5 .or. ixca.eq.50) then
        do 150 i=1,nmx
          do 160 j=1,2
            d(j)=rho(i,j)/pi4
            if(tloc) then
              dp(j)=c0
              dpp(j)=c0
            else
              dp(j)=rhop(i,j)/pi4
              dpp(j)=rhopp(i,j)/pi4
            endif
  160     continue
          call ggaxrad(2,r(i),d,dp,dpp,xpot,xen)
          if(ixca.eq.50) then
            vs(i,1)=xpot(1)
            vs(i,2)=xpot(2)
            dex(i)=xen
            dec(i)=c0
            excden(i)=xen
          else
            call ggacrad(3,r(i),d,dp,dpp,cpot,cen)
            vs(i,1)=xpot(1)+cpot(1)
            vs(i,2)=xpot(2)+cpot(2)
            vc_cs(i)=cpot(1)
            dex(i)=xen
            dec(i)=cen
            excden(i)=cen+xen
          endif
  150   continue
c
c perdew-wang 1991 gga
      else if(ixca.eq.6 .or. ixca.eq.60) then
        do i=1,nmx
          do j=1,2
            d(j)=rho(i,j)/pi4
            if(tloc) then
              dp(j)=c0
              dpp(j)=c0
            else
              dp(j)=rhop(i,j)/pi4
              dpp(j)=rhopp(i,j)/pi4
            endif
          enddo
          call ggaxrad(3,r(i),d,dp,dpp,xpot,xen)
          if(ixca.eq.60) then
            vs(i,1)=xpot(1)
            vs(i,2)=xpot(2)
            dex(i)=xen
            dec(i)=c0
            excden(i)=xen
          else
            call ggacrad(2,r(i),d,dp,dpp,cpot,cen)
            vs(i,1)=xpot(1)+cpot(1)
            vs(i,2)=xpot(2)+cpot(2)
            vc_cs(i)=cpot(1)
            dex(i)=xen
            dec(i)=cen
            excden(i)=cen+xen
          endif
        enddo
c
c perdew/burke/ernzerhof pbe 1996 gga
      else if(ixca .eq. 9 .or. ixca .eq. 90) then
        do i=1,nmx
          do j=1,2
            d(j)=rho(i,j)/pi4
            if(tloc) then
              dp(j)=c0
              dpp(j)=c0
            else
              dp(j)=rhop(i,j)/pi4
              dpp(j)=rhopp(i,j)/pi4
            endif
          enddo
          call ggaxrad(5,r(i),d,dp,dpp,xpot,xen)
          if(ixca.eq.60) then
            vs(i,1)=xpot(1)
            vs(i,2)=xpot(2)
            dex(i)=xen
            dec(i)=c0
            excden(i)=xen
          else
            call ggacrad(5,r(i),d,dp,dpp,cpot,cen)
            vs(i,1)=xpot(1)+cpot(1)
            vs(i,2)=xpot(2)+cpot(2)
            vc_cs(i)=cpot(1)
            dex(i)=xen
            dec(i)=cen
            excden(i)=cen+xen
          endif
        enddo
c@@edu>
c PBEint
      else if(ixca .eq. 21) then
        do i=1,nmx
          do j=1,2
            d(j)=rho(i,j)/pi4
            if(tloc) then
              dp(j)=c0
              dpp(j)=c0
            else
              dp(j)=rhop(i,j)/pi4
              dpp(j)=rhopp(i,j)/pi4
            endif
          enddo
          call ggaxrad(8,r(i),d,dp,dpp,xpot,xen)
          if(ixca.eq.60) then
            vs(i,1)=xpot(1)
            vs(i,2)=xpot(2)
            dex(i)=xen
            dec(i)=c0
            excden(i)=xen
          else
            call ggacrad(6,r(i),d,dp,dpp,cpot,cen)
            vs(i,1)=xpot(1)+cpot(1)
            vs(i,2)=xpot(2)+cpot(2)
            vc_cs(i)=cpot(1)
            dex(i)=xen
            dec(i)=cen
            excden(i)=cen+xen
          endif
        enddo
c@@edu<
c Adding RPBE and revPBE - MR(2009)
c First Hammer-Norskov RPBE
      else if(ixca .eq. 11) then
        do i=1,nmx
          do j=1,2
            d(j)=rho(i,j)/pi4
            if(tloc) then
              dp(j)=c0
              dpp(j)=c0
            else
              dp(j)=rhop(i,j)/pi4
              dpp(j)=rhopp(i,j)/pi4
            endif
          enddo
          call ggaxrad(6,r(i),d,dp,dpp,xpot,xen)
          if(ixca.eq.60) then
            vs(i,1)=xpot(1)
            vs(i,2)=xpot(2)
            dex(i)=xen
            dec(i)=c0
            excden(i)=xen
          else
            call ggacrad(5,r(i),d,dp,dpp,cpot,cen)
            vs(i,1)=xpot(1)+cpot(1)
            vs(i,2)=xpot(2)+cpot(2)
            vc_cs(i)=cpot(1)
            dex(i)=xen
            dec(i)=cen
            excden(i)=cen+xen
          endif
        enddo

c Then revPBE (Zhang-Yang)
      else if(ixca .eq. 12) then
        do i=1,nmx
          do j=1,2
            d(j)=rho(i,j)/pi4
            if(tloc) then
              dp(j)=c0
              dpp(j)=c0
            else
              dp(j)=rhop(i,j)/pi4
              dpp(j)=rhopp(i,j)/pi4
            endif
          enddo
          call ggaxrad(7,r(i),d,dp,dpp,xpot,xen)
          if(ixca.eq.60) then
            vs(i,1)=xpot(1)
            vs(i,2)=xpot(2)
            dex(i)=xen
            dec(i)=c0
            excden(i)=xen
          else
            call ggacrad(5,r(i),d,dp,dpp,cpot,cen)
            vs(i,1)=xpot(1)+cpot(1)
            vs(i,2)=xpot(2)+cpot(2)
            vc_cs(i)=cpot(1)
            dex(i)=xen
            dec(i)=cen
            excden(i)=cen+xen
          endif
        enddo
c
c

      else if(ixca .eq. 4) then
        do 250 i=1,nmx
          do 260 j=1,2
            dup=max(c0,rho(i,1)/pi4)
            ddn=max(c0,rho(i,2)/pi4)
  260     continue
cmf         call vxcfh(dup,ddn,vup,vdn,.true.)
cmf         call exc1fh(dup,ddn,xen,.true.)
          vs(i,1)=vup
          vs(i,2)=vdn
          dex(i)=xen
          dec(i)=c0
          excden(i)=xen
  250   continue
      else if(ixca .eq. 40) then
        do i=1,nmx
          do j=1,2
            dup=max(c0,rho(i,1)/pi4)
            ddn=max(c0,rho(i,2)/pi4)
          enddo
cmf          call vxcfh(dup,ddn,vup,vdn,.false.)
cmf          call exc1fh(dup,ddn,xen,.false.)
          vs(i,1)=vup
          vs(i,2)=vdn
          dex(i)=xen
cmf          call ggacrad(1,r(i),d,dp,dpp,cpot,cen)
          excden(i)=xen
        enddo
      endif
c colle-salvetti correlation
      if(ixc/10 .eq. 7) then
        do i=1,nmx
          do j=1,2
            d(j)=rho(i,j)/pi4
            dp(j)=rhop(i,j)/pi4
            dpp(j)=(2*rhop(i,j)/r(i)+rhopp(i,j))/pi4
          enddo
          statup=stat(i,1)/pi4
          statdn=stat(i,2)/pi4
cmf         call corrcs(.true.,r(i),d(1),d(2),dp(1),dp(2),dpp(1),dpp(2),
cmf     &        statup,statdn,
cmf     &        ec,cpot(1),cpot(2),vs1(i,1),vs1(i,2),vs2(i,1),vs2(i,2))
c implementation in post-LDA manner
c         ixcb=1
c         if(ixca .eq. 50) ixcb=3
c         if(ixca .eq. 60) ixcb=2
c         call ggacrad(ixcb,r(i),d,dp,dpp,cpot,cen)
          vc_cs(i)=cpot(1)
          vs(i,1)=vs(i,1)+cpot(1)
          vs(i,2)=vs(i,2)+cpot(2)
          dec(i)=ec
          excden(i)=excden(i)+ec
        enddo
c LYP correlation
      else if(ixc/10 .eq. 8) then
        do i=1,nmx
          do j=1,2
            d(j)=rho(i,j)/pi4
            dp(j)=rhop(i,j)/pi4
            dpp(j)=(2*rhop(i,j)/r(i)+rhopp(i,j))/pi4
          enddo
c ggacrad with ix=4 is always LYP and accounts for rho=0
          call ggacrad(4,r(i),d,dp,dpp,cpot,cen)
c          call corlyp(.true.,r(i),
c     &         d(1),d(2),dp(1),dp(2),dpp(1),dpp(2),
c     &         ec,cpot(1),cpot(2))
c implementation in post-LDA manner
c         ixcb=1
c         if(ixca .eq. 50) ixcb=3
c         if(ixca .eq. 60) ixcb=2
c         do j=1,2
c           dpp(j)=rhopp(i,j)/pi4
c         enddo
c         call ggacrad(ixcb,r(i),d,dp,dpp,cpot,cen)
          vc_cs(i)=cpot(1)
          vs(i,1)=vs(i,1)+cpot(1)
          vs(i,2)=vs(i,2)+cpot(2)
          dec(i)=ec
          excden(i)=excden(i)+ec
        enddo
      endif
c
      do i=nmx+1,mmax
        dex(i)=c0
        dec(i)=c0
        excden(i)=c0
        dexc(i)=c0
        vs(i,1)=c0
        vs(i,2)=c0
      enddo
c    
c ex   exchange energy 
c ec   correlation energy 
c eexc ex+ec - <n*vxc>
      if(tene) then 
        do i=1,mmax
          wksp1(i)=rho(i,1)+rho(i,2)
          wksp2(i)=rho(i,1)*vs(i,1)+rho(i,2)*vs(i,2)
          dexc(i)=excden(i)
        enddo
        ex=dmelm(mmax,al,r,wksp1,dex)
        ec=dmelm(mmax,al,r,wksp1,dec)
        eexc=ex+ec-fmom(0,mmax,al,c1,r,wksp2)
      endif
c
      return
      end
c
