c  Header:
cvb 23-02-05 interface modified for use with LocalOrb: 
cvb          now passing ves, vxc to the outside world for later use
cmf 27-06-03 added spin polarization for selected xc functionals
cmf 04-09-02 added kli exchange functional
c**********************************************************************
c self-consistent scalar-relativistic spherical all-electron atom 
c 1st & 2nd radial derivatives 
c frozen-core
c non-spinpolarized
c self-interaction correction
c
c Martin Fuchs, FHI der MPG, Berlin, 11-1996
c***********************************************************************
      subroutine sratom(
     1    it,irl,igr,iexc,tcore,z,nc,nv,n,l,m,f,e,ninu,de,
     1    mx,mev,ms,
     1    mmax,r,cut_atom,cut_type, r_cut,cut_scale,cut_width,
     1    rpk,vee,ves,vxc,d,dp,dpp,dcs,dcsp,dcspp,uu,
     1    tinpot,tsic,esic,vsic,tkli,svkli,svm,
     1    tmgga,uup,core_type)
c
      implicit none
      include  'default.h'
      include  'parameter.h'

cvb variable dimensions are now handed down:
c  mx : should be equal to n_max_grid, maximum number of logarithmic grid points
c  mev: maximum number of shells in free atom
c  ms: maximum number of different "channels" for effective potentials 
c      (e.g. spin channels if spin-polarised dft, orbitals if
c       if orbital-dependent effective potential, etc.)

cvb a bunch of comments:
c  vee(i_grid,i_channels) : effective potential tabulated on grid, for
c           different "channels" (LDA, no spin: only one channel;
c           with spin: two channels (spin up, spin down); 
C           HF: max_shells channels
C  vxc(i_grid, i_channels) : XC potential
C  ves(i_grid) : electrostatic potential

      logical  tconv,tcore,tkli,tspin
      integer  i,iexc,it,j,mch,mmax,nc,nv,nin,igr,irl,itmax
      integer  mx,mev,ms
      integer  ir,in,ispin,nspinmx,iexcs,mup
      integer  ninu(mev),n(mev),l(mev),m(mev),
     +         nstart(ms),ncore(ms),nend(ms)
      integer  ninmx(ms)
      real*8   al,als,amesh,bl,z,sf,eeel,eexc,ex,ec,et,uld,de,dmelm
      real*8   eeig,eeig_old,feps
      real*8   e(mev),f(mev),rpk(mev),r(mx),rho(mx),rhop(mx),
     +         rhopp(mx)
      real*8   dc(mx),dcp(mx),dcpp(mx),vi(mx)
      real*8   u(mx),up(mx),upp(mx),uu(mx,mev)
      real*8   ves(mx),wa(mx),wb(mx)
      real*8   fmom,sfval,eeel_val,ves_val(mx),rho_val(mx)
      real*8   d(mx,ms),dp(mx,ms),dpp(mx,ms)
      real*8   dcs(mx,ms),dcsp(mx,ms),dcspp(mx,ms)
      real*8   rhonow(mx),rhonowp(mx),rhonowpp(mx)
      real*8   vee(mx,ms),vxc(mx,ms)
      real*8   vi1(mx,ms),vn(mx,ms),vo(mx,ms),vo1(mx,ms)
      real*8   r1i(mx),r2i(mx),r3i(mx),r4i(mx)
c kli
      logical  tstart,tinpot
      integer  mp(mev),norb_kli,iukli,iexckli
      real*8   fp(mev),dens(mx),vckli(mx),vsl(mx),vxc_v(mx)
      real*8   vsl_v(mx),svkli(mev),svm(mev),svmdummy(mev)
      parameter (iukli=80)
c sic
      logical  tsic
      real*8   esic(mev),vsic(mx,mev),dorb(mx,2),dnil(mx,2)
      real*8   blcom
      real*8   vxc_spin(mx,2),vglob(mx),vsic_old(mx,ms)

c meta gga 
      logical  tmgga
      real*8   uup(mx,mev)
      external dmelm

c vb
C  r_outer: In subroutine dft_seq, specify outermost radial grid point (for onset of
C           boundary conditions). If this is zero, a Coulomb potential without cutoff
C           is assumed, and Martin Fuchs' default behavior is restored. 

      integer i_grid

cvb: Treatment of atomic cutoff potential
      logical :: cut_atom
      integer :: cut_type
      real*8  :: r_cut
      real*8 :: cut_scale,cut_width
cmf not used   real*8, external :: cutoff_pot
      real*8, external :: cutoff_pot
cmf f90   real*8, dimension(mx,ms) :: v_radial
      real*8 v_radial(mx,ms)
      real*8 r_outer

cvb: for rel_own
      logical :: core_type(mev)
      integer :: rela_treatment

Ctest
c       write(6,*) "In sratom."


c       write(6,*) "mx = ", mx
c       write(6,*) "mev= ", mev
c       write(6,*) "ms = ", ms
c       write(6,*) "it = ", it
c       write(6,*) "irl = ", irl
c       write(6,*) "igr = ", igr
c       write(6,*) "iexc = ", iexc
c       write(6,*) "tcore = ", tcore
c       write(6,*) "z = ", z
c       write(6,*) "nc = ", nc
c       write(6,*) "nv = ", nv
c       write(6,*) "n = ", n
c       write(6,*) "l = ", l
c       write(6,*) "m = ", m
c       write(6,*) "f = ", f
c       write(6,*) "e = ", e
c       write(6,*) "ninu = ", ninu
c       write(6,*) "de = ", de
c       write(6,*) "mmax = ", mmax
c       do i = 1, mmax, 1
c         write(6,*) "r(", i,") = ", r(i)
c       enddo
c       write(6,*) "rpk= ", rpk
c       do i = 1, mmax, 1
c         write(6,*) "vee(", i,") = ", vee(i,1)
c       enddo
c       do i = 1, mmax, 1
c         write(6,*) "d(", i,") = ", d(i,1)
c       enddo
c       write(6,*) "Missing stuff."
c       write(6,*) "tinpot= ", tinpot      
c       write(6,*) "tsic= ", tsic      
c       write(6,*) "esic= ", esic      
c       write(6,*) "vsic= ", vsic      
c       write(6,*) "tkli= ", tkli      
c       write(6,*) "svkli= ", svkli      
c       write(6,*) "svm= ", svm      
c       write(6,*) "tmgga= ", tmgga
Ctest end

Ctest
c       write(6,*) "Incoming: dcspp(8,1) = ", dcspp(8,1)
c       write(6,*) "Incoming: dcspp(9,1) = ", dcspp(9,1)
Ctest end

      feps=1.d-9
      itmax=it
      amesh=r(2)/r(1)
      al=log(amesh)

!     Zero out KS orbitals (very important !, otherwise breaks down) 
      do i=1,mev
        do j=1,mx
          uu(j,i) = 0.d0
        enddo
      enddo


C initialize storage array for mixing potential in previous iterations
cmf f90
      do i=1,ms
        do j=1,mx
          vo1(j,i) = 0.d0
          vi1(j,i) = 0.d0
          v_radial(j,i) = 0.d0
          u(j) = 0.d0
          dens(j) = 0.d0
        enddo
      enddo
      eeig_old = 0.d0

c spin initialize
      tspin=.false.
      nspinmx=1
      if(iexc.lt.0) nspinmx=2
      iexc=abs(iexc)
c kli initialize
      iexckli=iexc
      if(iexc.eq.12) iexckli=0
      do i=1,nc+nv
        mp(i) = n(i)
      enddo
      tstart=.true.
c sic initialize
      do i=1,mmax
        wa(i)=0.d0
        wb(i)=1.d0
        dorb(i,2)=0.d0
        dnil(i,1)=0.d0
        dnil(i,2)=0.d0
        r1i(i)=1.d0/r(i)
        r2i(i)=r1i(i)*r1i(i)
        r3i(i)=r1i(i)*r2i(i)
        r4i(i)=r2i(i)*r2i(i)
      enddo
c initialize frozen core
      if(.not. tcore) then
        do ispin=1,nspinmx
          do j=1,mmax
            dcs(j,ispin)=dcs(j,ispin)/r2i(j)
            dcsp(j,ispin)=dcsp(j,ispin)/r3i(j)
            dcspp(j,ispin)=dcspp(j,ispin)/r4i(j)
          enddo
        enddo
      endif

c
      sf=0.d0
      do ispin=1,nspinmx
        nstart(ispin)=(nc+nv)*(ispin-1)+1
        nend(ispin)=nstart(ispin)+nv+nc-1
        do i=nstart(ispin),nend(ispin)
          sf=sf+f(i)
          ninu(i)=1
        enddo
      enddo
c
c begin, bl is the mxing factor for the potential in the
c self-consistency loop
      bl=0.5d0
      blcom=1.d0-bl
      als=al*al
c
c return point for self-consistency loop
      do ispin=1,nspinmx
        nstart(ispin)=(nc+nv)*(ispin-1)+(nc+1)
        ncore(ispin)=0
        nend(ispin)=nstart(ispin)+nv-1
        if(tcore) then
           nstart(ispin)=(nc+nv)*(ispin-1)+1
           ncore(ispin)=nstart(ispin)+nc-1
           nend(ispin)=nstart(ispin)+nv+nc-1
        endif
cdebug
c       write(ie,*)
c    1   '%sratom - ispin,nstart,ncore,nend=',
c    1   ispin,nstart(ispin),ncore(ispin),nend(ispin) 

      enddo
      do it=1,itmax

ctest
c         if (it.eq.2) then
c           write(6,*) "sratom: Iter. ", it
c           do i_grid = 1, mmax, 1
c             write(6,*) "vee(i_grid,1): ", vee(i_grid,1)
c           enddo
c         end if
ctest end

c
c initialize for new iteration
        tconv=.true.
        mup=0
        do ispin=1,nspinmx
          do j=1,mmax
            rho(j)=0.0d0
            d(j,ispin)=0.0d0
            dp(j,ispin)=0.0d0
            dpp(j,ispin)=0.0d0
          enddo
        enddo
c
cvb - no sic        if(tsic) call dcpv(mx,1,mmax,vi,vglob)

c start spin loop ------------------------------------------------------
        do ispin=1,nspinmx
        ninmx(ispin)=0

        ! before solving the radial equation for a given electronic potential vee,
        ! add a cutoff potential to confine the radial functions, if desired.
          if (cut_atom) then
            ! Add cutoff radial cutoff potential to the potential
            ! that is actually used to construct the radial functions

            do j=1,mmax
              v_radial(j,ispin) = vee(j,ispin) 
     +        + cutoff_pot
     +          ( r(j), cut_type, r_cut, cut_width, cut_scale )
            enddo

            ! set outermost radius for radial function integration
            r_outer = r_cut + cut_width

          else
            ! simply use the electron-(electron+nucleus) potential
            ! as is, no cutoff ...

            do j=1,mmax
              v_radial(j,ispin) = vee(j,ispin) 
            enddo
            r_outer = 0.d0

          end if

ctest
c          write(6,*) "v_radial: ", v_radial
ctest end

        ! begin loop over all necessary radial functions
        do i=nstart(ispin),nend(ispin)

          if(f(i) .gt. feps .or. tconv) then
            et=e(i)
cvb - no sic           if(tsic) call at_dsuv(mx,1,mmax,vglob,vsic(1,i),vee(1,1))

ctest
c      write(6,*) "Before dftseq."
Ctest end
c TEST
c            write(ie,*) 'it state',it,i,n(i),l(i),m(i),e(i)
c            if (it==4) then
c              do j = 1, mx, 1
c                write(6,*) "vxc,dens:", v_radial(j,1), u(j) 
c              enddo 
c            endif
cvb wa is zero, wb is one - hardcoded.

            ! if we're doing the core scalar relativistically and this is a core state
            ! the pass rela_treatment = 10
            if ((irl.eq.10).and.(.not.(core_type(i)))) then
              ! valence state - regular ZORA
              rela_treatment = 9
            else
              rela_treatment = irl
            end if

c           if (.tkli.) then
c           call at_dftseq_kli(rela_treatment,z,mmax,r,n(i),l(i),1.d0,
c    +                  v_radial(1,ispin),
c    1                  wa,wb,nin,mch,uld,e(i),u,up,upp,r_outer) 
c           else 

c           TEST
c       write(6,*) "uld = ", uld
c       do j = 1, mmax, 1
c         write(6,*) "v_radial(", j,") = ", v_radial(j,ispin)
c       enddo
c           TEST


            call dftseq(rela_treatment,z,mmax,r,n(i),l(i),1.d0,
     +                  v_radial(1,ispin),
     1                  wa,wb,nin,mch,uld,e(i),u,up,upp,r_outer)
c           endif  

c debug
c            write(ie,*) 'it state',it,i,n(i),l(i),m(i),e(i)

            call dcpv(mx,1,nin,u,uu(1,i))
            if(e(i) .ne. et) tconv=.false.
            if(tconv) call dcpv(mx,1,nin,up,uup(1,i)) ! meta gga
            if(nin .gt. ninmx(ispin)) ninmx(ispin)=nin
            if(nin .gt. mup) mup=nin
            ninu(i)=nin
          endif

c accumulate charge & radial derivatives, rhop & rhopp are calculated as
c derivatives w.r.t. r, while up & upp are in terms of the logarithmic
c radial variable
c orbital density and update of spin densities
          do j=1,nin
            rhonow(j)=u(j)*u(j)
            d(j,ispin)=d(j,ispin)+f(i)*rhonow(j)
          enddo
          if(igr.eq.1) then
            do j=1,nin
              rhonowp(j)=2.d0*(up(j)/al-u(j))*u(j)
              rhonowpp(j)=2.d0*( (up(j)/al-u(j))*(up(j)/al
     1               -3.d0*u(j)) +(upp(j)-al*up(j))*u(j)/als )
              dp(j,ispin)=dp(j,ispin)+f(i)*rhonowp(j)
              dpp(j,ispin)=dpp(j,ispin)+f(i)*rhonowpp(j)
            enddo
          endif

c store core density
          if(i.eq.ncore(ispin)) then
            do j=1,ninmx(ispin)
              dcs(j,ispin)=d(j,ispin)
              dcsp(j,ispin)=dp(j,ispin)
              dcspp(j,ispin)=dpp(j,ispin)
              d(j,ispin)=0.d0
              dp(j,ispin)=0.d0
              dpp(j,ispin)=0.d0
            enddo
            do j=ninmx(ispin)+1,mmax
              dcs(j,ispin)=0.d0
              dcsp(j,ispin)=0.d0
              dcspp(j,ispin)=0.d0
              d(j,ispin)=0.d0
              dp(j,ispin)=0.d0
              dpp(j,ispin)=0.d0
            enddo
          endif

c
c find outermost peak of wavefunction (needed for peseudopotentials)
          if(tconv .or. it .ge. itmax) then
            do j=nin-1,2,-1
              if(up(j)*up(j+1) .lt. 0.0d0) goto 341
            enddo
  341       rpk(i)=r(j)
          endif
c
        enddo

c end spin loop --------------------------------------------------------
        enddo


Ctest
c        write(6,*) "Before total density: dpp(9,1)   = ", dpp(8,1)
c        write(6,*) "Before total density: dcspp(9,1) = ", dcspp(8,1)
c        write(6,*) "Before total density: dpp(9,1)   = ", dpp(9,1)
c        write(6,*) "Before total density: dcspp(9,1) = ", dcspp(9,1)
Ctest end

c
c make total density rho()
        do ispin=1,nspinmx
          do j=1,mmax
            d(j,ispin)=(d(j,ispin)+dcs(j,ispin))*r2i(j)
            rho(j)=rho(j)+d(j,ispin)
            dp(j,ispin)=(dp(j,ispin)+dcsp(j,ispin))*r3i(j)
            dpp(j,ispin)=(dpp(j,ispin)+dcspp(j,ispin))*r4i(j)
co          rho(j)=rho(j)/r(j)**2 +dc(j)
co          rhop(j)=rhop(j)/r(j)**3 +dcp(j)
co          rhopp(j)=rhopp(j)/r(j)**4 +dcpp(j)
          enddo
        enddo

c
c sic potentials
cvb - no sic        if(tsic) then
cvb - no sic          do i=nstart(ispin),nend(ispin)
cvb - no sic            do j=1,mmax
cvb - no sic              dorb(j,1)=(uu(j,i)/r(j))**2
cvb - no sic            enddo
cvb - no sic            call vestat(mmax,1.d0,eeel,r,dorb(1,1),ves,tconv)
cvb - no sicc only LSDA no GGA here, 3 corresponds to xc option 8 
cvb - no sicc           call vexcos(3,ninu(i),mmax,r,dorb,dnil,dnil,
cvb - no sicc    1        vxc_spin,eexc,ex,ec,tconv)
cvb - no sic            call dadv(mx,1,mmax,ves,vxc_spin(1,1),vsic(1,i))
cvb - no sic            if(it .gt. 1) then
cvb - no sic              do j=1,mmax
cvb - no sic                vsic(j,i)=bl*vsic(j,i)+blcom*vsic_old(j,i)
cvb - no sic              enddo
cvb - no sic            endif
cvb - no sic            call dcpv(mx,1,mmax,vsic(1,i),vsic_old(1,i))
cvb - no sic            esic(i)=-(0.5*eeel+ex+ec)
cvb - no sic          enddo
cvb - no sic        endif
c
c global effective potential
        call vestat(mmax,sf-z,eeel,r,rho,ves,.false.)
        if(.not. tinpot .and. tkli) then
          if(it < 5) then
!           write(ie,*) 'Calling vexcor ...' 
           call vexcor(8,mmax,r,d(1,1),dp(1,1),dpp(1,1),
     1                  vxc(1,1),eexc,ex,ec,.false.)
!           write(ie,*) 'Finished vexcor ...'

          else
            
            do in=1,nc+nv
              fp(in)=0.5d0*f(in)
            enddo
            do ir=1,mmax
              dens(ir)=0.5d0*rho(ir)
            enddo

            if(it == 5) then
             norb_kli=0
             do i=1,nc+nv
              if(fp(i) .gt. feps) norb_kli=norb_kli+1     
             enddo
            endif

!           write(ie,*) 'Calling at_vklix. norb =', norb_kli
!            if(.false.) then
!            call vklix(norb_kli,n,l,mp,fp,mmax,r,uu
!     +,                 dens,vxc(1,1),vsl,svkli,svm,1,tstart,tconv
!     +,                 .true.)
!            else
            call at_vklix(norb_kli,n,l,mp,fp,mmax,r,uu
     +,                 dens,vxc(1,1),vsl,svkli,svm,1,tstart,tconv
     +,                 .false.,mx,ms,mev)
!            endif
c          write(ie,*) 'Finished at_vklix ...'

            if(iexckli .ne. 0) then
co            call vexcor(iexckli,mmax,r,rho,rhop,rhopp,vckli,eexc,
              call vexcor(iexckli,mmax,r,d(1,1),dp(1,1),dpp(1,1),
     1                    vckli,eexc,ex,ec,tconv)
!             write(ie,*) '%sratom_n: ',
!     1        it, 'adding vckli',iexckli,vckli(1),vckli(100)
              do ir=1,mmax
                vxc(ir,1)=vxc(ir,1)+vckli(ir)
              enddo
            endif
            tstart=.false.

            if(tconv) then
c              write(iukli,'(a)') ' --- kli ae shifts (Ha) --- '
c              do i=1,norb_kli
c               write(iukli,'(a2,3i3,1x,e14.8)') '%',i,n(i),l(i),svkli(i)
c              enddo
c              write(iukli,*)
               
c              open(11,file='vkli_ae.dat')
c              write(11,'(a)') '#ae r(i) vxc(i) vslater(i)'
c              do i=1,mmax
c                write(11,*) r(i),vxc(i,1),vsl(i)
c              enddo
c              close(11)
c repeat klix but only for valence orbitals
              do ir=1,mmax
                dens(ir)=0.5d0*(rho(ir)-dc(ir))
              enddo
              i=norb_kli-nc
              call at_vklix(i,n(nc+1),l(nc+1),mp(nc+1),fp(nc+1),mmax,r
     1,          uu(1,nc+1),dens,vxc_v,vsl_v,svkli(nc+1),svmdummy
     1,          2,.false.,.false.,.false.,mx,ms,mev)

c              write(iukli,'(a)') ' --- kli valence shifts (Ha) --- '
c              do i=nc+1,norb_kli
c               write(iukli,'(a2,3i3,1x,e14.8)') '%',i,n(i),l(i),svkli(i)
c              enddo
c              write(iukli,*) 

c              open(11,file='vkli_valence_ae.dat')
c              write(11,'(a)') '#aevalence r(i) vxc(i) vslater(i)'
c              do i=1,mmax
c                write(11,*) r(i),vxc_v(i),vsl_v(i)
c              enddo
c              close(11)
            endif
          endif
        else
c lda or gga xc potential
          if(nspinmx.eq.1) then
c spin unpolarized
            call vexcor(iexc,mmax,r,d(1,1),dp(1,1),dpp(1,1),vxc(1,1),
     1                  eexc,ex,ec,.false.)
          else
c spin polarized
           call vexcos(iexc,mmax,mmax,mx,r,d,dp,dpp,vxc,eexc,ex,ec,
     +                 .false.)
          endif

        endif

Ctest
c        write(6,*) "After XC : d(9,1)   = ", d(9,1)
c        write(6,*) "After XC : dp(9,1)  = ", dp(9,1)
c        write(6,*) "After XC : dpp(9,1) = ", dpp(9,1)
c        write(6,*) "After XC : vxc(9,1) = ", vxc(9,1)
c        write(6,*) "After XC : ves(9)   = ", ves(9)
c        write(6,*) "After XC : vee(9,1) = ", vee(9,1)
Ctest end

c generate next iteration's potential using anderson method
        if(tsic) call dcpv(mx,1,mmax,vglob,vee(1,1))
        do ispin=1,nspinmx
          call dadv(mmax,1,mmax,ves,vxc(1,ispin),vo(1,ispin))
ctest
c          write(6,*) "Before anderson."
ctest end
          call anderson(bl,it,mmax,r,vo1(1,ispin),vo(1,ispin),
     1                  vi1(1,ispin),vee(1,ispin))
ctest
c          write(6,*) "After anderson."
ctest end
        enddo

c
c if iexc=0, overwrite with coulomb potential
        if(iexc.eq.0) then
          do ispin=1,nspinmx
            do j=1,mmax
              vee(j,ispin)=-z*r1i(j)
            enddo
          enddo
        endif

c total energy
        eeig=0.d0
        do ispin=1,nspinmx
          do i=nstart(ispin),nend(ispin)
            eeig=eeig+f(i)*e(i)
          enddo
        enddo
        de=abs((eeig-eeig_old)/eeig)
        eeig_old=eeig

      if(tconv) goto 24
      enddo
  24  continue

!      write (6,*) "sratom finished."

ctest
c     write charge density derivative
c      open (50,file='free_drho_dr.dat')
c      do j = 1, mmax
c        write(50,*) j, r(j), dp(j,1)
c      enddo
c      close(50)
ctest end

c proper core density for output
      do ispin=1,nspinmx
        do j=1,mmax
          dcs(j,ispin)=dcs(j,ispin)*r2i(j)
          dcsp(j,ispin)=dcsp(j,ispin)*r3i(j)
          dcspp(j,ispin)=dcspp(j,ispin)*r4i(j)
        enddo
      enddo

      return
      end
c

      subroutine at_dsuv(mx,ml,mh,a,b,c)
      implicit real*8 (a-h,o-z)
      integer mx,ml,mh
      dimension a(mx),b(mx),c(mx)
      integer i
      do i=ml,mh
        c(i)=a(i)-b(i)
      enddo
      return
      end

c Header: /f21/a/fuchs/Psp/Dkef/Ds/RCS/vklix.f,v 1.1 95/08/22 17:42:15
c fuchs
c Exp Locker: fuchs 
c
c compute Krieger-Li-Iafrate (kli) approximation to the optimized
ceffective
c potential of exchange-interaction, assuming spherical symmetry
c
c input
c iorb ........... number of orbitals in present spin channel
c n() ............ principal quantum number
c l() ............ angular momentum quantum number
c m() ............ spin projection 
c f() ............ orbital occupation numbers
c mmx ............ upper limit index on radial grid
c r() ............ radial mesh
c ud0() .......... radial kohn sham orbitals
c                  normalized to 4*pi
c                  ordered wrt to increasing eigenvalues
c dens() ......... spin density in present channel
c mode ........... 1 ........ evaluate orbital shifts
c                  2 ........ use input orbital shifts
c tene ........... .true. ... evaluate diagonal energy contributions
c                             i.e. exact exchange energy
c tstart ......... .true. ... initialize 
c
c output
c voep() ......... kli approximated oep
c sv()   ......... orbital shifts
c svx()  ......... orbital exchange energy
c
c Martin Fuchs 20-08-1995
c 
      subroutine at_vklix(iorb,n,l,m,f,mmx,r,ud0,dens,voep,vsl,sv,svx
     +,                   mode,tstart,tene,tprint,mx,ms,mev)

      implicit real*8 (a-h,o-z)
      integer iorb, n, l, m, mmx, mev
      include 'parameter.h'

      character*1 sni,sli,smi
      character*1 snj,slj,smj
      logical tene,tprint,tstart
      integer iukli,iuslater,mode,mx,ms, mup
      real *8 al,amesh,ahat
      integer la
      integer i, j, k, idm, ll

c     common/arr/r(mx),vi(mx),u(mx),up(mxr,),upp(mx)
c     common/par/z,amesh,al,mmax,nin,mch,iexc
c     common/vslater/vsl(mx)

      parameter(iukli=80,iuslater=81)
      parameter(zero=0.d0,one=1.d0,eps=1.d-12)

      integer indx
      dimension ud0(mx,mev),dens(mx),voep(mx)
     +,   dorb(mx,mev),vorb(mx,mev)
     +,   vxi(mx,mev),indx(mev),sm(mev,mev)
     +,   sv(mev),svx(mev),vav(mx)
     +,   n(mev),l(mev),m(mev),f(mev),svsave(mev)
     +,   smsave(mev,mev), uc(mx),wc(mx),rsqi(mx),r(mx),vsl(mx)
      integer icount
      data icount/0/
      save icount

c functions
      ahat(la)=2.d0*la+1

c files
      if(tprint) then
        open(unit=iuslater,file='vslater_accumul.dat',status='unknown')
        write(iukli,'(a,/)') 's List of Fock matrix elements ---'
      endif

c initializations and limiting radial range (overflow precaution)
      amesh=r(2)/r(1)
      al=log(amesh)
      if(tstart) icount = 0
      icount=icount+1
      mup=mmx
c      if(icount .eq. 1) then
c        write(iukli,*) '%at_vklix: check orbital parameters'
c        write(iukli,*) 'iorb=',iorb
c        do i=1,iorb
c          write(iukli,*) 'n l f', n(i),l(i),f(i)
c        enddo
c        write(iukli,*) '%at_vklix: end printout'
c      endif

      do j=1,mmx
        rsqi(j)=one/r(j)**2
      enddo


      do j=1,mmx
        voep(j)=zero
        vav(j)=zero
        if(dens(j) .lt. eps) mup=min(mup,j)
      enddo
      do j=1,iorb
        do i=1,iorb
          sm(j,i)=zero
      enddo
      enddo

! TEST
!      write(ie,*) 'mup:', mup

c compute (spherically averaged) orbital densities
      do i=1,iorb
        do j=1,mmx
          dorb(j,i)=ud0(j,i)*ud0(j,i)*rsqi(j)
          vorb(j,i)=zero
        enddo
      enddo

!TEST
!      write(ie,*) 'dorb:', dorb(:,1)
!      write(ie,*) 'dens:', dens(:)
!      write(ie,*) 'dens:', dens(), dens(1136),dens(1137)



c compute (spherical) exchange potentials
      do i=1,iorb
        if(tprint) then
          write(iukli,'(/,a)') 's Fock matrix elements FME:'
        endif

        do j=1,iorb
          if(i .ge. j .and. tprint) then
            call at_labels(n(i),l(i),m(i),sni,sli,smi)
            call at_labels(n(j),l(j),m(j),snj,slj,smj)
            write(iukli,'(a,4x,2i2,30x,a1,a1,1x,a1,a1)')
     +        's ij',i,j,sni,sli,snj,slj 
            write(iukli,'(a,1x,f8.4,2x,f8.4)')
     +        's f_i f_j    ',f(i),f(j)
          endif

!          write(ie,*) 'In at_vklix. i,j = ', i,j


          do ll=abs(l(i)-l(j)),l(i)+l(j)

            cj=-f(j)/ahat(ll)*at_acgc(l(j),l(i),ll)**2

!            write(ie,*) 'In at_vklix. ll = ', ll

c configuration averaging for open shells

            if(j .eq. i .and. f(j) .lt. ahat(l(j))) then
c monopole term
              if(ll .eq. 0) then
                cj=-one
c multipole terms
              else
                cj=-max(zero,(f(j)-one))
     +            *ahat(l(j))/(2*l(j)*ahat(ll))*at_acgc(l(j),l(j),ll)**2
              endif
            endif

cmf
c           cj=-f(j)*at_ack(l(i),l(j),ll)


            if(cj .lt. zero) then

c              if(icount .lt. 2) then
c                  print *,
c    + i,j,'  ',l(i),l(j),ll,l(i)+l(j)+ll,cj,at_acgc(l(j),l(i),ll)
c                endif

c             call arhf
c    +             (mup,r,l(i),l(j),ll,ud0(1,i),ud0(1,j),uc,wc,vxi(1,j))
              call at_arhf
     +             (mup,r,l(i),l(j),ll,ud0(1,i),ud0(1,j),uc,wc,vxi(1,j),
     +             mx)

              do k=1,mup
                vorb(k,i)=vorb(k,i)+cj*(vxi(k,j)*ud0(k,j)*rsqi(k))
              enddo

            endif

            if(i .ge. j .and. tprint) then
             do k=1,mup
              uc(k)=vxi(k,j)*ud0(k,j)*rsqi(k)
             enddo
             dfock=dmelm(mup,al,r,ud0(1,i),uc)
             write(iukli,'(a,28x,i2,i2,i2)')
     +         's l_i l_j L  ',l(i),l(j),ll
             write(iukli,'(a,f8.4,2x,f8.4)')
     +         's c_j c_j/f_j ',cj,cj/f(j)
             write(iukli,'(a,2(2x,e12.6))')    
     +         's FME FME*c/f ',dfock,dfock*cj/f(j)
             write(iukli,'(a)') 's'
            endif

          enddo

        enddo

c compute density weighted (Slater-like) potential
        do j=1,mup
          vav(j)=vav(j)+f(i)*ud0(j,i)*vorb(j,i)
        enddo
        if(tprint) then
         do j=1,mup
          vsl(j)=vav(j)/dens(j)
          write(iuslater,*) r(j),vsl(j) 
         enddo
         write(iuslater,*) '&',j
        endif

      enddo

      do j=1,mup
        vav(j)=vav(j)/dens(j)
        vsl(j)=vav(j)
      enddo

c average potential option
c     do j=1,mup
c       voep(j)=vav(j)
c     enddo
c     return

      if(mode .eq. 1) then
c potential constants
      if(tene) then
c        write(iukli,'(a,i2)') '%df diagonal fock elements',iorb
        do i=1,iorb
          do j=1,mup
            uc(j)=vorb(j,i)
          enddo
          do k=1,i
            sv(i)=dmelm(mup,al,r,ud0(1,k),uc)
c            write(iukli,'(a,5i3,f12.6,1x,f12.6)')
c     1        '%df',i,n(i),l(i),n(k),l(k),sv(i),sv(i)*27.2116
          enddo
        enddo
      endif
      do i=1,iorb
        do j=1,mup
          uc(j)=vav(j)*ud0(j,i)*rsqi(j)-vorb(j,i)
        enddo
        sv(i)=dmelm(mup,al,r,ud0(1,i),uc)
c       svsave(i)=sv(i)
      enddo

c linear equation setup 
      sv(iorb)=zero
      do i=1,iorb-1
        do j=1,mup
          uc(j)=dorb(j,i)/dens(j)
        enddo
        do j=1,i
          sm(i,j)=-f(j)*dmelm(mup,al,r,uc,dorb(1,j))
          sm(j,i)=sm(i,j)*f(i)/f(j)
          if(i .eq. j) sm(i,j)=sm(i,j)+1.d0
c         smsave(i,j)=sm(i,j)
c         smsave(j,i)=sm(j,i)
        enddo
      enddo

!     TEST
      if(tprint) then
      write(iukli,*) '--- vkli rhs vectors sv()---'
      do i=1,iorb
        write(iukli,*) i, sv(i)
      enddo
      write(iukli,*) '--- vkli corrector matrix sm() ---'
      do i=1,iorb-1
        write(iukli,'(i2,4(1x,e12.6))') i, (sm(i,j),j=1,min(4,iorb-1))
      enddo
      endif


c solve system --- sv(i) are the kli correctors
      idm=iorb-1
      if(idm.gt.0) then
        call at_ludcmp(sm,idm,mev,indx,xxx)
        call at_lubksb(sm,idm,mev,indx,sv)
      endif

c     call dgesv(idm,1,sm,ms,indx,sv,ms,ierr)

c     print *, '--- vkli average difference potential ---'
c     print *, 'iorb',iorb
c     do i=1,iorb
c       print *, i, n(i),l(i),m(i),f(i),sv(i)
c     enddo

c     print *, '--- vkli verify solution ---'
c     do i=1,iorb-1
c       test=zero
c       do j=1,iorb-1
c         test=test+smsave(i,j)*sv(j)
c       enddo
c       print *, i, test,'==?',svsave(i)
c     enddo
      else if (mode .ne. 2) then
        print *, '& - at_vklix - warning: unimplemented mode',mode
      endif

!     TEST
      if(tprint) then
      write(iukli,*) '--- vkli orbital correctors ---'
      do i=1,iorb-1
        write(iukli,*) i, sv(i)
      enddo
      endif


c compute kli potential
      do i=1,iorb

c       iunit=50+i
c       if(m(i).eq.1) then
c         open(iunit)
c         rewind(iunit)
c         do j=1,mup
c           write(iunit,*) r(j),f(i)*ud0(j,i)*vorb(j,i)/dens(j)
c         enddo
c         close(iunit)
c       endif

        do j=1,mup
          voep(j)=voep(j)+f(i)*sv(i)*dorb(j,i)
        enddo
      enddo
      do j=1,mup
        voep(j)=vav(j)+voep(j)/dens(j)
      enddo


      do j=mup+1,mmx
        voep(j)=-one/r(j)
      enddo
     
c orbital exchange energies
      if(tene) then
        do i=1,iorb
          svx(i)=dmelm(mup,al,r,ud0(1,i),vorb(1,i))
        enddo
      endif

c     print *, '--- vkli done ---'
c
      if(tprint) then
        close(iuslater)
        write(iukli,'(a,/)') 's ---------------------------'
      endif
      return
      end
c
c
c ancillary stuff
c
c
      real*8 function at_acgc(l1,l2,l3)

      implicit real*8 (a-h,o-z)
      integer l1, l2, l3
      integer mx
      parameter(mx=60)
      dimension fac(0:mx)
      logical ifc
      data    ifc/.true./
      save    ifc,fac
      integer i

      t(i)=fac(i/2)/sqrt(fac(i))
      hatl(i)=dble(2*i+1)

      if(ifc) then
        ifc=.false.
        fac(0)=1.d0
        do i=1,mx
          fac(i)=dble(i)*fac(i-1)
        enddo
      endif

      if(l1+l2+l3 .gt. mx) stop 'at_acgc --- l1+l2+l3'
      if(at_notri(l1,l2,l3) .gt.0) then
        at_acgc=sqrt(hatl(l3)/(l1+l2+l3+1))
     &         *t(l1+l2+l3)/
     &           (t(l1+l2-l3)*t(l1-l2+l3)*t(l2+l3-l1))
        if(mod((l1+l2-l3)/2,2) .ne. 0.) at_acgc=-at_acgc
      else
        at_acgc=0.d0
      endif

      return
      end

      subroutine at_arhf(mmax,r,la,lb,lc,ua,ub,uc,wc,yabc,mx)
c
      implicit real*8 (a-h,o-z)
      logical tfirst
      integer mmax, la, lb, lc, mx
c
      include 'parameter.h'
      integer ll, i, j
      integer llmx
      parameter(llmx=10)
      dimension r(mx),ua(mx),ub(mx),uc(mx),yabc(mx),wc(mx)
      dimension r_p(mx,llmx),r_i(mx,llmx),fi(mx),fo(mx)
      parameter(zero=0.d0)
      data tfirst/.true./
      save tfirst,al

c      if(tfirst) then
        al=0.1d0*log(r(11)/r(1))

      do ll=1,llmx
        do i=1,mmax
!         if(abs(r(i)).lt.1.e-10) then
!          write(6,*)"r(i):",i, r(i)
!         endif
          r_p(i,ll) = r(i)**(ll-1)
          r_i(i,ll) = r(i)**(1-ll)
        enddo
      enddo
c      endif

c evaluate radial matrix element w/ trapezoidal rule 
      do i=1,mmax
c       fi(i) = ua(i)*ub(i) * r(i)**(lc)        ! r**() goes into r_p(i,lc+1)
c       fo(i) = ua(i)*ub(i) * r(i)**(-(lc+1))   ! r**() goes into r_i(i,lc+2)
        fi(i) = ua(i)*ub(i) * r_p(i,lc+2)
        fo(i) = ua(i)*ub(i) * r_i(i,lc+1)
      enddo

      do i=1,mmax
        yabc(i)=zero
        do j=1,i-1
          yabc(i)=yabc(i)+fi(j)
        enddo
        yabc(i)=al*(.5d0*fi(i)+yabc(i))*r_i(i,lc+2)
        aux=zero
        do j=i+1,mmax-1
          aux=aux+fo(j)
        enddo
        aux=al*(aux+.5d0*(fo(i)+fo(mmax)))*r_p(i,lc+1)
        yabc(i)=yabc(i)+aux
      enddo

      return
      end

c
c***********************************************************************
c tags for angular momentum to literal conversion
c***********************************************************************
      subroutine at_labels(n,l,m,sn,sl,sm)
      integer n,l,m
      character*1 sn,sl,sm,cn(7),cl(6),cm(3)
      data cn/'1','2','3','4','5','6','.'/
      data cl/'s','p','d','f','g','.'/
      data cm/'+','-','.'/
      if(n .lt. 1 .or. n .gt. 6) n=7
      if(l .lt. 0 .or. l .gt. 4) l=5
      if(m .lt. 1 .or. m .gt. 2) m=3
      sn=cn(n)
      sl=cl(l+1)
      sm=cm(m)
      return
      end
c

      subroutine at_ludcmp(A,N,NP,INDX,D)
      implicit real*8 (a-h,o-z)
      integer N, NP, INDX
      integer I,J,K
      integer IMAX,NMAX
      real*8 TINY
      PARAMETER (NMAX=100,TINY=1.0d-20)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      D=1.d0
      DO 12 I=1,N
        AAMAX=0.d0
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1.d0/AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.d0
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1.d0/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END

      subroutine at_lubksb(A,N,NP,INDX,B)
      implicit real*8 (a-h,o-z)
      integer N, NP, INDX
      DIMENSION A(NP,NP),INDX(N),B(N)
      integer II, LL
      integer I, J

      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END

c verify triangle rule
      FUNCTION at_NOTRI(K,L,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer K, L, M
      at_NOTRI=-1
      IF(MOD((K+L+M),2).EQ.1) GOTO 10
      IF((K+L-M).LT.0) GOTO 10
      IF((K-L+M).LT.0) GOTO 10
      IF((M+L-K).LT.0) GOTO 10
      at_NOTRI=1
   10 RETURN
      END


