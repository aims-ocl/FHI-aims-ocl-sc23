!  This file is part of the ATOM_SPHERE code and is interfaced
!  by Rundong Zhao for use in FHI-aims.
!
!  ATOM_SPHERE was originally written by Stefan Goedecker
!  while he was at Cornell University, the Max Planck Institut Stuttgart
!  the CEA Grenoble and the university of Basel.
!  Santanu Saha at Basel university has interfaced it to LIBXC and
!  improved a few other things

! ATOM_SPHERE is free software: you can redistribute it and/or modify
! it under the terms of the version 3 of the license of the
! GNU Lesser General Public License as published by the Free
! Software Foundation.
!
! ATOM_SPHERE is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with ATOM_SPHERE. If not, see 
! https://www.gnu.org/licenses/lgpl-3.0.en.html
!
! ATOM_SPHERE reflects a substantial effort on the part of the
! developers, and we ask you to respect the spirit of the
! license that we chose and in particular to  keep any derivatives
! of ATOM_SPHERE under the same license that we chose for
! the original distribution, the GNU Lesser General Public License.

! We use atom_sphere in FHI-aims with permission from Stefan Goedecker 
! and his group. Please respect the license - thanks! VB




!      \int (r*R_l(r))**2 dr = \sum_j psi(j,...,l+1)**2 rw(j,2) = 1
!      Consequently a charge density is given by  R_l(r)**2/(4*pi)
!

 subroutine zorasolver(nameat,nprinx,nprin,nrad,lmaxp,lmax,norb,norbmax,idsx,nspin,nspol,occup,eval,chrg, &
  residue,ekin_orb,epot_orb,noae,no,il,lo,zo,so,ncovmax,ncov,potloc,pothart,potxc,pothf,pot0,ccleb,grad, &
  gradp,grads,gradhf,hhp,ssp,epslag,shift,rcov,rprb,znuc,psi,psid,kinetic,rr,rw,rr4,rw4,screened, &
  omega,hfmix,rhotot,rhoud,rhoij,screenexkernel,ppois,wrky,adiis,toten,epstres,pspcalcul,kinonly,time, &
  cut_atom,cutoff_type,r_cut,scale_cutoff,w_cutoff)

  use localorb_io, only: localorb_info
  use physics, only: finite_nuclear_radius
  use libxcModule
  implicit real*8 (a-h,o-z)
  !parameter(nprinx=7) !! Principal Quantum Number for each channel
  !parameter(nrad=232)
  !parameter(lmax=2,lmaxp=2*lmax)  !! What is lmax here Santanu??
  !parameter (nspol=2) ! just for better readability, any other value than 2 will not work at present
                      !! Spin Polarization 1  or 2
  !parameter(idsx=10)  ! length of the DIIS history
  logical screened,pspcalcul,kinonly
  character*2 nameat
  character*1 il(5)
  parameter(c=137.0359895d0)

  ! Array defined for storing the n,l,s of the AE and Pseudo Orbitals
  dimension no(norbmax),noae(norbmax),lo(norbmax),so(norbmax),  &
     zo(norbmax),nmod(nprinx,lmax+1,nspol)

  !! Radial Grid Variables
  dimension rr(nrad), rw(nrad,3), rr4(4*nrad),rw4(4*nrad,3),ww(4*nrad)
  ! rr, rr4 : Radial Grid Array, rr4 for the KE
  ! rw(*,1) : Pure wigths, 
  ! rw(*,2) : rw(*,1) mulitplied with r^2, 
  ! rw(*,3) : Inverse radial weight for GGA's


  !! Variables for describing Orbitals and density     
  dimension  nprin(lmax+1,nspol), &  ! Store n,l,s of pseudo orbitals
      occup(nprinx,nspol,lmax+1), &  ! occupation number 
      ene(nprinx,nspol,lmax+1), & ! Array to store energy
      psi(nrad,nspol,nprinx,lmax+1), &  ! radial wavefunctions (input and output of WF optimizationi)
      psid(nrad,nspol,nprinx,lmax+1,idsx), &   !WF history  used uniquely in the WF optimization part
      psidgood(nrad,nspol,nprinx,lmax+1,idsx), &   !WF history  used uniquely in the WF optimization part
      rhotot(nrad), &  ! total charge density
      rhoud(nrad,nspol), & ! spinup, down charge density
      rhocore(nrad), & ! core charge for nonlinear core correction
      rhoij(nrad) ! cross charge densities for Hartree Fock


  !! Array for storing potentials and HF related terms    
  dimension &
      ppois(nrad,-(lmaxp+1):lmaxp,2), &  !auxiliary array for solving Poisson's equation
      pothart(nrad), & ! Hartree potential
      Vnuc(nrad),& ! potential due to znuc (nucleus)
      potloc(nrad), &  ! Local potential of pseudopotential
      potxc(nrad,nspol), & ! Exchange Correlation potential
!!!   potrel(nrad), potrel4(4*nrad), & ! Relativistic Potential 
!!!   workrel(nrad), & ! Temporary work matrix calculating relativistic contribution in K.E.
      damprel(nrad), damprel4(4*nrad), & ! Relativistic dampening factor for kinetic energy
      screenexkernel(nrad,nrad,lmaxp+1), &  ! represents the action of the screened Coulomb operator
      pothf(nrad,lmaxp+1),   & ! Hartree Fock potential
      pot0(nrad,nprinx,lmax+1) ! potential used for preconditioning


  !! Array for storing different Gradients    
  dimension &
      grad(nrad,nspol,nprinx,lmax+1), &  ! total gradient
      gradp(nrad,nspol,nprinx,lmax+1,idsx), &  ! histopry list of preconditioner gradients for DIIS
      grads(nrad,nspol,nprinx,lmax+1), &  ! gradient of overlap
      gradhf(nrad,nspol,nprinx,lmax+1) ! gradient form Hartree Fock


  !! Array reqd for solving KS equation with orthonormalization constraint   
  dimension &
      hhp(2,nrad,nprinx,lmax+1), ssp(2,nrad), & ! Hamiltonian and overlap matrice for preconditioning
      epslag(nprinx,nprinx,nspol,lmax+1), &  ! Lagrange multiplier matrix
      residue(nprinx,nspol,lmax+1), &  ! residue of KS Orbitals
      shift(nprinx,lmax+1), & ! shifts used for preconditoning (CRITCAL FOR FAST CONVERGENCE!)
      zeffstore(nprinx,lmax+1), & ! shifts used for preconditoning (CRITCAL FOR FAST CONVERGENCE!)
      ccleb(lmaxp+1,lmax+1,lmax+1), & !contracted (overl m) Clebsch coefficients
      wrky(nrad,nprinx), &  !Work array 
      adiis(idsx+1,idsx+1,3),ipiv(idsx+1),rdiis(idsx+1) ! work arrays for DIIS


  !! Different energy, charge and related terms   
  dimension &
      ekin_orb(nprinx,nspol,lmax+1), &  ! Kinetic Energy of KS Orbitals
      epot_orb(nprinx,nspol,lmax+1), &  ! Potential Energy of KS Orbitals for znuc 
      eval(nprinx,nspol,lmax+1), & ! eigenvalues of the Kohn Sham Hamiltonian (Fock matrix)
      chrg(nprinx,nspol,lmax+1,ncovmax), & ! Charge of  the Kohn Sham Orbitals
    ! electric(nprinx,nspol,lmax+1,ncovmax), & ! Integration of r*psi**2 for the Electric field
      rcov(ncovmax),time(5) ! Array for cut off radii and Time

  logical, intent(in) :: cut_atom
  integer, intent(in) :: cutoff_type
  real*8, intent(in) :: r_cut
  real*8, intent(in) :: scale_cutoff
  real*8, intent(in) :: w_cutoff
  real*8, external :: cutoff_pot
  real*8, intent(out) :: kinetic(nrad,nspol,nprinx,lmax+1)
  character*300 :: info_str

  if (cut_atom) then 
    call detnp(nrad,rr,r_cut+w_cutoff,ncut)
    ncut=ncut-1
  endif

  !rhocore=0.d0
  pi=4.d0*atan(1.d0)
  c_sq = c*c
  itertot=0
  fourpi=16.d0*atan(1.d0)
  fourpiinv=1.d0/fourpi
  twopi=fourpi/2.d0
  ene=0.d0
  grad=0.d0
  epsdamp = 1.d-12
  damprel    = 0.5d0
  damprel4   = 0.5d0


  nmod=1
!! SETTING UP DIFFERENT SEGMENTS
  potloc=0.d0
  Vnuc=0.d0
! Finite nucleus Coulomb potential
  write(info_str,'(2X,A,F15.8,A)')'Using a finite nuclear model for the Coulombic potential in atom_sphere, with a radius of ',&
       finite_nuclear_radius, ' bohr'
  call localorb_info( info_str )
  size_nucleus=finite_nuclear_radius
  if (size_nucleus.gt.0.d0) then
    j=1
    if (rr(j).eq.0.d0) then
      potloc(j)=potloc(j)-2.d0*znuc/(size_nucleus*sqrt(pi))
      Vnuc(j)=potloc(j)
    else
      potloc(j)=potloc(j)-znuc*derf(rr(j)/size_nucleus)/rr(j)
      Vnuc(j)=potloc(j)
    endif
    do j=2,nrad
      potloc(j)=potloc(j)-znuc*derf(rr(j)/size_nucleus)/rr(j)
      Vnuc(j)=potloc(j)
    enddo
    !!znuc=0.d0 !the nucleonic potential is now contained in potloc and to avoid doublecounting we put znuc=0
  endif

!! confinement potential
  if (cut_atom) then
    do j=1,nrad
      potloc(j)=potloc(j)+cutoff_pot( rr(j), cutoff_type, r_cut, w_cutoff, scale_cutoff )
    end do
! Add parabolic confining potential
  else if (rprb.gt.0.d0) then
     do j=1,nrad
       potloc(j)=potloc(j)+.5d0*(rr(j)/rprb**2)**2
     enddo
  endif

  call crtssp(nrad,rr,ssp)

!!!!!FIXME: Building -1/2 d^2/dr^2 - Z/r + l(l+1)/2r^2 
!!        z=znuc
!!        call crthh(z,nrad,rr,rw,lmax,potloc,hh)

  call crtpois(nrad,lmaxp,rr,rw(1,1),ppois)
  !call crtpois(nrad,0,rr,rw(1,1),ppois)

  if (cut_atom) then
    !write(*,*)"Zerotail in the begin"
    call zerotail(lmax,nspol,nprinx,nprin,nrad,ncut,psi)
  end if

!GS Orthogonalization of psi
  call gsortho(nrad,nspol,nprinx,nprin,lmax,rw(1,2),wrky,psi,psid)
  call DCOPY(nrad*nspol*nprinx*(lmax+1),psid,1,psi,1)
  psidgood=psid !! STORE THE INITIAL GUESS WAVEFUNCTION

  !call checkorthogonalization(nrad,nspol,nprinx,nprin,lmax,rw(1,2),wrky,psi)

  call contract_clebsch(lmax,lmaxp,ccleb)
  write(6,*)'             '
  write(6,*)'             '
  write(6,*)'             '

!! SOLVING SCF CYCLE

! if there is no convergene within 50 steps of the 2000 loop the shift are presumably not optimal and they will 
! be set equal to the eigenvalues of the previous 3000 loop
  write(6,'(a)')"============================================================================================================"
  write(6,'(a)')" irel  ishift  ids    alpha         Total Energy          Energy Diff           Residue         Rel Pot Diff"
  write(6,'(a)')"============================================================================================================"

  call cpu_time(t)
  time(2)=time(2)-t

  tres=0.d0;tresnew=0.d0;tresold=0.d0;tresmin=1.d100
  diff=1.d10
  do 9201, irel = 1, 500
  do 3000, ishift=1,50 !! 5 
    !write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ',ishift
    eval=0.d0

    !DIIS Wavefunction minimization  loop begins?
    if (.not.cut_atom) then
    if (ishift.eq.11) then
       do isp=1,nspol
       do l=0,lmax
       do iprin=1,nprin(l+1,isp)
        !if (occup(iprin,isp,l+1).ge.1.d0) then
         do j=nmod(iprin,l+1,isp),nrad
            psi(j,isp,iprin,l+1)=0.d0
         end do
        !end if
       end do
       end do
       end do
    else if ((ishift.gt.11).and.(tres.gt.tresmin)) then
       do isp=1,nspol
       do l=0,lmax
       do iprin=1,nprin(l+1,isp)
       !if (occup(iprin,isp,l+1).ge.1.d0) then
         do j=nmod(iprin,l+1,isp),nrad
            psi(j,isp,iprin,l+1)=0.d0
         enddo
       !end if
       end do
       end do
       end do
    end if
    end if


    do l=0,lmax
    do isp=1,nspol
    do iprin=1,nprin(l+1,isp)
    do j=1,nrad
       psid(j,isp,iprin,l+1,1)=psi(j,isp,iprin,l+1)
    enddo
    enddo
    enddo
    enddo

    adiis=0.d0
    toten=0.d0
    alpha=0.75d0
    nitermax=25

    do 2000,ids=1,nitermax
      mids=mod(ids-1,idsx)+1
      if (idsx.eq.1 .or. (.not.pspcalcul)) mids=1

! calculate charge densities
      do isp=1,nspol
      do j=1,nrad
         rhoud(j,isp)=0.d0
         enddo
         do l=0,lmax
         do iprin=1,nprin(l+1,isp)
            do j=1,nrad
               rhoud(j,isp)=rhoud(j,isp)+occup(iprin,isp,l+1)*psid(j,isp,iprin,l+1,mids)**2
            enddo
         enddo
         enddo
      enddo

      if (nspol.eq.2) then
         do j=1,nrad
           rhoud(j,1)=rhoud(j,1)*fourpiinv
           rhoud(j,2)=rhoud(j,2)*fourpiinv
           rhotot(j)=rhoud(j,1)+rhoud(j,2)
         enddo
      else
         do j=1,nrad
           rhoud(j,1)=rhoud(j,1)*fourpiinv
           rhotot(j)=rhoud(j,1)
         enddo
      endif

! check normlization
      do isp=1,nspol
         sum=0.d0
      do j=1,nrad
         sum=sum+rhoud(j,isp)*rw(j,2)
      enddo
         sum=sum*fourpi
      enddo

! calculate gradient
        
      totenold=toten
      toten=0.d0
      z=0.d0  !! the nucleonic potential is now contained in potloc and to avoid doublecounting we put z=0
      kinonly=.false. 

! calculate unconstrained gradient
      do l=0,lmax
      do isp=1,nspol
      do iprin=1,nprin(l+1,isp)

! hh (kinetic energy + local potential) times psi
      !call mult(nrad,hh(1,1,l+1),psid(1,isp,iprin,l+1,mids),grad(1,isp,iprin,l+1))
!!   LOW RESIDUE RECIPE MODIFIED -SG
      call applyzorakinpot(nrad,rr,rw,rw4,kinonly,z,l,potloc,damprel,damprel4,psid(1,isp,iprin,l+1,mids),grad(1,isp,iprin,l+1))
      ene(iprin,isp,l+1)=DDOT(nrad,psid(1,isp,iprin,l+1,mids),1,grad(1,isp,iprin,l+1),1)
      !write(*,'(a,3(i3),1x,2(e12.5,2x))') 'energy expectation',l,iprin,isp,e, &
      !          occup(iprin,isp,l+1)
      enddo
      enddo
      enddo

      do l=0,lmax
      do isp=1,nspol
      do iprin=1,nprin(l+1,isp)
      toten=toten+ene(iprin,isp,l+1)*occup(iprin,isp,l+1)
      end do
      end do
      end do
      !write(*,'(a,e25.17)')"Kinetic + Local",toten

! Add hartree contribution
      !nlmp=(lmaxp+1)**2
      call poissoll(nrad,lmaxp,0,ccleb(1,1,1),ppois,rhotot,rw(1,2),pothart)
      ehart=0.d0
      do j=1,nrad
      ehart=ehart+pothart(j)*rhotot(j)!*rw(j,2)
      enddo
      ehart=.5d0*ehart*fourpi
      toten=toten+ehart

! Add XC terms
! Since the Hartree potential is calculated we can now add the core charge
      ! do j=1,nrad
      ! do isp=1,nspol
      !   rhoud(j,isp)=rhoud(j,isp)+.5d0*rhocore(j)
      ! enddo
      ! enddo

      call cpu_time(t)
      time(3)=time(3)-t
      enexc=0.d0;potxc=0.d0;vxc=0.d0

      if (hfmix.lt.1.d0) then
        call driveXCradial(nspol,nrad,rw(1,2),rw(1,3),rhoud,enexc,potxc)
        call cpu_time(t)
        enexc=enexc*fourpi
        !write(*,*) 'ehart,enexc',ehart,enexc
        !write(*,'(a,e25.17)')"enexc",enexc

        toten=toten+enexc  !! Modify this segment to include HF
        !! Weight Factor Already multiplied in driveXCradial
        !do j=1,nrad ! later on we will always need just the product of the potential times the radial weight
        !potxc(j,1)=rw(j,2)*potxc(j,1)
        !potxc(j,2)=rw(j,2)*potxc(j,2)
        !enddo

! This (vxc) is needed only for checking purposes
        vxc=0.d0
        do j=1,nrad
        !!!!!!vxc=vxc+potxc(j,1)*(rhoud(j,1)-.5d0*rhocore(j))
        vxc=vxc+potxc(j,1)*rhoud(j,1)
        !!!!vxc=vxc+potxc(j,2)*(rhoud(j,2)-.5d0*rhocore(j))
        vxc=vxc+potxc(j,2)*rhoud(j,2)
        enddo
        vxc=vxc*fourpi
      end if
      call cpu_time(t)
      time(3)=time(3)+t
      !write(*,'(a,e25.17)')"vxc",vxc

! Add Hartree Fock terms
      ehf=0.d0
      gradhf=0.d0
      call cpu_time(t)
      time(5)=time(5)-t

      if (hfmix.ne.0.d0) then
        do ksp=1,nspol
        do li=0,lmax
        do iprin=1,nprin(li+1,ksp)
           do lj=0,lmax
           do jprin=1,nprin(lj+1,ksp)
           if (occup(iprin,ksp,li+1).gt.1.d-20.or.occup(jprin,ksp,lj+1).gt.1.d-20) then
              do j=1,nrad
                 rhoij(j)=psid(j,ksp,iprin,li+1,mids)*psid(j,ksp,jprin,lj+1,mids)
              enddo

              if (screened) then
                 do lp=0,lmaxp
                    if (ccleb(lp+1,li+1,lj+1).ne.0.d0) then
                       call DGEMV('T', nrad, nrad, 1.d0, screenexkernel(1,1,lp+1),nrad,rhoij, 1, 0.d0, pothf(1,lp+1), 1)
                    endif
                 enddo
              else
                 call poissoll(nrad,lmaxp,lmaxp,ccleb(1,li+1,lj+1),ppois,rhoij,rw(1,2),pothf)
              endif
              
              do lp=0,lmaxp
                 if (ccleb(lp+1,li+1,lj+1).ne.0.d0) then
                    tt=0.d0
                    do j=1,nrad
                       tt=tt+pothf(j,lp+1)*rhoij(j)
                    enddo
                    ehf=ehf+tt*ccleb(lp+1,li+1,lj+1)*(occup(iprin,ksp,li+1)/(2*li+1))*(occup(jprin,ksp,lj+1)/(2*lj+1))
                 endif
              enddo

             do lp=0,lmaxp
                if (ccleb(lp+1,li+1,lj+1).ne.0.d0) then
                   fact=ccleb(lp+1,li+1,lj+1)*(occup(iprin,ksp,li+1)/(2*li+1))*(occup(jprin,ksp,lj+1)/(2*lj+1))
                do j=1,nrad
                   gradhf(j,ksp,iprin,li+1)=gradhf(j,ksp,iprin,li+1)+pothf(j,lp+1)*psid(j,ksp,jprin,lj+1,mids)*fact
                   gradhf(j,ksp,jprin,lj+1)=gradhf(j,ksp,jprin,lj+1)+pothf(j,lp+1)*psid(j,ksp,iprin,li+1,mids)*fact
                enddo
                endif
             enddo

           end if
           enddo
           enddo
        enddo
        enddo
        enddo
        ehf=-.5d0*ehf
        toten=toten+hfmix*ehf
      endif

      call cpu_time(t)
      time(5)=time(5)+t      


! Add hartree and XC contributions to gradient
      do l=0,lmax
      do isp=1,nspol
      do iprin=1,nprin(l+1,isp)
      do j=1,nrad
         grad(j,isp,iprin,l+1)=grad(j,isp,iprin,l+1)+(pothart(j)+potxc(j,isp))*psid(j,isp,iprin,l+1,mids) & 
                               - (.5d0*hfmix/occup(iprin,isp,l+1))*gradhf(j,isp,iprin,l+1)
      enddo
      enddo
      enddo
      enddo

! calculate Lagrange multiplier matrix
      evsum=0.d0
      do l=0,lmax
      do isp=1,nspol
      do iprin=1,nprin(l+1,isp)
         ! overlap matrix time psi
         do j=1,nrad
            grads(j,isp,iprin,l+1)=psid(j,isp,iprin,l+1,mids)*rw(j,2)
         enddo
         do jprin=1,nprin(l+1,isp)
            epslag(iprin,jprin,isp,l+1)=DDOT(nrad,psid(1,isp,iprin,l+1,mids),1,grad(1,isp,jprin,l+1),1)
         enddo
         evsum=evsum+occup(iprin,isp,l+1)*epslag(iprin,iprin,isp,l+1)
      enddo
      enddo
      enddo


!! THIS PART DECIDES THE VALUE OF NMOD FOR EACH ORBITAL
      if (ishift.eq.10.and.ids.eq.20) then
        rnum1=-1.d0*dlog(1.d-16)
        rnum2=dsqrt(-2.d0*dlog(1.d-16))
        do isp=1,nspol
        do l=0,lmax
        do iprin=1,nprin(l+1,isp)
        if (occup(iprin,isp,l+1).ge.1.d0) then
           !!write(6,*)"EPSLAG",iprin,l,isp,epslag(iprin,iprin,isp,l+1)
           if (epslag(iprin,iprin,isp,l+1).lt.0.d0) then
           rdenom1=(dsqrt(-2.d0*epslag(iprin,iprin,isp,l+1)))
           else
           rdenom1=(dsqrt(2.d0*epslag(iprin,iprin,isp,l+1)))
           end if
           r1=rnum1/rdenom1
           call detnp(nrad,rr,r1,nmod(iprin,l+1,isp))
           !! UNCOMMENT TO SET VIRTUAL ORBITALS
        else
           if (.not.cut_atom.and.rprb.gt.0.d0) then
              r2=rprb*rnum2
              call detnp(nrad,rr,r2,nmod(iprin,l+1,isp))
           else if (cut_atom) then
              !!r2=1.2d0*rnum2
              r2=r_cut
              !!r2=r_cut*rnum2
              call detnp(nrad,rr,r2,nmod(iprin,l+1,isp))
           else
              if (epslag(iprin,iprin,isp,l+1).lt.0.d0) then
              rdenom1=(dsqrt(-2.d0*epslag(iprin,iprin,isp,l+1)))
              else
              rdenom1=(dsqrt(2.d0*epslag(iprin,iprin,isp,l+1)))
              end if
              r1=rnum1/rdenom1
              call detnp(nrad,rr,r1,nmod(iprin,l+1,isp))
           end if
        end if
        end do
        end do
        end do
      end if

! transform psi, grad and grads to canonical ones
      call canonical(nrad,nspol,nprinx,nprin,lmax,epslag,psi,psid(1,1,1,1,mids),grad,grads,eval)  !Psi is used as work array

!     evtum=0.d0
!     do l=0,lmax
!     do isp=1,nspol
!     do iprin=1,nprin(l+1,isp)
!     evtum=evtum+occup(iprin,isp,l+1)*eval(iprin,isp,l+1)
!     enddo
!     enddo
!     enddo
!     write(*,*) 'evtum', evtum

! add constraints to gradient
      do l=0,lmax
      do isp=1,nspol
      do iprin=1,nprin(l+1,isp)
      do j=1,nrad-3,4
         grad(j+0,isp,iprin,l+1)=grad(j+0,isp,iprin,l+1)-eval(iprin,isp,l+1)*grads(j+0,isp,iprin,l+1) 
         grad(j+1,isp,iprin,l+1)=grad(j+1,isp,iprin,l+1)-eval(iprin,isp,l+1)*grads(j+1,isp,iprin,l+1) 
         grad(j+2,isp,iprin,l+1)=grad(j+2,isp,iprin,l+1)-eval(iprin,isp,l+1)*grads(j+2,isp,iprin,l+1) 
         grad(j+3,isp,iprin,l+1)=grad(j+3,isp,iprin,l+1)-eval(iprin,isp,l+1)*grads(j+3,isp,iprin,l+1) 
      enddo
      enddo
      enddo
      enddo

      if (cut_atom) then
         call zerotail(lmax,nspol,nprinx,nprin,nrad,ncut,grad)
         !write(6,*)"Zerotail in the tres calculation"
      endif

! residues
      tresold=tres
      tres=0.d0
      do l=0,lmax
      do isp=1,nspol
      do iprin=1,nprin(l+1,isp)
         t1=0.d0
         t2=0.d0
         t3=0.d0
         t4=0.d0
         do j=1,nrad-3,4
            t1=t1+grad(j+0,isp,iprin,l+1)**2
            t2=t2+grad(j+1,isp,iprin,l+1)**2
            t3=t3+grad(j+2,isp,iprin,l+1)**2
            t4=t4+grad(j+3,isp,iprin,l+1)**2
         enddo
         res=t1+t2+t3+t4
         residue(iprin,isp,l+1)=res
         tres=tres+res
      enddo
      enddo
      enddo
      tresnew=tres

! precondition gradient
      if (ids.eq.1) then ! Since preconditioner is fixed, it has to be calculated only in first iteration
        do l=0,lmax
        do iprin=1,max(nprin(l+1,1),nprin(l+1,2))
        do j=1,nrad
           pot0(j,iprin,l+1)=shift(iprin,l+1)  !use constant local potential
           if (cut_atom) then
             pot0(j,iprin,l+1)=pot0(j,iprin,l+1) +cutoff_pot( rr(j), cutoff_type, r_cut, w_cutoff, scale_cutoff )
           else if (rprb.gt.0.d0) then
             pot0(j,iprin,l+1)=pot0(j,iprin,l+1)+.5d0*(rr(j)/rprb**2)**2
           end if
        enddo
        enddo
        enddo
        z=0.d0 ! local potential entirely contained in potloc
        call crthhp(nspol,nprinx,nprin,nrad,lmax,rr,z,pot0,hhp) 
      endif

      do l=0,lmax
      do iprin=1,min(nprin(l+1,1),nprin(l+1,2))
         call tridag_two(nrad,hhp(1,1,iprin,l+1),grad(1,1,iprin,l+1),gradp(1,1,iprin,l+1,mids))
      enddo
      enddo

      do l=0,lmax
      do isp=1,nspol
      do iprin=min(nprin(l+1,1),nprin(l+1,2))+1,nprin(l+1,isp)
         call tridag(nrad,hhp(1,1,iprin,l+1),grad(1,isp,iprin,l+1),gradp(1,isp,iprin,l+1,mids))
      enddo
      enddo
      enddo

      !call DSCAL(nrad*nspol*nprinx*(lmax*1),.5d0,gradp(1,1,1,1,mids),1)

      if (idsx.eq.1 .or. (.not.pspcalcul)) then ! shortcut for steepest descent which is always used for all electron calculations

        ! determine stepsize alpha for wavefunction update 
         if (toten-totenold.gt.1.d-12) then
            alpha=alpha*.5d0
            alpha=max(alpha,1.d-2)
         else
            alpha=alpha*1.0905077d0
            alpha=min(alpha,0.75d0)
            !alpha=min(alpha,0.25d0)
         endif
! Update wavefunction
         do l=0,lmax
         do isp=1,nspol
         do iprin=1,nprin(l+1,isp)
         do j=1,nrad-3,4
            psi(j+0,isp,iprin,l+1)=psid(j+0,isp,iprin,l+1,1)-alpha*gradp(j+0,isp,iprin,l+1,1) 
            psi(j+1,isp,iprin,l+1)=psid(j+1,isp,iprin,l+1,1)-alpha*gradp(j+1,isp,iprin,l+1,1) 
            psi(j+2,isp,iprin,l+1)=psid(j+2,isp,iprin,l+1,1)-alpha*gradp(j+2,isp,iprin,l+1,1) 
            psi(j+3,isp,iprin,l+1)=psid(j+3,isp,iprin,l+1,1)-alpha*gradp(j+3,isp,iprin,l+1,1) 
         enddo
         enddo
         enddo
         enddo

      else    ! DIIS update scheme

! set up DIIS matrix (lower triangle)
         if (ids.gt.idsx) then
           ! shift left up matrix
            do i=1,idsx-1
            do j=1,i
             adiis(i,j,1)=adiis(i+1,j+1,1)
            enddo
            enddo
           ! reinitialize lowest line for summation
            do i=1,idsx
             adiis(idsx,i,1)=0.d0
            enddo
         endif

! calculate new line
         ist=max(1,ids-idsx+1)
         do i=ist,ids
            mi=mod(i-1,idsx)+1
            tt=0.d0
            do l=0,lmax
            do isp=1,nspol
            do iprin=1,nprin(l+1,isp)
               tt=tt+DDOT(nrad,gradp(1,isp,iprin,l+1,mids),1,gradp(1,isp,iprin,l+1,mi),1)
            enddo
            enddo
            enddo
            adiis(min(idsx,ids),i-ist+1,1)=adiis(min(idsx,ids),i-ist+1,1)+tt
         enddo

! copy to work array, right hand side, boundary elements
         do j=1,min(idsx,ids)
            adiis(min(idsx,ids)+1,j,2)=1.d0
            rdiis(j)=0.d0
            do i=j,min(idsx,ids)
               adiis(i,j,2)=adiis(i,j,1)
            enddo
         enddo
         adiis(min(idsx,ids)+1,min(idsx,ids)+1,2)=0.d0
         rdiis(min(idsx,ids)+1)=1.d0


! solve linear system
         call DSYSV('L',min(idsx,ids)+1,1,adiis(1,1,2),idsx+1,ipiv,rdiis,idsx+1,adiis(1,1,3),(idsx+1)**2,info)
         if (info.ne.0) print*, 'DGESV',info
         if (info.ne.0) stop 'DGESV'


! update wavefunction
         do l=0,lmax
         do isp=1,nspol
         do iprin=1,nprin(l+1,isp)
            call zero(nrad,psi(1,isp,iprin,l+1))
            jst=max(1,ids-idsx+1)
            jj=0
            do j=jst,ids
              jj=jj+1
              mj=mod(j-1,idsx)+1
              do k=1,nrad-3,4
              psi(k+0,isp,iprin,l+1)=psi(k+0,isp,iprin,l+1)+rdiis(jj)*(psid(k+0,isp,iprin,l+1,mj)-gradp(k+0,isp,iprin,l+1,mj))
              psi(k+1,isp,iprin,l+1)=psi(k+1,isp,iprin,l+1)+rdiis(jj)*(psid(k+1,isp,iprin,l+1,mj)-gradp(k+1,isp,iprin,l+1,mj))
              psi(k+2,isp,iprin,l+1)=psi(k+2,isp,iprin,l+1)+rdiis(jj)*(psid(k+2,isp,iprin,l+1,mj)-gradp(k+2,isp,iprin,l+1,mj))
              psi(k+3,isp,iprin,l+1)=psi(k+3,isp,iprin,l+1)+rdiis(jj)*(psid(k+3,isp,iprin,l+1,mj)-gradp(k+3,isp,iprin,l+1,mj))
              enddo
            enddo
         enddo
         enddo
         enddo

      endif

      if (cut_atom) then
        !write(6,*)"Zerotail after update"
        call zerotail(lmax,nspol,nprinx,nprin,nrad,ncut,psi)
      end if

      call cpu_time(t)
      time(4)=time(4)-t
      if (mod(ids,1).eq.0) & 
        write(6,'(1x,i3,3x,i2,2x,i6,2x,e11.4,1x,e20.10,(3(1x,e20.10)))') irel,ishift,ids,alpha,toten,toten-totenold,tres,diff
      call cpu_time(t)
      time(4)=time(4)+t

! orthogonalize wavefunctions
      if (idsx.eq.1 .or. (.not.pspcalcul)) then  ! same steepest descent shortcut condition
        midsp=1
      else
        midsp=mod(ids,idsx)+1
      endif
      call gsortho(nrad,nspol,nprinx,nprin,lmax,rw(1,2),wrky,psi,psid(1,1,1,1,midsp))
      !call checkorthogonalization(nrad,nspol,nprinx,nprin,lmax,rw(1,2),wrky,psid(1,1,1,1,midsp))
      !if (tres.lt.epstres) write(6,*)"Solution Converged"   !! Condition to Exit the SCF
      tresmin=min(tresmin,tres)
      if (tres.lt.epstres) goto 1100   !! Condition to Exit the SCF

2000 continue

    do l=0,lmax
    do isp=1,nspol
    do iprin=1,nprin(l+1,isp)
    do j=1,nrad-3,4
       psi(j+0,isp,iprin,l+1)=psid(j+0,isp,iprin,l+1,midsp)
       psi(j+1,isp,iprin,l+1)=psid(j+1,isp,iprin,l+1,midsp)
       psi(j+2,isp,iprin,l+1)=psid(j+2,isp,iprin,l+1,midsp)
       psi(j+3,isp,iprin,l+1)=psid(j+3,isp,iprin,l+1,midsp)
    enddo
    enddo
    enddo
    enddo

    !! Modifying the shift values
    do l=0,lmax
    do iprin=1,min(nprin(l+1,1),nprin(l+1,2))
       tt=-.5d0*(eval(iprin,1,l+1)+eval(iprin,2,l+1))
       tt=abs(tt)
       shift(iprin,l+1)=max(tt,1.d0)
       !write(6,*)"Before shift",iprin,l+1,shift(iprin,l+1)
       !if (tt.gt..25d0) shift(iprin,l+1)=max(tt,1.d0)
       !write(6,*)"After1  shift",iprin,l+1,shift(iprin,l+1)
    enddo
    do iprin=min(nprin(l+1,1),nprin(l+1,2))+1,max(nprin(l+1,1),nprin(l+1,2))
       tt=-(eval(iprin,1,l+1)+eval(iprin,2,l+1))
       tt=abs(tt)
       !if (tt.gt..25d0) shift(iprin,l+1)=max(tt,1.d0)
       shift(iprin,l+1)=max(tt,1.d0)   
    enddo
    enddo

    do l=0,lmax
       iprin=1
       do iprin=2,max(nprin(l+1,1),nprin(l+1,2))
          shift(iprin,l+1)=min(shift(iprin-1,l+1), shift(iprin,l+1))
          !shift(iprin,l+1)=max(shift(iprin,l+1),1.d0)
       enddo
    enddo

3000 continue
!! End of SCF

1100 continue
  call cpu_time(t)
  time(2)=time(2)+t

  itertot=itertot+(ishift-1)*nitermax+ids
  do l=0,lmax
  do isp=1,nspol
  do iprin=1,nprin(l+1,isp)
  do j=1,nrad-3,4
     psi(j+0,isp,iprin,l+1)=psid(j+0,isp,iprin,l+1,midsp)
     psi(j+1,isp,iprin,l+1)=psid(j+1,isp,iprin,l+1,midsp)
     psi(j+2,isp,iprin,l+1)=psid(j+2,isp,iprin,l+1,midsp)
     psi(j+3,isp,iprin,l+1)=psid(j+3,isp,iprin,l+1,midsp)
  enddo
  enddo
  enddo
  enddo

  diff   = 0.d0
  rmix   = 0.22d0

  if (rw(1,2)==0.d0) then
     damprel(1)=0.d0
     iradstart=2
  else
     iradstart=1
  endif

  do irad  = iradstart,nrad
    if (nspol == 1) then
      potnew =                       (potxc(irad,1) + pothart(irad))/rw(irad,2) + potloc(irad)
    else
      !potnew = 0.5d0*(potxc(irad,1) + potxc(irad,2) + pothart(irad))/rw(irad,2) + potloc(irad)  Bug discovered by Zhao
      potnew = (0.5d0*(potxc(irad,1) + potxc(irad,2)) + pothart(irad))/rw(irad,2) + potloc(irad) ! corrected line 
    end if
    
     damprel_new=c_sq/((2.d0*c_sq)-potnew)
     diff=diff+rw(irad,2)*(damprel_new-damprel(irad))**2
     damprel(irad) = rmix*damprel_new + (1.d0-rmix)*damprel(irad)
  end do
  diff=sqrt(diff)

  call intrapolate_nrad(damprel, damprel4, nrad)

9201 if (diff.lt.epsdamp) goto 9202

9202 continue

!
!     charge up to radius rcov or infinity
!
!!! Calculating Charge and the First Moment of Charge Density/ Electric Field
  if (lmax.gt.3) stop 'cannot calculate chrg'

  do l=0,lmax
  do isp=1,nspol
  do iprin=1,nprin(l+1,isp)
  do incov=1,ncovmax
     chrg(iprin,isp,l+1,incov)=0.d0
    !electric(iprin,isp,l+1,incov)=0.d0
  end do
  end do
  end do
  end do

  do isp=1,nspol
  do l=0,lmax
  do iprin=1,nprin(l+1,isp)
  do incov=1,ncov
     call detnp(nrad,rr,rcov(incov),np) !! Cutoff Value of Rcov 
     do j=1,np!-4
        chrg(iprin,isp,l+1,incov)=chrg(iprin,isp,l+1,incov)+ &
                   rw(j,2)*psi(j,isp,iprin,l+1)**2  !! INTEGRATION OF THE CHARGE  psi**2
       !electric(iprin,isp,l+1,incov)=electric(iprin,isp,l+1,incov)+ &
       !           rw(j,2)*rr(j)*psi(j,isp,iprin,l+1)**2 !! INTEGRATION OF r*psi**2
     end do
     !chrg(iprin,isp,l+1,incov)=chrg(iprin,isp,l+1,incov)+ &
     !        rw(4,2)*(psi(np-3,isp,iprin,l+1))**2
     !chrg(iprin,isp,l+1,incov)=chrg(iprin,isp,l+1,incov)+ &
     !        rw(3,2)*(psi(np-2,isp,iprin,l+1))**2
     !chrg(iprin,isp,l+1,incov)=chrg(iprin,isp,l+1,incov)+ &
     !        rw(2,2)*(psi(np-1,isp,iprin,l+1))**2
     !chrg(iprin,isp,l+1,incov)=chrg(iprin,isp,l+1,incov)+ &
     !        rw(1,2)*(psi(np,isp,iprin,l+1))**2
  end do
  end do
  end do
  end do

!! OPTIONAL
! calculate kinetic energy (just for checkin the virial theorem)
  potloc=0.d0  
  ekin=0.d0  !! Total KE
  eext=0.d0  !! PE due to el-ion interaction
  ekin_orb=0.d0 !! KE of individual orbital
  epot_orb=0.d0 !! PE of individual orbital due to el-ion interaction

  do l=0,lmax
  do isp=1,nspol
  do iprin=1,nprin(l+1,isp)
     kinonly=.true. !! Kinetic Energy Only
     call applyzorakinpot(nrad,rr,rw,rw4,kinonly,0.d0,l,potloc,damprel,damprel4,psi(1,isp,iprin,l+1),grad(1,isp,iprin,l+1))
     kinetic(:,isp,iprin,l+1) = grad(:,isp,iprin,l+1)
     ekin_orb(iprin,isp,l+1)=DDOT(nrad,psi(1,isp,iprin,l+1),1,grad(1,isp,iprin,l+1),1)
     ekin=ekin+occup(iprin,isp,l+1)*ekin_orb(iprin,isp,l+1)

     kinonly=.false. !! KE + Local Potential of the Nucleus
     call applyzorakinpot(nrad,rr,rw,rw4,kinonly,0.d0,l,Vnuc,damprel,damprel4,psi(1,isp,iprin,l+1),grad(1,isp,iprin,l+1))
     epot_orb(iprin,isp,l+1)=DDOT(nrad,psi(1,isp,iprin,l+1),1,grad(1,isp,iprin,l+1),1)
     epot_orb(iprin,isp,l+1)=epot_orb(iprin,isp,l+1)-ekin_orb(iprin,isp,l+1) !
     ! Only P.E. contribution

     eext=eext+occup(iprin,isp,l+1)*epot_orb(iprin,isp,l+1)
  enddo
  enddo
  enddo
! write(*,*) 'kinetic energy, ratio ', ekin,toten/ekin

  write(6,*)"     "
  write(6,*)"==========================================="
  write(6,*)"Zero and First Moment of Rho of Orbitals   "
  write(6,*)"Zero  Moment:Charge Density                "
  write(6,*)"First Moment:Contribution of Electric Field"
  write(6,*)"==========================================="
  write(6,*) 'charge  (rcov) = int_0^rcov   psi^2 r^2 dr'
 !write(6,*) 'electric(rcov) = int_0^rcov r*psi^2 r^2 dr'
  write(6,21)      
21 format(/,' nl   s    occ',5x,'rcov',4x,'charge(rcov)')
  do iorb=1,norb
    iprin=no(iorb)
    ll=lo(iorb)
    ss=so(iorb)
    isp=1
    if (ss.lt.0.1d0)isp=2
    write(6,31) noae(iorb),il(lo(iorb)+1),so(iorb),zo(iorb),rcov(1),chrg(iprin,isp,ll+1,1) !,electric(iprin,isp,ll+1,1)
    if (ncov.gt.1) then
       do incov=2,ncov
        write(6,'(16x,f8.3,e15.7,e15.7)')rcov(incov),chrg(iprin,isp,ll+1,incov) !,electric(iprin,isp,ll+1,incov)
       end do
    end if
  end do
31 format(1x,i1,a1,1x,f4.1,2(f8.3),e15.7,e15.7)

  write(6,*)"        "
  write(6,88) nameat
88 format(a3,25h output data for orbitals,/,1x,27('-'),//,17h nl    s      occ,9x,'eigenvalue',4x,14hkinetic energy,6x,'pot energy',6x,'residue',12x,'shift'/)
  do iorb=1,norb
     iprin=no(iorb)
     ll=lo(iorb)
     ss=so(iorb)
     isp=1
     if (ss.lt.0.1d0)isp=2
     zz=zo(iorb)
     ev=eval(iprin,isp,ll+1)
     ek=ekin_orb(iprin,isp,ll+1)
     ep=epot_orb(iprin,isp,ll+1)
     res1=residue(iprin,isp,ll+1)   
     write(6,95) noae(iorb),il(lo(iorb)+1),so(iorb),zo(iorb),ev,ek,ep,res1,shift(iprin,ll+1)
 95 format(1x,i1,a1,f6.1,f10.4,3f17.8,e17.8,e20.10)
  end do
  ekin2=evsum-eext-2.d0*ehart-vxc-hfmix*ehf ! KE from eigen value
  epot=eext+ehart+enexc+ehf*hfmix ! Potential energy
  toten2=evsum-ehart+(enexc-vxc)-ehf*hfmix ! Totel Ene from eigenvalue

  write(6,65) evsum,ekin,eext,ehart,vxc,enexc,hfmix*ehf, &
  ekin2,epot,toten2,toten
65 format(//,15h total energies,/,1x,14('-'),/,  &
    /,28h sum of eigenvalues        =,e25.15,  &
    /,28h kinetic energy from ek    =,e25.15,  &
    /,28h el-ion interaction energy =,e25.15,  &
    /,28h el-el  interaction energy =,e25.15,  &
    /,28h vxc    correction         =,e25.15,  &
    /,28h exchange + corr energy    =,e25.15,  &
    /,28h Hatree-Fock energy        =,e25.15,  &
    /,28h kinetic energy from ev    =,e25.15,  &
    /,28h potential energy          =,e25.15,/,1x,45('-'),  &
    /,28h total energy from ev      =,e25.15,  &
    /,28h total energy              =,e25.15)
   
  if (tres.gt.epstres) then
     write(*,*) "Total residue > epstres in zorasolver"
     stop
  end if

end subroutine


