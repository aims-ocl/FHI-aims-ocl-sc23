c  Header:
c**********************************************************************
c output electrostatic potential for spherical charge density
c 
c input
c    mmax        maximum grid index
c    r()         logarithmic radial grid
c    rho()       density array
c    zion        <control> == total charge for all-electron atom
c                          == 0  for neutral atom (Coulomb+Hartree potential)
c                          == number of electrons   (Hartree potential)
c    flag	     .true. : evaluate total energy contribution
c
c output (zion is control parameter)
c    eeel        electrostatic energy
c    vo()        electrostatic potential 
c
c original version by D.R. Hamann, gncpp
c**********************************************************************
c
      subroutine vestat(mmax,zion,eeel,r,rho,vo,flag)
c
      implicit none
      include  'parameter.h'
      logical  flag
      integer  mmax,i
      real*8   zion,eeel,tv,al,aii,dmelm
      real*8   r(mmax),rv(mmax),rvp(mmax),vo(mmax),rho(mmax)
      external aii,dmelm
c
      al=log(r(2)/r(1))
c
c integration for electrostatic potential
      do 60 i=1,mmax
   60   rvp(i)=rho(i)*al*r(i)**3
c
      rv(mmax)=zion
      rv(mmax-1)=zion
      rv(mmax-2)=zion
c
      do 70 i=mmax-2,2,-1
   70   rv(i-1)=rv(i)+aii(rvp,i,mmax)
c
      do 80 i=1,mmax
   80   rvp(i)=rho(i)*al*r(i)**2
c
      tv=0.0d0
      do 90 i=mmax-2,2,-1
        tv=tv+aii(rvp,i,mmax)
   90   rv(i-1)=rv(i-1)-r(i-1)*tv
c
      do 100 i=1,mmax
 100    vo(i)=rv(i)/r(i)
c
c electron-electron interaction for total energy
      if(flag) eeel=dmelm(mmax,al,r,vo,rho)
c
      return
      end
ce

