c $Header:$
c parameter file for GNCPP routines 
c rmin..........nuclear charge times minimum radius on grid
c rmax..........maximum radius of grid
c ainc..........mesh increment ainc=r(i+1)/r(i)
c
c number of gridpoints = log(rmax/rmin * z) / log(ainc)
c
      integer   itmx
      real*8    rmin,rmax,ainc,ry2,epsae,epspp
c coarse mesh not sufficient for GGA
c     parameter (rmin  = 0.01d0)
c     parameter (rmax  = 8.00d1)
c     parameter (ainc  = 1.05d0)

c fine mesh usually overdoing somwhat 
c     parameter (rmin  = 0.00625d0)
c     parameter (rmax  = 8.00d1)
c     parameter (ainc  = 1.0123d0)

c efficient one
cvb      parameter (rmin  = 0.00625d0)
cvb      parameter (rmax  = 8.00d1)
cvb      parameter (ainc  = 1.0247d0)

c a.u.
      parameter (ry2   = 27.2116d0)
c epsae.........relative accuracy to converge all-electron eigenvalues
c epspp.........relative accuracy to converge pseudo-wavefunction eigenvalue 
c               during construction of potentials
c itmx..........maximum number of self-consistency iterations
cvb orig      parameter (epsae = 1.d-12)
      parameter (epsae = 1.d-12)
      parameter (epspp = 1.d-12)
      parameter (itmx  = 100)
c
