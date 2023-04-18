c $Header:$
c parameter file for pseudopotential routines 
c mx .......... radial mesh index
c ms .......... states index

cvb all the fixed dimensions should be out, retained here only
c   in case I forgot anything.
c      include "../param.f"

c      integer   mx,ms,mzmx
cvb orig      parameter (mx=1200)
cvb orig      parameter (ms=80)
c      parameter (mx=max_free_grid)
c      parameter (ms=max_shells)
c      parameter (mzmx=120)
c ie .......... error output unit
c iogncpp ..... regular gncpp output unit 
c iopslp ...... regular pslp output unit 
      integer   ie,iofhipp,iopslp
      parameter (ie=6)
      parameter (iofhipp=23)
      parameter (iopslp=6)
