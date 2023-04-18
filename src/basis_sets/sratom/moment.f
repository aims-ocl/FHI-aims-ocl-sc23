c Header:
c***********************************************************************
c calculate moment m of function fr on logarithmic radial mesh r
c input
c	m	m-th moment 
c	mmx	maximum grid index
c	al	log(r(n+1)/r(n))
c	fnorm	normalization factor
c	r()	logarithmic mesh
c	fr()	function array
c***********************************************************************
c	
      real*8 function fmom(m,mmx,al,fnorm,r,fr)
c
      implicit none
      include  'parameter.h'
      integer  m,mmx,max,i
      real*8   al,fnorm,esum,osum,fm
      real*8   fr(mmx),r(mmx)
c
      fm(i)=fr(i)*r(i)*r(i)*r(i)*r(i)**m
c
c integration by simpson rule
      max=mmx
      if(mod(mmx,2) .eq. 0) then
c       write(ie,*) 
c    1     '& fmom - got even no. of mesh points max, assume mmx-1'
        max=max-1
      endif
      esum=0.d0
      do i=2,max-1,2
        esum=esum+fm(i)
      enddo
      osum=0.d0
      do i=3,max-2,2
        osum=osum+fm(i)
      enddo
      fmom=al*(4.0*esum+2.0*osum+fm(1)+fm(max))/3.d0
      fmom=(fmom+0.5*fm(1))/fnorm

      return
      end
c
c**********************************************************************
c output matrix element on radial logarithmic mesh
c mmx		maximum grid index
c al		mesh spacing
c r()		logarithmic mesh
c v()		potential
c usq()	density weight
c**********************************************************************
c
      real*8 function dmelm(mmx,al,r,v,usq)
c
      implicit none
      include  'parameter.h'
      integer  mmx,max,i
      real*8   al,esum,osum,fm
      real*8   v(mmx),r(mmx),usq(mmx)
c
      fm(i)=v(i)*usq(i)*r(i)*r(i)
c
c by simpson's rule
      max=mmx
      if(mod(mmx,2) .eq. 0) then
c       write(ie,*) 
c    1     '& dmelm - got even no. of mesh points max, assume mmx-1'
        max=max-1
      endif
      esum=0.d0
      do i=2,max-1,2
        esum=esum+r(i)*fm(i)
      enddo
      osum=0.d0
      do i=3,max-2,2
        osum=osum+r(i)*fm(i)
      enddo
      dmelm=(al*(4.0*esum+2.0*osum+r(max)*fm(max)+r(1)*fm(1)))/3.d0
      dmelm=dmelm+0.5*fm(1)*r(1)

      return
      end
c
c more matrix elements w/ Simpson rule
c logarithmic mesh 
c
      real*8 function gltfmv(mmx,al,r,ga,g,gb)
      implicit none
      include  'parameter.h'
      integer  mmx,i,max
      real*8   al,osum,esum
      real*8   r(mmx),ga(mmx),g(mmx),gb(mmx)

      max=mmx
      if(mod(mmx,2) .eq. 0) then
        write(ie,*) 
     1     '& gltfmv - got even no. of mesh points mmx, assume mmx-1'
        max=max-1
      endif
c
c simpson's rule
      esum=0.d0
      do i=2,max-1,2
        esum=esum+r(i)*ga(i)*g(i)*gb(i)
      enddo
      osum=0.d0
      do i=3,max-2,2
        osum=osum+r(i)*ga(i)*g(i)*gb(i)
      enddo
      gltfmv=(al*(4.0*esum+2.0*osum+
     1            ga(1)*g(1)*gb(1)*r(1)+ga(max)*g(max)*gb(max)*r(max)
     1           ))/3.d0
     1  +0.5*ga(1)*g(1)*gb(1)*r(1)
 
      return
      end
c
