c Header:
c array manipulations
c
c b = a
      subroutine dcpv(mx,ml,mh,a,b)
      real*8 a, b 
      integer mx, ml, mh
      dimension a(mx),b(mx)
      integer i
      do i=ml,mh
        b(i)=a(i)
      enddo
      return
      end

c
c c = a + b
      subroutine dadv(mx,ml,mh,a,b,c)
      real*8 a, b, c 
      integer mx, ml, mh
      dimension a(mx),b(mx),c(mx)
      integer i
      do i=ml,mh
        c(i)=a(i)+b(i)
      enddo
      return
      end
