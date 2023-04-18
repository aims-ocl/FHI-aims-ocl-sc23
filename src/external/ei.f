      
      Subroutine ei(x,eix)
c     use localorb_io
      implicit none
      INTEGER MAXIT
      REAL*8 eix,x,EPS,EULER,FPMIN
      PARAMETER (EPS=6.e-8,EULER=.57721566,MAXIT=100,FPMIN=1.e-30)
      INTEGER k
      Real*8 fact,prev,sum, term
      if(x.le.0.0d0) then 
         call E1(abs(x),eix)
         eix=-eix
         return
         write(6,*) "bad argument in exponetial integral function"
c     write(info_str,'(A)') "bad argument in exponetial integral function"
c     call localorb_info ( info_str )
         stop
      endif
      write(*,*)"hier sind wir falsch"
      if (x.lt.FPMIN) then
         eix=dlog(x) +Euler
      else if(x.le.-dlog(EPS))then
         sum=0.0
         fact=1.0
         do  k=1,MAXIT
            fact=fact*x/k
            term=fact/k
            sum=sum+term
            if(term.lt.EPS*sum)goto 1
         enddo
         write(6,*) "series faild in ei"
         stop
c     write(info_str,'(A)') "series faild in ei"
c     call localorb_info ( info_str )
         
 1       eix=sum+dlog(x)+EULER
      else
         sum=0.
         term=1.
         do k=1, MAXIT
            prev=term
            term=term*k/x
            if(term.lt.EPS) goto2
            if (term.lt.prev) then
               sum=sum+term
            else
               sum=sum-prev
               goto 2
            endif
         enddo
 2       eix=dexp(x)*(1.+sum)/x
      endif
      end
      
      
      subroutine E1(x,expint)
      implicit none
      integer n,maxit
      real*8 expint,x,eps,fpmin,euler
      parameter (maxit=1000,eps=1.e-7,fpmin=1.e-30,euler=.5772156649)
      integer i,ii,nm1
      real*8 a,b,c,d,del,fact,h,psi
      n=1
      nm1=n-1
      if(n.lt.0 .or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1))) then
         write(6,*) "bad argument in E1",x
         stop
c     pause 'bad argument in expint'
      else if(n.eq.0)then
         expint=dexp(-x)/x
      else if(x.eq.0.)then
         expint=dble(1.0d0/nm1)
      else if(x.gt.1.0)then
         b=x+n
         c=dble(1.0d0/fpmin)
         d=dble(1.0d0/b)
         h=d
         do i=1,maxit
            a=-i*(nm1+i)
            b=b+2.0
            d=1.0/(a*d+b)
            c=b+dble(a/c)
            del=c*d
            h=h*del
            if(abs(del-1.0d0).lt.eps)then
               expint=h*dexp(-x)
               return
            endif
         enddo
         write(6,*) "continued fraction failed in E1"      
         stop
      else

         if(nm1.ne.0)then
            expint=dble(1.0d0/nm1)
         else
            expint=-dlog(x)-euler
         endif

         fact=1.0d0

         do i=1,maxit,1
            
            fact=-fact*dble(x/i)

c            if(i.ne.nm1)then
c               write(*,*)"fact", fact
            del = dble(fact/(i-nm1))
c
c            else psi=-euler
c
c               do ii=1,nm1
c               
c                  psi=psi+dble(1.0d0/ii)
c
c               enddo
c               write(*,*)"hier sind wir falsch",i,nm1
c               del = fact*(-dlog(x)+psi)

c            endif
c
            expint = expint - del
c            write(*,*)"del", del
            if(abs(del).lt.(abs(expint)*eps)) return
         enddo
         write(6,*) "series failed in E1"      
         stop
      endif 
      return
      end
