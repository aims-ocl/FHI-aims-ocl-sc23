c Header:
c**********************************************************************c
c  MacDonald et al. relativistic correction to lda exchange
c  J Phys C 12 , 2977 (1979)
c  mode .... 0 retrun right away
c            1 calculate relativistic correction
c  d ....... density
c  dex ..... correction factor ex^rlda = dex * ex^lda 
c  dvx ..... correction factor vx^rlda = dvx * vx^lda
c**********************************************************************
c
      subroutine relxc(mode,d,dex,dvx)
c
      implicit real*8 (a-h,o-z)
      integer mode
      parameter(a=1.d0/44.3d0)
      parameter(thrd=0.333333333333333d0)
c
      dvx=1.d0
      dex=1.d0
      if(d .gt. 1d-15 .and. mode .ne. 0) then
        beta=a*d**thrd
        gamma=sqrt(1.d0+beta*beta)
        f1=log(beta+gamma)
        dvx=0.5d0*(3.d0*f1/(beta*gamma)-1.d0)
        dex=1.d0-1.5d0*(beta*gamma-f1)**2/beta**4
      endif
c
      return
      end
c
