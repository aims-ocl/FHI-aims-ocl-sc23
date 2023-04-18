      function pole_f (p,freq)
        implicit none
        complex*16 :: pole_f
        real*8, intent(in) :: freq
        real*8, intent(in) :: p
          pole_f =  1.d0/(p+(0.d0,1.d0)*freq)
      end function pole_f


      function pole_f_re (p,freq)
        implicit none
        real*8 :: pole_f_re
        real*8, intent(in) :: freq
        real*8, intent(in) :: p
          pole_f_re = (dabs(p))/(p*p+freq*freq)
      end function pole_f_re

      function pole_f_cc (p,freq)
        implicit none
        complex*16 :: pole_f_cc
        real*8, intent(in) :: freq
        real*8, intent(in) :: p
        pole_f_cc = 1.d0/(p-(0.d0,1.d0)*freq)
      end function pole_f_cc

      function exp_f (p,time)
        implicit none
        real*8 :: exp_f
        real*8, intent(in) :: time
        real*8, intent(in) :: p
          exp_f = exp (- p * time)
      end function exp_f


