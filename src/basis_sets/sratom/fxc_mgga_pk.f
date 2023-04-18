c enhancement factor over LDA exchnage of the meta gga of perdew & kurth
c rho .... spin density
c rhop ... |spin density gradient| 
c tkin ... kinetic energy spin density 
      real*8 function fx_xvar(rho,rhop,tkin)

      implicit none

      real*8   rho,rhop,tkin
      real*8   pvar,qvar,xvar
      real*8   c4i3,c5i3,c10i81,kappa,dpar,eps
      parameter (c4i3    = 4.d0/3.d0)
      parameter (c5i3    = 5.d0/3.d0)
      parameter (c10i81  = 1.d1/81.d0)
      parameter (dpar   = 0.113d0)
      parameter (kappa  = 0.804d0)
      parameter (eps    = 1.d-12)

      xvar = 0.d0
      if(rho < eps) return

      pvar = 0.261211729852336046d-1 * (rhop/rho**c4i3)**2  ! eq 10
    
      qvar = 0.156727037911401634 * tkin/rho**c5i3          ! eq 12
     1       - 0.45d0 - pvar/12.d0

      xvar = c10i81 * pvar + 146d0/2025d0 * qvar**2         ! eq 14
     1       - 73d0/405d0 * pvar*qvar
     1       + ( dpar + c10i81**2/kappa) * pvar**2

c effective scaled gradient
      fx_xvar = sqrt(xvar/0.21951d0)

      return
      end

      real*8 function fx_mgga_pk(rho,rhop,tkin)

      implicit none

      real*8   rho,rhop,tkin
      real*8   pvar,qvar,xvar
      real*8   c4i3,c5i3,c10i81,kappa,dpar,eps
      parameter (c4i3    = 4.d0/3.d0)
      parameter (c5i3    = 5.d0/3.d0)
      parameter (c10i81  = 1.d1/81.d0)
      parameter (dpar   = 0.113d0)
      parameter (kappa  = 0.804d0)
      parameter (eps    = 1.d-12)

      fx_mgga_pk = 1.d0 
      if(rho < eps) return

      pvar = 0.261211729852336046d-1 * (rhop/rho**c4i3)**2	! eq 10
     
      qvar = 0.156727037911401634 * tkin/rho**c5i3		! eq 12
     1       - 0.45d0 - pvar/12.d0 

      xvar = c10i81 * pvar + 146d0/2025d0 * qvar**2 		! eq 14
     1       - 73d0/405d0 * pvar*qvar
     1       + ( dpar + c10i81**2/kappa) * pvar**2

      fx_mgga_pk = 1.d0 + kappa - kappa/(1.d0 + xvar/kappa) 	! eq 13
 
      return
      end
      
 
c enhancement factor over LDA/GGA correlation of the meta gga of perdew & kurth
c rho_up ..... spin up density
c rho_dn ..... dto. down
c rhop_up .... |spin up density gradient| 
c rhop_dn .... dto. down
c tkin_up .... kinetic energy spin up density
c tkin_dn .... dto. down
c ec ......... correlation energy per electron spin-unpolarized
c ec_up ...... dto. spin up
c ec_dn ...... dto. spin down
c
      real*8 function fc_mgga_pk
     1   (rho_up,rhop_up,tkin_up,rho_dn,rhop_dn,tkin_dn,ec,ec_up,ec_dn)

      implicit none

      real*8   rho_up,rhop_up,tkin_up,tvw_up
      real*8   rho_dn,rhop_dn,tkin_dn,tvw_dn
      real*8   ec,ec_up,ec_dn
      real*8   cpar,eps
      parameter (cpar = 0.53d0)
      parameter (eps  = 1.d-12)
 
      fc_mgga_pk = 1.d0 
      if(rho_up < eps) return

c von Weizsaecker kinetic energy density
      tvw_up = 0.125d0 * rhop_up**2/rho_up 
      tvw_dn = 0.125d0 * rhop_dn**2/rho_dn 

c eq 15 correlation enhancement factor
      fc_mgga_pk = 1.d0 
     1       + cpar*( (tvw_up + tvw_dn)/(tkin_up + tkin_dn) )**2
     1       - (1.d0 + cpar) * (
     1         (tvw_up/tkin_up)**2 * rho_up * ec_up
     1       + (tvw_dn/tkin_dn)**2 * rho_dn * ec_dn
     1                         ) / ( (rho_up + rho_dn) *ec) 

      return
      end
