      subroutine xc_func(key, open, 
     >                   rhoc, rhoo, 
     >                   sigmacc, sigmaco, sigmaoo, tauc, tauo,
     >                   upsilonc, upsilono,
     >                   zk, vrhoc, vrhoo,
     >                   vsigmacc, vsigmaco, vsigmaoo, vtauc, vtauo,
     >                   vupsilonc, vupsilono)

      implicit real*8 (r,s,t,u,z,v)
      logical fderiv, open
      integer igrad, npt
      character(len=100) name
      character(len=*) key
      dimension rhoc(1), rhoo(1)
      dimension sigmacc(1), sigmaco(1), sigmaoo(1)
      dimension tauc(1), tauo(1), upsilonc(1), upsilono(1)
      dimension zk(1), vrhoc(1), vrhoo(1)
      dimension vsigmacc(1), vsigmaco(1), vsigmaoo(1)
      dimension vtauc(1), vtauo(1), vupsilonc(1), vupsilono(1)

      npt = 1
      fderiv = .true.

      if (.false.) then
c:DIRACXcallstart
        else if (key.eq.'DIRACX') then
        call dftacg_diracx(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,tauc,tauo,
     >                   upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,vtauc,vtauo,
     >                   vupsilonc,vupsilono)
c:DIRACXcallend
c:PBECcallstart
        else if (key.eq.'PBEC') then
        call dftacg_pbec(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,tauc,tauo,
     >                   upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,vtauc,vtauo,
     >                   vupsilonc,vupsilono)
c:PBECcallend
c:PBEXcallstart
        else if (key.eq.'PBEX') then
        call dftacg_pbex(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,tauc,tauo,
     >                   upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,vtauc,vtauo,
     >                   vupsilonc,vupsilono)
c:PBEXcallend
c:PW92Ccallstart
        else if (key.eq.'PW92C') then
        call dftacg_pw92c(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,tauc,tauo,
     >                   upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,vtauc,vtauo,
     >                   vupsilonc,vupsilono)
c:PW92Ccallend
c:SCANCcallstart
        else if (key.eq.'SCANC') then
        call dftacg_scanc(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,tauc,tauo,
     >                   upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,vtauc,vtauo,
     >                   vupsilonc,vupsilono)
c:SCANCcallend
c:SCANXcallstart
        else if (key.eq.'SCANX') then
        call dftacg_scanx(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,tauc,tauo,
     >                   upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,vtauc,vtauo,
     >                   vupsilonc,vupsilono)
c:SCANXcallend
c:TPSSCcallstart
        else if (key.eq.'TPSSC') then
        call dftacg_tpssc(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,tauc,tauo,
     >                   upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,vtauc,vtauo,
     >                   vupsilonc,vupsilono)
c:TPSSCcallend
c:TPSSXcallstart
        else if (key.eq.'TPSSX') then
        call dftacg_tpssx(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,tauc,tauo,
     >                   upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,vtauc,vtauo,
     >                   vupsilonc,vupsilono)
c:TPSSXcallend
c:dfauto
      end if

      end
c:DIRACXsubrstart

c    Generated: Mon Jun 20 16:05:12 CEST 2016

      subroutine dftacg_diracx
     > (name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   tauc,tauo,upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,
     >                   vtauc,vtauo,vupsilonc,vupsilono)
      implicit double precision (a-h,o-z)
      logical fderiv,open
      integer igrad,npt
      character*(*) name
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)
      double precision tauc(*),tauo(*)
      double precision upsilonc(*),upsilono(*)
      double precision zk(*),vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)
      double precision vtauc(*),vtauo(*)
      double precision vupsilonc(*),vupsilono(*)
      include "common/cdft"
      include "common/tapes"
      parameter(tol=1d-12)
      pi=acos(-1d0)
      name='Automatically generated DIRACX'
      igrad=0
       if(open) then
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      t1 = pi ** 2
      t3 = (t1 * rhob) ** (0.1D1 / 0.3D1)
      t5 = 0.1D1 / pi
      zk(i) = -0.136284044462410474417D1 * rhob * t3 * t5
      t10 = 0.681420222312052372084D0 * t3 * t5
      t11 = t3 ** 2
      t15 = 0.227140074104017457361D0 * rhob / t11 * pi
      vrhoc(i) = vrhoc(i) - t10 - t15
      vrhoo(i) = vrhoo(i) + t10 + t15

               elseif(rhob.lt.tol) then
             rho = rhoa
      t1 = pi ** 2
      t3 = (t1 * rhoa) ** (0.1D1 / 0.3D1)
      t5 = 0.1D1 / pi
      zk(i) = -0.136284044462410474417D1 * rhoa * t3 * t5
      t10 = 0.681420222312052372084D0 * t3 * t5
      t11 = t3 ** 2
      t15 = 0.227140074104017457361D0 * rhoa / t11 * pi
      vrhoc(i) = vrhoc(i) - t10 - t15
      vrhoo(i) = vrhoo(i) - t10 - t15

               else
             rho = rhoa + rhob
      t1 = pi ** 2
      t3 = (t1 * rhoa) ** (0.1D1 / 0.3D1)
      t5 = 0.1D1 / pi
      t9 = (t1 * rhob) ** (0.1D1 / 0.3D1)
      zk(i) = -0.136284044462410474417D1 * rhoa * t3 * t5 - 0.1362840444
     >62410474417D1 * rhob * t9 * t5
      t15 = 0.681420222312052372084D0 * t3 * t5
      t16 = t3 ** 2
      t20 = 0.227140074104017457361D0 * rhoa / t16 * pi
      t22 = 0.681420222312052372084D0 * t9 * t5
      t23 = t9 ** 2
      t27 = 0.227140074104017457361D0 * rhob / t23 * pi
      vrhoc(i) = vrhoc(i) - t15 - t20 - t22 - t27
      vrhoo(i) = vrhoo(i) - t15 - t20 + t22 + t27

               endif
             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      t1 = pi ** 2
      t3 = (t1 * rhob) ** (0.1D1 / 0.3D1)
      zk(i) = -0.136284044462410474417D1 * rhob * t3 / pi

               elseif(rhob.lt.tol) then
             rho = rhoa
      t1 = pi ** 2
      t3 = (t1 * rhoa) ** (0.1D1 / 0.3D1)
      zk(i) = -0.136284044462410474417D1 * rhoa * t3 / pi

               else
             rho = rhoa + rhob
      t1 = pi ** 2
      t3 = (t1 * rhoa) ** (0.1D1 / 0.3D1)
      t5 = 0.1D1 / pi
      t9 = (t1 * rhob) ** (0.1D1 / 0.3D1)
      zk(i) = -0.136284044462410474417D1 * rhoa * t3 * t5 - 0.1362840444
     >62410474417D1 * rhob * t9 * t5

               endif
             endif
           enddo
         endif
       else
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      t3 = pi ** 2
      t5 = (t3 * rhoa) ** (0.1D1 / 0.3D1)
      t7 = 0.1D1 / pi
      t11 = (t3 * rhob) ** (0.1D1 / 0.3D1)
      zk(i) = -0.136284044462410474417D1 * rhoa * t5 * t7 - 0.1362840444
     >62410474417D1 * rhob * t11 * t7
      t18 = t5 ** 2
      t25 = t11 ** 2
      vrhoc(i) = vrhoc(i) - 0.681420222312052372084D0 * t5 * t7 - 0.2271
     >40074104017457361D0 * rhoa / t18 * pi - 0.681420222312052372084D0
     >* t11 * t7 - 0.227140074104017457361D0 * rhob / t25 * pi

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      t3 = pi ** 2
      t5 = (t3 * rhoa) ** (0.1D1 / 0.3D1)
      t7 = 0.1D1 / pi
      t11 = (t3 * rhob) ** (0.1D1 / 0.3D1)
      zk(i) = -0.136284044462410474417D1 * rhoa * t5 * t7 - 0.1362840444
     >62410474417D1 * rhob * t11 * t7

             endif
           enddo
         endif
       endif

      return
      end

c:DIRACXsubrend
c:PBECsubrstart

c    Generated: Mon Jun 20 16:05:13 CEST 2016

      subroutine dftacg_pbec
     > (name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   tauc,tauo,upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,
     >                   vtauc,vtauo,vupsilonc,vupsilono)
      implicit double precision (a-h,o-z)
      logical fderiv,open
      integer igrad,npt
      character*(*) name
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)
      double precision tauc(*),tauo(*)
      double precision upsilonc(*),upsilono(*)
      double precision zk(*),vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)
      double precision vtauc(*),vtauo(*)
      double precision vupsilonc(*),vupsilono(*)
      include "common/cdft"
      include "common/tapes"
      parameter(tol=1d-12)
      pi=acos(-1d0)
      name='Automatically generated PBEC'
      igrad=1
       if(open) then
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      sigmabb = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) - 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmabb
      t8 = 0.1D1 / pi
      t10 = t8 / rhob
      t11 = t10 ** (0.1D1 / 0.3D1)
      t13 = 0.1D1 + 0.186690969707574028554D0 * t11
      t14 = t10 ** (0.1D1 / 0.6D1)
      t17 = sqrt(t10)
      t19 = t11 ** 2
      t21 = 0.134579137143944477912D2 * t14 + 0.563098414909787598194D1 
     >* t11 + 0.291521471421917737271D1 * t17 + 0.516066464547863440989D
     >0 * t19
      t24 = 0.1D1 + 0.321639589973850701335D2 / t21
      t25 = log(t24)
      t26 = t13 * t25
      t28 = pi ** 2
      t29 = 0.1D1 / t28
      t30 = t28 * pi
      t32 = t28 * rhob
      t33 = t32 ** (0.1D1 / 0.3D1)
      t34 = 0.1D1 / t33
      t35 = t30 * sigmabb * t34
      t36 = rhob ** 2
      t37 = 0.1D1 / t36
      t40 = exp(0.202642426794280972788D0 * t26 * t28)
      t41 = t40 - 0.1D1
      t42 = 0.1D1 / t41
      t43 = t30 * t42
      t44 = sigmabb * t34
      t47 = 0.149583857013095144140D-1 * t43 * t44 * t37
      t48 = 0.1D1 + t47
      t49 = t37 * t48
      t50 = t28 ** 2
      t51 = t50 * t28
      t52 = t41 ** 2
      t53 = 0.1D1 / t52
      t54 = t51 * t53
      t55 = sigmabb ** 2
      t56 = t33 ** 2
      t57 = 0.1D1 / t56
      t58 = t55 * t57
      t59 = t36 ** 2
      t60 = 0.1D1 / t59
      t64 = 0.1D1 + t47 + 0.223753302789140933370D-3 * t54 * t58 * t60
      t65 = 0.1D1 / t64
      t66 = t49 * t65
      t69 = 0.1D1 + 0.149583857013095144140D-1 * t35 * t66
      t70 = log(t69)
      t71 = t29 * t70
      zk(i) = rhob * (-0.3109070D-1 * t26 + 0.153426409720027345292D0 * 
     >t71)
      t75 = 0.155453500000000000000D-1 * t26
      t76 = 0.767132048600136726458D-1 * t71
      t77 = 0.1D1 / t19
      t78 = t77 * t8
      t79 = t37 * t25
      t82 = t21 ** 2
      t84 = t13 / t82
      t85 = t14 ** 2
      t86 = t85 ** 2
      t104 = (-0.224298561906574129855D1 / t86 / t14 * t8 * t37 - 0.1876
     >99471636595866065D1 * t78 * t37 - 0.145760735710958868637D1 / t17
     >* t8 * t37 - 0.344044309698575627326D0 / t11 * t8 * t37) / t24
      t107 = t50 * pi
      t110 = 0.1D1 / t33 / t32
      t115 = 0.1D1 / t36 / rhob
      t122 = t34 * t37
      t130 = (-0.126105037207067989169D-1 * t77 * pi * t79 - 0.651778270
     >654185890920D1 * t84 * t104 * t28) * t40
      t133 = 0.149583857013095144140D-1 * t30 * t53 * sigmabb * t122 * t
     >130
      t138 = 0.498612856710317147135D-2 * t107 * t42 * sigmabb * t110 * 
     >t37
      t141 = 0.299167714026190288281D-1 * t43 * t44 * t115
      t147 = t64 ** 2
      t148 = 0.1D1 / t147
      t157 = t50 ** 2
      t177 = 0.1D1 / t69
      t182 = 0.500000000000000000000D0 * rhob * (0.193478431062909061652
     >D-2 * t78 * t79 + 0.100000000000000000000D1 * t84 * t104 + 0.15342
     >6409720027345292D0 * t29 * (-0.498612856710317147135D-2 * t107 * s
     >igmabb * t110 * t66 - 0.299167714026190288281D-1 * t35 * t115 * t4
     >8 * t65 + 0.149583857013095144140D-1 * t35 * t37 * (-t133 - t138 -
     > t141) * t65 - 0.149583857013095144140D-1 * t35 * t49 * t148 * (-t
     >133 - t138 - t141 - 0.447506605578281866741D-3 * t51 / t52 / t41 *
     > t55 * t57 * t60 * t130 - 0.149168868526093955580D-3 * t157 * t53
     >* t55 / t56 / t32 * t60 - 0.895013211156563733482D-3 * t54 * t58 /
     > t59 / rhob)) * t177)
      vrhoc(i) = vrhoc(i) - t75 + t76 + t182
      vrhoo(i) = vrhoo(i) + t75 - t76 - t182
      t208 = rhob * t29 * (0.149583857013095144140D-1 * t30 * t34 * t66 
     >+ 0.223753302789140933370D-3 * t51 * sigmabb * t57 * t60 * t42 * t
     >65 - 0.149583857013095144140D-1 * t35 * t49 * t148 * (0.1495838570
     >13095144140D-1 * t43 * t122 + 0.447506605578281866741D-3 * t54 * s
     >igmabb * t57 * t60)) * t177
      t209 = 0.383566024300068363229D-1 * t208
      vsigmacc(i) = vsigmacc(i) + t209
      vsigmaco(i) = vsigmaco(i) - 0.767132048600136726458D-1 * t208
      vsigmaoo(i) = vsigmaoo(i) + t209

               elseif(rhob.lt.tol) then
             rho = rhoa
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) + 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmaaa
      t8 = 0.1D1 / pi
      t10 = t8 / rhoa
      t11 = t10 ** (0.1D1 / 0.3D1)
      t13 = 0.1D1 + 0.186690969707574028554D0 * t11
      t14 = t10 ** (0.1D1 / 0.6D1)
      t17 = sqrt(t10)
      t19 = t11 ** 2
      t21 = 0.134579137143944477912D2 * t14 + 0.563098414909787598194D1 
     >* t11 + 0.291521471421917737271D1 * t17 + 0.516066464547863440989D
     >0 * t19
      t24 = 0.1D1 + 0.321639589973850701335D2 / t21
      t25 = log(t24)
      t26 = t13 * t25
      t28 = pi ** 2
      t29 = 0.1D1 / t28
      t30 = t28 * pi
      t32 = t28 * rhoa
      t33 = t32 ** (0.1D1 / 0.3D1)
      t34 = 0.1D1 / t33
      t35 = t30 * sigmaaa * t34
      t36 = rhoa ** 2
      t37 = 0.1D1 / t36
      t40 = exp(0.202642426794280972788D0 * t26 * t28)
      t41 = t40 - 0.1D1
      t42 = 0.1D1 / t41
      t43 = t30 * t42
      t44 = sigmaaa * t34
      t47 = 0.149583857013095144140D-1 * t43 * t44 * t37
      t48 = 0.1D1 + t47
      t49 = t37 * t48
      t50 = t28 ** 2
      t51 = t50 * t28
      t52 = t41 ** 2
      t53 = 0.1D1 / t52
      t54 = t51 * t53
      t55 = sigmaaa ** 2
      t56 = t33 ** 2
      t57 = 0.1D1 / t56
      t58 = t55 * t57
      t59 = t36 ** 2
      t60 = 0.1D1 / t59
      t64 = 0.1D1 + t47 + 0.223753302789140933370D-3 * t54 * t58 * t60
      t65 = 0.1D1 / t64
      t66 = t49 * t65
      t69 = 0.1D1 + 0.149583857013095144140D-1 * t35 * t66
      t70 = log(t69)
      t71 = t29 * t70
      zk(i) = rhoa * (-0.3109070D-1 * t26 + 0.153426409720027345292D0 * 
     >t71)
      t75 = 0.155453500000000000000D-1 * t26
      t76 = 0.767132048600136726458D-1 * t71
      t77 = 0.1D1 / t19
      t78 = t77 * t8
      t79 = t37 * t25
      t82 = t21 ** 2
      t84 = t13 / t82
      t85 = t14 ** 2
      t86 = t85 ** 2
      t104 = (-0.224298561906574129855D1 / t86 / t14 * t8 * t37 - 0.1876
     >99471636595866065D1 * t78 * t37 - 0.145760735710958868637D1 / t17
     >* t8 * t37 - 0.344044309698575627326D0 / t11 * t8 * t37) / t24
      t107 = t50 * pi
      t110 = 0.1D1 / t33 / t32
      t115 = 0.1D1 / t36 / rhoa
      t122 = t34 * t37
      t130 = (-0.126105037207067989169D-1 * t77 * pi * t79 - 0.651778270
     >654185890920D1 * t84 * t104 * t28) * t40
      t133 = 0.149583857013095144140D-1 * t30 * t53 * sigmaaa * t122 * t
     >130
      t138 = 0.498612856710317147135D-2 * t107 * t42 * sigmaaa * t110 * 
     >t37
      t141 = 0.299167714026190288281D-1 * t43 * t44 * t115
      t147 = t64 ** 2
      t148 = 0.1D1 / t147
      t157 = t50 ** 2
      t177 = 0.1D1 / t69
      t182 = 0.500000000000000000000D0 * rhoa * (0.193478431062909061652
     >D-2 * t78 * t79 + 0.100000000000000000000D1 * t84 * t104 + 0.15342
     >6409720027345292D0 * t29 * (-0.498612856710317147135D-2 * t107 * s
     >igmaaa * t110 * t66 - 0.299167714026190288281D-1 * t35 * t115 * t4
     >8 * t65 + 0.149583857013095144140D-1 * t35 * t37 * (-t133 - t138 -
     > t141) * t65 - 0.149583857013095144140D-1 * t35 * t49 * t148 * (-t
     >133 - t138 - t141 - 0.447506605578281866741D-3 * t51 / t52 / t41 *
     > t55 * t57 * t60 * t130 - 0.149168868526093955580D-3 * t157 * t53
     >* t55 / t56 / t32 * t60 - 0.895013211156563733482D-3 * t54 * t58 /
     > t59 / rhoa)) * t177)
      vrhoc(i) = vrhoc(i) - t75 + t76 + t182
      vrhoo(i) = vrhoo(i) - t75 + t76 + t182
      t208 = rhoa * t29 * (0.149583857013095144140D-1 * t30 * t34 * t66 
     >+ 0.223753302789140933370D-3 * t51 * sigmaaa * t57 * t60 * t42 * t
     >65 - 0.149583857013095144140D-1 * t35 * t49 * t148 * (0.1495838570
     >13095144140D-1 * t43 * t122 + 0.447506605578281866741D-3 * t54 * s
     >igmaaa * t57 * t60)) * t177
      t209 = 0.383566024300068363229D-1 * t208
      vsigmacc(i) = vsigmacc(i) + t209
      vsigmaco(i) = vsigmaco(i) + 0.767132048600136726458D-1 * t208
      vsigmaoo(i) = vsigmaoo(i) + t209

               else
             rho = rhoa + rhob
      t2 = 0.250000000000000000000D0 * sigmacc(i)
      t4 = 0.250000000000000000000D0 * sigmaoo(i)
      t6 = 0.500000000000000000000D0 * sigmaco(i)
      sigmaaa = max(0.0D0, t2 + t4 + t6)
      sigmaab = t2 - t4
      sigmabb = max(0.0D0, t2 + t4 - t6)
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t10 = 0.1D1 / pi
      t11 = 0.1D1 / rho
      t12 = t10 * t11
      t13 = t12 ** (0.1D1 / 0.3D1)
      t15 = 0.1D1 + 0.194159335344114122552D0 * t13
      t16 = t12 ** (0.1D1 / 0.6D1)
      t19 = sqrt(t12)
      t21 = t13 ** 2
      t23 = 0.724010193431683113327D1 * t16 + 0.325955091942229212011D1 
     >* t13 + 0.141872281647966739112D1 * t19 + 0.406913004517529319387D
     >0 * t21
      t26 = 0.1D1 + 0.160819794986925350668D2 / t23
      t27 = log(t26)
      t29 = 0.621814D-1 * t15 * t27
      t31 = 0.1D1 + 0.101077332976287768525D0 * t13
      t36 = 0.987212972256927209438D1 * t16 + 0.329180480994506259905D1 
     >* t13 + 0.762327521935289963194D0 * t19 + 0.410025070949612505036D
     >0 * t21
      t39 = 0.1D1 + 0.296087499777934375166D2 / t36
      t40 = log(t39)
      t41 = t31 * t40
      t43 = rhoa - 0.1D1 * rhob
      t44 = t43 * t11
      t45 = 0.1D1 + t44
      t46 = t45 ** (0.1D1 / 0.3D1)
      t49 = 0.1D1 - 0.1D1 * t44
      t50 = t49 ** (0.1D1 / 0.3D1)
      t52 = t46 * t45 + t50 * t49 - 0.2D1
      t53 = t43 ** 2
      t54 = t53 ** 2
      t55 = rho ** 2
      t56 = t55 ** 2
      t57 = 0.1D1 / t56
      t58 = t54 * t57
      t60 = 0.1D1 - 0.1D1 * t58
      t63 = 0.379955235370239451738D-1 * t41 * t52 * t60
      t65 = 0.1D1 + 0.186690969707574028554D0 * t13
      t70 = 0.134579137143944477912D2 * t16 + 0.563098414909787598194D1 
     >* t13 + 0.291521471421917737271D1 * t19 + 0.516066464547863440989D
     >0 * t21
      t73 = 0.1D1 + 0.321639589973850701335D2 / t70
      t74 = log(t73)
      t77 = -0.3109070D-1 * t65 * t74 + t29
      t78 = t77 * t52
      t80 = 0.192366105093153631974D1 * t78 * t58
      t81 = pi ** 2
      t82 = 0.1D1 / t81
      t83 = t46 ** 2
      t85 = t50 ** 2
      t87 = 0.500000000000000000000D0 * t83 + 0.500000000000000000000D0 
     >* t85
      t88 = t87 ** 2
      t89 = t88 * t87
      t90 = t82 * t89
      t91 = t81 * pi
      t92 = t91 * sigma
      t93 = 0.1D1 / t88
      t94 = t92 * t93
      t95 = t81 * rho
      t96 = t95 ** (0.1D1 / 0.3D1)
      t97 = 0.1D1 / t96
      t98 = 0.1D1 / t55
      t99 = t97 * t98
      t101 = (-t29 + t63 + t80) * t81
      t102 = 0.1D1 / t89
      t105 = exp(-0.325889135327092945460D1 * t101 * t102)
      t106 = t105 - 0.1D1
      t107 = 0.1D1 / t106
      t108 = t91 * t107
      t109 = t108 * sigma
      t110 = t93 * t97
      t111 = t110 * t98
      t113 = 0.942319250876317101329D-2 * t109 * t111
      t114 = 0.1D1 + t113
      t115 = t81 ** 2
      t116 = t115 * t81
      t117 = t106 ** 2
      t118 = 0.1D1 / t117
      t119 = t116 * t118
      t120 = sigma ** 2
      t121 = t119 * t120
      t122 = t88 ** 2
      t123 = 0.1D1 / t122
      t124 = t96 ** 2
      t125 = 0.1D1 / t124
      t126 = t123 * t125
      t127 = t126 * t57
      t130 = 0.1D1 + t113 + 0.887965570572103448132D-4 * t121 * t127
      t131 = 0.1D1 / t130
      t132 = t114 * t131
      t136 = 0.1D1 + 0.942319250876317101329D-2 * t94 * t99 * t132
      t137 = log(t136)
      t139 = 0.306852819440054690583D0 * t90 * t137
      zk(i) = rho * (-t29 + t63 + t80 + t139)
      t146 = 0.133333333333333333333D1 * t46 * t11 - 0.13333333333333333
     >3333D1 * t50 * t11
      t149 = 0.379955235370239451738D-1 * t41 * t146 * t60
      t150 = t53 * t43
      t154 = 0.151982094148095780695D0 * t41 * t52 * t150 * t57
      t157 = 0.192366105093153631974D1 * t77 * t146 * t58
      t160 = 0.769464420372614527896D1 * t78 * t150 * t57
      t161 = t82 * t88
      t162 = 0.1D1 / t46
      t165 = 0.1D1 / t50
      t168 = 0.333333333333333333333D0 * t162 * t11 - 0.3333333333333333
     >33333D0 * t165 * t11
      t172 = t102 * t97
      t173 = t92 * t172
      t174 = t98 * t114
      t181 = t91 * t118 * sigma * t93
      t190 = (-0.325889135327092945460D1 * (t149 - t154 + t157 + t160) *
     > t81 * t102 + 0.977667405981278836380D1 * t101 * t123 * t168) * t1
     >05
      t193 = 0.942319250876317101329D-2 * t181 * t99 * t190
      t197 = 0.188463850175263420266D-1 * t109 * t172 * t98 * t168
      t203 = t92 * t110
      t204 = t130 ** 2
      t205 = 0.1D1 / t204
      t210 = t116 / t117 / t106 * t120 * t123
      t211 = t125 * t57
      t217 = 0.1D1 / t122 / t87 * t125
      t228 = 0.1D1 / t136
      t234 = 0.500000000000000000000D0 * rho * (t149 - t154 + t157 + t16
     >0 + 0.920558458320164071749D0 * t161 * t137 * t168 + 0.30685281944
     >0054690583D0 * t90 * (-0.188463850175263420266D-1 * t173 * t174 *
     >t131 * t168 + 0.942319250876317101329D-2 * t94 * t99 * (-t193 - t1
     >97) * t131 - 0.942319250876317101329D-2 * t203 * t174 * t205 * (-t
     >193 - t197 - 0.177593114114420689627D-3 * t210 * t211 * t190 - 0.3
     >55186228228841379252D-3 * t121 * t217 * t57 * t168)) * t228)
      t235 = -t146
      t238 = 0.379955235370239451738D-1 * t41 * t235 * t60
      t241 = 0.192366105093153631974D1 * t77 * t235 * t58
      t242 = -t168
      t258 = (-0.325889135327092945460D1 * (t238 + t154 + t241 - t160) *
     > t81 * t102 + 0.977667405981278836380D1 * t101 * t123 * t242) * t1
     >05
      t261 = 0.942319250876317101329D-2 * t181 * t99 * t258
      t265 = 0.188463850175263420266D-1 * t109 * t172 * t98 * t242
      t289 = 0.500000000000000000000D0 * rho * (t238 + t154 + t241 - t16
     >0 + 0.920558458320164071749D0 * t161 * t137 * t242 + 0.30685281944
     >0054690583D0 * t90 * (-0.188463850175263420266D-1 * t173 * t174 *
     >t131 * t242 + 0.942319250876317101329D-2 * t94 * t99 * (-t261 - t2
     >65) * t131 - 0.942319250876317101329D-2 * t203 * t174 * t205 * (-t
     >261 - t265 - 0.177593114114420689627D-3 * t210 * t211 * t258 - 0.3
     >55186228228841379252D-3 * t121 * t217 * t57 * t242)) * t228)
      t291 = 0.1D1 / t21 * t10
      t294 = 0.402436643158883263334D-2 * t291 * t98 * t27
      t295 = t23 ** 2
      t298 = t16 ** 2
      t299 = t298 ** 2
      t303 = 0.1D1 / t299 / t16 * t10 * t98
      t305 = t291 * t98
      t309 = 0.1D1 / t19 * t10 * t98
      t313 = 0.1D1 / t13 * t10 * t98
      t319 = 0.100000000000000000000D1 * t15 / t295 * (-0.12066836557194
     >7185555D1 * t303 - 0.108651697314076404004D1 * t305 - 0.7093614082
     >39833695563D0 * t309 - 0.271275336345019546258D0 * t313) / t26
      t323 = 0.128016206138671616206D-2 * t305 * t40 * t52 * t60
      t324 = t36 ** 2
      t337 = 0.112499995668310776915D1 * t31 / t324 * (-0.16453549537615
     >4534908D1 * t303 - 0.109726826998168753302D1 * t305 - 0.3811637609
     >67644981600D0 * t309 - 0.273350047299741670024D0 * t313) / t39 * t
     >52 * t60
      t344 = -0.133333333333333333333D1 * t46 * t43 * t98 + 0.1333333333
     >33333333333D1 * t50 * t43 * t98
      t347 = 0.379955235370239451738D-1 * t41 * t344 * t60
      t350 = 0.1D1 / t56 / rho
      t353 = 0.151982094148095780695D0 * t41 * t52 * t54 * t350
      t357 = t70 ** 2
      t372 = 0.192366105093153631974D1 * (0.193478431062909061652D-2 * t
     >291 * t98 * t74 + 0.100000000000000000000D1 * t65 / t357 * (-0.224
     >298561906574129855D1 * t303 - 0.187699471636595866065D1 * t305 - 0
     >.145760735710958868637D1 * t309 - 0.344044309698575627326D0 * t313
     >) / t73 - t294 - t319) * t52 * t58
      t375 = 0.192366105093153631974D1 * t77 * t344 * t58
      t378 = 0.769464420372614527896D1 * t78 * t54 * t350
      t385 = -0.333333333333333333333D0 * t162 * t43 * t98 + 0.333333333
     >333333333333D0 * t165 * t43 * t98
      t393 = t115 * pi
      t397 = 0.1D1 / t96 / t95
      t403 = 0.1D1 / t55 / rho
      t416 = (-0.325889135327092945460D1 * (t294 + t319 - t323 - t337 + 
     >t347 + t353 + t372 + t375 - t378) * t81 * t102 + 0.977667405981278
     >836380D1 * t101 * t123 * t385) * t105
      t419 = 0.942319250876317101329D-2 * t181 * t99 * t416
      t423 = 0.188463850175263420266D-1 * t109 * t172 * t98 * t385
      t429 = 0.314106416958772367110D-2 * t393 * t107 * sigma * t93 * t3
     >97 * t98
      t432 = 0.188463850175263420266D-1 * t109 * t110 * t403
      t445 = t115 ** 2
      t466 = t294 + t319 - t323 - t337 + t347 + t353 + t372 + t375 - t37
     >8 + 0.920558458320164071749D0 * t161 * t137 * t385 + 0.30685281944
     >0054690583D0 * t90 * (-0.188463850175263420266D-1 * t173 * t174 *
     >t131 * t385 - 0.314106416958772367110D-2 * t393 * sigma * t93 * t3
     >97 * t98 * t132 - 0.188463850175263420266D-1 * t94 * t97 * t403 *
     >t132 + 0.942319250876317101329D-2 * t94 * t99 * (-t419 - t423 - t4
     >29 - t432) * t131 - 0.942319250876317101329D-2 * t203 * t174 * t20
     >5 * (-t419 - t423 - t429 - t432 - 0.177593114114420689627D-3 * t21
     >0 * t211 * t416 - 0.355186228228841379252D-3 * t121 * t217 * t57 *
     > t385 - 0.591977047048068965421D-4 * t445 * t118 * t120 * t123 / t
     >124 / t95 * t57 - 0.355186228228841379252D-3 * t121 * t126 * t350)
     >) * t228
      vrhoc(i) = vrhoc(i) + t234 + t289 - t29 + t63 + t80 + t139 + rho *
     > t466
      vrhoo(i) = vrhoo(i) + t234 - t289
      vsigmacc(i) = vsigmacc(i) + 0.306852819440054690583D0 * rho * t82 
     >* t89 * (0.942319250876317101329D-2 * t91 * t93 * t97 * t174 * t13
     >1 + 0.887965570572103448136D-4 * t116 * sigma * t123 * t211 * t107
     > * t131 - 0.942319250876317101329D-2 * t203 * t174 * t205 * (0.942
     >319250876317101329D-2 * t108 * t111 + 0.177593114114420689627D-3 *
     > t119 * sigma * t127)) * t228
      vsigmaco(i) = vsigmaco(i)
      vsigmaoo(i) = vsigmaoo(i)

               endif
             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      sigmabb = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) - 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmabb
      t10 = 0.1D1 / pi / rhob
      t11 = t10 ** (0.1D1 / 0.3D1)
      t14 = t10 ** (0.1D1 / 0.6D1)
      t17 = sqrt(t10)
      t19 = t11 ** 2
      t25 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t14 + 0.563098414909787598194D1 * t11 + 0.291521471421917
     >737271D1 * t17 + 0.516066464547863440989D0 * t19))
      t26 = (0.1D1 + 0.186690969707574028554D0 * t11) * t25
      t28 = pi ** 2
      t30 = t28 * pi
      t33 = (t28 * rhob) ** (0.1D1 / 0.3D1)
      t34 = 0.1D1 / t33
      t36 = rhob ** 2
      t37 = 0.1D1 / t36
      t40 = exp(0.202642426794280972788D0 * t26 * t28)
      t41 = t40 - 0.1D1
      t47 = 0.149583857013095144140D-1 * t30 / t41 * sigmabb * t34 * t37
      t50 = t28 ** 2
      t52 = t41 ** 2
      t55 = sigmabb ** 2
      t56 = t33 ** 2
      t59 = t36 ** 2
      t70 = log(0.1D1 + 0.149583857013095144140D-1 * t30 * sigmabb * t34
     > * t37 * (0.1D1 + t47) / (0.1D1 + t47 + 0.223753302789140933370D-3
     > * t50 * t28 / t52 * t55 / t56 / t59))
      zk(i) = rhob * (-0.3109070D-1 * t26 + 0.153426409720027345292D0 / 
     >t28 * t70)

               elseif(rhob.lt.tol) then
             rho = rhoa
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) + 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmaaa
      t10 = 0.1D1 / pi / rhoa
      t11 = t10 ** (0.1D1 / 0.3D1)
      t14 = t10 ** (0.1D1 / 0.6D1)
      t17 = sqrt(t10)
      t19 = t11 ** 2
      t25 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t14 + 0.563098414909787598194D1 * t11 + 0.291521471421917
     >737271D1 * t17 + 0.516066464547863440989D0 * t19))
      t26 = (0.1D1 + 0.186690969707574028554D0 * t11) * t25
      t28 = pi ** 2
      t30 = t28 * pi
      t33 = (t28 * rhoa) ** (0.1D1 / 0.3D1)
      t34 = 0.1D1 / t33
      t36 = rhoa ** 2
      t37 = 0.1D1 / t36
      t40 = exp(0.202642426794280972788D0 * t26 * t28)
      t41 = t40 - 0.1D1
      t47 = 0.149583857013095144140D-1 * t30 / t41 * sigmaaa * t34 * t37
      t50 = t28 ** 2
      t52 = t41 ** 2
      t55 = sigmaaa ** 2
      t56 = t33 ** 2
      t59 = t36 ** 2
      t70 = log(0.1D1 + 0.149583857013095144140D-1 * t30 * sigmaaa * t34
     > * t37 * (0.1D1 + t47) / (0.1D1 + t47 + 0.223753302789140933370D-3
     > * t50 * t28 / t52 * t55 / t56 / t59))
      zk(i) = rhoa * (-0.3109070D-1 * t26 + 0.153426409720027345292D0 / 
     >t28 * t70)

               else
             rho = rhoa + rhob
      t2 = 0.250000000000000000000D0 * sigmacc(i)
      t4 = 0.250000000000000000000D0 * sigmaoo(i)
      t6 = 0.500000000000000000000D0 * sigmaco(i)
      sigmaaa = max(0.0D0, t2 + t4 + t6)
      sigmaab = t2 - t4
      sigmabb = max(0.0D0, t2 + t4 - t6)
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t11 = 0.1D1 / rho
      t12 = 0.1D1 / pi * t11
      t13 = t12 ** (0.1D1 / 0.3D1)
      t16 = t12 ** (0.1D1 / 0.6D1)
      t19 = sqrt(t12)
      t21 = t13 ** 2
      t27 = log(0.1D1 + 0.160819794986925350668D2 / (0.72401019343168311
     >3327D1 * t16 + 0.325955091942229212011D1 * t13 + 0.141872281647966
     >739112D1 * t19 + 0.406913004517529319387D0 * t21))
      t29 = 0.621814D-1 * (0.1D1 + 0.194159335344114122552D0 * t13) * t2
     >7
      t40 = log(0.1D1 + 0.296087499777934375166D2 / (0.98721297225692720
     >9438D1 * t16 + 0.329180480994506259905D1 * t13 + 0.762327521935289
     >963194D0 * t19 + 0.410025070949612505036D0 * t21))
      t43 = rhoa - 0.1D1 * rhob
      t44 = t43 * t11
      t45 = 0.1D1 + t44
      t46 = t45 ** (0.1D1 / 0.3D1)
      t49 = 0.1D1 - 0.1D1 * t44
      t50 = t49 ** (0.1D1 / 0.3D1)
      t52 = t46 * t45 + t50 * t49 - 0.2D1
      t53 = t43 ** 2
      t54 = t53 ** 2
      t55 = rho ** 2
      t56 = t55 ** 2
      t57 = 0.1D1 / t56
      t58 = t54 * t57
      t63 = 0.379955235370239451738D-1 * (0.1D1 + 0.10107733297628776852
     >5D0 * t13) * t40 * t52 * (0.1D1 - 0.1D1 * t58)
      t74 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t16 + 0.563098414909787598194D1 * t13 + 0.291521471421917
     >737271D1 * t19 + 0.516066464547863440989D0 * t21))
      t80 = 0.192366105093153631974D1 * (-0.3109070D-1 * (0.1D1 + 0.1866
     >90969707574028554D0 * t13) * t74 + t29) * t52 * t58
      t81 = pi ** 2
      t83 = t46 ** 2
      t85 = t50 ** 2
      t87 = 0.500000000000000000000D0 * t83 + 0.500000000000000000000D0 
     >* t85
      t88 = t87 ** 2
      t89 = t88 * t87
      t91 = t81 * pi
      t93 = 0.1D1 / t88
      t96 = (t81 * rho) ** (0.1D1 / 0.3D1)
      t97 = 0.1D1 / t96
      t98 = 0.1D1 / t55
      t105 = exp(-0.325889135327092945460D1 * (-t29 + t63 + t80) * t81 /
     > t89)
      t106 = t105 - 0.1D1
      t113 = 0.942319250876317101329D-2 * t91 / t106 * sigma * t93 * t97
     > * t98
      t115 = t81 ** 2
      t117 = t106 ** 2
      t120 = sigma ** 2
      t122 = t88 ** 2
      t124 = t96 ** 2
      t137 = log(0.1D1 + 0.942319250876317101329D-2 * t91 * sigma * t93 
     >* t97 * t98 * (0.1D1 + t113) / (0.1D1 + t113 + 0.88796557057210344
     >8132D-4 * t115 * t81 / t117 * t120 / t122 / t124 * t57))
      zk(i) = rho * (-t29 + t63 + t80 + 0.306852819440054690583D0 / t81 
     >* t89 * t137)

               endif
             endif
           enddo
         endif
       else
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i))
      sigmaab = sigmaaa
      sigmabb = sigmaab
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t6 = 0.1D1 / pi
      t7 = 0.1D1 / rho
      t8 = t6 * t7
      t9 = t8 ** (0.1D1 / 0.3D1)
      t11 = 0.1D1 + 0.194159335344114122552D0 * t9
      t12 = t8 ** (0.1D1 / 0.6D1)
      t15 = sqrt(t8)
      t17 = t9 ** 2
      t19 = 0.724010193431683113327D1 * t12 + 0.325955091942229212011D1 
     >* t9 + 0.141872281647966739112D1 * t15 + 0.406913004517529319387D0
     > * t17
      t22 = 0.1D1 + 0.160819794986925350668D2 / t19
      t23 = log(t22)
      t25 = 0.621814D-1 * t11 * t23
      t27 = 0.1D1 + 0.101077332976287768525D0 * t9
      t32 = 0.987212972256927209438D1 * t12 + 0.329180480994506259905D1 
     >* t9 + 0.762327521935289963194D0 * t15 + 0.410025070949612505036D0
     > * t17
      t35 = 0.1D1 + 0.296087499777934375166D2 / t32
      t36 = log(t35)
      t37 = t27 * t36
      t39 = rhoa - 0.1D1 * rhob
      t40 = t39 * t7
      t41 = 0.1D1 + t40
      t42 = t41 ** (0.1D1 / 0.3D1)
      t45 = 0.1D1 - 0.1D1 * t40
      t46 = t45 ** (0.1D1 / 0.3D1)
      t48 = t42 * t41 + t46 * t45 - 0.2D1
      t49 = t39 ** 2
      t50 = t49 ** 2
      t51 = rho ** 2
      t52 = t51 ** 2
      t53 = 0.1D1 / t52
      t54 = t50 * t53
      t56 = 0.1D1 - 0.1D1 * t54
      t59 = 0.379955235370239451738D-1 * t37 * t48 * t56
      t61 = 0.1D1 + 0.186690969707574028554D0 * t9
      t66 = 0.134579137143944477912D2 * t12 + 0.563098414909787598194D1 
     >* t9 + 0.291521471421917737271D1 * t15 + 0.516066464547863440989D0
     > * t17
      t69 = 0.1D1 + 0.321639589973850701335D2 / t66
      t70 = log(t69)
      t73 = -0.3109070D-1 * t61 * t70 + t25
      t74 = t73 * t48
      t76 = 0.192366105093153631974D1 * t74 * t54
      t77 = pi ** 2
      t78 = 0.1D1 / t77
      t79 = t42 ** 2
      t81 = t46 ** 2
      t83 = 0.500000000000000000000D0 * t79 + 0.500000000000000000000D0 
     >* t81
      t84 = t83 ** 2
      t85 = t84 * t83
      t86 = t78 * t85
      t87 = t77 * pi
      t88 = t87 * sigma
      t89 = 0.1D1 / t84
      t90 = t88 * t89
      t91 = t77 * rho
      t92 = t91 ** (0.1D1 / 0.3D1)
      t93 = 0.1D1 / t92
      t94 = 0.1D1 / t51
      t95 = t93 * t94
      t97 = (-t25 + t59 + t76) * t77
      t98 = 0.1D1 / t85
      t101 = exp(-0.325889135327092945460D1 * t97 * t98)
      t102 = t101 - 0.1D1
      t103 = 0.1D1 / t102
      t104 = t87 * t103
      t105 = t104 * sigma
      t106 = t89 * t93
      t107 = t106 * t94
      t109 = 0.942319250876317101329D-2 * t105 * t107
      t110 = 0.1D1 + t109
      t111 = t77 ** 2
      t112 = t111 * t77
      t113 = t102 ** 2
      t114 = 0.1D1 / t113
      t115 = t112 * t114
      t116 = sigma ** 2
      t117 = t115 * t116
      t118 = t84 ** 2
      t119 = 0.1D1 / t118
      t120 = t92 ** 2
      t121 = 0.1D1 / t120
      t122 = t119 * t121
      t123 = t122 * t53
      t126 = 0.1D1 + t109 + 0.887965570572103448132D-4 * t117 * t123
      t127 = 0.1D1 / t126
      t128 = t110 * t127
      t132 = 0.1D1 + 0.942319250876317101329D-2 * t90 * t95 * t128
      t133 = log(t132)
      t135 = 0.306852819440054690583D0 * t86 * t133
      zk(i) = rho * (-t25 + t59 + t76 + t135)
      t142 = 0.133333333333333333333D1 * t42 * t7 - 0.133333333333333333
     >333D1 * t46 * t7
      t145 = 0.379955235370239451738D-1 * t37 * t142 * t56
      t146 = t49 * t39
      t150 = 0.151982094148095780695D0 * t37 * t48 * t146 * t53
      t153 = 0.192366105093153631974D1 * t73 * t142 * t54
      t156 = 0.769464420372614527896D1 * t74 * t146 * t53
      t157 = t78 * t84
      t158 = 0.1D1 / t42
      t161 = 0.1D1 / t46
      t164 = 0.333333333333333333333D0 * t158 * t7 - 0.33333333333333333
     >3333D0 * t161 * t7
      t168 = t98 * t93
      t169 = t88 * t168
      t170 = t94 * t110
      t177 = t87 * t114 * sigma * t89
      t186 = (-0.325889135327092945460D1 * (t145 - t150 + t153 + t156) *
     > t77 * t98 + 0.977667405981278836380D1 * t97 * t119 * t164) * t101
      t189 = 0.942319250876317101329D-2 * t177 * t95 * t186
      t193 = 0.188463850175263420266D-1 * t105 * t168 * t94 * t164
      t199 = t88 * t106
      t200 = t126 ** 2
      t201 = 0.1D1 / t200
      t206 = t112 / t113 / t102 * t116 * t119
      t207 = t121 * t53
      t213 = 0.1D1 / t118 / t83 * t121
      t224 = 0.1D1 / t132
      t231 = -t142
      t234 = 0.379955235370239451738D-1 * t37 * t231 * t56
      t237 = 0.192366105093153631974D1 * t73 * t231 * t54
      t238 = -t164
      t254 = (-0.325889135327092945460D1 * (t234 + t150 + t237 - t156) *
     > t77 * t98 + 0.977667405981278836380D1 * t97 * t119 * t238) * t101
      t257 = 0.942319250876317101329D-2 * t177 * t95 * t254
      t261 = 0.188463850175263420266D-1 * t105 * t168 * t94 * t238
      t287 = 0.1D1 / t17 * t6
      t290 = 0.402436643158883263334D-2 * t287 * t94 * t23
      t291 = t19 ** 2
      t294 = t12 ** 2
      t295 = t294 ** 2
      t299 = 0.1D1 / t295 / t12 * t6 * t94
      t301 = t287 * t94
      t305 = 0.1D1 / t15 * t6 * t94
      t309 = 0.1D1 / t9 * t6 * t94
      t315 = 0.100000000000000000000D1 * t11 / t291 * (-0.12066836557194
     >7185555D1 * t299 - 0.108651697314076404004D1 * t301 - 0.7093614082
     >39833695563D0 * t305 - 0.271275336345019546258D0 * t309) / t22
      t319 = 0.128016206138671616206D-2 * t301 * t36 * t48 * t56
      t320 = t32 ** 2
      t333 = 0.112499995668310776915D1 * t27 / t320 * (-0.16453549537615
     >4534908D1 * t299 - 0.109726826998168753302D1 * t301 - 0.3811637609
     >67644981600D0 * t305 - 0.273350047299741670024D0 * t309) / t35 * t
     >48 * t56
      t340 = -0.133333333333333333333D1 * t42 * t39 * t94 + 0.1333333333
     >33333333333D1 * t46 * t39 * t94
      t343 = 0.379955235370239451738D-1 * t37 * t340 * t56
      t346 = 0.1D1 / t52 / rho
      t349 = 0.151982094148095780695D0 * t37 * t48 * t50 * t346
      t353 = t66 ** 2
      t368 = 0.192366105093153631974D1 * (0.193478431062909061652D-2 * t
     >287 * t94 * t70 + 0.100000000000000000000D1 * t61 / t353 * (-0.224
     >298561906574129855D1 * t299 - 0.187699471636595866065D1 * t301 - 0
     >.145760735710958868637D1 * t305 - 0.344044309698575627326D0 * t309
     >) / t69 - t290 - t315) * t48 * t54
      t371 = 0.192366105093153631974D1 * t73 * t340 * t54
      t374 = 0.769464420372614527896D1 * t74 * t50 * t346
      t381 = -0.333333333333333333333D0 * t158 * t39 * t94 + 0.333333333
     >333333333333D0 * t161 * t39 * t94
      t389 = t111 * pi
      t393 = 0.1D1 / t92 / t91
      t399 = 0.1D1 / t51 / rho
      t412 = (-0.325889135327092945460D1 * (t290 + t315 - t319 - t333 + 
     >t343 + t349 + t368 + t371 - t374) * t77 * t98 + 0.9776674059812788
     >36380D1 * t97 * t119 * t381) * t101
      t415 = 0.942319250876317101329D-2 * t177 * t95 * t412
      t419 = 0.188463850175263420266D-1 * t105 * t168 * t94 * t381
      t425 = 0.314106416958772367110D-2 * t389 * t103 * sigma * t89 * t3
     >93 * t94
      t428 = 0.188463850175263420266D-1 * t105 * t106 * t399
      t441 = t111 ** 2
      t462 = t290 + t315 - t319 - t333 + t343 + t349 + t368 + t371 - t37
     >4 + 0.920558458320164071749D0 * t157 * t133 * t381 + 0.30685281944
     >0054690583D0 * t86 * (-0.188463850175263420266D-1 * t169 * t170 *
     >t127 * t381 - 0.314106416958772367110D-2 * t389 * sigma * t89 * t3
     >93 * t94 * t128 - 0.188463850175263420266D-1 * t90 * t93 * t399 *
     >t128 + 0.942319250876317101329D-2 * t90 * t95 * (-t415 - t419 - t4
     >25 - t428) * t127 - 0.942319250876317101329D-2 * t199 * t170 * t20
     >1 * (-t415 - t419 - t425 - t428 - 0.177593114114420689627D-3 * t20
     >6 * t207 * t412 - 0.355186228228841379252D-3 * t117 * t213 * t53 *
     > t381 - 0.591977047048068965421D-4 * t441 * t114 * t116 * t119 / t
     >120 / t91 * t53 - 0.355186228228841379252D-3 * t117 * t122 * t346)
     >) * t224
      vrhoc(i) = vrhoc(i) + 0.5D0 * rho * (t145 - t150 + t153 + t156 + 0
     >.920558458320164071749D0 * t157 * t133 * t164 + 0.3068528194400546
     >90583D0 * t86 * (-0.188463850175263420266D-1 * t169 * t170 * t127
     >* t164 + 0.942319250876317101329D-2 * t90 * t95 * (-t189 - t193) *
     > t127 - 0.942319250876317101329D-2 * t199 * t170 * t201 * (-t189 -
     > t193 - 0.177593114114420689627D-3 * t206 * t207 * t186 - 0.355186
     >228228841379252D-3 * t117 * t213 * t53 * t164)) * t224) + 0.5D0 *
     >rho * (t234 + t150 + t237 - t156 + 0.920558458320164071749D0 * t15
     >7 * t133 * t238 + 0.306852819440054690583D0 * t86 * (-0.1884638501
     >75263420266D-1 * t169 * t170 * t127 * t238 + 0.9423192508763171013
     >29D-2 * t90 * t95 * (-t257 - t261) * t127 - 0.94231925087631710132
     >9D-2 * t199 * t170 * t201 * (-t257 - t261 - 0.17759311411442068962
     >7D-3 * t206 * t207 * t254 - 0.355186228228841379252D-3 * t117 * t2
     >13 * t53 * t238)) * t224) - t25 + t59 + t76 + t135 + rho * t462
      vsigmacc(i) = vsigmacc(i) + 0.306852819440054690583D0 * rho * t78 
     >* t85 * (0.942319250876317101329D-2 * t87 * t89 * t93 * t170 * t12
     >7 + 0.887965570572103448136D-4 * t112 * sigma * t119 * t207 * t103
     > * t127 - 0.942319250876317101329D-2 * t199 * t170 * t201 * (0.942
     >319250876317101329D-2 * t104 * t107 + 0.177593114114420689627D-3 *
     > t115 * sigma * t123)) * t224

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i))
      sigmaab = sigmaaa
      sigmabb = sigmaab
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t7 = 0.1D1 / rho
      t8 = 0.1D1 / pi * t7
      t9 = t8 ** (0.1D1 / 0.3D1)
      t12 = t8 ** (0.1D1 / 0.6D1)
      t15 = sqrt(t8)
      t17 = t9 ** 2
      t23 = log(0.1D1 + 0.160819794986925350668D2 / (0.72401019343168311
     >3327D1 * t12 + 0.325955091942229212011D1 * t9 + 0.1418722816479667
     >39112D1 * t15 + 0.406913004517529319387D0 * t17))
      t25 = 0.621814D-1 * (0.1D1 + 0.194159335344114122552D0 * t9) * t23
      t36 = log(0.1D1 + 0.296087499777934375166D2 / (0.98721297225692720
     >9438D1 * t12 + 0.329180480994506259905D1 * t9 + 0.7623275219352899
     >63194D0 * t15 + 0.410025070949612505036D0 * t17))
      t39 = rhoa - 0.1D1 * rhob
      t40 = t39 * t7
      t41 = 0.1D1 + t40
      t42 = t41 ** (0.1D1 / 0.3D1)
      t45 = 0.1D1 - 0.1D1 * t40
      t46 = t45 ** (0.1D1 / 0.3D1)
      t48 = t42 * t41 + t46 * t45 - 0.2D1
      t49 = t39 ** 2
      t50 = t49 ** 2
      t51 = rho ** 2
      t52 = t51 ** 2
      t53 = 0.1D1 / t52
      t54 = t50 * t53
      t59 = 0.379955235370239451738D-1 * (0.1D1 + 0.10107733297628776852
     >5D0 * t9) * t36 * t48 * (0.1D1 - 0.1D1 * t54)
      t70 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t12 + 0.563098414909787598194D1 * t9 + 0.2915214714219177
     >37271D1 * t15 + 0.516066464547863440989D0 * t17))
      t76 = 0.192366105093153631974D1 * (-0.3109070D-1 * (0.1D1 + 0.1866
     >90969707574028554D0 * t9) * t70 + t25) * t48 * t54
      t77 = pi ** 2
      t79 = t42 ** 2
      t81 = t46 ** 2
      t83 = 0.500000000000000000000D0 * t79 + 0.500000000000000000000D0 
     >* t81
      t84 = t83 ** 2
      t85 = t84 * t83
      t87 = t77 * pi
      t89 = 0.1D1 / t84
      t92 = (t77 * rho) ** (0.1D1 / 0.3D1)
      t93 = 0.1D1 / t92
      t94 = 0.1D1 / t51
      t101 = exp(-0.325889135327092945460D1 * (-t25 + t59 + t76) * t77 /
     > t85)
      t102 = t101 - 0.1D1
      t109 = 0.942319250876317101329D-2 * t87 / t102 * sigma * t89 * t93
     > * t94
      t111 = t77 ** 2
      t113 = t102 ** 2
      t116 = sigma ** 2
      t118 = t84 ** 2
      t120 = t92 ** 2
      t133 = log(0.1D1 + 0.942319250876317101329D-2 * t87 * sigma * t89 
     >* t93 * t94 * (0.1D1 + t109) / (0.1D1 + t109 + 0.88796557057210344
     >8132D-4 * t111 * t77 / t113 * t116 / t118 / t120 * t53))
      zk(i) = rho * (-t25 + t59 + t76 + 0.306852819440054690583D0 / t77 
     >* t85 * t133)

             endif
           enddo
         endif
       endif

      return
      end

c:PBECsubrend
c:PBEXsubrstart

c    Generated: Mon Jun 20 16:05:15 CEST 2016

      subroutine dftacg_pbex
     > (name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   tauc,tauo,upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,
     >                   vtauc,vtauo,vupsilonc,vupsilono)
      implicit double precision (a-h,o-z)
      logical fderiv,open
      integer igrad,npt
      character*(*) name
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)
      double precision tauc(*),tauo(*)
      double precision upsilonc(*),upsilono(*)
      double precision zk(*),vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)
      double precision vtauc(*),vtauo(*)
      double precision vupsilonc(*),vupsilono(*)
      include "common/cdft"
      include "common/tapes"
      parameter(tol=1d-12)
      pi=acos(-1d0)
      name='Automatically generated PBEX'
      igrad=1
       if(open) then
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      sigmabb = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) - 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmabb
      t8 = pi ** 2
      t9 = t8 * rhob
      t10 = t9 ** (0.1D1 / 0.3D1)
      t11 = rhob * t10
      t12 = 0.1D1 / pi
      t13 = t8 * sigmabb
      t14 = t10 ** 2
      t15 = 0.1D1 / t14
      t16 = rhob ** 2
      t17 = 0.1D1 / t16
      t21 = 0.1D1 + 0.209451650699150979027D-2 * t13 * t15 * t17
      t24 = 0.1804D1 - 0.804D0 / t21
      zk(i) = -0.136284044462410474417D1 * t11 * t12 * t24
      t31 = 0.681420222312052372084D0 * t10 * t12 * t24
      t35 = 0.227140074104017457361D0 * rhob * t15 * pi * t24
      t36 = t21 ** 2
      t37 = 0.1D1 / t36
      t39 = t8 ** 2
      t54 = 0.547861858738890107155D0 * t11 * t12 * t37 * (-0.1396344337
     >99433986018D-2 * t39 * sigmabb / t14 / t9 * t17 - 0.41890330139830
     >1958055D-2 * t13 * t15 / t16 / rhob)
      vrhoc(i) = vrhoc(i) - t31 - t35 - t54
      vrhoo(i) = vrhoo(i) + t31 + t35 + t54
      t61 = 0.1D1 / rhob / t10 * pi * t37
      t62 = 0.573752853339828035106D-3 * t61
      vsigmacc(i) = vsigmacc(i) - t62
      vsigmaco(i) = vsigmaco(i) + 0.114750570667965607021D-2 * t61
      vsigmaoo(i) = vsigmaoo(i) - t62

               elseif(rhob.lt.tol) then
             rho = rhoa
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) + 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmaaa
      t8 = pi ** 2
      t9 = t8 * rhoa
      t10 = t9 ** (0.1D1 / 0.3D1)
      t11 = rhoa * t10
      t12 = 0.1D1 / pi
      t13 = t8 * sigmaaa
      t14 = t10 ** 2
      t15 = 0.1D1 / t14
      t16 = rhoa ** 2
      t17 = 0.1D1 / t16
      t21 = 0.1D1 + 0.209451650699150979027D-2 * t13 * t15 * t17
      t24 = 0.1804D1 - 0.804D0 / t21
      zk(i) = -0.136284044462410474417D1 * t11 * t12 * t24
      t31 = 0.681420222312052372084D0 * t10 * t12 * t24
      t35 = 0.227140074104017457361D0 * rhoa * t15 * pi * t24
      t36 = t21 ** 2
      t37 = 0.1D1 / t36
      t39 = t8 ** 2
      t54 = 0.547861858738890107155D0 * t11 * t12 * t37 * (-0.1396344337
     >99433986018D-2 * t39 * sigmaaa / t14 / t9 * t17 - 0.41890330139830
     >1958055D-2 * t13 * t15 / t16 / rhoa)
      vrhoc(i) = vrhoc(i) - t31 - t35 - t54
      vrhoo(i) = vrhoo(i) - t31 - t35 - t54
      t61 = 0.1D1 / rhoa / t10 * pi * t37
      t62 = 0.573752853339828035106D-3 * t61
      vsigmacc(i) = vsigmacc(i) - t62
      vsigmaco(i) = vsigmaco(i) - 0.114750570667965607021D-2 * t61
      vsigmaoo(i) = vsigmaoo(i) - t62

               else
             rho = rhoa + rhob
      t2 = 0.250000000000000000000D0 * sigmacc(i)
      t4 = 0.250000000000000000000D0 * sigmaoo(i)
      t6 = 0.500000000000000000000D0 * sigmaco(i)
      sigmaaa = max(0.0D0, t2 + t4 + t6)
      sigmaab = t2 - t4
      sigmabb = max(0.0D0, t2 + t4 - t6)
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t10 = pi ** 2
      t11 = rhoa * t10
      t12 = t11 ** (0.1D1 / 0.3D1)
      t13 = rhoa * t12
      t14 = 0.1D1 / pi
      t15 = t10 * sigmaaa
      t16 = t12 ** 2
      t17 = 0.1D1 / t16
      t18 = rhoa ** 2
      t19 = 0.1D1 / t18
      t23 = 0.1D1 + 0.209451650699150979027D-2 * t15 * t17 * t19
      t26 = 0.1804D1 - 0.804D0 / t23
      t30 = rhob * t10
      t31 = t30 ** (0.1D1 / 0.3D1)
      t32 = rhob * t31
      t33 = t10 * sigmabb
      t34 = t31 ** 2
      t35 = 0.1D1 / t34
      t36 = rhob ** 2
      t37 = 0.1D1 / t36
      t41 = 0.1D1 + 0.209451650699150979027D-2 * t33 * t35 * t37
      t44 = 0.1804D1 - 0.804D0 / t41
      zk(i) = -0.136284044462410474417D1 * t13 * t14 * t26 - 0.136284044
     >462410474417D1 * t32 * t14 * t44
      t51 = 0.681420222312052372084D0 * t12 * t14 * t26
      t55 = 0.227140074104017457361D0 * rhoa * t17 * pi * t26
      t56 = t23 ** 2
      t57 = 0.1D1 / t56
      t59 = t10 ** 2
      t74 = 0.547861858738890107155D0 * t13 * t14 * t57 * (-0.1396344337
     >99433986018D-2 * t59 * sigmaaa / t16 / t11 * t19 - 0.4189033013983
     >01958055D-2 * t15 * t17 / t18 / rhoa)
      t77 = 0.681420222312052372084D0 * t31 * t14 * t44
      t81 = 0.227140074104017457361D0 * rhob * t35 * pi * t44
      t82 = t41 ** 2
      t83 = 0.1D1 / t82
      t99 = 0.547861858738890107155D0 * t32 * t14 * t83 * (-0.1396344337
     >99433986018D-2 * t59 * sigmabb / t34 / t30 * t37 - 0.4189033013983
     >01958055D-2 * t33 * t35 / t36 / rhob)
      vrhoc(i) = vrhoc(i) - t51 - t55 - t74 - t77 - t81 - t99
      vrhoo(i) = vrhoo(i) - t51 - t55 - t74 + t77 + t81 + t99
      t106 = 0.1D1 / rhoa / t12 * pi * t57
      t107 = 0.573752853339828035106D-3 * t106
      t112 = 0.1D1 / rhob / t31 * pi * t83
      t113 = 0.573752853339828035106D-3 * t112
      vsigmacc(i) = vsigmacc(i) - t107 - t113
      vsigmaco(i) = vsigmaco(i) - 0.114750570667965607021D-2 * t106 + 0.
     >114750570667965607021D-2 * t112
      vsigmaoo(i) = vsigmaoo(i) - t107 - t113

               endif
             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      sigmabb = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) - 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmabb
      t8 = pi ** 2
      t10 = (t8 * rhob) ** (0.1D1 / 0.3D1)
      t14 = t10 ** 2
      t16 = rhob ** 2
      zk(i) = -0.136284044462410474417D1 * rhob * t10 / pi * (0.1804D1 -
     > 0.804D0 / (0.1D1 + 0.209451650699150979027D-2 * t8 * sigmabb / t1
     >4 / t16))

               elseif(rhob.lt.tol) then
             rho = rhoa
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) + 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmaaa
      t8 = pi ** 2
      t10 = (t8 * rhoa) ** (0.1D1 / 0.3D1)
      t14 = t10 ** 2
      t16 = rhoa ** 2
      zk(i) = -0.136284044462410474417D1 * rhoa * t10 / pi * (0.1804D1 -
     > 0.804D0 / (0.1D1 + 0.209451650699150979027D-2 * t8 * sigmaaa / t1
     >4 / t16))

               else
             rho = rhoa + rhob
      t2 = 0.250000000000000000000D0 * sigmacc(i)
      t4 = 0.250000000000000000000D0 * sigmaoo(i)
      t6 = 0.500000000000000000000D0 * sigmaco(i)
      sigmaaa = max(0.0D0, t2 + t4 + t6)
      sigmaab = t2 - t4
      sigmabb = max(0.0D0, t2 + t4 - t6)
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t10 = pi ** 2
      t12 = (rhoa * t10) ** (0.1D1 / 0.3D1)
      t14 = 0.1D1 / pi
      t16 = t12 ** 2
      t18 = rhoa ** 2
      t31 = (rhob * t10) ** (0.1D1 / 0.3D1)
      t34 = t31 ** 2
      t36 = rhob ** 2
      zk(i) = -0.136284044462410474417D1 * rhoa * t12 * t14 * (0.1804D1 
     >- 0.804D0 / (0.1D1 + 0.209451650699150979027D-2 * t10 * sigmaaa /
     >t16 / t18)) - 0.136284044462410474417D1 * rhob * t31 * t14 * (0.18
     >04D1 - 0.804D0 / (0.1D1 + 0.209451650699150979027D-2 * t10 * sigma
     >bb / t34 / t36))

               endif
             endif
           enddo
         endif
       else
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i))
      sigmaab = sigmaaa
      sigmabb = sigmaab
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t6 = pi ** 2
      t7 = t6 * rhoa
      t8 = t7 ** (0.1D1 / 0.3D1)
      t9 = t8 * rhoa
      t10 = 0.1D1 / pi
      t11 = t6 * sigmaaa
      t12 = t8 ** 2
      t13 = 0.1D1 / t12
      t14 = rhoa ** 2
      t15 = 0.1D1 / t14
      t19 = 0.1D1 + 0.209451650699150979027D-2 * t11 * t13 * t15
      t22 = 0.1804D1 - 0.804D0 / t19
      t26 = t6 * rhob
      t27 = t26 ** (0.1D1 / 0.3D1)
      t28 = rhob * t27
      t29 = t6 * sigmabb
      t30 = t27 ** 2
      t31 = 0.1D1 / t30
      t32 = rhob ** 2
      t33 = 0.1D1 / t32
      t37 = 0.1D1 + 0.209451650699150979027D-2 * t29 * t31 * t33
      t40 = 0.1804D1 - 0.804D0 / t37
      zk(i) = -0.136284044462410474417D1 * t9 * t10 * t22 - 0.1362840444
     >62410474417D1 * t28 * t10 * t40
      t52 = t19 ** 2
      t53 = 0.1D1 / t52
      t55 = t6 ** 2
      t78 = t37 ** 2
      t79 = 0.1D1 / t78
      vrhoc(i) = vrhoc(i) - 0.681420222312052372084D0 * t8 * t10 * t22 -
     > 0.227140074104017457361D0 * rhoa * t13 * pi * t22 - 0.54786185873
     >8890107155D0 * t9 * t10 * t53 * (-0.139634433799433986018D-2 * t55
     > * sigmaaa / t12 / t7 * t15 - 0.418903301398301958055D-2 * t11 * t
     >13 / t14 / rhoa) - 0.681420222312052372084D0 * t27 * t10 * t40 - 0
     >.227140074104017457361D0 * rhob * t31 * pi * t40 - 0.5478618587388
     >90107155D0 * t28 * t10 * t79 * (-0.139634433799433986018D-2 * t55
     >* sigmabb / t30 / t26 * t33 - 0.418903301398301958055D-2 * t29 * t
     >31 / t32 / rhob)
      vsigmacc(i) = vsigmacc(i) - 0.573752853339828035106D-3 / rhoa / t8
     > * pi * t53 - 0.573752853339828035106D-3 / rhob / t27 * pi * t79

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i))
      sigmaab = sigmaaa
      sigmabb = sigmaab
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t6 = pi ** 2
      t8 = (t6 * rhoa) ** (0.1D1 / 0.3D1)
      t10 = 0.1D1 / pi
      t12 = t8 ** 2
      t14 = rhoa ** 2
      t27 = (t6 * rhob) ** (0.1D1 / 0.3D1)
      t30 = t27 ** 2
      t32 = rhob ** 2
      zk(i) = -0.136284044462410474417D1 * rhoa * t8 * t10 * (0.1804D1 -
     > 0.804D0 / (0.1D1 + 0.209451650699150979027D-2 * t6 * sigmaaa / t1
     >2 / t14)) - 0.136284044462410474417D1 * rhob * t27 * t10 * (0.1804
     >D1 - 0.804D0 / (0.1D1 + 0.209451650699150979027D-2 * t6 * sigmabb
     >/ t30 / t32))

             endif
           enddo
         endif
       endif

      return
      end

c:PBEXsubrend
c:PW92Csubrstart

c    Generated: Mon Jun 20 16:05:17 CEST 2016

      subroutine dftacg_pw92c
     > (name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   tauc,tauo,upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,
     >                   vtauc,vtauo,vupsilonc,vupsilono)
      implicit double precision (a-h,o-z)
      logical fderiv,open
      integer igrad,npt
      character*(*) name
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)
      double precision tauc(*),tauo(*)
      double precision upsilonc(*),upsilono(*)
      double precision zk(*),vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)
      double precision vtauc(*),vtauo(*)
      double precision vupsilonc(*),vupsilono(*)
      include "common/cdft"
      include "common/tapes"
      parameter(tol=1d-12)
      pi=acos(-1d0)
      name='Automatically generated PW92C'
      igrad=0
       if(open) then
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

             rho = rhoa + rhob
      t1 = 0.1D1 / pi
      t2 = 0.1D1 / rho
      t3 = t1 * t2
      t4 = t3 ** (0.1D1 / 0.3D1)
      t6 = 0.1D1 + 0.194159335344114122552D0 * t4
      t7 = t3 ** (0.1D1 / 0.6D1)
      t10 = sqrt(t3)
      t12 = t4 ** 2
      t14 = 0.724010193431683113327D1 * t7 + 0.325955091942229212011D1 *
     > t4 + 0.141872281647966739112D1 * t10 + 0.406913004517529319387D0
     >* t12
      t17 = 0.1D1 + 0.160819794986925350668D2 / t14
      t18 = log(t17)
      t20 = 0.621814D-1 * t6 * t18
      t22 = 0.1D1 + 0.101077332976287768525D0 * t4
      t27 = 0.987212972256927209438D1 * t7 + 0.329180480994506259905D1 *
     > t4 + 0.762327521935289963194D0 * t10 + 0.410025070949612505036D0
     >* t12
      t30 = 0.1D1 + 0.296087499777934375166D2 / t27
      t31 = log(t30)
      t32 = t22 * t31
      t34 = rhoa - 0.1D1 * rhob
      t35 = t34 * t2
      t36 = 0.1D1 + t35
      t37 = t36 ** (0.1D1 / 0.3D1)
      t40 = 0.1D1 - 0.1D1 * t35
      t41 = t40 ** (0.1D1 / 0.3D1)
      t43 = t37 * t36 + t41 * t40 - 0.2D1
      t44 = t34 ** 2
      t45 = t44 ** 2
      t46 = rho ** 2
      t47 = t46 ** 2
      t48 = 0.1D1 / t47
      t49 = t45 * t48
      t51 = 0.1D1 - 0.1D1 * t49
      t54 = 0.379955235370239451738D-1 * t32 * t43 * t51
      t56 = 0.1D1 + 0.186690969707574028554D0 * t4
      t61 = 0.134579137143944477912D2 * t7 + 0.563098414909787598194D1 *
     > t4 + 0.291521471421917737271D1 * t10 + 0.516066464547863440989D0
     >* t12
      t64 = 0.1D1 + 0.321639589973850701335D2 / t61
      t65 = log(t64)
      t68 = -0.3109070D-1 * t56 * t65 + t20
      t69 = t68 * t43
      t71 = 0.192366105093153631974D1 * t69 * t49
      zk(i) = rho * (-t20 + t54 + t71)
      t78 = 0.133333333333333333333D1 * t37 * t2 - 0.1333333333333333333
     >33D1 * t41 * t2
      t82 = t44 * t34
      t86 = 0.151982094148095780695D0 * t32 * t43 * t82 * t48
      t92 = 0.769464420372614527896D1 * t69 * t82 * t48
      t95 = 0.500000000000000000000D0 * rho * (0.379955235370239451738D-
     >1 * t32 * t78 * t51 - t86 + 0.192366105093153631974D1 * t68 * t78
     >* t49 + t92)
      t96 = -t78
      t105 = 0.500000000000000000000D0 * rho * (0.379955235370239451738D
     >-1 * t32 * t96 * t51 + t86 + 0.192366105093153631974D1 * t68 * t96
     > * t49 - t92)
      t107 = 0.1D1 / t12 * t1
      t108 = 0.1D1 / t46
      t111 = 0.402436643158883263334D-2 * t107 * t108 * t18
      t112 = t14 ** 2
      t115 = t7 ** 2
      t116 = t115 ** 2
      t120 = 0.1D1 / t116 / t7 * t1 * t108
      t122 = t107 * t108
      t126 = 0.1D1 / t10 * t1 * t108
      t130 = 0.1D1 / t4 * t1 * t108
      t136 = 0.100000000000000000000D1 * t6 / t112 * (-0.120668365571947
     >185555D1 * t120 - 0.108651697314076404004D1 * t122 - 0.70936140823
     >9833695563D0 * t126 - 0.271275336345019546258D0 * t130) / t17
      t141 = t27 ** 2
      t161 = -0.133333333333333333333D1 * t37 * t34 * t108 + 0.133333333
     >333333333333D1 * t41 * t34 * t108
      t167 = 0.1D1 / t47 / rho
      t174 = t61 ** 2
      vrhoc(i) = vrhoc(i) + t95 + t105 - t20 + t54 + t71 + rho * (t111 +
     > t136 - 0.128016206138671616206D-2 * t122 * t31 * t43 * t51 - 0.11
     >2499995668310776915D1 * t22 / t141 * (-0.164535495376154534908D1 *
     > t120 - 0.109726826998168753302D1 * t122 - 0.381163760967644981600
     >D0 * t126 - 0.273350047299741670024D0 * t130) / t30 * t43 * t51 +
     >0.379955235370239451738D-1 * t32 * t161 * t51 + 0.1519820941480957
     >80695D0 * t32 * t43 * t45 * t167 + 0.192366105093153631974D1 * (0.
     >193478431062909061652D-2 * t107 * t108 * t65 + 0.10000000000000000
     >0000D1 * t56 / t174 * (-0.224298561906574129855D1 * t120 - 0.18769
     >9471636595866065D1 * t122 - 0.145760735710958868637D1 * t126 - 0.3
     >44044309698575627326D0 * t130) / t64 - t111 - t136) * t43 * t49 +
     >0.192366105093153631974D1 * t68 * t161 * t49 - 0.76946442037261452
     >7896D1 * t69 * t45 * t167)
      vrhoo(i) = vrhoo(i) + t95 - t105

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

             rho = rhoa + rhob
      t2 = 0.1D1 / rho
      t3 = 0.1D1 / pi * t2
      t4 = t3 ** (0.1D1 / 0.3D1)
      t7 = t3 ** (0.1D1 / 0.6D1)
      t10 = sqrt(t3)
      t12 = t4 ** 2
      t18 = log(0.1D1 + 0.160819794986925350668D2 / (0.72401019343168311
     >3327D1 * t7 + 0.325955091942229212011D1 * t4 + 0.14187228164796673
     >9112D1 * t10 + 0.406913004517529319387D0 * t12))
      t20 = 0.621814D-1 * (0.1D1 + 0.194159335344114122552D0 * t4) * t18
      t31 = log(0.1D1 + 0.296087499777934375166D2 / (0.98721297225692720
     >9438D1 * t7 + 0.329180480994506259905D1 * t4 + 0.76232752193528996
     >3194D0 * t10 + 0.410025070949612505036D0 * t12))
      t34 = rhoa - 0.1D1 * rhob
      t35 = t34 * t2
      t36 = 0.1D1 + t35
      t37 = t36 ** (0.1D1 / 0.3D1)
      t40 = 0.1D1 - 0.1D1 * t35
      t41 = t40 ** (0.1D1 / 0.3D1)
      t43 = t37 * t36 + t41 * t40 - 0.2D1
      t44 = t34 ** 2
      t45 = t44 ** 2
      t46 = rho ** 2
      t47 = t46 ** 2
      t49 = t45 / t47
      t65 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t7 + 0.563098414909787598194D1 * t4 + 0.29152147142191773
     >7271D1 * t10 + 0.516066464547863440989D0 * t12))
      zk(i) = rho * (-t20 + 0.379955235370239451738D-1 * (0.1D1 + 0.1010
     >77332976287768525D0 * t4) * t31 * t43 * (0.1D1 - 0.1D1 * t49) + 0.
     >192366105093153631974D1 * (-0.3109070D-1 * (0.1D1 + 0.186690969707
     >574028554D0 * t4) * t65 + t20) * t43 * t49)

             endif
           enddo
         endif
       else
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      t3 = 0.1D1 / pi
      t4 = 0.1D1 / rho
      t5 = t3 * t4
      t6 = t5 ** (0.1D1 / 0.3D1)
      t8 = 0.1D1 + 0.194159335344114122552D0 * t6
      t9 = t5 ** (0.1D1 / 0.6D1)
      t12 = sqrt(t5)
      t14 = t6 ** 2
      t16 = 0.724010193431683113327D1 * t9 + 0.325955091942229212011D1 *
     > t6 + 0.141872281647966739112D1 * t12 + 0.406913004517529319387D0
     >* t14
      t19 = 0.1D1 + 0.160819794986925350668D2 / t16
      t20 = log(t19)
      t22 = 0.621814D-1 * t8 * t20
      t24 = 0.1D1 + 0.101077332976287768525D0 * t6
      t29 = 0.987212972256927209438D1 * t9 + 0.329180480994506259905D1 *
     > t6 + 0.762327521935289963194D0 * t12 + 0.410025070949612505036D0
     >* t14
      t32 = 0.1D1 + 0.296087499777934375166D2 / t29
      t33 = log(t32)
      t34 = t24 * t33
      t36 = rhoa - 0.1D1 * rhob
      t37 = t36 * t4
      t38 = 0.1D1 + t37
      t39 = t38 ** (0.1D1 / 0.3D1)
      t42 = 0.1D1 - 0.1D1 * t37
      t43 = t42 ** (0.1D1 / 0.3D1)
      t45 = t39 * t38 + t43 * t42 - 0.2D1
      t46 = t36 ** 2
      t47 = t46 ** 2
      t48 = rho ** 2
      t49 = t48 ** 2
      t50 = 0.1D1 / t49
      t51 = t47 * t50
      t53 = 0.1D1 - 0.1D1 * t51
      t56 = 0.379955235370239451738D-1 * t34 * t45 * t53
      t58 = 0.1D1 + 0.186690969707574028554D0 * t6
      t63 = 0.134579137143944477912D2 * t9 + 0.563098414909787598194D1 *
     > t6 + 0.291521471421917737271D1 * t12 + 0.516066464547863440989D0
     >* t14
      t66 = 0.1D1 + 0.321639589973850701335D2 / t63
      t67 = log(t66)
      t70 = -0.3109070D-1 * t58 * t67 + t22
      t71 = t70 * t45
      t73 = 0.192366105093153631974D1 * t71 * t51
      zk(i) = rho * (-t22 + t56 + t73)
      t80 = 0.133333333333333333333D1 * t39 * t4 - 0.1333333333333333333
     >33D1 * t43 * t4
      t84 = t46 * t36
      t88 = 0.151982094148095780695D0 * t34 * t45 * t84 * t50
      t94 = 0.769464420372614527896D1 * t71 * t84 * t50
      t98 = -t80
      t109 = 0.1D1 / t14 * t3
      t110 = 0.1D1 / t48
      t113 = 0.402436643158883263334D-2 * t109 * t110 * t20
      t114 = t16 ** 2
      t117 = t9 ** 2
      t118 = t117 ** 2
      t122 = 0.1D1 / t118 / t9 * t3 * t110
      t124 = t109 * t110
      t128 = 0.1D1 / t12 * t3 * t110
      t132 = 0.1D1 / t6 * t3 * t110
      t138 = 0.100000000000000000000D1 * t8 / t114 * (-0.120668365571947
     >185555D1 * t122 - 0.108651697314076404004D1 * t124 - 0.70936140823
     >9833695563D0 * t128 - 0.271275336345019546258D0 * t132) / t19
      t143 = t29 ** 2
      t163 = -0.133333333333333333333D1 * t39 * t36 * t110 + 0.133333333
     >333333333333D1 * t43 * t36 * t110
      t169 = 0.1D1 / t49 / rho
      t176 = t63 ** 2
      vrhoc(i) = vrhoc(i) + 0.5D0 * rho * (0.379955235370239451738D-1 * 
     >t34 * t80 * t53 - t88 + 0.192366105093153631974D1 * t70 * t80 * t5
     >1 + t94) + 0.5D0 * rho * (0.379955235370239451738D-1 * t34 * t98 *
     > t53 + t88 + 0.192366105093153631974D1 * t70 * t98 * t51 - t94) -
     >t22 + t56 + t73 + rho * (t113 + t138 - 0.128016206138671616206D-2
     >* t124 * t33 * t45 * t53 - 0.112499995668310776915D1 * t24 / t143
     >* (-0.164535495376154534908D1 * t122 - 0.109726826998168753302D1 *
     > t124 - 0.381163760967644981600D0 * t128 - 0.273350047299741670024
     >D0 * t132) / t32 * t45 * t53 + 0.379955235370239451738D-1 * t34 *
     >t163 * t53 + 0.151982094148095780695D0 * t34 * t45 * t47 * t169 +
     >0.192366105093153631974D1 * (0.193478431062909061652D-2 * t109 * t
     >110 * t67 + 0.100000000000000000000D1 * t58 / t176 * (-0.224298561
     >906574129855D1 * t122 - 0.187699471636595866065D1 * t124 - 0.14576
     >0735710958868637D1 * t128 - 0.344044309698575627326D0 * t132) / t6
     >6 - t113 - t138) * t45 * t51 + 0.192366105093153631974D1 * t70 * t
     >163 * t51 - 0.769464420372614527896D1 * t71 * t47 * t169)

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      t4 = 0.1D1 / rho
      t5 = 0.1D1 / pi * t4
      t6 = t5 ** (0.1D1 / 0.3D1)
      t9 = t5 ** (0.1D1 / 0.6D1)
      t12 = sqrt(t5)
      t14 = t6 ** 2
      t20 = log(0.1D1 + 0.160819794986925350668D2 / (0.72401019343168311
     >3327D1 * t9 + 0.325955091942229212011D1 * t6 + 0.14187228164796673
     >9112D1 * t12 + 0.406913004517529319387D0 * t14))
      t22 = 0.621814D-1 * (0.1D1 + 0.194159335344114122552D0 * t6) * t20
      t33 = log(0.1D1 + 0.296087499777934375166D2 / (0.98721297225692720
     >9438D1 * t9 + 0.329180480994506259905D1 * t6 + 0.76232752193528996
     >3194D0 * t12 + 0.410025070949612505036D0 * t14))
      t36 = rhoa - 0.1D1 * rhob
      t37 = t36 * t4
      t38 = 0.1D1 + t37
      t39 = t38 ** (0.1D1 / 0.3D1)
      t42 = 0.1D1 - 0.1D1 * t37
      t43 = t42 ** (0.1D1 / 0.3D1)
      t45 = t39 * t38 + t43 * t42 - 0.2D1
      t46 = t36 ** 2
      t47 = t46 ** 2
      t48 = rho ** 2
      t49 = t48 ** 2
      t51 = t47 / t49
      t67 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t9 + 0.563098414909787598194D1 * t6 + 0.29152147142191773
     >7271D1 * t12 + 0.516066464547863440989D0 * t14))
      zk(i) = rho * (-t22 + 0.379955235370239451738D-1 * (0.1D1 + 0.1010
     >77332976287768525D0 * t6) * t33 * t45 * (0.1D1 - 0.1D1 * t51) + 0.
     >192366105093153631974D1 * (-0.3109070D-1 * (0.1D1 + 0.186690969707
     >574028554D0 * t6) * t67 + t22) * t45 * t51)

             endif
           enddo
         endif
       endif

      return
      end

c:PW92Csubrend
c:SCANCsubrstart

c    Generated: Mon Jun 20 16:05:18 CEST 2016

      subroutine dftacg_scanc
     > (name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   tauc,tauo,upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,
     >                   vtauc,vtauo,vupsilonc,vupsilono)
      implicit double precision (a-h,o-z)
      logical fderiv,open
      integer igrad,npt
      character*(*) name
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)
      double precision tauc(*),tauo(*)
      double precision upsilonc(*),upsilono(*)
      double precision zk(*),vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)
      double precision vtauc(*),vtauo(*)
      double precision vupsilonc(*),vupsilono(*)
      include "common/cdft"
      include "common/tapes"
      parameter(tol=1d-12)
      pi=acos(-1d0)
      name='Automatically generated SCANC'
      igrad=2
       if(open) then
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      sigmabb = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) - 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmabb
      taub = max(0.0D0, 0.500000000000000000000D0 * tauc(i) - 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taub
      t13 = 0.1D1 / pi
      t14 = 0.1D1 / rhob
      t15 = t13 * t14
      t16 = t15 ** (0.1D1 / 0.3D1)
      t18 = 0.1D1 + 0.186690969707574028554D0 * t16
      t19 = t15 ** (0.1D1 / 0.6D1)
      t22 = sqrt(t15)
      t24 = t16 ** 2
      t26 = 0.134579137143944477912D2 * t19 + 0.563098414909787598194D1 
     >* t16 + 0.291521471421917737271D1 * t22 + 0.516066464547863440989D
     >0 * t24
      t29 = 0.1D1 + 0.321639589973850701335D2 / t26
      t30 = log(t29)
      t31 = t18 * t30
      t32 = 0.3109070D-1 * t31
      t34 = exp(0.199998070181081341868D1 * t31)
      t35 = t34 - 0.1D1
      t37 = 0.1D1 + 0.908560296416069829442D-1 * t16
      t39 = 0.1D1 + 0.161542020702777215675D0 * t16
      t40 = 0.1D1 / t39
      t41 = t37 * t40
      t42 = 0.1D1 / t35
      t43 = pi ** 2
      t44 = t43 ** (0.1D1 / 0.3D1)
      t45 = t44 ** 2
      t47 = t41 * t42 * t45
      t48 = t43 * rhob
      t49 = t48 ** (0.1D1 / 0.3D1)
      t50 = t49 ** 2
      t51 = 0.1D1 / t50
      t52 = sigmabb * t51
      t53 = rhob ** 2
      t54 = 0.1D1 / t53
      t55 = 0.1D1 / t16
      t56 = t54 * t55
      t60 = 0.1D1 + 0.590527525871618703916D0 * t47 * t52 * t56
      t61 = t60 ** (0.1D1 / 0.4D1)
      t64 = 0.1D1 - 0.1D1 / t61
      t66 = 0.1D1 + t35 * t64
      t67 = log(t66)
      t68 = 0.155455000000000000000D-1 * t67
      t72 = 0.500000000000000000000D0 * taub - 0.125000000000000000000D0
     > * sigmabb * t14
      t73 = 0.1D1 / t45
      t74 = t72 * t73
      t75 = rhob ** (0.1D1 / 0.3D1)
      t76 = t75 ** 2
      t78 = 0.1D1 / t76 / rhob
      t79 = t74 * t78
      t81 = 0.1D1 - 0.100951144046229981050D1 * t79
      t82 = 0.1D1 / t81
      t83 = t78 * t82
      t86 = tanh(0.100951144046229981050D31 * t74 * t83)
      t88 = 0.500000000000000000000D0 + 0.500000000000000000000D0 * t86
      t92 = exp(-0.646087321895871878717D0 * t74 * t83 * t88)
      t95 = tanh(0.1D31 - 0.100951144046229981050D31 * t79)
      t96 = 0.500000000000000000000D0 * t95
      t97 = 0.500000000000000000000D0 + t96
      t100 = tanh(0.1D31 * t82)
      t102 = 0.500000000000000000000D0 - 0.500000000000000000000D0 * t10
     >0
      t105 = exp(0.15D1 * t82 * t102)
      t106 = 0.500000000000000000000D0 - t96
      t109 = t92 * t97 - 0.7D0 * t105 * t106
      t110 = t32 - t68
      t111 = t109 * t110
      zk(i) = rhob * (-t32 + t68 + t111)
      t114 = 0.155453500000000000000D-1 * t31
      t115 = 0.777275000000000000000D-2 * t67
      t116 = 0.500000000000000000000D0 * t111
      t118 = 0.1D1 / t24 * t13
      t120 = t118 * t54 * t30
      t121 = 0.193478431062909061652D-2 * t120
      t122 = t26 ** 2
      t125 = t19 ** 2
      t126 = t125 ** 2
      t144 = t18 / t122 * (-0.224298561906574129855D1 / t126 / t19 * t13
     > * t54 - 0.187699471636595866065D1 * t118 * t54 - 0.14576073571095
     >8868637D1 / t22 * t13 * t54 - 0.344044309698575627326D0 * t55 * t1
     >3 * t54) / t29
      t145 = 0.100000000000000000000D1 * t144
      t148 = -0.124459445539165071340D0 * t120 - 0.643272972886044192855
     >D2 * t144
      t152 = 0.1D1 / t61 / t60
      t154 = t53 * rhob
      t155 = 0.1D1 / t154
      t158 = t45 * sigmabb
      t162 = t39 ** 2
      t170 = t35 ** 2
      t175 = t51 * t54
      t192 = t53 ** 2
      t205 = 0.1D1 / t66
      t207 = 0.155455000000000000000D-1 * (t148 * t34 * t64 + 0.25000000
     >0000000000000D0 * t35 * t152 * (-0.178843287982588744867D-1 * t155
     > * t40 * t42 * t158 * t51 + 0.317983366033042788373D-1 * t37 / t16
     >2 * t42 * t158 * t51 * t155 - 0.590527525871618703916D0 * t41 / t1
     >70 * t45 * sigmabb * t175 * t55 * t148 * t34 - 0.39368501724774580
     >2611D0 * t47 * sigmabb / t50 / t48 * t56 * t43 - 0.118105505174323
     >740783D1 * t47 * t52 * t155 * t55 + 0.196842508623872901305D0 * t4
     >7 * t52 / t192 / t16 / t15 * t13)) * t205
      t210 = sigmabb / t76 / t154
      t211 = t73 * t82
      t216 = 0.1D1 / t76 / t53
      t217 = t216 * t82
      t221 = t81 ** 2
      t222 = 0.1D1 / t221
      t224 = t210 * t73
      t226 = t74 * t216
      t228 = -0.126188930057787476312D0 * t224 + 0.168251906743716635082
     >D1 * t226
      t232 = t86 ** 2
      t235 = t82 * (0.1D1 - 0.1D1 * t232)
      t251 = t95 ** 2
      t253 = 0.1D1 - 0.1D1 * t251
      t254 = t92 * t253
      t257 = -0.126188930057787476312D30 * t224 + 0.16825190674371663508
     >3D31 * t226
      t260 = t222 * t102
      t265 = t100 ** 2
      t268 = 0.1D1 / t221 / t81 * (0.1D1 - 0.1D1 * t265)
      t275 = t105 * t253
      t284 = 0.500000000000000000000D0 * rhob * (t121 + t145 + t207 + ((
     >-0.807609152369839848397D-1 * t210 * t211 * t88 + 0.10768122031597
     >8646453D1 * t74 * t217 * t88 + 0.646087321895871878717D0 * t79 * t
     >222 * t88 * t228 - 0.323043660947935939359D0 * t79 * t235 * (0.126
     >188930057787476312D30 * t210 * t211 - 0.168251906743716635083D31 *
     > t74 * t217 - 0.100951144046229981050D31 * t74 * t78 * t222 * t228
     >)) * t92 * t97 + 0.500000000000000000000D0 * t254 * t257 - 0.7D0 *
     > (-0.15D1 * t260 * t228 + 0.750000000000000000000D30 * t268 * t228
     >) * t105 * t106 + 0.350000000000000000000D0 * t275 * t257) * t110
     >+ t109 * (-t121 - t145 - t207))
      vrhoc(i) = vrhoc(i) - t114 + t115 + t116 + t284
      vrhoo(i) = vrhoo(i) + t114 - t115 - t116 - t284
      t294 = t216 * t73
      t295 = t82 * t88
      t300 = t72 / t44 / t43
      t303 = 0.1D1 / t75 / t192 * t222
      t340 = rhob * (0.229501141335931214043D-2 * t152 * t37 * t40 * t45
     > * t175 * t55 * t205 + ((0.807609152369839848397D-1 * t294 * t295
     >+ 0.815290678739413996018D-1 * t300 * t303 * t88 - 0.3230436609479
     >35939359D0 * t79 * t235 * (-0.126188930057787476312D30 * t294 * t8
     >2 - 0.127389168553033436878D30 * t300 * t303)) * t92 * t97 + 0.630
     >944650288937381559D29 * t254 * t294 - 0.7D0 * (-0.1892833950866812
     >14468D0 * t260 * t294 + 0.946416975433406072338D29 * t268 * t294)
     >* t105 * t106 + 0.441661255202256167090D29 * t275 * t294) * t110 -
     > 0.229501141335931214043D-2 * t109 * t152 * t41 * t45 * t51 * t56
     >* t205)
      t341 = 0.25D0 * t340
      vsigmacc(i) = vsigmacc(i) + t341
      vsigmaco(i) = vsigmaco(i) - 0.5D0 * t340
      vsigmaoo(i) = vsigmaoo(i) + t341
      t346 = t73 * t78
      t351 = 0.1D1 / t75 / t154 * t222
      t381 = 0.500000000000000000000D0 * rhob * ((-0.3230436609479359393
     >59D0 * t346 * t295 - 0.326116271495765598407D0 * t300 * t351 * t88
     > - 0.323043660947935939359D0 * t79 * t235 * (0.5047557202311499052
     >48D30 * t346 * t82 + 0.509556674212133747510D30 * t300 * t351)) *
     >t92 * t97 - 0.252377860115574952624D30 * t254 * t346 - 0.7D0 * (0.
     >757133580346724857871D0 * t260 * t346 - 0.378566790173362428937D30
     > * t268 * t346) * t105 * t106 - 0.176664502080902466837D30 * t275
     >* t346) * t110
      vtauc(i) = vtauc(i) + t381
      vtauo(i) = vtauo(i) - t381

               elseif(rhob.lt.tol) then
             rho = rhoa
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) + 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmaaa
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i) + 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taua
      t13 = 0.1D1 / pi
      t14 = 0.1D1 / rhoa
      t15 = t13 * t14
      t16 = t15 ** (0.1D1 / 0.3D1)
      t18 = 0.1D1 + 0.186690969707574028554D0 * t16
      t19 = t15 ** (0.1D1 / 0.6D1)
      t22 = sqrt(t15)
      t24 = t16 ** 2
      t26 = 0.134579137143944477912D2 * t19 + 0.563098414909787598194D1 
     >* t16 + 0.291521471421917737271D1 * t22 + 0.516066464547863440989D
     >0 * t24
      t29 = 0.1D1 + 0.321639589973850701335D2 / t26
      t30 = log(t29)
      t31 = t18 * t30
      t32 = 0.3109070D-1 * t31
      t34 = exp(0.199998070181081341868D1 * t31)
      t35 = t34 - 0.1D1
      t37 = 0.1D1 + 0.908560296416069829442D-1 * t16
      t39 = 0.1D1 + 0.161542020702777215675D0 * t16
      t40 = 0.1D1 / t39
      t41 = t37 * t40
      t42 = 0.1D1 / t35
      t43 = pi ** 2
      t44 = t43 ** (0.1D1 / 0.3D1)
      t45 = t44 ** 2
      t47 = t41 * t42 * t45
      t48 = t43 * rhoa
      t49 = t48 ** (0.1D1 / 0.3D1)
      t50 = t49 ** 2
      t51 = 0.1D1 / t50
      t52 = sigmaaa * t51
      t53 = rhoa ** 2
      t54 = 0.1D1 / t53
      t55 = 0.1D1 / t16
      t56 = t54 * t55
      t60 = 0.1D1 + 0.590527525871618703916D0 * t47 * t52 * t56
      t61 = t60 ** (0.1D1 / 0.4D1)
      t64 = 0.1D1 - 0.1D1 / t61
      t66 = 0.1D1 + t35 * t64
      t67 = log(t66)
      t68 = 0.155455000000000000000D-1 * t67
      t72 = 0.500000000000000000000D0 * taua - 0.125000000000000000000D0
     > * sigmaaa * t14
      t73 = 0.1D1 / t45
      t74 = t72 * t73
      t75 = rhoa ** (0.1D1 / 0.3D1)
      t76 = t75 ** 2
      t78 = 0.1D1 / t76 / rhoa
      t79 = t74 * t78
      t81 = 0.1D1 - 0.100951144046229981050D1 * t79
      t82 = 0.1D1 / t81
      t83 = t78 * t82
      t86 = tanh(0.100951144046229981050D31 * t74 * t83)
      t88 = 0.500000000000000000000D0 + 0.500000000000000000000D0 * t86
      t92 = exp(-0.646087321895871878717D0 * t74 * t83 * t88)
      t95 = tanh(0.1D31 - 0.100951144046229981050D31 * t79)
      t96 = 0.500000000000000000000D0 * t95
      t97 = 0.500000000000000000000D0 + t96
      t100 = tanh(0.1D31 * t82)
      t102 = 0.500000000000000000000D0 - 0.500000000000000000000D0 * t10
     >0
      t105 = exp(0.15D1 * t82 * t102)
      t106 = 0.500000000000000000000D0 - t96
      t109 = t92 * t97 - 0.7D0 * t105 * t106
      t110 = t32 - t68
      t111 = t109 * t110
      zk(i) = rhoa * (-t32 + t68 + t111)
      t114 = 0.155453500000000000000D-1 * t31
      t115 = 0.777275000000000000000D-2 * t67
      t116 = 0.500000000000000000000D0 * t111
      t118 = 0.1D1 / t24 * t13
      t120 = t118 * t54 * t30
      t121 = 0.193478431062909061652D-2 * t120
      t122 = t26 ** 2
      t125 = t19 ** 2
      t126 = t125 ** 2
      t144 = t18 / t122 * (-0.224298561906574129855D1 / t126 / t19 * t13
     > * t54 - 0.187699471636595866065D1 * t118 * t54 - 0.14576073571095
     >8868637D1 / t22 * t13 * t54 - 0.344044309698575627326D0 * t55 * t1
     >3 * t54) / t29
      t145 = 0.100000000000000000000D1 * t144
      t148 = -0.124459445539165071340D0 * t120 - 0.643272972886044192855
     >D2 * t144
      t152 = 0.1D1 / t61 / t60
      t154 = t53 * rhoa
      t155 = 0.1D1 / t154
      t158 = t45 * sigmaaa
      t162 = t39 ** 2
      t170 = t35 ** 2
      t175 = t51 * t54
      t192 = t53 ** 2
      t205 = 0.1D1 / t66
      t207 = 0.155455000000000000000D-1 * (t148 * t34 * t64 + 0.25000000
     >0000000000000D0 * t35 * t152 * (-0.178843287982588744867D-1 * t155
     > * t40 * t42 * t158 * t51 + 0.317983366033042788373D-1 * t37 / t16
     >2 * t42 * t158 * t51 * t155 - 0.590527525871618703916D0 * t41 / t1
     >70 * t45 * sigmaaa * t175 * t55 * t148 * t34 - 0.39368501724774580
     >2611D0 * t47 * sigmaaa / t50 / t48 * t56 * t43 - 0.118105505174323
     >740783D1 * t47 * t52 * t155 * t55 + 0.196842508623872901305D0 * t4
     >7 * t52 / t192 / t16 / t15 * t13)) * t205
      t210 = sigmaaa / t76 / t154
      t211 = t73 * t82
      t216 = 0.1D1 / t76 / t53
      t217 = t216 * t82
      t221 = t81 ** 2
      t222 = 0.1D1 / t221
      t224 = t210 * t73
      t226 = t74 * t216
      t228 = -0.126188930057787476312D0 * t224 + 0.168251906743716635082
     >D1 * t226
      t232 = t86 ** 2
      t235 = t82 * (0.1D1 - 0.1D1 * t232)
      t251 = t95 ** 2
      t253 = 0.1D1 - 0.1D1 * t251
      t254 = t92 * t253
      t257 = -0.126188930057787476312D30 * t224 + 0.16825190674371663508
     >3D31 * t226
      t260 = t222 * t102
      t265 = t100 ** 2
      t268 = 0.1D1 / t221 / t81 * (0.1D1 - 0.1D1 * t265)
      t275 = t105 * t253
      t284 = 0.500000000000000000000D0 * rhoa * (t121 + t145 + t207 + ((
     >-0.807609152369839848397D-1 * t210 * t211 * t88 + 0.10768122031597
     >8646453D1 * t74 * t217 * t88 + 0.646087321895871878717D0 * t79 * t
     >222 * t88 * t228 - 0.323043660947935939359D0 * t79 * t235 * (0.126
     >188930057787476312D30 * t210 * t211 - 0.168251906743716635083D31 *
     > t74 * t217 - 0.100951144046229981050D31 * t74 * t78 * t222 * t228
     >)) * t92 * t97 + 0.500000000000000000000D0 * t254 * t257 - 0.7D0 *
     > (-0.15D1 * t260 * t228 + 0.750000000000000000000D30 * t268 * t228
     >) * t105 * t106 + 0.350000000000000000000D0 * t275 * t257) * t110
     >+ t109 * (-t121 - t145 - t207))
      vrhoc(i) = vrhoc(i) - t114 + t115 + t116 + t284
      vrhoo(i) = vrhoo(i) - t114 + t115 + t116 + t284
      t294 = t216 * t73
      t295 = t82 * t88
      t300 = t72 / t44 / t43
      t303 = 0.1D1 / t75 / t192 * t222
      t340 = rhoa * (0.229501141335931214043D-2 * t152 * t37 * t40 * t45
     > * t175 * t55 * t205 + ((0.807609152369839848397D-1 * t294 * t295
     >+ 0.815290678739413996018D-1 * t300 * t303 * t88 - 0.3230436609479
     >35939359D0 * t79 * t235 * (-0.126188930057787476312D30 * t294 * t8
     >2 - 0.127389168553033436878D30 * t300 * t303)) * t92 * t97 + 0.630
     >944650288937381559D29 * t254 * t294 - 0.7D0 * (-0.1892833950866812
     >14468D0 * t260 * t294 + 0.946416975433406072338D29 * t268 * t294)
     >* t105 * t106 + 0.441661255202256167090D29 * t275 * t294) * t110 -
     > 0.229501141335931214043D-2 * t109 * t152 * t41 * t45 * t51 * t56
     >* t205)
      t341 = 0.25D0 * t340
      vsigmacc(i) = vsigmacc(i) + t341
      vsigmaco(i) = vsigmaco(i) + 0.5D0 * t340
      vsigmaoo(i) = vsigmaoo(i) + t341
      t346 = t73 * t78
      t351 = 0.1D1 / t75 / t154 * t222
      t381 = 0.500000000000000000000D0 * rhoa * ((-0.3230436609479359393
     >59D0 * t346 * t295 - 0.326116271495765598407D0 * t300 * t351 * t88
     > - 0.323043660947935939359D0 * t79 * t235 * (0.5047557202311499052
     >48D30 * t346 * t82 + 0.509556674212133747510D30 * t300 * t351)) *
     >t92 * t97 - 0.252377860115574952624D30 * t254 * t346 - 0.7D0 * (0.
     >757133580346724857871D0 * t260 * t346 - 0.378566790173362428937D30
     > * t268 * t346) * t105 * t106 - 0.176664502080902466837D30 * t275
     >* t346) * t110
      vtauc(i) = vtauc(i) + t381
      vtauo(i) = vtauo(i) + t381

               else
             rho = rhoa + rhob
      t2 = 0.250000000000000000000D0 * sigmacc(i)
      t4 = 0.250000000000000000000D0 * sigmaoo(i)
      t6 = 0.500000000000000000000D0 * sigmaco(i)
      sigmaaa = max(0.0D0, t2 + t4 + t6)
      sigmaab = t2 - t4
      sigmabb = max(0.0D0, t2 + t4 - t6)
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t11 = 0.500000000000000000000D0 * tauc(i)
      t13 = 0.500000000000000000000D0 * tauo(i)
      taua = max(0.0D0, t11 + t13)
      taub = max(0.0D0, t11 - t13)
      tau = taua + taub
      t16 = 0.1D1 / pi
      t17 = 0.1D1 / rho
      t18 = t16 * t17
      t19 = t18 ** (0.1D1 / 0.3D1)
      t21 = 0.1D1 + 0.194159335344114122552D0 * t19
      t22 = t18 ** (0.1D1 / 0.6D1)
      t25 = sqrt(t18)
      t27 = t19 ** 2
      t29 = 0.724010193431683113327D1 * t22 + 0.325955091942229212011D1 
     >* t19 + 0.141872281647966739112D1 * t25 + 0.406913004517529319387D
     >0 * t27
      t32 = 0.1D1 + 0.160819794986925350668D2 / t29
      t33 = log(t32)
      t35 = 0.621814D-1 * t21 * t33
      t37 = 0.1D1 + 0.101077332976287768525D0 * t19
      t42 = 0.987212972256927209438D1 * t22 + 0.329180480994506259905D1 
     >* t19 + 0.762327521935289963194D0 * t25 + 0.410025070949612505036D
     >0 * t27
      t45 = 0.1D1 + 0.296087499777934375166D2 / t42
      t46 = log(t45)
      t47 = t37 * t46
      t49 = rhoa - 0.1D1 * rhob
      t50 = t49 * t17
      t51 = 0.1D1 + t50
      t52 = t51 ** (0.1D1 / 0.3D1)
      t53 = t52 * t51
      t55 = 0.1D1 - 0.1D1 * t50
      t56 = t55 ** (0.1D1 / 0.3D1)
      t57 = t56 * t55
      t58 = t53 + t57 - 0.2D1
      t59 = t49 ** 2
      t60 = t59 ** 2
      t61 = rho ** 2
      t62 = t61 ** 2
      t63 = 0.1D1 / t62
      t64 = t60 * t63
      t66 = 0.1D1 - 0.1D1 * t64
      t69 = 0.379955235370239451738D-1 * t47 * t58 * t66
      t71 = 0.1D1 + 0.186690969707574028554D0 * t19
      t76 = 0.134579137143944477912D2 * t22 + 0.563098414909787598194D1 
     >* t19 + 0.291521471421917737271D1 * t25 + 0.516066464547863440989D
     >0 * t27
      t79 = 0.1D1 + 0.321639589973850701335D2 / t76
      t80 = log(t79)
      t83 = -0.3109070D-1 * t71 * t80 + t35
      t84 = t83 * t58
      t86 = 0.192366105093153631974D1 * t84 * t64
      t87 = t52 ** 2
      t89 = t56 ** 2
      t91 = 0.500000000000000000000D0 * t87 + 0.500000000000000000000D0 
     >* t89
      t92 = t91 ** 2
      t93 = t92 * t91
      t94 = -t35 + t69 + t86
      t95 = 0.1D1 / t93
      t98 = exp(-0.321636486443022096427D2 * t94 * t95)
      t99 = t98 - 0.1D1
      t101 = 0.1D1 + 0.908560296416069829442D-1 * t19
      t103 = 0.1D1 + 0.161542020702777215675D0 * t19
      t104 = 0.1D1 / t103
      t105 = t101 * t104
      t106 = 0.1D1 / t99
      t107 = pi ** 2
      t108 = t107 ** (0.1D1 / 0.3D1)
      t109 = t108 ** 2
      t110 = t106 * t109
      t111 = t105 * t110
      t112 = t107 * rho
      t113 = t112 ** (0.1D1 / 0.3D1)
      t114 = t113 ** 2
      t115 = 0.1D1 / t114
      t116 = sigma * t115
      t117 = 0.1D1 / t61
      t118 = 0.1D1 / t92
      t120 = 0.1D1 / t19
      t125 = 0.1D1 + 0.372009030193995856361D0 * t111 * t116 * t117 * t1
     >18 * t120
      t126 = t125 ** (0.1D1 / 0.4D1)
      t129 = 0.1D1 - 0.1D1 / t126
      t131 = 0.1D1 + t99 * t129
      t132 = log(t131)
      t134 = 0.31091D-1 * t93 * t132
      t138 = 0.500000000000000000000D0 * tau - 0.125000000000000000000D0
     > * sigma * t17
      t139 = 0.1D1 / t109
      t140 = t138 * t139
      t141 = rho ** (0.1D1 / 0.3D1)
      t142 = t141 ** 2
      t144 = 0.1D1 / t142 / rho
      t145 = t140 * t144
      t150 = 0.500000000000000000000D0 * t87 * t51 + 0.50000000000000000
     >0000D0 * t89 * t55
      t151 = 0.1D1 / t150
      t152 = t144 * t151
      t153 = t140 * t152
      t155 = 0.1D1 - 0.160249952256378709147D1 * t153
      t156 = 0.1D1 / t155
      t157 = t151 * t156
      t161 = tanh(0.160249952256378709147D31 * t140 * t152 * t156)
      t163 = 0.500000000000000000000D0 + 0.500000000000000000000D0 * t16
     >1
      t164 = t157 * t163
      t167 = exp(-0.102559969444082373854D1 * t145 * t164)
      t170 = tanh(0.1D31 - 0.160249952256378709147D31 * t153)
      t171 = 0.500000000000000000000D0 * t170
      t172 = 0.500000000000000000000D0 + t171
      t175 = tanh(0.1D31 * t156)
      t177 = 0.500000000000000000000D0 - 0.500000000000000000000D0 * t17
     >5
      t180 = exp(0.15D1 * t156 * t177)
      t181 = 0.500000000000000000000D0 - t171
      t184 = t167 * t172 - 0.7D0 * t180 * t181
      t187 = 0.1D1 + 0.847380836474276614068D-1 * t22 + 0.11406156817236
     >9822458D0 * t19
      t188 = 0.1D1 / t187
      t191 = exp(0.100000000000000000000D1 * t188)
      t192 = t191 - 0.1D1
      t195 = 0.1D1 + 0.615484811627254218514D-1 * t116 * t117
      t196 = t195 ** (0.1D1 / 0.4D1)
      t199 = 0.1D1 - 0.1D1 / t196
      t201 = 0.1D1 + t192 * t199
      t202 = log(t201)
      t204 = -0.285764D-1 * t188 + 0.285764D-1 * t202
      t207 = 0.33631D1 - 0.118155000000000000000D1 * t53 - 0.11815500000
     >0000000000D1 * t57
      t208 = t204 * t207
      t209 = t60 ** 2
      t210 = t209 * t60
      t211 = t62 ** 2
      t213 = 0.1D1 / t211 / t62
      t216 = 0.1D1 - 0.1D1 * t210 * t213
      t218 = t208 * t216 + t35 - t69 - t86 - t134
      t219 = t184 * t218
      zk(i) = rho * (-t35 + t69 + t86 + t134 + t219)
      t222 = t52 * t17
      t224 = t56 * t17
      t226 = 0.133333333333333333333D1 * t222 - 0.133333333333333333333D
     >1 * t224
      t229 = 0.379955235370239451738D-1 * t47 * t226 * t66
      t230 = t59 * t49
      t234 = 0.151982094148095780695D0 * t47 * t58 * t230 * t63
      t237 = 0.192366105093153631974D1 * t83 * t226 * t64
      t240 = 0.769464420372614527896D1 * t84 * t230 * t63
      t241 = t92 * t132
      t242 = 0.1D1 / t52
      t245 = 0.1D1 / t56
      t248 = 0.333333333333333333333D0 * t242 * t17 - 0.3333333333333333
     >33333D0 * t245 * t17
      t250 = 0.93273D-1 * t241 * t248
      t254 = t92 ** 2
      t256 = t94 / t254
      t259 = -0.321636486443022096427D2 * (t229 - t234 + t237 + t240) * 
     >t95 + 0.964909459329066289281D2 * t256 * t248
      t263 = 0.1D1 / t126 / t125
      t264 = t99 * t263
      t265 = t99 ** 2
      t269 = t105 / t265 * t109 * sigma
      t270 = t115 * t117
      t271 = t270 * t118
      t278 = t105 * t110 * sigma
      t279 = t95 * t120
      t289 = 0.1D1 / t131
      t291 = 0.31091D-1 * t93 * (t259 * t98 * t129 + 0.25000000000000000
     >0000D0 * t264 * (-0.372009030193995856361D0 * t269 * t271 * t120 *
     > t259 * t98 - 0.744018060387991712722D0 * t278 * t270 * t279 * t24
     >8)) * t289
      t292 = t150 ** 2
      t293 = 0.1D1 / t292
      t294 = t293 * t156
      t299 = 0.833333333333333333333D0 * t87 * t17 - 0.83333333333333333
     >3333D0 * t89 * t17
      t300 = t163 * t299
      t304 = t138 ** 2
      t306 = 0.1D1 / t108 / t107
      t308 = t61 * rho
      t310 = 0.1D1 / t141 / t308
      t311 = t304 * t306 * t310
      t314 = t155 ** 2
      t315 = 0.1D1 / t314
      t316 = 0.1D1 / t292 / t150 * t315
      t320 = t161 ** 2
      t322 = 0.1D1 - 0.1D1 * t320
      t337 = t170 ** 2
      t339 = 0.1D1 - 0.1D1 * t337
      t340 = t167 * t339
      t341 = t340 * t138
      t342 = t139 * t144
      t344 = t342 * t293 * t299
      t347 = t315 * t177
      t348 = t347 * t138
      t353 = t175 ** 2
      t356 = 0.1D1 / t314 / t155 * (0.1D1 - 0.1D1 * t353)
      t357 = t356 * t138
      t364 = t180 * t339
      t365 = t364 * t138
      t372 = -0.157540000000000000000D1 * t222 + 0.157540000000000000000
     >D1 * t224
      t378 = 0.12D2 * t208 * t209 * t230 * t213
      t383 = 0.500000000000000000000D0 * rho * (t229 - t234 + t237 + t24
     >0 + t250 + t291 + ((0.102559969444082373854D1 * t145 * t294 * t300
     > + 0.164352302068298596704D1 * t311 * t316 * t300 - 0.512799847220
     >411869270D0 * t145 * t157 * t322 * (-0.160249952256378709147D31 *
     >t145 * t294 * t299 - 0.256800471981716557348D31 * t311 * t316 * t2
     >99)) * t167 * t172 + 0.801249761281893545733D30 * t341 * t344 - 0.
     >7D0 * (-0.240374928384568063720D1 * t348 * t344 + 0.12018746419228
     >4031860D31 * t357 * t344) * t180 * t181 + 0.560874832897325482012D
     >30 * t365 * t344) * t218 + t184 * (t204 * t372 * t216 - t378 - t22
     >9 + t234 - t237 - t240 - t250 - t291))
      t384 = -t226
      t387 = 0.379955235370239451738D-1 * t47 * t384 * t66
      t390 = 0.192366105093153631974D1 * t83 * t384 * t64
      t391 = -t248
      t393 = 0.93273D-1 * t241 * t391
      t399 = -0.321636486443022096427D2 * (t387 + t234 + t390 - t240) * 
     >t95 + 0.964909459329066289281D2 * t256 * t391
      t417 = 0.31091D-1 * t93 * (t399 * t98 * t129 + 0.25000000000000000
     >0000D0 * t264 * (-0.372009030193995856361D0 * t269 * t271 * t120 *
     > t399 * t98 - 0.744018060387991712722D0 * t278 * t270 * t279 * t39
     >1)) * t289
      t418 = -t299
      t419 = t163 * t418
      t441 = t342 * t293 * t418
      t463 = 0.500000000000000000000D0 * rho * (t387 + t234 + t390 - t24
     >0 + t393 + t417 + ((0.102559969444082373854D1 * t145 * t294 * t419
     > + 0.164352302068298596704D1 * t311 * t316 * t419 - 0.512799847220
     >411869270D0 * t145 * t157 * t322 * (-0.160249952256378709147D31 *
     >t145 * t294 * t418 - 0.256800471981716557348D31 * t311 * t316 * t4
     >18)) * t167 * t172 + 0.801249761281893545733D30 * t341 * t441 - 0.
     >7D0 * (-0.240374928384568063720D1 * t348 * t441 + 0.12018746419228
     >4031860D31 * t357 * t441) * t180 * t181 + 0.560874832897325482012D
     >30 * t365 * t441) * t218 + t184 * (-t204 * t372 * t216 + t378 - t3
     >87 - t234 - t390 + t240 - t393 - t417))
      t465 = 0.1D1 / t27 * t16
      t468 = 0.402436643158883263334D-2 * t465 * t117 * t33
      t469 = t29 ** 2
      t472 = t22 ** 2
      t473 = t472 ** 2
      t477 = 0.1D1 / t473 / t22 * t16 * t117
      t479 = t465 * t117
      t483 = 0.1D1 / t25 * t16 * t117
      t486 = t120 * t16 * t117
      t492 = 0.100000000000000000000D1 * t21 / t469 * (-0.12066836557194
     >7185555D1 * t477 - 0.108651697314076404004D1 * t479 - 0.7093614082
     >39833695563D0 * t483 - 0.271275336345019546258D0 * t486) / t32
      t496 = 0.128016206138671616206D-2 * t479 * t46 * t58 * t66
      t497 = t42 ** 2
      t510 = 0.112499995668310776915D1 * t37 / t497 * (-0.16453549537615
     >4534908D1 * t477 - 0.109726826998168753302D1 * t479 - 0.3811637609
     >67644981600D0 * t483 - 0.273350047299741670024D0 * t486) / t45 * t
     >58 * t66
      t512 = t52 * t49 * t117
      t515 = t56 * t49 * t117
      t517 = -0.133333333333333333333D1 * t512 + 0.133333333333333333333
     >D1 * t515
      t520 = 0.379955235370239451738D-1 * t47 * t517 * t66
      t522 = t62 * rho
      t523 = 0.1D1 / t522
      t526 = 0.151982094148095780695D0 * t47 * t58 * t60 * t523
      t530 = t76 ** 2
      t545 = 0.192366105093153631974D1 * (0.193478431062909061652D-2 * t
     >465 * t117 * t80 + 0.100000000000000000000D1 * t71 / t530 * (-0.22
     >4298561906574129855D1 * t477 - 0.187699471636595866065D1 * t479 -
     >0.145760735710958868637D1 * t483 - 0.344044309698575627326D0 * t48
     >6) / t79 - t468 - t492) * t58 * t64
      t548 = 0.192366105093153631974D1 * t83 * t517 * t64
      t551 = 0.769464420372614527896D1 * t84 * t60 * t523
      t558 = -0.333333333333333333333D0 * t242 * t49 * t117 + 0.33333333
     >3333333333333D0 * t245 * t49 * t117
      t560 = 0.93273D-1 * t241 * t558
      t566 = -0.321636486443022096427D2 * (t468 + t492 - t496 - t510 + t
     >520 + t526 + t545 + t548 - t551) * t95 + 0.964909459329066289281D2
     > * t256 * t558
      t569 = 0.1D1 / t308
      t577 = t103 ** 2
      t581 = t569 * t118
      t591 = 0.1D1 / t114 / t112
      t620 = 0.31091D-1 * t93 * (t566 * t98 * t129 + 0.25000000000000000
     >0000D0 * t264 * (-0.112664211580837182141D-1 * t569 * t104 * t106
     >* t109 * sigma * t115 * t118 + 0.200316968190728509847D-1 * t101 /
     > t577 * t110 * t116 * t581 - 0.372009030193995856361D0 * t269 * t2
     >71 * t120 * t566 * t98 - 0.248006020129330570907D0 * t278 * t591 *
     > t117 * t118 * t120 * t107 - 0.744018060387991712722D0 * t111 * t1
     >16 * t581 * t120 - 0.744018060387991712722D0 * t278 * t270 * t279
     >* t558 + 0.124003010064665285454D0 * t278 * t115 * t63 * t118 / t1
     >9 / t18 * t16)) * t289
      t623 = sigma / t142 / t308
      t628 = 0.1D1 / t142 / t61
      t638 = -0.833333333333333333333D0 * t87 * t49 * t117 + 0.833333333
     >333333333333D0 * t89 * t49 * t117
      t643 = t151 * t315
      t644 = t139 * t151
      t645 = t623 * t644
      t647 = t628 * t151
      t648 = t140 * t647
      t652 = t140 * t144 * t293 * t638
      t654 = -0.200312440320473386433D0 * t645 + 0.267083253760631181911
     >D1 * t648 + 0.160249952256378709147D1 * t652
      t682 = -0.200312440320473386433D30 * t645 + 0.26708325376063118191
     >1D31 * t648 + 0.160249952256378709147D31 * t652
      t697 = t187 ** 2
      t702 = 0.1D1 / t697 * (-0.141230139412379435678D-1 * t477 - 0.3802
     >05227241232741527D-1 * t479)
      t709 = t192 / t196 / t195
      t720 = 0.1D1 / t201
      t736 = (0.285764D-1 * t702 + 0.285764D-1 * (-0.1000000000000000000
     >00D1 * t702 * t191 * t199 + 0.250000000000000000000D0 * t709 * (-0
     >.410323207751502812342D-1 * sigma * t591 * t117 * t107 - 0.1230969
     >62325450843703D0 * t116 * t569)) * t720) * t207 * t216 + t204 * (0
     >.157540000000000000000D1 * t512 - 0.157540000000000000000D1 * t515
     >) * t216 + 0.12D2 * t208 * t210 / t211 / t522 - t468 - t492 + t496
     > + t510 - t520 - t526 - t545 - t548 + t551 - t560 - t620
      t738 = t468 + t492 - t496 - t510 + t520 + t526 + t545 + t548 - t55
     >1 + t560 + t620 + ((-0.128199961805102967317D0 * t623 * t139 * t16
     >4 + 0.170933282406803956424D1 * t140 * t628 * t164 + 0.10255996944
     >4082373854D1 * t145 * t294 * t163 * t638 + 0.102559969444082373854
     >D1 * t145 * t643 * t163 * t654 - 0.512799847220411869270D0 * t145
     >* t157 * t322 * (0.200312440320473386433D30 * t623 * t644 * t156 -
     > 0.267083253760631181911D31 * t140 * t647 * t156 - 0.1602499522563
     >78709147D31 * t145 * t294 * t638 - 0.160249952256378709147D31 * t1
     >45 * t643 * t654)) * t167 * t172 + 0.500000000000000000000D0 * t34
     >0 * t682 - 0.7D0 * (-0.15D1 * t347 * t654 + 0.75000000000000000000
     >0D30 * t356 * t654) * t180 * t181 + 0.350000000000000000000D0 * t3
     >64 * t682) * t218 + t184 * t736
      vrhoc(i) = vrhoc(i) + t383 + t463 - t35 + t69 + t86 + t134 + t219 
     >+ rho * t738
      vrhoo(i) = vrhoo(i) + t383 - t463
      t749 = 0.289153318944038129253D-2 * t91 * t263 * t105 * t109 * t11
     >5 * t117 * t120 * t289
      t750 = t628 * t139
      t753 = t138 * t306
      t755 = 0.1D1 / t141 / t62
      t758 = t293 * t315 * t163
      t775 = t750 * t151
      vsigmacc(i) = vsigmacc(i) + rho * (t749 + ((0.12819996180510296731
     >7D0 * t750 * t164 + 0.205440377585373245880D0 * t753 * t755 * t758
     > - 0.512799847220411869270D0 * t145 * t157 * t322 * (-0.2003124403
     >20473386433D30 * t750 * t157 - 0.321000589977145696686D30 * t753 *
     > t755 * t293 * t315)) * t167 * t172 + 0.100156220160236693217D30 *
     > t340 * t775 - 0.7D0 * (-0.300468660480710079650D0 * t347 * t775 +
     > 0.150234330240355039825D30 * t356 * t775) * t180 * t181 + 0.70109
     >3541121656852518D29 * t364 * t775) * t218 + t184 * (0.439708504274
     >626686249D-3 * t709 * t115 * t117 * t720 * t207 * t216 - t749))
      vsigmaco(i) = vsigmaco(i)
      vsigmaoo(i) = vsigmaoo(i)
      t820 = t342 * t151
      vtauc(i) = vtauc(i) + rho * ((-0.512799847220411869270D0 * t342 * 
     >t164 - 0.821761510341492983518D0 * t753 * t310 * t758 - 0.51279984
     >7220411869270D0 * t145 * t157 * t322 * (0.801249761281893545733D30
     > * t342 * t157 + 0.128400235990858278675D31 * t753 * t310 * t293 *
     > t315)) * t167 * t172 - 0.400624880640946772867D30 * t340 * t820 -
     > 0.7D0 * (0.120187464192284031860D1 * t347 * t820 - 0.600937320961
     >420159300D30 * t356 * t820) * t180 * t181 - 0.28043741644866274100
     >7D30 * t364 * t820) * t218
      vtauo(i) = vtauo(i)

               endif
             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      sigmabb = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) - 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmabb
      taub = max(0.0D0, 0.500000000000000000000D0 * tauc(i) - 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taub
      t14 = 0.1D1 / rhob
      t15 = 0.1D1 / pi * t14
      t16 = t15 ** (0.1D1 / 0.3D1)
      t19 = t15 ** (0.1D1 / 0.6D1)
      t22 = sqrt(t15)
      t24 = t16 ** 2
      t30 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t19 + 0.563098414909787598194D1 * t16 + 0.291521471421917
     >737271D1 * t22 + 0.516066464547863440989D0 * t24))
      t31 = (0.1D1 + 0.186690969707574028554D0 * t16) * t30
      t32 = 0.3109070D-1 * t31
      t34 = exp(0.199998070181081341868D1 * t31)
      t35 = t34 - 0.1D1
      t43 = pi ** 2
      t44 = t43 ** (0.1D1 / 0.3D1)
      t45 = t44 ** 2
      t49 = (t43 * rhob) ** (0.1D1 / 0.3D1)
      t50 = t49 ** 2
      t53 = rhob ** 2
      t61 = (0.1D1 + 0.590527525871618703916D0 * (0.1D1 + 0.908560296416
     >069829442D-1 * t16) / (0.1D1 + 0.161542020702777215675D0 * t16) /
     >t35 * t45 * sigmabb / t50 / t53 / t16) ** (0.1D1 / 0.4D1)
      t67 = log(0.1D1 + t35 * (0.1D1 - 0.1D1 / t61))
      t68 = 0.155455000000000000000D-1 * t67
      t74 = (0.500000000000000000000D0 * taub - 0.125000000000000000000D
     >0 * sigmabb * t14) / t45
      t75 = rhob ** (0.1D1 / 0.3D1)
      t76 = t75 ** 2
      t78 = 0.1D1 / t76 / rhob
      t79 = t74 * t78
      t82 = 0.1D1 / (0.1D1 - 0.100951144046229981050D1 * t79)
      t83 = t78 * t82
      t86 = tanh(0.100951144046229981050D31 * t74 * t83)
      t92 = exp(-0.646087321895871878717D0 * t74 * t83 * (0.500000000000
     >000000000D0 + 0.500000000000000000000D0 * t86))
      t95 = tanh(0.1D31 - 0.100951144046229981050D31 * t79)
      t96 = 0.500000000000000000000D0 * t95
      t100 = tanh(0.1D31 * t82)
      t105 = exp(0.15D1 * t82 * (0.500000000000000000000D0 - 0.500000000
     >000000000000D0 * t100))
      zk(i) = rhob * (-t32 + t68 + (t92 * (0.500000000000000000000D0 + t
     >96) - 0.7D0 * t105 * (0.500000000000000000000D0 - t96)) * (t32 - t
     >68))

               elseif(rhob.lt.tol) then
             rho = rhoa
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) + 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmaaa
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i) + 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taua
      t14 = 0.1D1 / rhoa
      t15 = 0.1D1 / pi * t14
      t16 = t15 ** (0.1D1 / 0.3D1)
      t19 = t15 ** (0.1D1 / 0.6D1)
      t22 = sqrt(t15)
      t24 = t16 ** 2
      t30 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t19 + 0.563098414909787598194D1 * t16 + 0.291521471421917
     >737271D1 * t22 + 0.516066464547863440989D0 * t24))
      t31 = (0.1D1 + 0.186690969707574028554D0 * t16) * t30
      t32 = 0.3109070D-1 * t31
      t34 = exp(0.199998070181081341868D1 * t31)
      t35 = t34 - 0.1D1
      t43 = pi ** 2
      t44 = t43 ** (0.1D1 / 0.3D1)
      t45 = t44 ** 2
      t49 = (t43 * rhoa) ** (0.1D1 / 0.3D1)
      t50 = t49 ** 2
      t53 = rhoa ** 2
      t61 = (0.1D1 + 0.590527525871618703916D0 * (0.1D1 + 0.908560296416
     >069829442D-1 * t16) / (0.1D1 + 0.161542020702777215675D0 * t16) /
     >t35 * t45 * sigmaaa / t50 / t53 / t16) ** (0.1D1 / 0.4D1)
      t67 = log(0.1D1 + t35 * (0.1D1 - 0.1D1 / t61))
      t68 = 0.155455000000000000000D-1 * t67
      t74 = (0.500000000000000000000D0 * taua - 0.125000000000000000000D
     >0 * sigmaaa * t14) / t45
      t75 = rhoa ** (0.1D1 / 0.3D1)
      t76 = t75 ** 2
      t78 = 0.1D1 / t76 / rhoa
      t79 = t74 * t78
      t82 = 0.1D1 / (0.1D1 - 0.100951144046229981050D1 * t79)
      t83 = t78 * t82
      t86 = tanh(0.100951144046229981050D31 * t74 * t83)
      t92 = exp(-0.646087321895871878717D0 * t74 * t83 * (0.500000000000
     >000000000D0 + 0.500000000000000000000D0 * t86))
      t95 = tanh(0.1D31 - 0.100951144046229981050D31 * t79)
      t96 = 0.500000000000000000000D0 * t95
      t100 = tanh(0.1D31 * t82)
      t105 = exp(0.15D1 * t82 * (0.500000000000000000000D0 - 0.500000000
     >000000000000D0 * t100))
      zk(i) = rhoa * (-t32 + t68 + (t92 * (0.500000000000000000000D0 + t
     >96) - 0.7D0 * t105 * (0.500000000000000000000D0 - t96)) * (t32 - t
     >68))

               else
             rho = rhoa + rhob
      t2 = 0.250000000000000000000D0 * sigmacc(i)
      t4 = 0.250000000000000000000D0 * sigmaoo(i)
      t6 = 0.500000000000000000000D0 * sigmaco(i)
      sigmaaa = max(0.0D0, t2 + t4 + t6)
      sigmaab = t2 - t4
      sigmabb = max(0.0D0, t2 + t4 - t6)
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t11 = 0.500000000000000000000D0 * tauc(i)
      t13 = 0.500000000000000000000D0 * tauo(i)
      taua = max(0.0D0, t11 + t13)
      taub = max(0.0D0, t11 - t13)
      tau = taua + taub
      t17 = 0.1D1 / rho
      t18 = 0.1D1 / pi * t17
      t19 = t18 ** (0.1D1 / 0.3D1)
      t22 = t18 ** (0.1D1 / 0.6D1)
      t25 = sqrt(t18)
      t27 = t19 ** 2
      t33 = log(0.1D1 + 0.160819794986925350668D2 / (0.72401019343168311
     >3327D1 * t22 + 0.325955091942229212011D1 * t19 + 0.141872281647966
     >739112D1 * t25 + 0.406913004517529319387D0 * t27))
      t35 = 0.621814D-1 * (0.1D1 + 0.194159335344114122552D0 * t19) * t3
     >3
      t46 = log(0.1D1 + 0.296087499777934375166D2 / (0.98721297225692720
     >9438D1 * t22 + 0.329180480994506259905D1 * t19 + 0.762327521935289
     >963194D0 * t25 + 0.410025070949612505036D0 * t27))
      t49 = rhoa - 0.1D1 * rhob
      t50 = t49 * t17
      t51 = 0.1D1 + t50
      t52 = t51 ** (0.1D1 / 0.3D1)
      t53 = t52 * t51
      t55 = 0.1D1 - 0.1D1 * t50
      t56 = t55 ** (0.1D1 / 0.3D1)
      t57 = t56 * t55
      t58 = t53 + t57 - 0.2D1
      t59 = t49 ** 2
      t60 = t59 ** 2
      t61 = rho ** 2
      t62 = t61 ** 2
      t64 = t60 / t62
      t69 = 0.379955235370239451738D-1 * (0.1D1 + 0.10107733297628776852
     >5D0 * t19) * t46 * t58 * (0.1D1 - 0.1D1 * t64)
      t80 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t22 + 0.563098414909787598194D1 * t19 + 0.291521471421917
     >737271D1 * t25 + 0.516066464547863440989D0 * t27))
      t86 = 0.192366105093153631974D1 * (-0.3109070D-1 * (0.1D1 + 0.1866
     >90969707574028554D0 * t19) * t80 + t35) * t58 * t64
      t87 = t52 ** 2
      t89 = t56 ** 2
      t91 = 0.500000000000000000000D0 * t87 + 0.500000000000000000000D0 
     >* t89
      t92 = t91 ** 2
      t93 = t92 * t91
      t98 = exp(-0.321636486443022096427D2 * (-t35 + t69 + t86) / t93)
      t99 = t98 - 0.1D1
      t107 = pi ** 2
      t108 = t107 ** (0.1D1 / 0.3D1)
      t109 = t108 ** 2
      t113 = (t107 * rho) ** (0.1D1 / 0.3D1)
      t114 = t113 ** 2
      t116 = sigma / t114
      t117 = 0.1D1 / t61
      t126 = (0.1D1 + 0.372009030193995856361D0 * (0.1D1 + 0.90856029641
     >6069829442D-1 * t19) / (0.1D1 + 0.161542020702777215675D0 * t19) /
     > t99 * t109 * t116 * t117 / t92 / t19) ** (0.1D1 / 0.4D1)
      t132 = log(0.1D1 + t99 * (0.1D1 - 0.1D1 / t126))
      t134 = 0.31091D-1 * t93 * t132
      t140 = (0.500000000000000000000D0 * tau - 0.125000000000000000000D
     >0 * sigma * t17) / t109
      t141 = rho ** (0.1D1 / 0.3D1)
      t142 = t141 ** 2
      t144 = 0.1D1 / t142 / rho
      t151 = 0.1D1 / (0.500000000000000000000D0 * t87 * t51 + 0.50000000
     >0000000000000D0 * t89 * t55)
      t152 = t144 * t151
      t153 = t140 * t152
      t156 = 0.1D1 / (0.1D1 - 0.160249952256378709147D1 * t153)
      t161 = tanh(0.160249952256378709147D31 * t140 * t152 * t156)
      t167 = exp(-0.102559969444082373854D1 * t140 * t144 * t151 * t156 
     >* (0.500000000000000000000D0 + 0.500000000000000000000D0 * t161))
      t170 = tanh(0.1D31 - 0.160249952256378709147D31 * t153)
      t171 = 0.500000000000000000000D0 * t170
      t175 = tanh(0.1D31 * t156)
      t180 = exp(0.15D1 * t156 * (0.500000000000000000000D0 - 0.50000000
     >0000000000000D0 * t175))
      t188 = 0.1D1 / (0.1D1 + 0.847380836474276614068D-1 * t22 + 0.11406
     >1568172369822458D0 * t19)
      t191 = exp(0.100000000000000000000D1 * t188)
      t196 = (0.1D1 + 0.615484811627254218514D-1 * t116 * t117) ** (0.1D
     >1 / 0.4D1)
      t202 = log(0.1D1 + (t191 - 0.1D1) * (0.1D1 - 0.1D1 / t196))
      t209 = t60 ** 2
      t211 = t62 ** 2
      zk(i) = rho * (-t35 + t69 + t86 + t134 + (t167 * (0.50000000000000
     >0000000D0 + t171) - 0.7D0 * t180 * (0.500000000000000000000D0 - t1
     >71)) * ((-0.285764D-1 * t188 + 0.285764D-1 * t202) * (0.33631D1 -
     >0.118155000000000000000D1 * t53 - 0.118155000000000000000D1 * t57)
     > * (0.1D1 - 0.1D1 * t209 * t60 / t211 / t62) + t35 - t69 - t86 - t
     >134))

               endif
             endif
           enddo
         endif
       else
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i))
      sigmaab = sigmaaa
      sigmabb = sigmaab
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i))
      taub = taua
      tau = taua + taub
      t8 = 0.1D1 / pi
      t9 = 0.1D1 / rho
      t10 = t8 * t9
      t11 = t10 ** (0.1D1 / 0.3D1)
      t13 = 0.1D1 + 0.194159335344114122552D0 * t11
      t14 = t10 ** (0.1D1 / 0.6D1)
      t17 = sqrt(t10)
      t19 = t11 ** 2
      t21 = 0.724010193431683113327D1 * t14 + 0.325955091942229212011D1 
     >* t11 + 0.141872281647966739112D1 * t17 + 0.406913004517529319387D
     >0 * t19
      t24 = 0.1D1 + 0.160819794986925350668D2 / t21
      t25 = log(t24)
      t27 = 0.621814D-1 * t13 * t25
      t29 = 0.1D1 + 0.101077332976287768525D0 * t11
      t34 = 0.987212972256927209438D1 * t14 + 0.329180480994506259905D1 
     >* t11 + 0.762327521935289963194D0 * t17 + 0.410025070949612505036D
     >0 * t19
      t37 = 0.1D1 + 0.296087499777934375166D2 / t34
      t38 = log(t37)
      t39 = t29 * t38
      t41 = rhoa - 0.1D1 * rhob
      t42 = t41 * t9
      t43 = 0.1D1 + t42
      t44 = t43 ** (0.1D1 / 0.3D1)
      t45 = t44 * t43
      t47 = 0.1D1 - 0.1D1 * t42
      t48 = t47 ** (0.1D1 / 0.3D1)
      t49 = t48 * t47
      t50 = t45 + t49 - 0.2D1
      t51 = t41 ** 2
      t52 = t51 ** 2
      t53 = rho ** 2
      t54 = t53 ** 2
      t55 = 0.1D1 / t54
      t56 = t52 * t55
      t58 = 0.1D1 - 0.1D1 * t56
      t61 = 0.379955235370239451738D-1 * t39 * t50 * t58
      t63 = 0.1D1 + 0.186690969707574028554D0 * t11
      t68 = 0.134579137143944477912D2 * t14 + 0.563098414909787598194D1 
     >* t11 + 0.291521471421917737271D1 * t17 + 0.516066464547863440989D
     >0 * t19
      t71 = 0.1D1 + 0.321639589973850701335D2 / t68
      t72 = log(t71)
      t75 = -0.3109070D-1 * t63 * t72 + t27
      t76 = t75 * t50
      t78 = 0.192366105093153631974D1 * t76 * t56
      t79 = t44 ** 2
      t81 = t48 ** 2
      t83 = 0.500000000000000000000D0 * t79 + 0.500000000000000000000D0 
     >* t81
      t84 = t83 ** 2
      t85 = t84 * t83
      t86 = -t27 + t61 + t78
      t87 = 0.1D1 / t85
      t90 = exp(-0.321636486443022096427D2 * t86 * t87)
      t91 = t90 - 0.1D1
      t93 = 0.1D1 + 0.908560296416069829442D-1 * t11
      t95 = 0.1D1 + 0.161542020702777215675D0 * t11
      t96 = 0.1D1 / t95
      t97 = t93 * t96
      t98 = 0.1D1 / t91
      t99 = pi ** 2
      t100 = t99 ** (0.1D1 / 0.3D1)
      t101 = t100 ** 2
      t102 = t98 * t101
      t103 = t97 * t102
      t104 = t99 * rho
      t105 = t104 ** (0.1D1 / 0.3D1)
      t106 = t105 ** 2
      t107 = 0.1D1 / t106
      t108 = sigma * t107
      t109 = 0.1D1 / t53
      t110 = 0.1D1 / t84
      t112 = 0.1D1 / t11
      t117 = 0.1D1 + 0.372009030193995856361D0 * t103 * t108 * t109 * t1
     >10 * t112
      t118 = t117 ** (0.1D1 / 0.4D1)
      t121 = 0.1D1 - 0.1D1 / t118
      t123 = 0.1D1 + t91 * t121
      t124 = log(t123)
      t126 = 0.31091D-1 * t85 * t124
      t130 = 0.500000000000000000000D0 * tau - 0.125000000000000000000D0
     > * sigma * t9
      t131 = 0.1D1 / t101
      t132 = t130 * t131
      t133 = rho ** (0.1D1 / 0.3D1)
      t134 = t133 ** 2
      t136 = 0.1D1 / t134 / rho
      t137 = t132 * t136
      t142 = 0.500000000000000000000D0 * t79 * t43 + 0.50000000000000000
     >0000D0 * t81 * t47
      t143 = 0.1D1 / t142
      t144 = t136 * t143
      t145 = t132 * t144
      t147 = 0.1D1 - 0.160249952256378709147D1 * t145
      t148 = 0.1D1 / t147
      t149 = t143 * t148
      t153 = tanh(0.160249952256378709147D31 * t132 * t144 * t148)
      t155 = 0.500000000000000000000D0 + 0.500000000000000000000D0 * t15
     >3
      t156 = t149 * t155
      t159 = exp(-0.102559969444082373854D1 * t137 * t156)
      t162 = tanh(0.1D31 - 0.160249952256378709147D31 * t145)
      t163 = 0.500000000000000000000D0 * t162
      t164 = 0.500000000000000000000D0 + t163
      t167 = tanh(0.1D31 * t148)
      t169 = 0.500000000000000000000D0 - 0.500000000000000000000D0 * t16
     >7
      t172 = exp(0.15D1 * t148 * t169)
      t173 = 0.500000000000000000000D0 - t163
      t176 = t159 * t164 - 0.7D0 * t172 * t173
      t179 = 0.1D1 + 0.847380836474276614068D-1 * t14 + 0.11406156817236
     >9822458D0 * t11
      t180 = 0.1D1 / t179
      t183 = exp(0.100000000000000000000D1 * t180)
      t184 = t183 - 0.1D1
      t187 = 0.1D1 + 0.615484811627254218514D-1 * t108 * t109
      t188 = t187 ** (0.1D1 / 0.4D1)
      t191 = 0.1D1 - 0.1D1 / t188
      t193 = 0.1D1 + t184 * t191
      t194 = log(t193)
      t196 = -0.285764D-1 * t180 + 0.285764D-1 * t194
      t199 = 0.33631D1 - 0.118155000000000000000D1 * t45 - 0.11815500000
     >0000000000D1 * t49
      t200 = t196 * t199
      t201 = t52 ** 2
      t202 = t201 * t52
      t203 = t54 ** 2
      t205 = 0.1D1 / t203 / t54
      t208 = 0.1D1 - 0.1D1 * t202 * t205
      t210 = t200 * t208 + t27 - t61 - t78 - t126
      t211 = t176 * t210
      zk(i) = rho * (-t27 + t61 + t78 + t126 + t211)
      t214 = t44 * t9
      t216 = t48 * t9
      t218 = 0.133333333333333333333D1 * t214 - 0.133333333333333333333D
     >1 * t216
      t221 = 0.379955235370239451738D-1 * t39 * t218 * t58
      t222 = t51 * t41
      t226 = 0.151982094148095780695D0 * t39 * t50 * t222 * t55
      t229 = 0.192366105093153631974D1 * t75 * t218 * t56
      t232 = 0.769464420372614527896D1 * t76 * t222 * t55
      t233 = t84 * t124
      t234 = 0.1D1 / t44
      t237 = 0.1D1 / t48
      t240 = 0.333333333333333333333D0 * t234 * t9 - 0.33333333333333333
     >3333D0 * t237 * t9
      t242 = 0.93273D-1 * t233 * t240
      t246 = t84 ** 2
      t248 = t86 / t246
      t251 = -0.321636486443022096427D2 * (t221 - t226 + t229 + t232) * 
     >t87 + 0.964909459329066289281D2 * t248 * t240
      t255 = 0.1D1 / t118 / t117
      t256 = t91 * t255
      t257 = t91 ** 2
      t261 = t97 / t257 * t101 * sigma
      t262 = t107 * t109
      t263 = t262 * t110
      t270 = t97 * t102 * sigma
      t271 = t87 * t112
      t281 = 0.1D1 / t123
      t283 = 0.31091D-1 * t85 * (t251 * t90 * t121 + 0.25000000000000000
     >0000D0 * t256 * (-0.372009030193995856361D0 * t261 * t263 * t112 *
     > t251 * t90 - 0.744018060387991712722D0 * t270 * t262 * t271 * t24
     >0)) * t281
      t284 = t142 ** 2
      t285 = 0.1D1 / t284
      t286 = t285 * t148
      t291 = 0.833333333333333333333D0 * t79 * t9 - 0.833333333333333333
     >333D0 * t81 * t9
      t292 = t155 * t291
      t296 = t130 ** 2
      t298 = 0.1D1 / t100 / t99
      t300 = t53 * rho
      t302 = 0.1D1 / t133 / t300
      t303 = t296 * t298 * t302
      t306 = t147 ** 2
      t307 = 0.1D1 / t306
      t308 = 0.1D1 / t284 / t142 * t307
      t312 = t153 ** 2
      t314 = 0.1D1 - 0.1D1 * t312
      t329 = t162 ** 2
      t331 = 0.1D1 - 0.1D1 * t329
      t332 = t159 * t331
      t333 = t332 * t130
      t334 = t131 * t136
      t336 = t334 * t285 * t291
      t339 = t307 * t169
      t340 = t339 * t130
      t345 = t167 ** 2
      t348 = 0.1D1 / t306 / t147 * (0.1D1 - 0.1D1 * t345)
      t349 = t348 * t130
      t356 = t172 * t331
      t357 = t356 * t130
      t364 = -0.157540000000000000000D1 * t214 + 0.157540000000000000000
     >D1 * t216
      t370 = 0.12D2 * t200 * t201 * t222 * t205
      t376 = -t218
      t379 = 0.379955235370239451738D-1 * t39 * t376 * t58
      t382 = 0.192366105093153631974D1 * t75 * t376 * t56
      t383 = -t240
      t385 = 0.93273D-1 * t233 * t383
      t391 = -0.321636486443022096427D2 * (t379 + t226 + t382 - t232) * 
     >t87 + 0.964909459329066289281D2 * t248 * t383
      t409 = 0.31091D-1 * t85 * (t391 * t90 * t121 + 0.25000000000000000
     >0000D0 * t256 * (-0.372009030193995856361D0 * t261 * t263 * t112 *
     > t391 * t90 - 0.744018060387991712722D0 * t270 * t262 * t271 * t38
     >3)) * t281
      t410 = -t291
      t411 = t155 * t410
      t433 = t334 * t285 * t410
      t457 = 0.1D1 / t19 * t8
      t460 = 0.402436643158883263334D-2 * t457 * t109 * t25
      t461 = t21 ** 2
      t464 = t14 ** 2
      t465 = t464 ** 2
      t469 = 0.1D1 / t465 / t14 * t8 * t109
      t471 = t457 * t109
      t475 = 0.1D1 / t17 * t8 * t109
      t478 = t112 * t8 * t109
      t484 = 0.100000000000000000000D1 * t13 / t461 * (-0.12066836557194
     >7185555D1 * t469 - 0.108651697314076404004D1 * t471 - 0.7093614082
     >39833695563D0 * t475 - 0.271275336345019546258D0 * t478) / t24
      t488 = 0.128016206138671616206D-2 * t471 * t38 * t50 * t58
      t489 = t34 ** 2
      t502 = 0.112499995668310776915D1 * t29 / t489 * (-0.16453549537615
     >4534908D1 * t469 - 0.109726826998168753302D1 * t471 - 0.3811637609
     >67644981600D0 * t475 - 0.273350047299741670024D0 * t478) / t37 * t
     >50 * t58
      t504 = t44 * t41 * t109
      t507 = t48 * t41 * t109
      t509 = -0.133333333333333333333D1 * t504 + 0.133333333333333333333
     >D1 * t507
      t512 = 0.379955235370239451738D-1 * t39 * t509 * t58
      t514 = t54 * rho
      t515 = 0.1D1 / t514
      t518 = 0.151982094148095780695D0 * t39 * t50 * t52 * t515
      t522 = t68 ** 2
      t537 = 0.192366105093153631974D1 * (0.193478431062909061652D-2 * t
     >457 * t109 * t72 + 0.100000000000000000000D1 * t63 / t522 * (-0.22
     >4298561906574129855D1 * t469 - 0.187699471636595866065D1 * t471 -
     >0.145760735710958868637D1 * t475 - 0.344044309698575627326D0 * t47
     >8) / t71 - t460 - t484) * t50 * t56
      t540 = 0.192366105093153631974D1 * t75 * t509 * t56
      t543 = 0.769464420372614527896D1 * t76 * t52 * t515
      t550 = -0.333333333333333333333D0 * t234 * t41 * t109 + 0.33333333
     >3333333333333D0 * t237 * t41 * t109
      t552 = 0.93273D-1 * t233 * t550
      t558 = -0.321636486443022096427D2 * (t460 + t484 - t488 - t502 + t
     >512 + t518 + t537 + t540 - t543) * t87 + 0.964909459329066289281D2
     > * t248 * t550
      t561 = 0.1D1 / t300
      t569 = t95 ** 2
      t573 = t561 * t110
      t583 = 0.1D1 / t106 / t104
      t612 = 0.31091D-1 * t85 * (t558 * t90 * t121 + 0.25000000000000000
     >0000D0 * t256 * (-0.112664211580837182141D-1 * t561 * t96 * t98 *
     >t101 * sigma * t107 * t110 + 0.200316968190728509847D-1 * t93 / t5
     >69 * t102 * t108 * t573 - 0.372009030193995856361D0 * t261 * t263
     >* t112 * t558 * t90 - 0.248006020129330570907D0 * t270 * t583 * t1
     >09 * t110 * t112 * t99 - 0.744018060387991712722D0 * t103 * t108 *
     > t573 * t112 - 0.744018060387991712722D0 * t270 * t262 * t271 * t5
     >50 + 0.124003010064665285454D0 * t270 * t107 * t55 * t110 / t11 /
     >t10 * t8)) * t281
      t615 = sigma / t134 / t300
      t620 = 0.1D1 / t134 / t53
      t630 = -0.833333333333333333333D0 * t79 * t41 * t109 + 0.833333333
     >333333333333D0 * t81 * t41 * t109
      t635 = t143 * t307
      t636 = t131 * t143
      t637 = t615 * t636
      t639 = t620 * t143
      t640 = t132 * t639
      t644 = t132 * t136 * t285 * t630
      t646 = -0.200312440320473386433D0 * t637 + 0.267083253760631181911
     >D1 * t640 + 0.160249952256378709147D1 * t644
      t674 = -0.200312440320473386433D30 * t637 + 0.26708325376063118191
     >1D31 * t640 + 0.160249952256378709147D31 * t644
      t689 = t179 ** 2
      t694 = 0.1D1 / t689 * (-0.141230139412379435678D-1 * t469 - 0.3802
     >05227241232741527D-1 * t471)
      t701 = t184 / t188 / t187
      t712 = 0.1D1 / t193
      t728 = (0.285764D-1 * t694 + 0.285764D-1 * (-0.1000000000000000000
     >00D1 * t694 * t183 * t191 + 0.250000000000000000000D0 * t701 * (-0
     >.410323207751502812342D-1 * sigma * t583 * t109 * t99 - 0.12309696
     >2325450843703D0 * t108 * t561)) * t712) * t199 * t208 + t196 * (0.
     >157540000000000000000D1 * t504 - 0.157540000000000000000D1 * t507)
     > * t208 + 0.12D2 * t200 * t202 / t203 / t514 - t460 - t484 + t488
     >+ t502 - t512 - t518 - t537 - t540 + t543 - t552 - t612
      t730 = t460 + t484 - t488 - t502 + t512 + t518 + t537 + t540 - t54
     >3 + t552 + t612 + ((-0.128199961805102967317D0 * t615 * t131 * t15
     >6 + 0.170933282406803956424D1 * t132 * t620 * t156 + 0.10255996944
     >4082373854D1 * t137 * t286 * t155 * t630 + 0.102559969444082373854
     >D1 * t137 * t635 * t155 * t646 - 0.512799847220411869270D0 * t137
     >* t149 * t314 * (0.200312440320473386433D30 * t615 * t636 * t148 -
     > 0.267083253760631181911D31 * t132 * t639 * t148 - 0.1602499522563
     >78709147D31 * t137 * t286 * t630 - 0.160249952256378709147D31 * t1
     >37 * t635 * t646)) * t159 * t164 + 0.500000000000000000000D0 * t33
     >2 * t674 - 0.7D0 * (-0.15D1 * t339 * t646 + 0.75000000000000000000
     >0D30 * t348 * t646) * t172 * t173 + 0.350000000000000000000D0 * t3
     >56 * t674) * t210 + t176 * t728
      vrhoc(i) = vrhoc(i) + 0.5D0 * rho * (t221 - t226 + t229 + t232 + t
     >242 + t283 + ((0.102559969444082373854D1 * t137 * t286 * t292 + 0.
     >164352302068298596704D1 * t303 * t308 * t292 - 0.51279984722041186
     >9270D0 * t137 * t149 * t314 * (-0.160249952256378709147D31 * t137
     >* t286 * t291 - 0.256800471981716557348D31 * t303 * t308 * t291))
     >* t159 * t164 + 0.801249761281893545733D30 * t333 * t336 - 0.7D0 *
     > (-0.240374928384568063720D1 * t340 * t336 + 0.1201874641922840318
     >60D31 * t349 * t336) * t172 * t173 + 0.560874832897325482012D30 *
     >t357 * t336) * t210 + t176 * (t196 * t364 * t208 - t370 - t221 + t
     >226 - t229 - t232 - t242 - t283)) + 0.5D0 * rho * (t379 + t226 + t
     >382 - t232 + t385 + t409 + ((0.102559969444082373854D1 * t137 * t2
     >86 * t411 + 0.164352302068298596704D1 * t303 * t308 * t411 - 0.512
     >799847220411869270D0 * t137 * t149 * t314 * (-0.160249952256378709
     >147D31 * t137 * t286 * t410 - 0.256800471981716557348D31 * t303 *
     >t308 * t410)) * t159 * t164 + 0.801249761281893545733D30 * t333 *
     >t433 - 0.7D0 * (-0.240374928384568063720D1 * t340 * t433 + 0.12018
     >7464192284031860D31 * t349 * t433) * t172 * t173 + 0.5608748328973
     >25482012D30 * t357 * t433) * t210 + t176 * (-t196 * t364 * t208 +
     >t370 - t379 - t226 - t382 + t232 - t385 - t409)) - t27 + t61 + t78
     > + t126 + t211 + rho * t730
      t740 = 0.289153318944038129253D-2 * t83 * t255 * t97 * t101 * t107
     > * t109 * t112 * t281
      t741 = t620 * t131
      t744 = t130 * t298
      t746 = 0.1D1 / t133 / t54
      t749 = t285 * t307 * t155
      t766 = t741 * t143
      vsigmacc(i) = vsigmacc(i) + rho * (t740 + ((0.12819996180510296731
     >7D0 * t741 * t156 + 0.205440377585373245880D0 * t744 * t746 * t749
     > - 0.512799847220411869270D0 * t137 * t149 * t314 * (-0.2003124403
     >20473386433D30 * t741 * t149 - 0.321000589977145696686D30 * t744 *
     > t746 * t285 * t307)) * t159 * t164 + 0.100156220160236693217D30 *
     > t332 * t766 - 0.7D0 * (-0.300468660480710079650D0 * t339 * t766 +
     > 0.150234330240355039825D30 * t348 * t766) * t172 * t173 + 0.70109
     >3541121656852518D29 * t356 * t766) * t210 + t176 * (0.439708504274
     >626686249D-3 * t701 * t107 * t109 * t712 * t199 * t208 - t740))
      t811 = t334 * t143
      vtauc(i) = vtauc(i) + rho * ((-0.512799847220411869270D0 * t334 * 
     >t156 - 0.821761510341492983518D0 * t744 * t302 * t749 - 0.51279984
     >7220411869270D0 * t137 * t149 * t314 * (0.801249761281893545733D30
     > * t334 * t149 + 0.128400235990858278675D31 * t744 * t302 * t285 *
     > t307)) * t159 * t164 - 0.400624880640946772867D30 * t332 * t811 -
     > 0.7D0 * (0.120187464192284031860D1 * t339 * t811 - 0.600937320961
     >420159300D30 * t348 * t811) * t172 * t173 - 0.28043741644866274100
     >7D30 * t356 * t811) * t210

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i))
      sigmaab = sigmaaa
      sigmabb = sigmaab
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i))
      taub = taua
      tau = taua + taub
      t9 = 0.1D1 / rho
      t10 = 0.1D1 / pi * t9
      t11 = t10 ** (0.1D1 / 0.3D1)
      t14 = t10 ** (0.1D1 / 0.6D1)
      t17 = sqrt(t10)
      t19 = t11 ** 2
      t25 = log(0.1D1 + 0.160819794986925350668D2 / (0.72401019343168311
     >3327D1 * t14 + 0.325955091942229212011D1 * t11 + 0.141872281647966
     >739112D1 * t17 + 0.406913004517529319387D0 * t19))
      t27 = 0.621814D-1 * (0.1D1 + 0.194159335344114122552D0 * t11) * t2
     >5
      t38 = log(0.1D1 + 0.296087499777934375166D2 / (0.98721297225692720
     >9438D1 * t14 + 0.329180480994506259905D1 * t11 + 0.762327521935289
     >963194D0 * t17 + 0.410025070949612505036D0 * t19))
      t41 = rhoa - 0.1D1 * rhob
      t42 = t41 * t9
      t43 = 0.1D1 + t42
      t44 = t43 ** (0.1D1 / 0.3D1)
      t45 = t44 * t43
      t47 = 0.1D1 - 0.1D1 * t42
      t48 = t47 ** (0.1D1 / 0.3D1)
      t49 = t48 * t47
      t50 = t45 + t49 - 0.2D1
      t51 = t41 ** 2
      t52 = t51 ** 2
      t53 = rho ** 2
      t54 = t53 ** 2
      t56 = t52 / t54
      t61 = 0.379955235370239451738D-1 * (0.1D1 + 0.10107733297628776852
     >5D0 * t11) * t38 * t50 * (0.1D1 - 0.1D1 * t56)
      t72 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t14 + 0.563098414909787598194D1 * t11 + 0.291521471421917
     >737271D1 * t17 + 0.516066464547863440989D0 * t19))
      t78 = 0.192366105093153631974D1 * (-0.3109070D-1 * (0.1D1 + 0.1866
     >90969707574028554D0 * t11) * t72 + t27) * t50 * t56
      t79 = t44 ** 2
      t81 = t48 ** 2
      t83 = 0.500000000000000000000D0 * t79 + 0.500000000000000000000D0 
     >* t81
      t84 = t83 ** 2
      t85 = t84 * t83
      t90 = exp(-0.321636486443022096427D2 * (-t27 + t61 + t78) / t85)
      t91 = t90 - 0.1D1
      t99 = pi ** 2
      t100 = t99 ** (0.1D1 / 0.3D1)
      t101 = t100 ** 2
      t105 = (t99 * rho) ** (0.1D1 / 0.3D1)
      t106 = t105 ** 2
      t108 = sigma / t106
      t109 = 0.1D1 / t53
      t118 = (0.1D1 + 0.372009030193995856361D0 * (0.1D1 + 0.90856029641
     >6069829442D-1 * t11) / (0.1D1 + 0.161542020702777215675D0 * t11) /
     > t91 * t101 * t108 * t109 / t84 / t11) ** (0.1D1 / 0.4D1)
      t124 = log(0.1D1 + t91 * (0.1D1 - 0.1D1 / t118))
      t126 = 0.31091D-1 * t85 * t124
      t132 = (0.500000000000000000000D0 * tau - 0.125000000000000000000D
     >0 * sigma * t9) / t101
      t133 = rho ** (0.1D1 / 0.3D1)
      t134 = t133 ** 2
      t136 = 0.1D1 / t134 / rho
      t143 = 0.1D1 / (0.500000000000000000000D0 * t79 * t43 + 0.50000000
     >0000000000000D0 * t81 * t47)
      t144 = t136 * t143
      t145 = t132 * t144
      t148 = 0.1D1 / (0.1D1 - 0.160249952256378709147D1 * t145)
      t153 = tanh(0.160249952256378709147D31 * t132 * t144 * t148)
      t159 = exp(-0.102559969444082373854D1 * t132 * t136 * t143 * t148 
     >* (0.500000000000000000000D0 + 0.500000000000000000000D0 * t153))
      t162 = tanh(0.1D31 - 0.160249952256378709147D31 * t145)
      t163 = 0.500000000000000000000D0 * t162
      t167 = tanh(0.1D31 * t148)
      t172 = exp(0.15D1 * t148 * (0.500000000000000000000D0 - 0.50000000
     >0000000000000D0 * t167))
      t180 = 0.1D1 / (0.1D1 + 0.847380836474276614068D-1 * t14 + 0.11406
     >1568172369822458D0 * t11)
      t183 = exp(0.100000000000000000000D1 * t180)
      t188 = (0.1D1 + 0.615484811627254218514D-1 * t108 * t109) ** (0.1D
     >1 / 0.4D1)
      t194 = log(0.1D1 + (t183 - 0.1D1) * (0.1D1 - 0.1D1 / t188))
      t201 = t52 ** 2
      t203 = t54 ** 2
      zk(i) = rho * (-t27 + t61 + t78 + t126 + (t159 * (0.50000000000000
     >0000000D0 + t163) - 0.7D0 * t172 * (0.500000000000000000000D0 - t1
     >63)) * ((-0.285764D-1 * t180 + 0.285764D-1 * t194) * (0.33631D1 -
     >0.118155000000000000000D1 * t45 - 0.118155000000000000000D1 * t49)
     > * (0.1D1 - 0.1D1 * t201 * t52 / t203 / t54) + t27 - t61 - t78 - t
     >126))

             endif
           enddo
         endif
       endif

      return
      end

c:SCANCsubrend
c:SCANXsubrstart

c    Generated: Mon Jun 20 16:05:21 CEST 2016

      subroutine dftacg_scanx
     > (name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   tauc,tauo,upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,
     >                   vtauc,vtauo,vupsilonc,vupsilono)
      implicit double precision (a-h,o-z)
      logical fderiv,open
      integer igrad,npt
      character*(*) name
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)
      double precision tauc(*),tauo(*)
      double precision upsilonc(*),upsilono(*)
      double precision zk(*),vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)
      double precision vtauc(*),vtauo(*)
      double precision vupsilonc(*),vupsilono(*)
      include "common/cdft"
      include "common/tapes"
      parameter(tol=1d-12)
      pi=acos(-1d0)
      name='Automatically generated SCANX'
      igrad=2
       if(open) then
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      sigmabb = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) - 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmabb
      taub = max(0.0D0, 0.500000000000000000000D0 * tauc(i) - 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taub
      t13 = pi ** 2
      t14 = t13 * rhob
      t15 = t14 ** (0.1D1 / 0.3D1)
      t16 = rhob * t15
      t17 = 0.1D1 / pi
      t18 = t15 ** 2
      t19 = 0.1D1 / t18
      t20 = sigmabb * t19
      t21 = rhob ** 2
      t22 = 0.1D1 / t21
      t23 = t20 * t22
      t25 = exp(-0.747166092922644636268D-1 * t23)
      t26 = t22 * t25
      t29 = 0.1D1 + 0.747166092922644636268D-1 * t20 * t26
      t30 = t22 * t29
      t34 = 0.1D1 / rhob
      t37 = taub - 0.250000000000000000000D0 * sigmabb * t34
      t38 = t13 ** (0.1D1 / 0.3D1)
      t39 = t38 ** 2
      t40 = 0.1D1 / t39
      t41 = t37 * t40
      t42 = rhob ** (0.1D1 / 0.3D1)
      t43 = t42 ** 2
      t45 = 0.1D1 / t43 / rhob
      t46 = t41 * t45
      t48 = 0.1D1 - 0.504755720231149905248D0 * t46
      t49 = t48 ** 2
      t51 = exp(-0.5D0 * t49)
      t54 = 0.118591405585874358362D-1 * t23 + 0.120830459735945720683D0
     > * t48 * t51
      t55 = t54 ** 2
      t57 = 0.1D1 + 0.143805048498903106908D0 * t20 * t30 + 0.1538461538
     >46153846154D2 * t55
      t59 = 0.65D-1 / t57
      t60 = 0.1D1 / t48
      t61 = t45 * t60
      t64 = tanh(0.504755720231149905248D30 * t41 * t61)
      t66 = 0.500000000000000000000D0 + 0.500000000000000000000D0 * t64
      t70 = exp(-0.336672065394176986800D0 * t41 * t61 * t66)
      t73 = tanh(0.1D31 - 0.504755720231149905248D30 * t46)
      t74 = 0.500000000000000000000D0 * t73
      t75 = 0.500000000000000000000D0 + t74
      t78 = tanh(0.1D31 * t60)
      t80 = 0.500000000000000000000D0 - 0.500000000000000000000D0 * t78
      t83 = exp(0.8D0 * t60 * t80)
      t84 = 0.500000000000000000000D0 - t74
      t87 = t70 * t75 - 0.124D1 * t83 * t84
      t88 = 0.109D0 + t59
      t90 = 0.1065D1 - t59 + t87 * t88
      t91 = t17 * t90
      t92 = sqrt(sigmabb)
      t94 = t92 / t15
      t95 = t94 * t34
      t96 = sqrt(t95)
      t99 = exp(-0.943252112663908494661D1 / t96)
      t101 = 0.1D1 - 0.1D1 * t99
      zk(i) = -0.136284044462410474417D1 * t16 * t91 * t101
      t109 = 0.681420222312052372084D0 * t15 * t17 * t90 * t101
      t114 = 0.227140074104017457361D0 * rhob * t19 * pi * t90 * t101
      t115 = t57 ** 2
      t116 = 0.1D1 / t115
      t119 = sigmabb / t18 / t14
      t123 = t21 * rhob
      t124 = 0.1D1 / t123
      t135 = t119 * t22 * t13
      t137 = t20 * t124
      t152 = sigmabb / t43 / t123
      t153 = t152 * t40
      t156 = 0.1D1 / t43 / t21
      t157 = t41 * t156
      t159 = -0.126188930057787476312D0 * t153 + 0.841259533718583175412
     >D0 * t157
      t168 = -0.958700323326020712721D-1 * t119 * t30 * t13 - 0.28761009
     >6997806213817D0 * t20 * t124 * t29 + 0.143805048498903106908D0 * t
     >20 * t22 * (-0.498110728615096424178D-1 * t119 * t26 * t13 - 0.149
     >433218584528927254D0 * t20 * t124 * t25 + 0.747166092922644636268D
     >-1 * t20 * t22 * (0.498110728615096424178D-1 * t135 + 0.1494332185
     >84528927254D0 * t137) * t25) + 0.307692307692307692308D2 * t54 * (
     >-0.790609370572495722413D-2 * t135 - 0.237182811171748716724D-1 *
     >t137 + 0.120830459735945720683D0 * t159 * t51 - 0.1208304597359457
     >20683D0 * t49 * t159 * t51)
      t171 = t40 * t60
      t175 = t156 * t60
      t179 = 0.1D1 / t49
      t184 = t64 ** 2
      t187 = t60 * (0.1D1 - 0.1D1 * t184)
      t203 = t73 ** 2
      t205 = 0.1D1 - 0.1D1 * t203
      t206 = t70 * t205
      t209 = -0.126188930057787476312D30 * t153 + 0.84125953371858317541
     >2D30 * t157
      t212 = t179 * t80
      t217 = t78 ** 2
      t220 = 0.1D1 / t49 / t48 * (0.1D1 - 0.1D1 * t217)
      t227 = t83 * t205
      t232 = t87 * t116
      t239 = 0.681420222312052372084D0 * t16 * t17 * (0.65D-1 * t116 * t
     >168 + ((-0.841680163485442467001D-1 * t152 * t171 * t66 + 0.561120
     >108990294978001D0 * t41 * t175 * t66 + 0.336672065394176986800D0 *
     > t46 * t179 * t66 * t159 - 0.168336032697088493400D0 * t46 * t187
     >* (0.126188930057787476312D30 * t152 * t171 - 0.841259533718583175
     >412D30 * t41 * t175 - 0.504755720231149905248D30 * t41 * t45 * t17
     >9 * t159)) * t70 * t75 + 0.500000000000000000000D0 * t206 * t209 -
     > 0.124D1 * (-0.8D0 * t212 * t159 + 0.400000000000000000000D30 * t2
     >20 * t159) * t83 * t84 + 0.620000000000000000000D0 * t227 * t209)
     >* t88 - 0.65D-1 * t232 * t168) * t101
      t242 = 0.1D1 / t96 / t95
      t245 = 0.1D1 / t15 / t14
      t256 = 0.973296829181994948224D0 * t16 * t17 * t90 * t242 * (-0.11
     >0064241629820889462D1 * t92 * t245 * t34 * t13 - 0.330192724889462
     >668387D1 * t94 * t22) * t99
      vrhoc(i) = vrhoc(i) - t109 - t114 - t239 + t256
      vrhoo(i) = vrhoo(i) + t109 + t114 + t239 - t256
      t259 = t19 * t22
      t265 = t21 ** 2
      t275 = t156 * t40
      t285 = 0.143805048498903106908D0 * t259 * t29 + 0.1438050484989031
     >06908D0 * t20 * t22 * (0.747166092922644636268D-1 * t259 * t25 - 0
     >.558257170413290039222D-2 * sigmabb * t245 / t265 * t25) + 0.30769
     >2307692307692308D2 * t54 * (0.118591405585874358362D-1 * t259 + 0.
     >152474664324695603608D-1 * t275 * t51 - 0.152474664324695603608D-1
     > * t49 * t156 * t40 * t51)
      t288 = t60 * t66
      t293 = t37 / t38 / t13
      t296 = 0.1D1 / t42 / t265 * t179
      t330 = t16 * t17 * (0.65D-1 * t116 * t285 + ((0.841680163485442467
     >001D-1 * t275 * t288 + 0.424842877124366511988D-1 * t293 * t296 *
     >t66 - 0.168336032697088493400D0 * t46 * t187 * (-0.126188930057787
     >476312D30 * t275 * t60 - 0.636945842765167184389D29 * t293 * t296)
     >) * t70 * t75 + 0.630944650288937381559D29 * t206 * t275 - 0.124D1
     > * (-0.100951144046229981050D0 * t212 * t275 + 0.50475572023114990
     >5248D29 * t220 * t275) * t83 * t84 + 0.782371366358282353132D29 *
     >t227 * t275) * t88 - 0.65D-1 * t232 * t285) * t101
      t331 = 0.340710111156026186042D0 * t330
      t335 = t91 * t242 / t92 * t99
      t336 = 0.803438830384691996314D0 * t335
      vsigmacc(i) = vsigmacc(i) - t331 + t336
      vsigmaco(i) = vsigmaco(i) + 0.681420222312052372084D0 * t330 - 0.1
     >60687766076938399263D1 * t335
      vsigmaoo(i) = vsigmaoo(i) - t331 + t336
      t343 = t40 * t45
      t350 = -0.609898657298782414434D-1 * t343 * t51 + 0.60989865729878
     >2414434D-1 * t49 * t40 * t45 * t51
      t357 = 0.1D1 / t42 / t123 * t179
      t393 = 0.681420222312052372084D0 * t16 * t17 * (0.2000000000000000
     >00000D1 * t116 * t54 * t350 + ((-0.336672065394176986800D0 * t343
     >* t288 - 0.169937150849746604795D0 * t293 * t357 * t66 - 0.1683360
     >32697088493400D0 * t46 * t187 * (0.504755720231149905248D30 * t343
     > * t60 + 0.254778337106066873756D30 * t293 * t357)) * t70 * t75 -
     >0.252377860115574952624D30 * t206 * t343 - 0.124D1 * (0.4038045761
     >84919924197D0 * t212 * t343 - 0.201902288092459962099D30 * t220 *
     >t343) * t83 * t84 - 0.312948546543312941253D30 * t227 * t343) * t8
     >8 - 0.200000000000000000000D1 * t232 * t54 * t350) * t101
      vtauc(i) = vtauc(i) - t393
      vtauo(i) = vtauo(i) + t393

               elseif(rhob.lt.tol) then
             rho = rhoa
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) + 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmaaa
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i) + 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taua
      t13 = pi ** 2
      t14 = t13 * rhoa
      t15 = t14 ** (0.1D1 / 0.3D1)
      t16 = rhoa * t15
      t17 = 0.1D1 / pi
      t18 = t15 ** 2
      t19 = 0.1D1 / t18
      t20 = sigmaaa * t19
      t21 = rhoa ** 2
      t22 = 0.1D1 / t21
      t23 = t20 * t22
      t25 = exp(-0.747166092922644636268D-1 * t23)
      t26 = t22 * t25
      t29 = 0.1D1 + 0.747166092922644636268D-1 * t20 * t26
      t30 = t22 * t29
      t34 = 0.1D1 / rhoa
      t37 = taua - 0.250000000000000000000D0 * sigmaaa * t34
      t38 = t13 ** (0.1D1 / 0.3D1)
      t39 = t38 ** 2
      t40 = 0.1D1 / t39
      t41 = t37 * t40
      t42 = rhoa ** (0.1D1 / 0.3D1)
      t43 = t42 ** 2
      t45 = 0.1D1 / t43 / rhoa
      t46 = t41 * t45
      t48 = 0.1D1 - 0.504755720231149905248D0 * t46
      t49 = t48 ** 2
      t51 = exp(-0.5D0 * t49)
      t54 = 0.118591405585874358362D-1 * t23 + 0.120830459735945720683D0
     > * t48 * t51
      t55 = t54 ** 2
      t57 = 0.1D1 + 0.143805048498903106908D0 * t20 * t30 + 0.1538461538
     >46153846154D2 * t55
      t59 = 0.65D-1 / t57
      t60 = 0.1D1 / t48
      t61 = t45 * t60
      t64 = tanh(0.504755720231149905248D30 * t41 * t61)
      t66 = 0.500000000000000000000D0 + 0.500000000000000000000D0 * t64
      t70 = exp(-0.336672065394176986800D0 * t41 * t61 * t66)
      t73 = tanh(0.1D31 - 0.504755720231149905248D30 * t46)
      t74 = 0.500000000000000000000D0 * t73
      t75 = 0.500000000000000000000D0 + t74
      t78 = tanh(0.1D31 * t60)
      t80 = 0.500000000000000000000D0 - 0.500000000000000000000D0 * t78
      t83 = exp(0.8D0 * t60 * t80)
      t84 = 0.500000000000000000000D0 - t74
      t87 = t70 * t75 - 0.124D1 * t83 * t84
      t88 = 0.109D0 + t59
      t90 = 0.1065D1 - t59 + t87 * t88
      t91 = t17 * t90
      t92 = sqrt(sigmaaa)
      t94 = t92 / t15
      t95 = t94 * t34
      t96 = sqrt(t95)
      t99 = exp(-0.943252112663908494661D1 / t96)
      t101 = 0.1D1 - 0.1D1 * t99
      zk(i) = -0.136284044462410474417D1 * t16 * t91 * t101
      t109 = 0.681420222312052372084D0 * t15 * t17 * t90 * t101
      t114 = 0.227140074104017457361D0 * rhoa * t19 * pi * t90 * t101
      t115 = t57 ** 2
      t116 = 0.1D1 / t115
      t119 = sigmaaa / t18 / t14
      t123 = t21 * rhoa
      t124 = 0.1D1 / t123
      t135 = t119 * t22 * t13
      t137 = t20 * t124
      t152 = sigmaaa / t43 / t123
      t153 = t152 * t40
      t156 = 0.1D1 / t43 / t21
      t157 = t41 * t156
      t159 = -0.126188930057787476312D0 * t153 + 0.841259533718583175412
     >D0 * t157
      t168 = -0.958700323326020712721D-1 * t119 * t30 * t13 - 0.28761009
     >6997806213817D0 * t20 * t124 * t29 + 0.143805048498903106908D0 * t
     >20 * t22 * (-0.498110728615096424178D-1 * t119 * t26 * t13 - 0.149
     >433218584528927254D0 * t20 * t124 * t25 + 0.747166092922644636268D
     >-1 * t20 * t22 * (0.498110728615096424178D-1 * t135 + 0.1494332185
     >84528927254D0 * t137) * t25) + 0.307692307692307692308D2 * t54 * (
     >-0.790609370572495722413D-2 * t135 - 0.237182811171748716724D-1 *
     >t137 + 0.120830459735945720683D0 * t159 * t51 - 0.1208304597359457
     >20683D0 * t49 * t159 * t51)
      t171 = t40 * t60
      t175 = t156 * t60
      t179 = 0.1D1 / t49
      t184 = t64 ** 2
      t187 = t60 * (0.1D1 - 0.1D1 * t184)
      t203 = t73 ** 2
      t205 = 0.1D1 - 0.1D1 * t203
      t206 = t70 * t205
      t209 = -0.126188930057787476312D30 * t153 + 0.84125953371858317541
     >2D30 * t157
      t212 = t179 * t80
      t217 = t78 ** 2
      t220 = 0.1D1 / t49 / t48 * (0.1D1 - 0.1D1 * t217)
      t227 = t83 * t205
      t232 = t87 * t116
      t239 = 0.681420222312052372084D0 * t16 * t17 * (0.65D-1 * t116 * t
     >168 + ((-0.841680163485442467001D-1 * t152 * t171 * t66 + 0.561120
     >108990294978001D0 * t41 * t175 * t66 + 0.336672065394176986800D0 *
     > t46 * t179 * t66 * t159 - 0.168336032697088493400D0 * t46 * t187
     >* (0.126188930057787476312D30 * t152 * t171 - 0.841259533718583175
     >412D30 * t41 * t175 - 0.504755720231149905248D30 * t41 * t45 * t17
     >9 * t159)) * t70 * t75 + 0.500000000000000000000D0 * t206 * t209 -
     > 0.124D1 * (-0.8D0 * t212 * t159 + 0.400000000000000000000D30 * t2
     >20 * t159) * t83 * t84 + 0.620000000000000000000D0 * t227 * t209)
     >* t88 - 0.65D-1 * t232 * t168) * t101
      t242 = 0.1D1 / t96 / t95
      t245 = 0.1D1 / t15 / t14
      t256 = 0.973296829181994948224D0 * t16 * t17 * t90 * t242 * (-0.11
     >0064241629820889462D1 * t92 * t245 * t34 * t13 - 0.330192724889462
     >668387D1 * t94 * t22) * t99
      vrhoc(i) = vrhoc(i) - t109 - t114 - t239 + t256
      vrhoo(i) = vrhoo(i) - t109 - t114 - t239 + t256
      t259 = t19 * t22
      t265 = t21 ** 2
      t275 = t156 * t40
      t285 = 0.143805048498903106908D0 * t259 * t29 + 0.1438050484989031
     >06908D0 * t20 * t22 * (0.747166092922644636268D-1 * t259 * t25 - 0
     >.558257170413290039222D-2 * sigmaaa * t245 / t265 * t25) + 0.30769
     >2307692307692308D2 * t54 * (0.118591405585874358362D-1 * t259 + 0.
     >152474664324695603608D-1 * t275 * t51 - 0.152474664324695603608D-1
     > * t49 * t156 * t40 * t51)
      t288 = t60 * t66
      t293 = t37 / t38 / t13
      t296 = 0.1D1 / t42 / t265 * t179
      t330 = t16 * t17 * (0.65D-1 * t116 * t285 + ((0.841680163485442467
     >001D-1 * t275 * t288 + 0.424842877124366511988D-1 * t293 * t296 *
     >t66 - 0.168336032697088493400D0 * t46 * t187 * (-0.126188930057787
     >476312D30 * t275 * t60 - 0.636945842765167184389D29 * t293 * t296)
     >) * t70 * t75 + 0.630944650288937381559D29 * t206 * t275 - 0.124D1
     > * (-0.100951144046229981050D0 * t212 * t275 + 0.50475572023114990
     >5248D29 * t220 * t275) * t83 * t84 + 0.782371366358282353132D29 *
     >t227 * t275) * t88 - 0.65D-1 * t232 * t285) * t101
      t331 = 0.340710111156026186042D0 * t330
      t335 = t91 * t242 / t92 * t99
      t336 = 0.803438830384691996314D0 * t335
      vsigmacc(i) = vsigmacc(i) - t331 + t336
      vsigmaco(i) = vsigmaco(i) - 0.681420222312052372084D0 * t330 + 0.1
     >60687766076938399263D1 * t335
      vsigmaoo(i) = vsigmaoo(i) - t331 + t336
      t343 = t40 * t45
      t350 = -0.609898657298782414434D-1 * t343 * t51 + 0.60989865729878
     >2414434D-1 * t49 * t40 * t45 * t51
      t357 = 0.1D1 / t42 / t123 * t179
      t393 = 0.681420222312052372084D0 * t16 * t17 * (0.2000000000000000
     >00000D1 * t116 * t54 * t350 + ((-0.336672065394176986800D0 * t343
     >* t288 - 0.169937150849746604795D0 * t293 * t357 * t66 - 0.1683360
     >32697088493400D0 * t46 * t187 * (0.504755720231149905248D30 * t343
     > * t60 + 0.254778337106066873756D30 * t293 * t357)) * t70 * t75 -
     >0.252377860115574952624D30 * t206 * t343 - 0.124D1 * (0.4038045761
     >84919924197D0 * t212 * t343 - 0.201902288092459962099D30 * t220 *
     >t343) * t83 * t84 - 0.312948546543312941253D30 * t227 * t343) * t8
     >8 - 0.200000000000000000000D1 * t232 * t54 * t350) * t101
      vtauc(i) = vtauc(i) - t393
      vtauo(i) = vtauo(i) - t393

               else
             rho = rhoa + rhob
      t2 = 0.250000000000000000000D0 * sigmacc(i)
      t4 = 0.250000000000000000000D0 * sigmaoo(i)
      t6 = 0.500000000000000000000D0 * sigmaco(i)
      sigmaaa = max(0.0D0, t2 + t4 + t6)
      sigmaab = t2 - t4
      sigmabb = max(0.0D0, t2 + t4 - t6)
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t11 = 0.500000000000000000000D0 * tauc(i)
      t13 = 0.500000000000000000000D0 * tauo(i)
      taua = max(0.0D0, t11 + t13)
      taub = max(0.0D0, t11 - t13)
      tau = taua + taub
      t16 = pi ** 2
      t17 = t16 * rhoa
      t18 = t17 ** (0.1D1 / 0.3D1)
      t19 = rhoa * t18
      t20 = 0.1D1 / pi
      t21 = t18 ** 2
      t22 = 0.1D1 / t21
      t23 = sigmaaa * t22
      t24 = rhoa ** 2
      t25 = 0.1D1 / t24
      t26 = t23 * t25
      t28 = exp(-0.747166092922644636268D-1 * t26)
      t29 = t25 * t28
      t32 = 0.1D1 + 0.747166092922644636268D-1 * t23 * t29
      t33 = t25 * t32
      t37 = 0.1D1 / rhoa
      t40 = taua - 0.250000000000000000000D0 * sigmaaa * t37
      t41 = t16 ** (0.1D1 / 0.3D1)
      t42 = t41 ** 2
      t43 = 0.1D1 / t42
      t44 = t40 * t43
      t45 = rhoa ** (0.1D1 / 0.3D1)
      t46 = t45 ** 2
      t48 = 0.1D1 / t46 / rhoa
      t49 = t44 * t48
      t51 = 0.1D1 - 0.504755720231149905248D0 * t49
      t52 = t51 ** 2
      t54 = exp(-0.5D0 * t52)
      t57 = 0.118591405585874358362D-1 * t26 + 0.120830459735945720683D0
     > * t51 * t54
      t58 = t57 ** 2
      t60 = 0.1D1 + 0.143805048498903106908D0 * t23 * t33 + 0.1538461538
     >46153846154D2 * t58
      t62 = 0.65D-1 / t60
      t63 = 0.1D1 / t51
      t64 = t48 * t63
      t67 = tanh(0.504755720231149905248D30 * t44 * t64)
      t69 = 0.500000000000000000000D0 + 0.500000000000000000000D0 * t67
      t73 = exp(-0.336672065394176986800D0 * t44 * t64 * t69)
      t76 = tanh(0.1D31 - 0.504755720231149905248D30 * t49)
      t77 = 0.500000000000000000000D0 * t76
      t78 = 0.500000000000000000000D0 + t77
      t81 = tanh(0.1D31 * t63)
      t83 = 0.500000000000000000000D0 - 0.500000000000000000000D0 * t81
      t86 = exp(0.8D0 * t63 * t83)
      t87 = 0.500000000000000000000D0 - t77
      t90 = t73 * t78 - 0.124D1 * t86 * t87
      t91 = 0.109D0 + t62
      t93 = 0.1065D1 - t62 + t90 * t91
      t94 = t20 * t93
      t95 = sqrt(sigmaaa)
      t97 = t95 / t18
      t98 = t97 * t37
      t99 = sqrt(t98)
      t102 = exp(-0.943252112663908494661D1 / t99)
      t104 = 0.1D1 - 0.1D1 * t102
      t108 = t16 * rhob
      t109 = t108 ** (0.1D1 / 0.3D1)
      t110 = rhob * t109
      t111 = t109 ** 2
      t112 = 0.1D1 / t111
      t113 = sigmabb * t112
      t114 = rhob ** 2
      t115 = 0.1D1 / t114
      t116 = t113 * t115
      t118 = exp(-0.747166092922644636268D-1 * t116)
      t119 = t115 * t118
      t122 = 0.1D1 + 0.747166092922644636268D-1 * t113 * t119
      t123 = t115 * t122
      t127 = 0.1D1 / rhob
      t130 = taub - 0.250000000000000000000D0 * sigmabb * t127
      t131 = t130 * t43
      t132 = rhob ** (0.1D1 / 0.3D1)
      t133 = t132 ** 2
      t135 = 0.1D1 / t133 / rhob
      t136 = t131 * t135
      t138 = 0.1D1 - 0.504755720231149905248D0 * t136
      t139 = t138 ** 2
      t141 = exp(-0.5D0 * t139)
      t144 = 0.118591405585874358362D-1 * t116 + 0.120830459735945720683
     >D0 * t138 * t141
      t145 = t144 ** 2
      t147 = 0.1D1 + 0.143805048498903106908D0 * t113 * t123 + 0.1538461
     >53846153846154D2 * t145
      t149 = 0.65D-1 / t147
      t150 = 0.1D1 / t138
      t151 = t135 * t150
      t154 = tanh(0.504755720231149905248D30 * t131 * t151)
      t156 = 0.500000000000000000000D0 + 0.500000000000000000000D0 * t15
     >4
      t160 = exp(-0.336672065394176986800D0 * t131 * t151 * t156)
      t163 = tanh(0.1D31 - 0.504755720231149905248D30 * t136)
      t164 = 0.500000000000000000000D0 * t163
      t165 = 0.500000000000000000000D0 + t164
      t168 = tanh(0.1D31 * t150)
      t170 = 0.500000000000000000000D0 - 0.500000000000000000000D0 * t16
     >8
      t173 = exp(0.8D0 * t150 * t170)
      t174 = 0.500000000000000000000D0 - t164
      t177 = t160 * t165 - 0.124D1 * t173 * t174
      t178 = 0.109D0 + t149
      t180 = 0.1065D1 - t149 + t177 * t178
      t181 = t20 * t180
      t182 = sqrt(sigmabb)
      t184 = t182 / t109
      t185 = t184 * t127
      t186 = sqrt(t185)
      t189 = exp(-0.943252112663908494661D1 / t186)
      t191 = 0.1D1 - 0.1D1 * t189
      zk(i) = -0.136284044462410474417D1 * t19 * t94 * t104 - 0.13628404
     >4462410474417D1 * t110 * t181 * t191
      t199 = 0.681420222312052372084D0 * t18 * t20 * t93 * t104
      t204 = 0.227140074104017457361D0 * rhoa * t22 * pi * t93 * t104
      t205 = t60 ** 2
      t206 = 0.1D1 / t205
      t209 = sigmaaa / t21 / t17
      t213 = t24 * rhoa
      t214 = 0.1D1 / t213
      t225 = t209 * t25 * t16
      t227 = t23 * t214
      t242 = sigmaaa / t46 / t213
      t243 = t242 * t43
      t246 = 0.1D1 / t46 / t24
      t247 = t44 * t246
      t249 = -0.126188930057787476312D0 * t243 + 0.841259533718583175412
     >D0 * t247
      t258 = -0.958700323326020712721D-1 * t209 * t33 * t16 - 0.28761009
     >6997806213817D0 * t23 * t214 * t32 + 0.143805048498903106908D0 * t
     >23 * t25 * (-0.498110728615096424178D-1 * t209 * t29 * t16 - 0.149
     >433218584528927254D0 * t23 * t214 * t28 + 0.747166092922644636268D
     >-1 * t23 * t25 * (0.498110728615096424178D-1 * t225 + 0.1494332185
     >84528927254D0 * t227) * t28) + 0.307692307692307692308D2 * t57 * (
     >-0.790609370572495722413D-2 * t225 - 0.237182811171748716724D-1 *
     >t227 + 0.120830459735945720683D0 * t249 * t54 - 0.1208304597359457
     >20683D0 * t52 * t249 * t54)
      t261 = t43 * t63
      t265 = t246 * t63
      t269 = 0.1D1 / t52
      t274 = t67 ** 2
      t277 = t63 * (0.1D1 - 0.1D1 * t274)
      t293 = t76 ** 2
      t295 = 0.1D1 - 0.1D1 * t293
      t296 = t73 * t295
      t299 = -0.126188930057787476312D30 * t243 + 0.84125953371858317541
     >2D30 * t247
      t302 = t269 * t83
      t307 = t81 ** 2
      t310 = 0.1D1 / t52 / t51 * (0.1D1 - 0.1D1 * t307)
      t317 = t86 * t295
      t322 = t90 * t206
      t329 = 0.681420222312052372084D0 * t19 * t20 * (0.65D-1 * t206 * t
     >258 + ((-0.841680163485442467001D-1 * t242 * t261 * t69 + 0.561120
     >108990294978001D0 * t44 * t265 * t69 + 0.336672065394176986800D0 *
     > t49 * t269 * t69 * t249 - 0.168336032697088493400D0 * t49 * t277
     >* (0.126188930057787476312D30 * t242 * t261 - 0.841259533718583175
     >412D30 * t44 * t265 - 0.504755720231149905248D30 * t44 * t48 * t26
     >9 * t249)) * t73 * t78 + 0.500000000000000000000D0 * t296 * t299 -
     > 0.124D1 * (-0.8D0 * t302 * t249 + 0.400000000000000000000D30 * t3
     >10 * t249) * t86 * t87 + 0.620000000000000000000D0 * t317 * t299)
     >* t91 - 0.65D-1 * t322 * t258) * t104
      t332 = 0.1D1 / t99 / t98
      t335 = 0.1D1 / t18 / t17
      t346 = 0.973296829181994948224D0 * t19 * t20 * t93 * t332 * (-0.11
     >0064241629820889462D1 * t95 * t335 * t37 * t16 - 0.330192724889462
     >668387D1 * t97 * t25) * t102
      t350 = 0.681420222312052372084D0 * t109 * t20 * t180 * t191
      t355 = 0.227140074104017457361D0 * rhob * t112 * pi * t180 * t191
      t356 = t147 ** 2
      t357 = 0.1D1 / t356
      t360 = sigmabb / t111 / t108
      t364 = t114 * rhob
      t365 = 0.1D1 / t364
      t376 = t360 * t115 * t16
      t378 = t113 * t365
      t393 = sigmabb / t133 / t364
      t394 = t393 * t43
      t397 = 0.1D1 / t133 / t114
      t398 = t131 * t397
      t400 = -0.126188930057787476312D0 * t394 + 0.841259533718583175412
     >D0 * t398
      t409 = -0.958700323326020712721D-1 * t360 * t123 * t16 - 0.2876100
     >96997806213817D0 * t113 * t365 * t122 + 0.143805048498903106908D0
     >* t113 * t115 * (-0.498110728615096424178D-1 * t360 * t119 * t16 -
     > 0.149433218584528927254D0 * t113 * t365 * t118 + 0.74716609292264
     >4636268D-1 * t113 * t115 * (0.498110728615096424178D-1 * t376 + 0.
     >149433218584528927254D0 * t378) * t118) + 0.307692307692307692308D
     >2 * t144 * (-0.790609370572495722413D-2 * t376 - 0.237182811171748
     >716724D-1 * t378 + 0.120830459735945720683D0 * t400 * t141 - 0.120
     >830459735945720683D0 * t139 * t400 * t141)
      t412 = t43 * t150
      t416 = t397 * t150
      t420 = 0.1D1 / t139
      t425 = t154 ** 2
      t428 = t150 * (0.1D1 - 0.1D1 * t425)
      t444 = t163 ** 2
      t446 = 0.1D1 - 0.1D1 * t444
      t447 = t160 * t446
      t450 = -0.126188930057787476312D30 * t394 + 0.84125953371858317541
     >2D30 * t398
      t453 = t420 * t170
      t458 = t168 ** 2
      t461 = 0.1D1 / t139 / t138 * (0.1D1 - 0.1D1 * t458)
      t468 = t173 * t446
      t473 = t177 * t357
      t480 = 0.681420222312052372084D0 * t110 * t20 * (0.65D-1 * t357 * 
     >t409 + ((-0.841680163485442467001D-1 * t393 * t412 * t156 + 0.5611
     >20108990294978001D0 * t131 * t416 * t156 + 0.336672065394176986800
     >D0 * t136 * t420 * t156 * t400 - 0.168336032697088493400D0 * t136
     >* t428 * (0.126188930057787476312D30 * t393 * t412 - 0.84125953371
     >8583175412D30 * t131 * t416 - 0.504755720231149905248D30 * t131 *
     >t135 * t420 * t400)) * t160 * t165 + 0.500000000000000000000D0 * t
     >447 * t450 - 0.124D1 * (-0.8D0 * t453 * t400 + 0.40000000000000000
     >0000D30 * t461 * t400) * t173 * t174 + 0.620000000000000000000D0 *
     > t468 * t450) * t178 - 0.65D-1 * t473 * t409) * t191
      t483 = 0.1D1 / t186 / t185
      t486 = 0.1D1 / t109 / t108
      t497 = 0.973296829181994948224D0 * t110 * t20 * t180 * t483 * (-0.
     >110064241629820889462D1 * t182 * t486 * t127 * t16 - 0.33019272488
     >9462668387D1 * t184 * t115) * t189
      vrhoc(i) = vrhoc(i) - t199 - t204 - t329 + t346 - t350 - t355 - t4
     >80 + t497
      vrhoo(i) = vrhoo(i) - t199 - t204 - t329 + t346 + t350 + t355 + t4
     >80 - t497
      t500 = t22 * t25
      t506 = t24 ** 2
      t516 = t246 * t43
      t526 = 0.143805048498903106908D0 * t500 * t32 + 0.1438050484989031
     >06908D0 * t23 * t25 * (0.747166092922644636268D-1 * t500 * t28 - 0
     >.558257170413290039222D-2 * sigmaaa * t335 / t506 * t28) + 0.30769
     >2307692307692308D2 * t57 * (0.118591405585874358362D-1 * t500 + 0.
     >152474664324695603608D-1 * t516 * t54 - 0.152474664324695603608D-1
     > * t52 * t246 * t43 * t54)
      t529 = t63 * t69
      t533 = 0.1D1 / t41 / t16
      t534 = t40 * t533
      t537 = 0.1D1 / t45 / t506 * t269
      t571 = t19 * t20 * (0.65D-1 * t206 * t526 + ((0.841680163485442467
     >001D-1 * t516 * t529 + 0.424842877124366511988D-1 * t534 * t537 *
     >t69 - 0.168336032697088493400D0 * t49 * t277 * (-0.126188930057787
     >476312D30 * t516 * t63 - 0.636945842765167184389D29 * t534 * t537)
     >) * t73 * t78 + 0.630944650288937381559D29 * t296 * t516 - 0.124D1
     > * (-0.100951144046229981050D0 * t302 * t516 + 0.50475572023114990
     >5248D29 * t310 * t516) * t86 * t87 + 0.782371366358282353132D29 *
     >t317 * t516) * t91 - 0.65D-1 * t322 * t526) * t104
      t572 = 0.340710111156026186042D0 * t571
      t576 = t94 * t332 / t95 * t102
      t577 = 0.803438830384691996314D0 * t576
      t578 = t112 * t115
      t584 = t114 ** 2
      t594 = t397 * t43
      t604 = 0.143805048498903106908D0 * t578 * t122 + 0.143805048498903
     >106908D0 * t113 * t115 * (0.747166092922644636268D-1 * t578 * t118
     > - 0.558257170413290039222D-2 * sigmabb * t486 / t584 * t118) + 0.
     >307692307692307692308D2 * t144 * (0.118591405585874358362D-1 * t57
     >8 + 0.152474664324695603608D-1 * t594 * t141 - 0.15247466432469560
     >3608D-1 * t139 * t397 * t43 * t141)
      t607 = t150 * t156
      t610 = t130 * t533
      t613 = 0.1D1 / t132 / t584 * t420
      t647 = t110 * t20 * (0.65D-1 * t357 * t604 + ((0.84168016348544246
     >7001D-1 * t594 * t607 + 0.424842877124366511988D-1 * t610 * t613 *
     > t156 - 0.168336032697088493400D0 * t136 * t428 * (-0.126188930057
     >787476312D30 * t594 * t150 - 0.636945842765167184389D29 * t610 * t
     >613)) * t160 * t165 + 0.630944650288937381559D29 * t447 * t594 - 0
     >.124D1 * (-0.100951144046229981050D0 * t453 * t594 + 0.50475572023
     >1149905248D29 * t461 * t594) * t173 * t174 + 0.7823713663582823531
     >32D29 * t468 * t594) * t178 - 0.65D-1 * t473 * t604) * t191
      t648 = 0.340710111156026186042D0 * t647
      t652 = t181 * t483 / t182 * t189
      t653 = 0.803438830384691996314D0 * t652
      vsigmacc(i) = vsigmacc(i) - t572 + t577 - t648 + t653
      vsigmaco(i) = vsigmaco(i) - 0.681420222312052372084D0 * t571 + 0.1
     >60687766076938399263D1 * t576 + 0.681420222312052372084D0 * t647 -
     > 0.160687766076938399263D1 * t652
      vsigmaoo(i) = vsigmaoo(i) - t572 + t577 - t648 + t653
      t662 = t43 * t48
      t669 = -0.609898657298782414434D-1 * t662 * t54 + 0.60989865729878
     >2414434D-1 * t52 * t43 * t48 * t54
      t676 = 0.1D1 / t45 / t213 * t269
      t712 = 0.681420222312052372084D0 * t19 * t20 * (0.2000000000000000
     >00000D1 * t206 * t57 * t669 + ((-0.336672065394176986800D0 * t662
     >* t529 - 0.169937150849746604795D0 * t534 * t676 * t69 - 0.1683360
     >32697088493400D0 * t49 * t277 * (0.504755720231149905248D30 * t662
     > * t63 + 0.254778337106066873756D30 * t534 * t676)) * t73 * t78 -
     >0.252377860115574952624D30 * t296 * t662 - 0.124D1 * (0.4038045761
     >84919924197D0 * t302 * t662 - 0.201902288092459962099D30 * t310 *
     >t662) * t86 * t87 - 0.312948546543312941253D30 * t317 * t662) * t9
     >1 - 0.200000000000000000000D1 * t322 * t57 * t669) * t104
      t714 = t43 * t135
      t721 = -0.609898657298782414434D-1 * t714 * t141 + 0.6098986572987
     >82414434D-1 * t139 * t43 * t135 * t141
      t728 = 0.1D1 / t132 / t364 * t420
      t764 = 0.681420222312052372084D0 * t110 * t20 * (0.200000000000000
     >000000D1 * t357 * t144 * t721 + ((-0.336672065394176986800D0 * t71
     >4 * t607 - 0.169937150849746604795D0 * t610 * t728 * t156 - 0.1683
     >36032697088493400D0 * t136 * t428 * (0.504755720231149905248D30 *
     >t714 * t150 + 0.254778337106066873756D30 * t610 * t728)) * t160 *
     >t165 - 0.252377860115574952624D30 * t447 * t714 - 0.124D1 * (0.403
     >804576184919924197D0 * t453 * t714 - 0.201902288092459962099D30 *
     >t461 * t714) * t173 * t174 - 0.312948546543312941253D30 * t468 * t
     >714) * t178 - 0.200000000000000000000D1 * t473 * t144 * t721) * t1
     >91
      vtauc(i) = vtauc(i) - t712 - t764
      vtauo(i) = vtauo(i) - t712 + t764

               endif
             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      sigmabb = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) - 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmabb
      taub = max(0.0D0, 0.500000000000000000000D0 * tauc(i) - 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taub
      t13 = pi ** 2
      t15 = (t13 * rhob) ** (0.1D1 / 0.3D1)
      t18 = t15 ** 2
      t20 = sigmabb / t18
      t21 = rhob ** 2
      t22 = 0.1D1 / t21
      t23 = t20 * t22
      t25 = exp(-0.747166092922644636268D-1 * t23)
      t34 = 0.1D1 / rhob
      t38 = t13 ** (0.1D1 / 0.3D1)
      t39 = t38 ** 2
      t41 = (taub - 0.250000000000000000000D0 * sigmabb * t34) / t39
      t42 = rhob ** (0.1D1 / 0.3D1)
      t43 = t42 ** 2
      t45 = 0.1D1 / t43 / rhob
      t46 = t41 * t45
      t48 = 0.1D1 - 0.504755720231149905248D0 * t46
      t49 = t48 ** 2
      t51 = exp(-0.5D0 * t49)
      t55 = (0.118591405585874358362D-1 * t23 + 0.120830459735945720683D
     >0 * t48 * t51) ** 2
      t59 = 0.65D-1 / (0.1D1 + 0.143805048498903106908D0 * t20 * t22 * (
     >0.1D1 + 0.747166092922644636268D-1 * t20 * t22 * t25) + 0.15384615
     >3846153846154D2 * t55)
      t60 = 0.1D1 / t48
      t61 = t45 * t60
      t64 = tanh(0.504755720231149905248D30 * t41 * t61)
      t70 = exp(-0.336672065394176986800D0 * t41 * t61 * (0.500000000000
     >000000000D0 + 0.500000000000000000000D0 * t64))
      t73 = tanh(0.1D31 - 0.504755720231149905248D30 * t46)
      t74 = 0.500000000000000000000D0 * t73
      t78 = tanh(0.1D31 * t60)
      t83 = exp(0.8D0 * t60 * (0.500000000000000000000D0 - 0.50000000000
     >0000000000D0 * t78))
      t92 = sqrt(sigmabb)
      t96 = sqrt(t92 / t15 * t34)
      t99 = exp(-0.943252112663908494661D1 / t96)
      zk(i) = -0.136284044462410474417D1 * rhob * t15 / pi * (0.1065D1 -
     > t59 + (t70 * (0.500000000000000000000D0 + t74) - 0.124D1 * t83 *
     >(0.500000000000000000000D0 - t74)) * (0.109D0 + t59)) * (0.1D1 - 0
     >.1D1 * t99)

               elseif(rhob.lt.tol) then
             rho = rhoa
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) + 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmaaa
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i) + 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taua
      t13 = pi ** 2
      t15 = (t13 * rhoa) ** (0.1D1 / 0.3D1)
      t18 = t15 ** 2
      t20 = sigmaaa / t18
      t21 = rhoa ** 2
      t22 = 0.1D1 / t21
      t23 = t20 * t22
      t25 = exp(-0.747166092922644636268D-1 * t23)
      t34 = 0.1D1 / rhoa
      t38 = t13 ** (0.1D1 / 0.3D1)
      t39 = t38 ** 2
      t41 = (taua - 0.250000000000000000000D0 * sigmaaa * t34) / t39
      t42 = rhoa ** (0.1D1 / 0.3D1)
      t43 = t42 ** 2
      t45 = 0.1D1 / t43 / rhoa
      t46 = t41 * t45
      t48 = 0.1D1 - 0.504755720231149905248D0 * t46
      t49 = t48 ** 2
      t51 = exp(-0.5D0 * t49)
      t55 = (0.118591405585874358362D-1 * t23 + 0.120830459735945720683D
     >0 * t48 * t51) ** 2
      t59 = 0.65D-1 / (0.1D1 + 0.143805048498903106908D0 * t20 * t22 * (
     >0.1D1 + 0.747166092922644636268D-1 * t20 * t22 * t25) + 0.15384615
     >3846153846154D2 * t55)
      t60 = 0.1D1 / t48
      t61 = t45 * t60
      t64 = tanh(0.504755720231149905248D30 * t41 * t61)
      t70 = exp(-0.336672065394176986800D0 * t41 * t61 * (0.500000000000
     >000000000D0 + 0.500000000000000000000D0 * t64))
      t73 = tanh(0.1D31 - 0.504755720231149905248D30 * t46)
      t74 = 0.500000000000000000000D0 * t73
      t78 = tanh(0.1D31 * t60)
      t83 = exp(0.8D0 * t60 * (0.500000000000000000000D0 - 0.50000000000
     >0000000000D0 * t78))
      t92 = sqrt(sigmaaa)
      t96 = sqrt(t92 / t15 * t34)
      t99 = exp(-0.943252112663908494661D1 / t96)
      zk(i) = -0.136284044462410474417D1 * rhoa * t15 / pi * (0.1065D1 -
     > t59 + (t70 * (0.500000000000000000000D0 + t74) - 0.124D1 * t83 *
     >(0.500000000000000000000D0 - t74)) * (0.109D0 + t59)) * (0.1D1 - 0
     >.1D1 * t99)

               else
             rho = rhoa + rhob
      t2 = 0.250000000000000000000D0 * sigmacc(i)
      t4 = 0.250000000000000000000D0 * sigmaoo(i)
      t6 = 0.500000000000000000000D0 * sigmaco(i)
      sigmaaa = max(0.0D0, t2 + t4 + t6)
      sigmaab = t2 - t4
      sigmabb = max(0.0D0, t2 + t4 - t6)
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t11 = 0.500000000000000000000D0 * tauc(i)
      t13 = 0.500000000000000000000D0 * tauo(i)
      taua = max(0.0D0, t11 + t13)
      taub = max(0.0D0, t11 - t13)
      tau = taua + taub
      t16 = pi ** 2
      t18 = (t16 * rhoa) ** (0.1D1 / 0.3D1)
      t20 = 0.1D1 / pi
      t21 = t18 ** 2
      t23 = sigmaaa / t21
      t24 = rhoa ** 2
      t25 = 0.1D1 / t24
      t26 = t23 * t25
      t28 = exp(-0.747166092922644636268D-1 * t26)
      t37 = 0.1D1 / rhoa
      t41 = t16 ** (0.1D1 / 0.3D1)
      t42 = t41 ** 2
      t43 = 0.1D1 / t42
      t44 = (taua - 0.250000000000000000000D0 * sigmaaa * t37) * t43
      t45 = rhoa ** (0.1D1 / 0.3D1)
      t46 = t45 ** 2
      t48 = 0.1D1 / t46 / rhoa
      t49 = t44 * t48
      t51 = 0.1D1 - 0.504755720231149905248D0 * t49
      t52 = t51 ** 2
      t54 = exp(-0.5D0 * t52)
      t58 = (0.118591405585874358362D-1 * t26 + 0.120830459735945720683D
     >0 * t51 * t54) ** 2
      t62 = 0.65D-1 / (0.1D1 + 0.143805048498903106908D0 * t23 * t25 * (
     >0.1D1 + 0.747166092922644636268D-1 * t23 * t25 * t28) + 0.15384615
     >3846153846154D2 * t58)
      t63 = 0.1D1 / t51
      t64 = t48 * t63
      t67 = tanh(0.504755720231149905248D30 * t44 * t64)
      t73 = exp(-0.336672065394176986800D0 * t44 * t64 * (0.500000000000
     >000000000D0 + 0.500000000000000000000D0 * t67))
      t76 = tanh(0.1D31 - 0.504755720231149905248D30 * t49)
      t77 = 0.500000000000000000000D0 * t76
      t81 = tanh(0.1D31 * t63)
      t86 = exp(0.8D0 * t63 * (0.500000000000000000000D0 - 0.50000000000
     >0000000000D0 * t81))
      t95 = sqrt(sigmaaa)
      t99 = sqrt(t95 / t18 * t37)
      t102 = exp(-0.943252112663908494661D1 / t99)
      t109 = (t16 * rhob) ** (0.1D1 / 0.3D1)
      t111 = t109 ** 2
      t113 = sigmabb / t111
      t114 = rhob ** 2
      t115 = 0.1D1 / t114
      t116 = t113 * t115
      t118 = exp(-0.747166092922644636268D-1 * t116)
      t127 = 0.1D1 / rhob
      t131 = (taub - 0.250000000000000000000D0 * sigmabb * t127) * t43
      t132 = rhob ** (0.1D1 / 0.3D1)
      t133 = t132 ** 2
      t135 = 0.1D1 / t133 / rhob
      t136 = t131 * t135
      t138 = 0.1D1 - 0.504755720231149905248D0 * t136
      t139 = t138 ** 2
      t141 = exp(-0.5D0 * t139)
      t145 = (0.118591405585874358362D-1 * t116 + 0.12083045973594572068
     >3D0 * t138 * t141) ** 2
      t149 = 0.65D-1 / (0.1D1 + 0.143805048498903106908D0 * t113 * t115 
     >* (0.1D1 + 0.747166092922644636268D-1 * t113 * t115 * t118) + 0.15
     >3846153846153846154D2 * t145)
      t150 = 0.1D1 / t138
      t151 = t135 * t150
      t154 = tanh(0.504755720231149905248D30 * t131 * t151)
      t160 = exp(-0.336672065394176986800D0 * t131 * t151 * (0.500000000
     >000000000000D0 + 0.500000000000000000000D0 * t154))
      t163 = tanh(0.1D31 - 0.504755720231149905248D30 * t136)
      t164 = 0.500000000000000000000D0 * t163
      t168 = tanh(0.1D31 * t150)
      t173 = exp(0.8D0 * t150 * (0.500000000000000000000D0 - 0.500000000
     >000000000000D0 * t168))
      t182 = sqrt(sigmabb)
      t186 = sqrt(t182 / t109 * t127)
      t189 = exp(-0.943252112663908494661D1 / t186)
      zk(i) = -0.136284044462410474417D1 * rhoa * t18 * t20 * (0.1065D1 
     >- t62 + (t73 * (0.500000000000000000000D0 + t77) - 0.124D1 * t86 *
     > (0.500000000000000000000D0 - t77)) * (0.109D0 + t62)) * (0.1D1 -
     >0.1D1 * t102) - 0.136284044462410474417D1 * rhob * t109 * t20 * (0
     >.1065D1 - t149 + (t160 * (0.500000000000000000000D0 + t164) - 0.12
     >4D1 * t173 * (0.500000000000000000000D0 - t164)) * (0.109D0 + t149
     >)) * (0.1D1 - 0.1D1 * t189)

               endif
             endif
           enddo
         endif
       else
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i))
      sigmaab = sigmaaa
      sigmabb = sigmaab
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i))
      taub = taua
      tau = taua + taub
      t8 = pi ** 2
      t9 = t8 * rhoa
      t10 = t9 ** (0.1D1 / 0.3D1)
      t11 = rhoa * t10
      t12 = 0.1D1 / pi
      t13 = t10 ** 2
      t14 = 0.1D1 / t13
      t15 = sigmaaa * t14
      t16 = rhoa ** 2
      t17 = 0.1D1 / t16
      t18 = t15 * t17
      t20 = exp(-0.747166092922644636268D-1 * t18)
      t21 = t17 * t20
      t24 = 0.1D1 + 0.747166092922644636268D-1 * t15 * t21
      t25 = t17 * t24
      t29 = 0.1D1 / rhoa
      t32 = taua - 0.250000000000000000000D0 * sigmaaa * t29
      t33 = t8 ** (0.1D1 / 0.3D1)
      t34 = t33 ** 2
      t35 = 0.1D1 / t34
      t36 = t32 * t35
      t37 = rhoa ** (0.1D1 / 0.3D1)
      t38 = t37 ** 2
      t40 = 0.1D1 / t38 / rhoa
      t41 = t36 * t40
      t43 = 0.1D1 - 0.504755720231149905248D0 * t41
      t44 = t43 ** 2
      t46 = exp(-0.5D0 * t44)
      t49 = 0.118591405585874358362D-1 * t18 + 0.120830459735945720683D0
     > * t43 * t46
      t50 = t49 ** 2
      t52 = 0.1D1 + 0.143805048498903106908D0 * t15 * t25 + 0.1538461538
     >46153846154D2 * t50
      t54 = 0.65D-1 / t52
      t55 = 0.1D1 / t43
      t56 = t40 * t55
      t59 = tanh(0.504755720231149905248D30 * t36 * t56)
      t61 = 0.500000000000000000000D0 + 0.500000000000000000000D0 * t59
      t65 = exp(-0.336672065394176986800D0 * t36 * t56 * t61)
      t68 = tanh(0.1D31 - 0.504755720231149905248D30 * t41)
      t69 = 0.500000000000000000000D0 * t68
      t70 = 0.500000000000000000000D0 + t69
      t73 = tanh(0.1D31 * t55)
      t75 = 0.500000000000000000000D0 - 0.500000000000000000000D0 * t73
      t78 = exp(0.8D0 * t55 * t75)
      t79 = 0.500000000000000000000D0 - t69
      t82 = t65 * t70 - 0.124D1 * t78 * t79
      t83 = 0.109D0 + t54
      t85 = 0.1065D1 - t54 + t82 * t83
      t86 = t12 * t85
      t87 = sqrt(sigmaaa)
      t89 = t87 / t10
      t90 = t89 * t29
      t91 = sqrt(t90)
      t94 = exp(-0.943252112663908494661D1 / t91)
      t96 = 0.1D1 - 0.1D1 * t94
      t100 = t8 * rhob
      t101 = t100 ** (0.1D1 / 0.3D1)
      t102 = rhob * t101
      t103 = t101 ** 2
      t104 = 0.1D1 / t103
      t105 = sigmabb * t104
      t106 = rhob ** 2
      t107 = 0.1D1 / t106
      t108 = t105 * t107
      t110 = exp(-0.747166092922644636268D-1 * t108)
      t111 = t107 * t110
      t114 = 0.1D1 + 0.747166092922644636268D-1 * t105 * t111
      t115 = t107 * t114
      t119 = 0.1D1 / rhob
      t122 = taub - 0.250000000000000000000D0 * sigmabb * t119
      t123 = t122 * t35
      t124 = rhob ** (0.1D1 / 0.3D1)
      t125 = t124 ** 2
      t127 = 0.1D1 / t125 / rhob
      t128 = t123 * t127
      t130 = 0.1D1 - 0.504755720231149905248D0 * t128
      t131 = t130 ** 2
      t133 = exp(-0.5D0 * t131)
      t136 = 0.118591405585874358362D-1 * t108 + 0.120830459735945720683
     >D0 * t130 * t133
      t137 = t136 ** 2
      t139 = 0.1D1 + 0.143805048498903106908D0 * t105 * t115 + 0.1538461
     >53846153846154D2 * t137
      t141 = 0.65D-1 / t139
      t142 = 0.1D1 / t130
      t143 = t127 * t142
      t146 = tanh(0.504755720231149905248D30 * t123 * t143)
      t148 = 0.500000000000000000000D0 + 0.500000000000000000000D0 * t14
     >6
      t152 = exp(-0.336672065394176986800D0 * t123 * t143 * t148)
      t155 = tanh(0.1D31 - 0.504755720231149905248D30 * t128)
      t156 = 0.500000000000000000000D0 * t155
      t157 = 0.500000000000000000000D0 + t156
      t160 = tanh(0.1D31 * t142)
      t162 = 0.500000000000000000000D0 - 0.500000000000000000000D0 * t16
     >0
      t165 = exp(0.8D0 * t142 * t162)
      t166 = 0.500000000000000000000D0 - t156
      t169 = t152 * t157 - 0.124D1 * t165 * t166
      t170 = 0.109D0 + t141
      t172 = 0.1065D1 - t141 + t169 * t170
      t173 = t12 * t172
      t174 = sqrt(sigmabb)
      t176 = t174 / t101
      t177 = t176 * t119
      t178 = sqrt(t177)
      t181 = exp(-0.943252112663908494661D1 / t178)
      t183 = 0.1D1 - 0.1D1 * t181
      zk(i) = -0.136284044462410474417D1 * t11 * t86 * t96 - 0.136284044
     >462410474417D1 * t102 * t173 * t183
      t197 = t52 ** 2
      t198 = 0.1D1 / t197
      t201 = sigmaaa / t13 / t9
      t205 = t16 * rhoa
      t206 = 0.1D1 / t205
      t217 = t201 * t17 * t8
      t219 = t15 * t206
      t234 = sigmaaa / t38 / t205
      t235 = t234 * t35
      t238 = 0.1D1 / t38 / t16
      t239 = t36 * t238
      t241 = -0.126188930057787476312D0 * t235 + 0.841259533718583175412
     >D0 * t239
      t250 = -0.958700323326020712721D-1 * t201 * t25 * t8 - 0.287610096
     >997806213817D0 * t15 * t206 * t24 + 0.143805048498903106908D0 * t1
     >5 * t17 * (-0.498110728615096424178D-1 * t201 * t21 * t8 - 0.14943
     >3218584528927254D0 * t15 * t206 * t20 + 0.747166092922644636268D-1
     > * t15 * t17 * (0.498110728615096424178D-1 * t217 + 0.149433218584
     >528927254D0 * t219) * t20) + 0.307692307692307692308D2 * t49 * (-0
     >.790609370572495722413D-2 * t217 - 0.237182811171748716724D-1 * t2
     >19 + 0.120830459735945720683D0 * t241 * t46 - 0.120830459735945720
     >683D0 * t44 * t241 * t46)
      t253 = t35 * t55
      t257 = t238 * t55
      t261 = 0.1D1 / t44
      t266 = t59 ** 2
      t269 = t55 * (0.1D1 - 0.1D1 * t266)
      t285 = t68 ** 2
      t287 = 0.1D1 - 0.1D1 * t285
      t288 = t65 * t287
      t291 = -0.126188930057787476312D30 * t235 + 0.84125953371858317541
     >2D30 * t239
      t294 = t261 * t75
      t299 = t73 ** 2
      t302 = 0.1D1 / t44 / t43 * (0.1D1 - 0.1D1 * t299)
      t309 = t78 * t287
      t314 = t82 * t198
      t324 = 0.1D1 / t91 / t90
      t327 = 0.1D1 / t10 / t9
      t348 = t139 ** 2
      t349 = 0.1D1 / t348
      t352 = sigmabb / t103 / t100
      t356 = t106 * rhob
      t357 = 0.1D1 / t356
      t368 = t352 * t107 * t8
      t370 = t105 * t357
      t385 = sigmabb / t125 / t356
      t386 = t385 * t35
      t389 = 0.1D1 / t125 / t106
      t390 = t123 * t389
      t392 = -0.126188930057787476312D0 * t386 + 0.841259533718583175412
     >D0 * t390
      t401 = -0.958700323326020712721D-1 * t352 * t115 * t8 - 0.28761009
     >6997806213817D0 * t105 * t357 * t114 + 0.143805048498903106908D0 *
     > t105 * t107 * (-0.498110728615096424178D-1 * t352 * t111 * t8 - 0
     >.149433218584528927254D0 * t105 * t357 * t110 + 0.7471660929226446
     >36268D-1 * t105 * t107 * (0.498110728615096424178D-1 * t368 + 0.14
     >9433218584528927254D0 * t370) * t110) + 0.307692307692307692308D2
     >* t136 * (-0.790609370572495722413D-2 * t368 - 0.23718281117174871
     >6724D-1 * t370 + 0.120830459735945720683D0 * t392 * t133 - 0.12083
     >0459735945720683D0 * t131 * t392 * t133)
      t404 = t35 * t142
      t408 = t389 * t142
      t412 = 0.1D1 / t131
      t417 = t146 ** 2
      t420 = t142 * (0.1D1 - 0.1D1 * t417)
      t436 = t155 ** 2
      t438 = 0.1D1 - 0.1D1 * t436
      t439 = t152 * t438
      t442 = -0.126188930057787476312D30 * t386 + 0.84125953371858317541
     >2D30 * t390
      t445 = t412 * t162
      t450 = t160 ** 2
      t453 = 0.1D1 / t131 / t130 * (0.1D1 - 0.1D1 * t450)
      t460 = t165 * t438
      t465 = t169 * t349
      t475 = 0.1D1 / t178 / t177
      t478 = 0.1D1 / t101 / t100
      vrhoc(i) = vrhoc(i) - 0.681420222312052372084D0 * t10 * t12 * t85 
     >* t96 - 0.227140074104017457361D0 * rhoa * t14 * pi * t85 * t96 -
     >0.681420222312052372084D0 * t11 * t12 * (0.65D-1 * t198 * t250 + (
     >(-0.841680163485442467001D-1 * t234 * t253 * t61 + 0.5611201089902
     >94978001D0 * t36 * t257 * t61 + 0.336672065394176986800D0 * t41 *
     >t261 * t61 * t241 - 0.168336032697088493400D0 * t41 * t269 * (0.12
     >6188930057787476312D30 * t234 * t253 - 0.841259533718583175412D30
     >* t36 * t257 - 0.504755720231149905248D30 * t36 * t40 * t261 * t24
     >1)) * t65 * t70 + 0.500000000000000000000D0 * t288 * t291 - 0.124D
     >1 * (-0.8D0 * t294 * t241 + 0.400000000000000000000D30 * t302 * t2
     >41) * t78 * t79 + 0.620000000000000000000D0 * t309 * t291) * t83 -
     > 0.65D-1 * t314 * t250) * t96 + 0.973296829181994948224D0 * t11 *
     >t12 * t85 * t324 * (-0.110064241629820889462D1 * t87 * t327 * t29
     >* t8 - 0.330192724889462668387D1 * t89 * t17) * t94 - 0.6814202223
     >12052372084D0 * t101 * t12 * t172 * t183 - 0.227140074104017457361
     >D0 * rhob * t104 * pi * t172 * t183 - 0.681420222312052372084D0 *
     >t102 * t12 * (0.65D-1 * t349 * t401 + ((-0.841680163485442467001D-
     >1 * t385 * t404 * t148 + 0.561120108990294978001D0 * t123 * t408 *
     > t148 + 0.336672065394176986800D0 * t128 * t412 * t148 * t392 - 0.
     >168336032697088493400D0 * t128 * t420 * (0.126188930057787476312D3
     >0 * t385 * t404 - 0.841259533718583175412D30 * t123 * t408 - 0.504
     >755720231149905248D30 * t123 * t127 * t412 * t392)) * t152 * t157
     >+ 0.500000000000000000000D0 * t439 * t442 - 0.124D1 * (-0.8D0 * t4
     >45 * t392 + 0.400000000000000000000D30 * t453 * t392) * t165 * t16
     >6 + 0.620000000000000000000D0 * t460 * t442) * t170 - 0.65D-1 * t4
     >65 * t401) * t183 + 0.973296829181994948224D0 * t102 * t12 * t172
     >* t475 * (-0.110064241629820889462D1 * t174 * t478 * t119 * t8 - 0
     >.330192724889462668387D1 * t176 * t107) * t181
      t491 = t14 * t17
      t497 = t16 ** 2
      t507 = t238 * t35
      t517 = 0.143805048498903106908D0 * t491 * t24 + 0.1438050484989031
     >06908D0 * t15 * t17 * (0.747166092922644636268D-1 * t491 * t20 - 0
     >.558257170413290039222D-2 * sigmaaa * t327 / t497 * t20) + 0.30769
     >2307692307692308D2 * t49 * (0.118591405585874358362D-1 * t491 + 0.
     >152474664324695603608D-1 * t507 * t46 - 0.152474664324695603608D-1
     > * t44 * t238 * t35 * t46)
      t520 = t55 * t61
      t524 = 0.1D1 / t33 / t8
      t525 = t32 * t524
      t528 = 0.1D1 / t37 / t497 * t261
      t569 = t104 * t107
      t575 = t106 ** 2
      t585 = t389 * t35
      t595 = 0.143805048498903106908D0 * t569 * t114 + 0.143805048498903
     >106908D0 * t105 * t107 * (0.747166092922644636268D-1 * t569 * t110
     > - 0.558257170413290039222D-2 * sigmabb * t478 / t575 * t110) + 0.
     >307692307692307692308D2 * t136 * (0.118591405585874358362D-1 * t56
     >9 + 0.152474664324695603608D-1 * t585 * t133 - 0.15247466432469560
     >3608D-1 * t131 * t389 * t35 * t133)
      t598 = t142 * t148
      t601 = t122 * t524
      t604 = 0.1D1 / t124 / t575 * t412
      vsigmacc(i) = vsigmacc(i) - 0.340710111156026186042D0 * t11 * t12 
     >* (0.65D-1 * t198 * t517 + ((0.841680163485442467001D-1 * t507 * t
     >520 + 0.424842877124366511988D-1 * t525 * t528 * t61 - 0.168336032
     >697088493400D0 * t41 * t269 * (-0.126188930057787476312D30 * t507
     >* t55 - 0.636945842765167184389D29 * t525 * t528)) * t65 * t70 + 0
     >.630944650288937381559D29 * t288 * t507 - 0.124D1 * (-0.1009511440
     >46229981050D0 * t294 * t507 + 0.504755720231149905248D29 * t302 *
     >t507) * t78 * t79 + 0.782371366358282353132D29 * t309 * t507) * t8
     >3 - 0.65D-1 * t314 * t517) * t96 + 0.803438830384691996314D0 * t86
     > * t324 / t87 * t94 - 0.340710111156026186042D0 * t102 * t12 * (0.
     >65D-1 * t349 * t595 + ((0.841680163485442467001D-1 * t585 * t598 +
     > 0.424842877124366511988D-1 * t601 * t604 * t148 - 0.1683360326970
     >88493400D0 * t128 * t420 * (-0.126188930057787476312D30 * t585 * t
     >142 - 0.636945842765167184389D29 * t601 * t604)) * t152 * t157 + 0
     >.630944650288937381559D29 * t439 * t585 - 0.124D1 * (-0.1009511440
     >46229981050D0 * t445 * t585 + 0.504755720231149905248D29 * t453 *
     >t585) * t165 * t166 + 0.782371366358282353132D29 * t460 * t585) *
     >t170 - 0.65D-1 * t465 * t595) * t183 + 0.803438830384691996314D0 *
     > t173 * t475 / t174 * t181
      t647 = t35 * t40
      t654 = -0.609898657298782414434D-1 * t647 * t46 + 0.60989865729878
     >2414434D-1 * t44 * t35 * t40 * t46
      t661 = 0.1D1 / t37 / t205 * t261
      t699 = t35 * t127
      t706 = -0.609898657298782414434D-1 * t699 * t133 + 0.6098986572987
     >82414434D-1 * t131 * t35 * t127 * t133
      t713 = 0.1D1 / t124 / t356 * t412
      vtauc(i) = vtauc(i) - 0.681420222312052372084D0 * t11 * t12 * (0.2
     >00000000000000000000D1 * t198 * t49 * t654 + ((-0.3366720653941769
     >86800D0 * t647 * t520 - 0.169937150849746604795D0 * t525 * t661 *
     >t61 - 0.168336032697088493400D0 * t41 * t269 * (0.5047557202311499
     >05248D30 * t647 * t55 + 0.254778337106066873756D30 * t525 * t661))
     > * t65 * t70 - 0.252377860115574952624D30 * t288 * t647 - 0.124D1
     >* (0.403804576184919924197D0 * t294 * t647 - 0.2019022880924599620
     >99D30 * t302 * t647) * t78 * t79 - 0.312948546543312941253D30 * t3
     >09 * t647) * t83 - 0.200000000000000000000D1 * t314 * t49 * t654)
     >* t96 - 0.681420222312052372084D0 * t102 * t12 * (0.20000000000000
     >0000000D1 * t349 * t136 * t706 + ((-0.336672065394176986800D0 * t6
     >99 * t598 - 0.169937150849746604795D0 * t601 * t713 * t148 - 0.168
     >336032697088493400D0 * t128 * t420 * (0.504755720231149905248D30 *
     > t699 * t142 + 0.254778337106066873756D30 * t601 * t713)) * t152 *
     > t157 - 0.252377860115574952624D30 * t439 * t699 - 0.124D1 * (0.40
     >3804576184919924197D0 * t445 * t699 - 0.201902288092459962099D30 *
     > t453 * t699) * t165 * t166 - 0.312948546543312941253D30 * t460 *
     >t699) * t170 - 0.200000000000000000000D1 * t465 * t136 * t706) * t
     >183

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i))
      sigmaab = sigmaaa
      sigmabb = sigmaab
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i))
      taub = taua
      tau = taua + taub
      t8 = pi ** 2
      t10 = (t8 * rhoa) ** (0.1D1 / 0.3D1)
      t12 = 0.1D1 / pi
      t13 = t10 ** 2
      t15 = sigmaaa / t13
      t16 = rhoa ** 2
      t17 = 0.1D1 / t16
      t18 = t15 * t17
      t20 = exp(-0.747166092922644636268D-1 * t18)
      t29 = 0.1D1 / rhoa
      t33 = t8 ** (0.1D1 / 0.3D1)
      t34 = t33 ** 2
      t35 = 0.1D1 / t34
      t36 = (taua - 0.250000000000000000000D0 * sigmaaa * t29) * t35
      t37 = rhoa ** (0.1D1 / 0.3D1)
      t38 = t37 ** 2
      t40 = 0.1D1 / t38 / rhoa
      t41 = t36 * t40
      t43 = 0.1D1 - 0.504755720231149905248D0 * t41
      t44 = t43 ** 2
      t46 = exp(-0.5D0 * t44)
      t50 = (0.118591405585874358362D-1 * t18 + 0.120830459735945720683D
     >0 * t43 * t46) ** 2
      t54 = 0.65D-1 / (0.1D1 + 0.143805048498903106908D0 * t15 * t17 * (
     >0.1D1 + 0.747166092922644636268D-1 * t15 * t17 * t20) + 0.15384615
     >3846153846154D2 * t50)
      t55 = 0.1D1 / t43
      t56 = t40 * t55
      t59 = tanh(0.504755720231149905248D30 * t36 * t56)
      t65 = exp(-0.336672065394176986800D0 * t36 * t56 * (0.500000000000
     >000000000D0 + 0.500000000000000000000D0 * t59))
      t68 = tanh(0.1D31 - 0.504755720231149905248D30 * t41)
      t69 = 0.500000000000000000000D0 * t68
      t73 = tanh(0.1D31 * t55)
      t78 = exp(0.8D0 * t55 * (0.500000000000000000000D0 - 0.50000000000
     >0000000000D0 * t73))
      t87 = sqrt(sigmaaa)
      t91 = sqrt(t87 / t10 * t29)
      t94 = exp(-0.943252112663908494661D1 / t91)
      t101 = (t8 * rhob) ** (0.1D1 / 0.3D1)
      t103 = t101 ** 2
      t105 = sigmabb / t103
      t106 = rhob ** 2
      t107 = 0.1D1 / t106
      t108 = t105 * t107
      t110 = exp(-0.747166092922644636268D-1 * t108)
      t119 = 0.1D1 / rhob
      t123 = (taub - 0.250000000000000000000D0 * sigmabb * t119) * t35
      t124 = rhob ** (0.1D1 / 0.3D1)
      t125 = t124 ** 2
      t127 = 0.1D1 / t125 / rhob
      t128 = t123 * t127
      t130 = 0.1D1 - 0.504755720231149905248D0 * t128
      t131 = t130 ** 2
      t133 = exp(-0.5D0 * t131)
      t137 = (0.118591405585874358362D-1 * t108 + 0.12083045973594572068
     >3D0 * t130 * t133) ** 2
      t141 = 0.65D-1 / (0.1D1 + 0.143805048498903106908D0 * t105 * t107 
     >* (0.1D1 + 0.747166092922644636268D-1 * t105 * t107 * t110) + 0.15
     >3846153846153846154D2 * t137)
      t142 = 0.1D1 / t130
      t143 = t127 * t142
      t146 = tanh(0.504755720231149905248D30 * t123 * t143)
      t152 = exp(-0.336672065394176986800D0 * t123 * t143 * (0.500000000
     >000000000000D0 + 0.500000000000000000000D0 * t146))
      t155 = tanh(0.1D31 - 0.504755720231149905248D30 * t128)
      t156 = 0.500000000000000000000D0 * t155
      t160 = tanh(0.1D31 * t142)
      t165 = exp(0.8D0 * t142 * (0.500000000000000000000D0 - 0.500000000
     >000000000000D0 * t160))
      t174 = sqrt(sigmabb)
      t178 = sqrt(t174 / t101 * t119)
      t181 = exp(-0.943252112663908494661D1 / t178)
      zk(i) = -0.136284044462410474417D1 * rhoa * t10 * t12 * (0.1065D1 
     >- t54 + (t65 * (0.500000000000000000000D0 + t69) - 0.124D1 * t78 *
     > (0.500000000000000000000D0 - t69)) * (0.109D0 + t54)) * (0.1D1 -
     >0.1D1 * t94) - 0.136284044462410474417D1 * rhob * t101 * t12 * (0.
     >1065D1 - t141 + (t152 * (0.500000000000000000000D0 + t156) - 0.124
     >D1 * t165 * (0.500000000000000000000D0 - t156)) * (0.109D0 + t141)
     >) * (0.1D1 - 0.1D1 * t181)

             endif
           enddo
         endif
       endif

      return
      end

c:SCANXsubrend
c:TPSSCsubrstart

c    Generated: Mon Jun 20 16:05:25 CEST 2016

      subroutine dftacg_tpssc
     > (name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   tauc,tauo,upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,
     >                   vtauc,vtauo,vupsilonc,vupsilono)
      implicit double precision (a-h,o-z)
      logical fderiv,open
      integer igrad,npt
      character*(*) name
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)
      double precision tauc(*),tauo(*)
      double precision upsilonc(*),upsilono(*)
      double precision zk(*),vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)
      double precision vtauc(*),vtauo(*)
      double precision vupsilonc(*),vupsilono(*)
      include "common/cdft"
      include "common/tapes"
      parameter(tol=1d-12)
      pi=acos(-1d0)
      name='Automatically generated TPSSC'
      igrad=2
       if(open) then
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      sigmabb = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) - 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmabb
      taub = max(0.0D0, 0.500000000000000000000D0 * tauc(i) - 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taub
      t13 = 0.1D1 / pi
      t14 = 0.1D1 / rhob
      t15 = t13 * t14
      t16 = t15 ** (0.1D1 / 0.3D1)
      t18 = 0.1D1 + 0.186690969707574028554D0 * t16
      t19 = t15 ** (0.1D1 / 0.6D1)
      t22 = sqrt(t15)
      t24 = t16 ** 2
      t26 = 0.134579137143944477912D2 * t19 + 0.563098414909787598194D1 
     >* t16 + 0.291521471421917737271D1 * t22 + 0.516066464547863440989D
     >0 * t24
      t29 = 0.1D1 + 0.321639589973850701335D2 / t26
      t30 = log(t29)
      t31 = t18 * t30
      t33 = pi ** 2
      t34 = 0.1D1 / t33
      t35 = t33 * pi
      t37 = t33 * rhob
      t38 = t37 ** (0.1D1 / 0.3D1)
      t39 = 0.1D1 / t38
      t40 = t35 * sigmabb * t39
      t41 = rhob ** 2
      t42 = 0.1D1 / t41
      t45 = exp(0.202642426794280972788D0 * t31 * t33)
      t46 = t45 - 0.1D1
      t47 = 0.1D1 / t46
      t48 = t35 * t47
      t49 = sigmabb * t39
      t52 = 0.149583857013095144140D-1 * t48 * t49 * t42
      t53 = 0.1D1 + t52
      t54 = t42 * t53
      t55 = t33 ** 2
      t56 = t55 * t33
      t57 = t46 ** 2
      t58 = 0.1D1 / t57
      t59 = t56 * t58
      t60 = sigmabb ** 2
      t61 = t38 ** 2
      t62 = 0.1D1 / t61
      t63 = t60 * t62
      t64 = t41 ** 2
      t65 = 0.1D1 / t64
      t69 = 0.1D1 + t52 + 0.223753302789140933370D-3 * t59 * t63 * t65
      t70 = 0.1D1 / t69
      t71 = t54 * t70
      t74 = 0.1D1 + 0.149583857013095144140D-1 * t40 * t71
      t75 = log(t74)
      t78 = -0.3109070D-1 * t31 + 0.153426409720027345292D0 * t34 * t75
      t79 = t60 * t42
      t80 = taub ** 2
      t81 = 0.1D1 / t80
      t82 = t79 * t81
      t84 = 0.1D1 + 0.260000000000000000000D0 * t82
      t89 = t78 * t84 - 0.322500000000000000000D0 * t79 * t81 * t78
      t90 = rhob * t89
      t91 = t60 * sigmabb
      t92 = t89 * t91
      t94 = 0.1D1 / t41 / rhob
      t96 = 0.1D1 / t80 / taub
      t97 = t94 * t96
      t100 = 0.1D1 + 0.437500000000000000000D-1 * t92 * t97
      zk(i) = t90 * t100
      t103 = 0.500000000000000000000D0 * t89 * t100
      t104 = 0.1D1 / t24
      t105 = t104 * t13
      t106 = t42 * t30
      t109 = t26 ** 2
      t111 = t18 / t109
      t112 = t19 ** 2
      t113 = t112 ** 2
      t131 = (-0.224298561906574129855D1 / t113 / t19 * t13 * t42 - 0.18
     >7699471636595866065D1 * t105 * t42 - 0.145760735710958868637D1 / t
     >22 * t13 * t42 - 0.344044309698575627326D0 / t16 * t13 * t42) / t2
     >9
      t134 = t55 * pi
      t137 = 0.1D1 / t38 / t37
      t147 = t39 * t42
      t155 = (-0.126105037207067989169D-1 * t104 * pi * t106 - 0.6517782
     >70654185890920D1 * t111 * t131 * t33) * t45
      t158 = 0.149583857013095144140D-1 * t35 * t58 * sigmabb * t147 * t
     >155
      t163 = 0.498612856710317147135D-2 * t134 * t47 * sigmabb * t137 * 
     >t42
      t166 = 0.299167714026190288281D-1 * t48 * t49 * t94
      t172 = t69 ** 2
      t173 = 0.1D1 / t172
      t182 = t55 ** 2
      t191 = 0.1D1 / t64 / rhob
      t202 = 0.1D1 / t74
      t205 = 0.193478431062909061652D-2 * t105 * t106 + 0.10000000000000
     >0000000D1 * t111 * t131 + 0.153426409720027345292D0 * t34 * (-0.49
     >8612856710317147135D-2 * t134 * sigmabb * t137 * t71 - 0.299167714
     >026190288281D-1 * t40 * t94 * t53 * t70 + 0.149583857013095144140D
     >-1 * t40 * t42 * (-t158 - t163 - t166) * t70 - 0.14958385701309514
     >4140D-1 * t40 * t54 * t173 * (-t158 - t163 - t166 - 0.447506605578
     >281866741D-3 * t56 / t57 / t46 * t60 * t62 * t65 * t155 - 0.149168
     >868526093955580D-3 * t182 * t58 * t60 / t61 / t37 * t65 - 0.895013
     >211156563733482D-3 * t59 * t63 * t191)) * t202
      t214 = t205 * t84 + 0.125000000000000000000D0 * t78 * t60 * t94 * 
     >t81 - 0.322500000000000000000D0 * t79 * t81 * t205
      t217 = 0.500000000000000000000D0 * rhob * t214 * t100
      t226 = 0.500000000000000000000D0 * t90 * (0.437500000000000000000D
     >-1 * t214 * t91 * t97 - 0.131250000000000000000D0 * t92 * t65 * t9
     >6)
      vrhoc(i) = vrhoc(i) + t103 + t217 + t226
      vrhoo(i) = vrhoo(i) - t103 - t217 - t226
      t250 = t34 * (0.149583857013095144140D-1 * t35 * t39 * t71 + 0.223
     >753302789140933370D-3 * t56 * sigmabb * t62 * t65 * t47 * t70 - 0.
     >149583857013095144140D-1 * t40 * t54 * t173 * (0.14958385701309514
     >4140D-1 * t48 * t147 + 0.447506605578281866741D-3 * t59 * sigmabb
     >* t62 * t65))
      t261 = 0.153426409720027345292D0 * t250 * t202 * t84 - 0.125000000
     >000000000000D0 * t78 * sigmabb * t42 * t81 - 0.4948001713470881885
     >65D-1 * t82 * t250 * t202
      t263 = rhob * t261 * t100
      t264 = 0.25D0 * t263
      t272 = t90 * (0.437500000000000000000D-1 * t261 * t91 * t97 + 0.13
     >1250000000000000000D0 * t89 * t60 * t97)
      t273 = 0.25D0 * t272
      vsigmacc(i) = vsigmacc(i) + t264 + t273
      vsigmaco(i) = vsigmaco(i) - 0.5D0 * t263 - 0.5D0 * t272
      vsigmaoo(i) = vsigmaoo(i) + t264 + t273
      t283 = 0.625000000000000000000D-1 * t14 * t78 * t60 * t96 * t100
      t284 = t60 ** 2
      t287 = t80 ** 2
      t299 = 0.500000000000000000000D0 * t90 * (0.546875000000000000000D
     >-2 * t78 * t284 * sigmabb * t191 / t287 / t80 - 0.1312500000000000
     >00000D0 * t92 * t94 / t287)
      vtauc(i) = vtauc(i) + t283 + t299
      vtauo(i) = vtauo(i) - t283 - t299

               elseif(rhob.lt.tol) then
             rho = rhoa
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) + 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmaaa
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i) + 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taua
      t13 = 0.1D1 / pi
      t14 = 0.1D1 / rhoa
      t15 = t13 * t14
      t16 = t15 ** (0.1D1 / 0.3D1)
      t18 = 0.1D1 + 0.186690969707574028554D0 * t16
      t19 = t15 ** (0.1D1 / 0.6D1)
      t22 = sqrt(t15)
      t24 = t16 ** 2
      t26 = 0.134579137143944477912D2 * t19 + 0.563098414909787598194D1 
     >* t16 + 0.291521471421917737271D1 * t22 + 0.516066464547863440989D
     >0 * t24
      t29 = 0.1D1 + 0.321639589973850701335D2 / t26
      t30 = log(t29)
      t31 = t18 * t30
      t33 = pi ** 2
      t34 = 0.1D1 / t33
      t35 = t33 * pi
      t37 = t33 * rhoa
      t38 = t37 ** (0.1D1 / 0.3D1)
      t39 = 0.1D1 / t38
      t40 = t35 * sigmaaa * t39
      t41 = rhoa ** 2
      t42 = 0.1D1 / t41
      t45 = exp(0.202642426794280972788D0 * t31 * t33)
      t46 = t45 - 0.1D1
      t47 = 0.1D1 / t46
      t48 = t35 * t47
      t49 = sigmaaa * t39
      t52 = 0.149583857013095144140D-1 * t48 * t49 * t42
      t53 = 0.1D1 + t52
      t54 = t42 * t53
      t55 = t33 ** 2
      t56 = t55 * t33
      t57 = t46 ** 2
      t58 = 0.1D1 / t57
      t59 = t56 * t58
      t60 = sigmaaa ** 2
      t61 = t38 ** 2
      t62 = 0.1D1 / t61
      t63 = t60 * t62
      t64 = t41 ** 2
      t65 = 0.1D1 / t64
      t69 = 0.1D1 + t52 + 0.223753302789140933370D-3 * t59 * t63 * t65
      t70 = 0.1D1 / t69
      t71 = t54 * t70
      t74 = 0.1D1 + 0.149583857013095144140D-1 * t40 * t71
      t75 = log(t74)
      t78 = -0.3109070D-1 * t31 + 0.153426409720027345292D0 * t34 * t75
      t79 = t60 * t42
      t80 = taua ** 2
      t81 = 0.1D1 / t80
      t82 = t79 * t81
      t84 = 0.1D1 + 0.260000000000000000000D0 * t82
      t89 = t78 * t84 - 0.322500000000000000000D0 * t79 * t81 * t78
      t90 = rhoa * t89
      t91 = t60 * sigmaaa
      t92 = t89 * t91
      t94 = 0.1D1 / t41 / rhoa
      t96 = 0.1D1 / t80 / taua
      t97 = t94 * t96
      t100 = 0.1D1 + 0.437500000000000000000D-1 * t92 * t97
      zk(i) = t90 * t100
      t103 = 0.500000000000000000000D0 * t89 * t100
      t104 = 0.1D1 / t24
      t105 = t104 * t13
      t106 = t42 * t30
      t109 = t26 ** 2
      t111 = t18 / t109
      t112 = t19 ** 2
      t113 = t112 ** 2
      t131 = (-0.224298561906574129855D1 / t113 / t19 * t13 * t42 - 0.18
     >7699471636595866065D1 * t105 * t42 - 0.145760735710958868637D1 / t
     >22 * t13 * t42 - 0.344044309698575627326D0 / t16 * t13 * t42) / t2
     >9
      t134 = t55 * pi
      t137 = 0.1D1 / t38 / t37
      t147 = t39 * t42
      t155 = (-0.126105037207067989169D-1 * t104 * pi * t106 - 0.6517782
     >70654185890920D1 * t111 * t131 * t33) * t45
      t158 = 0.149583857013095144140D-1 * t35 * t58 * sigmaaa * t147 * t
     >155
      t163 = 0.498612856710317147135D-2 * t134 * t47 * sigmaaa * t137 * 
     >t42
      t166 = 0.299167714026190288281D-1 * t48 * t49 * t94
      t172 = t69 ** 2
      t173 = 0.1D1 / t172
      t182 = t55 ** 2
      t191 = 0.1D1 / t64 / rhoa
      t202 = 0.1D1 / t74
      t205 = 0.193478431062909061652D-2 * t105 * t106 + 0.10000000000000
     >0000000D1 * t111 * t131 + 0.153426409720027345292D0 * t34 * (-0.49
     >8612856710317147135D-2 * t134 * sigmaaa * t137 * t71 - 0.299167714
     >026190288281D-1 * t40 * t94 * t53 * t70 + 0.149583857013095144140D
     >-1 * t40 * t42 * (-t158 - t163 - t166) * t70 - 0.14958385701309514
     >4140D-1 * t40 * t54 * t173 * (-t158 - t163 - t166 - 0.447506605578
     >281866741D-3 * t56 / t57 / t46 * t60 * t62 * t65 * t155 - 0.149168
     >868526093955580D-3 * t182 * t58 * t60 / t61 / t37 * t65 - 0.895013
     >211156563733482D-3 * t59 * t63 * t191)) * t202
      t214 = t205 * t84 + 0.125000000000000000000D0 * t78 * t60 * t94 * 
     >t81 - 0.322500000000000000000D0 * t79 * t81 * t205
      t217 = 0.500000000000000000000D0 * rhoa * t214 * t100
      t226 = 0.500000000000000000000D0 * t90 * (0.437500000000000000000D
     >-1 * t214 * t91 * t97 - 0.131250000000000000000D0 * t92 * t65 * t9
     >6)
      vrhoc(i) = vrhoc(i) + t103 + t217 + t226
      vrhoo(i) = vrhoo(i) + t103 + t217 + t226
      t250 = t34 * (0.149583857013095144140D-1 * t35 * t39 * t71 + 0.223
     >753302789140933370D-3 * t56 * sigmaaa * t62 * t65 * t47 * t70 - 0.
     >149583857013095144140D-1 * t40 * t54 * t173 * (0.14958385701309514
     >4140D-1 * t48 * t147 + 0.447506605578281866741D-3 * t59 * sigmaaa
     >* t62 * t65))
      t261 = 0.153426409720027345292D0 * t250 * t202 * t84 - 0.125000000
     >000000000000D0 * t78 * sigmaaa * t42 * t81 - 0.4948001713470881885
     >65D-1 * t82 * t250 * t202
      t263 = rhoa * t261 * t100
      t264 = 0.25D0 * t263
      t272 = t90 * (0.437500000000000000000D-1 * t261 * t91 * t97 + 0.13
     >1250000000000000000D0 * t89 * t60 * t97)
      t273 = 0.25D0 * t272
      vsigmacc(i) = vsigmacc(i) + t264 + t273
      vsigmaco(i) = vsigmaco(i) + 0.5D0 * t263 + 0.5D0 * t272
      vsigmaoo(i) = vsigmaoo(i) + t264 + t273
      t283 = 0.625000000000000000000D-1 * t14 * t78 * t60 * t96 * t100
      t284 = t60 ** 2
      t287 = t80 ** 2
      t299 = 0.500000000000000000000D0 * t90 * (0.546875000000000000000D
     >-2 * t78 * t284 * sigmaaa * t191 / t287 / t80 - 0.1312500000000000
     >00000D0 * t92 * t94 / t287)
      vtauc(i) = vtauc(i) + t283 + t299
      vtauo(i) = vtauo(i) + t283 + t299

               else
             rho = rhoa + rhob
      t2 = 0.250000000000000000000D0 * sigmacc(i)
      t4 = 0.250000000000000000000D0 * sigmaoo(i)
      t6 = 0.500000000000000000000D0 * sigmaco(i)
      sigmaaa = max(0.0D0, t2 + t4 + t6)
      sigmaab = t2 - t4
      sigmabb = max(0.0D0, t2 + t4 - t6)
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t11 = 0.500000000000000000000D0 * tauc(i)
      t13 = 0.500000000000000000000D0 * tauo(i)
      taua = max(0.0D0, t11 + t13)
      taub = max(0.0D0, t11 - t13)
      tau = taua + taub
      t16 = 0.1D1 / pi
      t17 = 0.1D1 / rho
      t18 = t16 * t17
      t19 = t18 ** (0.1D1 / 0.3D1)
      t21 = 0.1D1 + 0.194159335344114122552D0 * t19
      t22 = t18 ** (0.1D1 / 0.6D1)
      t25 = sqrt(t18)
      t27 = t19 ** 2
      t29 = 0.724010193431683113327D1 * t22 + 0.325955091942229212011D1 
     >* t19 + 0.141872281647966739112D1 * t25 + 0.406913004517529319387D
     >0 * t27
      t32 = 0.1D1 + 0.160819794986925350668D2 / t29
      t33 = log(t32)
      t34 = t21 * t33
      t35 = 0.621814D-1 * t34
      t37 = 0.1D1 + 0.101077332976287768525D0 * t19
      t42 = 0.987212972256927209438D1 * t22 + 0.329180480994506259905D1 
     >* t19 + 0.762327521935289963194D0 * t25 + 0.410025070949612505036D
     >0 * t27
      t45 = 0.1D1 + 0.296087499777934375166D2 / t42
      t46 = log(t45)
      t47 = t37 * t46
      t49 = rhoa - 0.1D1 * rhob
      t50 = t49 * t17
      t51 = 0.1D1 + t50
      t52 = t51 ** (0.1D1 / 0.3D1)
      t53 = t52 * t51
      t55 = 0.1D1 - 0.1D1 * t50
      t56 = t55 ** (0.1D1 / 0.3D1)
      t57 = t56 * t55
      t58 = t53 + t57 - 0.2D1
      t59 = t49 ** 2
      t60 = t59 ** 2
      t61 = rho ** 2
      t62 = t61 ** 2
      t63 = 0.1D1 / t62
      t64 = t60 * t63
      t66 = 0.1D1 - 0.1D1 * t64
      t68 = t47 * t58 * t66
      t69 = 0.379955235370239451738D-1 * t68
      t71 = 0.1D1 + 0.186690969707574028554D0 * t19
      t76 = 0.134579137143944477912D2 * t22 + 0.563098414909787598194D1 
     >* t19 + 0.291521471421917737271D1 * t25 + 0.516066464547863440989D
     >0 * t27
      t79 = 0.1D1 + 0.321639589973850701335D2 / t76
      t80 = log(t79)
      t83 = -0.3109070D-1 * t71 * t80 + t35
      t84 = t83 * t58
      t85 = t84 * t64
      t86 = 0.192366105093153631974D1 * t85
      t87 = pi ** 2
      t88 = 0.1D1 / t87
      t89 = t52 ** 2
      t91 = t56 ** 2
      t93 = 0.500000000000000000000D0 * t89 + 0.500000000000000000000D0 
     >* t91
      t94 = t93 ** 2
      t95 = t94 * t93
      t96 = t88 * t95
      t97 = t87 * pi
      t98 = t97 * sigma
      t99 = 0.1D1 / t94
      t100 = t98 * t99
      t101 = t87 * rho
      t102 = t101 ** (0.1D1 / 0.3D1)
      t103 = 0.1D1 / t102
      t104 = 0.1D1 / t61
      t105 = t103 * t104
      t107 = (-t35 + t69 + t86) * t87
      t108 = 0.1D1 / t95
      t111 = exp(-0.325889135327092945460D1 * t107 * t108)
      t112 = t111 - 0.1D1
      t113 = 0.1D1 / t112
      t114 = t97 * t113
      t115 = t114 * sigma
      t116 = t99 * t103
      t117 = t116 * t104
      t119 = 0.942319250876317101329D-2 * t115 * t117
      t120 = 0.1D1 + t119
      t121 = t87 ** 2
      t122 = t121 * t87
      t123 = t112 ** 2
      t124 = 0.1D1 / t123
      t125 = t122 * t124
      t126 = sigma ** 2
      t127 = t125 * t126
      t128 = t94 ** 2
      t129 = 0.1D1 / t128
      t130 = t102 ** 2
      t131 = 0.1D1 / t130
      t132 = t129 * t131
      t133 = t132 * t63
      t136 = 0.1D1 + t119 + 0.887965570572103448132D-4 * t127 * t133
      t137 = 0.1D1 / t136
      t138 = t120 * t137
      t142 = 0.1D1 + 0.942319250876317101329D-2 * t100 * t105 * t138
      t143 = log(t142)
      t144 = t96 * t143
      t146 = -t35 + t69 + t86 + 0.306852819440054690583D0 * t144
      t150 = t60 * t59
      t152 = 0.1D1 / t62 / t61
      t155 = 0.53D0 + 0.87D0 * t59 * t104 + 0.50D0 * t64 + 0.226D1 * t15
     >0 * t152
      t156 = t55 ** 2
      t158 = t51 ** 2
      t160 = t55 * t51
      t163 = t156 * sigmaaa + t158 * sigmabb - 0.2D1 * t160 * sigmaab
      t164 = t163 * t104
      t167 = 0.1D1 / t53 + 0.1D1 / t57
      t168 = t131 * t167
      t171 = 0.1D1 + 0.151426716069344971574D0 * t164 * t168
      t172 = t171 ** 2
      t173 = t172 ** 2
      t174 = 0.1D1 / t173
      t175 = t155 * t174
      t176 = t126 * t104
      t177 = tau ** 2
      t178 = 0.1D1 / t177
      t179 = t176 * t178
      t182 = 0.1D1 + 0.625000000000000000000D-1 * t175 * t179
      t184 = 0.1D1 + t175
      t185 = t184 * t126
      t186 = t104 * t178
      t187 = rhoa * t17
      t188 = 0.621814D29 * t34
      t189 = 0.379955235370239451738D29 * t68
      t190 = 0.192366105093153631974D31 * t85
      t191 = 0.306852819440054690583D30 * t144
      t193 = t16 / rhoa
      t194 = t193 ** (0.1D1 / 0.3D1)
      t196 = 0.1D1 + 0.186690969707574028554D0 * t194
      t197 = t193 ** (0.1D1 / 0.6D1)
      t200 = sqrt(t193)
      t202 = t194 ** 2
      t204 = 0.134579137143944477912D2 * t197 + 0.563098414909787598194D
     >1 * t194 + 0.291521471421917737271D1 * t200 + 0.516066464547863440
     >989D0 * t202
      t207 = 0.1D1 + 0.321639589973850701335D2 / t204
      t208 = log(t207)
      t209 = t196 * t208
      t212 = t87 * rhoa
      t213 = t212 ** (0.1D1 / 0.3D1)
      t214 = 0.1D1 / t213
      t215 = t97 * sigmaaa * t214
      t216 = rhoa ** 2
      t217 = 0.1D1 / t216
      t220 = exp(0.202642426794280972788D0 * t209 * t87)
      t221 = t220 - 0.1D1
      t222 = 0.1D1 / t221
      t223 = t97 * t222
      t224 = sigmaaa * t214
      t227 = 0.149583857013095144140D-1 * t223 * t224 * t217
      t228 = 0.1D1 + t227
      t229 = t217 * t228
      t230 = t221 ** 2
      t231 = 0.1D1 / t230
      t232 = t122 * t231
      t233 = sigmaaa ** 2
      t234 = t213 ** 2
      t235 = 0.1D1 / t234
      t236 = t233 * t235
      t237 = t216 ** 2
      t238 = 0.1D1 / t237
      t242 = 0.1D1 + t227 + 0.223753302789140933370D-3 * t232 * t236 * t
     >238
      t243 = 0.1D1 / t242
      t244 = t229 * t243
      t247 = 0.1D1 + 0.149583857013095144140D-1 * t215 * t244
      t248 = log(t247)
      t249 = t88 * t248
      t251 = t188 - t189 - t190 - t191 - 0.3109070D29 * t209 + 0.1534264
     >09720027345292D30 * t249
      t252 = tanh(t251)
      t253 = 0.500000000000000000000D0 * t252
      t254 = 0.500000000000000000000D0 - t253
      t256 = 0.500000000000000000000D0 + t253
      t259 = -0.3109070D-1 * t209 + 0.153426409720027345292D0 * t249
      t261 = t254 * t146 + t256 * t259
      t263 = rhob * t17
      t265 = t16 / rhob
      t266 = t265 ** (0.1D1 / 0.3D1)
      t268 = 0.1D1 + 0.186690969707574028554D0 * t266
      t269 = t265 ** (0.1D1 / 0.6D1)
      t272 = sqrt(t265)
      t274 = t266 ** 2
      t276 = 0.134579137143944477912D2 * t269 + 0.563098414909787598194D
     >1 * t266 + 0.291521471421917737271D1 * t272 + 0.516066464547863440
     >989D0 * t274
      t279 = 0.1D1 + 0.321639589973850701335D2 / t276
      t280 = log(t279)
      t281 = t268 * t280
      t284 = t87 * rhob
      t285 = t284 ** (0.1D1 / 0.3D1)
      t286 = 0.1D1 / t285
      t287 = t97 * sigmabb * t286
      t288 = rhob ** 2
      t289 = 0.1D1 / t288
      t292 = exp(0.202642426794280972788D0 * t281 * t87)
      t293 = t292 - 0.1D1
      t294 = 0.1D1 / t293
      t295 = t97 * t294
      t296 = sigmabb * t286
      t299 = 0.149583857013095144140D-1 * t295 * t296 * t289
      t300 = 0.1D1 + t299
      t301 = t289 * t300
      t302 = t293 ** 2
      t303 = 0.1D1 / t302
      t304 = t122 * t303
      t305 = sigmabb ** 2
      t306 = t285 ** 2
      t307 = 0.1D1 / t306
      t308 = t305 * t307
      t309 = t288 ** 2
      t310 = 0.1D1 / t309
      t314 = 0.1D1 + t299 + 0.223753302789140933370D-3 * t304 * t308 * t
     >310
      t315 = 0.1D1 / t314
      t316 = t301 * t315
      t319 = 0.1D1 + 0.149583857013095144140D-1 * t287 * t316
      t320 = log(t319)
      t321 = t88 * t320
      t323 = t188 - t189 - t190 - t191 - 0.3109070D29 * t281 + 0.1534264
     >09720027345292D30 * t321
      t324 = tanh(t323)
      t325 = 0.500000000000000000000D0 * t324
      t326 = 0.500000000000000000000D0 - t325
      t328 = 0.500000000000000000000D0 + t325
      t331 = -0.3109070D-1 * t281 + 0.153426409720027345292D0 * t321
      t333 = t326 * t146 + t328 * t331
      t335 = t187 * t261 + t263 * t333
      t336 = t186 * t335
      t339 = t146 * t182 - 0.625000000000000000000D-1 * t185 * t336
      t340 = rho * t339
      t341 = t126 * sigma
      t342 = t339 * t341
      t343 = t61 * rho
      t344 = 0.1D1 / t343
      t346 = 0.1D1 / t177 / tau
      t347 = t344 * t346
      t350 = 0.1D1 + 0.437500000000000000000D-1 * t342 * t347
      zk(i) = t340 * t350
      t353 = 0.1D1 / t27 * t16
      t355 = t353 * t104 * t33
      t356 = 0.402436643158883263334D-2 * t355
      t357 = t29 ** 2
      t360 = t22 ** 2
      t361 = t360 ** 2
      t365 = 0.1D1 / t361 / t22 * t16 * t104
      t367 = t353 * t104
      t371 = 0.1D1 / t25 * t16 * t104
      t375 = 0.1D1 / t19 * t16 * t104
      t380 = t21 / t357 * (-0.120668365571947185555D1 * t365 - 0.1086516
     >97314076404004D1 * t367 - 0.709361408239833695563D0 * t371 - 0.271
     >275336345019546258D0 * t375) / t32
      t381 = 0.100000000000000000000D1 * t380
      t384 = t367 * t46 * t58 * t66
      t385 = 0.128016206138671616206D-2 * t384
      t386 = t42 ** 2
      t398 = t37 / t386 * (-0.164535495376154534908D1 * t365 - 0.1097268
     >26998168753302D1 * t367 - 0.381163760967644981600D0 * t371 - 0.273
     >350047299741670024D0 * t375) / t45 * t58 * t66
      t399 = 0.112499995668310776915D1 * t398
      t400 = t49 * t104
      t401 = 0.1D1 * t400
      t402 = t17 - t401
      t405 = 0.1D1 * t17
      t406 = -t405 + t400
      t409 = 0.133333333333333333333D1 * t52 * t402 + 0.1333333333333333
     >33333D1 * t56 * t406
      t411 = t47 * t409 * t66
      t412 = 0.379955235370239451738D-1 * t411
      t414 = t59 * t49 * t63
      t415 = 0.4D1 * t414
      t417 = 0.1D1 / t62 / rho
      t418 = t60 * t417
      t419 = 0.4D1 * t418
      t422 = t47 * t58 * (-t415 + t419)
      t423 = 0.379955235370239451738D-1 * t422
      t427 = t76 ** 2
      t441 = (0.193478431062909061652D-2 * t353 * t104 * t80 + 0.1000000
     >00000000000000D1 * t71 / t427 * (-0.224298561906574129855D1 * t365
     > - 0.187699471636595866065D1 * t367 - 0.145760735710958868637D1 *
     >t371 - 0.344044309698575627326D0 * t375) / t79 - t356 - t381) * t5
     >8 * t64
      t442 = 0.192366105093153631974D1 * t441
      t444 = t83 * t409 * t64
      t445 = 0.192366105093153631974D1 * t444
      t446 = t84 * t414
      t447 = 0.769464420372614527896D1 * t446
      t448 = t84 * t418
      t449 = 0.769464420372614527896D1 * t448
      t450 = t88 * t94
      t451 = 0.1D1 / t52
      t454 = 0.1D1 / t56
      t457 = 0.333333333333333333333D0 * t451 * t402 + 0.333333333333333
     >333333D0 * t454 * t406
      t459 = t450 * t143 * t457
      t461 = t108 * t103
      t462 = t98 * t461
      t463 = t104 * t120
      t468 = t121 * pi
      t472 = 0.1D1 / t102 / t101
      t476 = 0.314106416958772367110D-2 * t468 * sigma * t99 * t472 * t1
     >04 * t138
      t480 = 0.188463850175263420266D-1 * t100 * t103 * t344 * t138
      t483 = t97 * t124 * sigma * t99
      t492 = (-0.325889135327092945460D1 * (t356 + t381 - t385 - t399 + 
     >t412 + t423 + t442 + t445 + t447 - t449) * t87 * t108 + 0.97766740
     >5981278836380D1 * t107 * t129 * t457) * t111
      t495 = 0.942319250876317101329D-2 * t483 * t105 * t492
      t499 = 0.188463850175263420266D-1 * t115 * t461 * t104 * t457
      t505 = 0.314106416958772367110D-2 * t468 * t113 * sigma * t99 * t4
     >72 * t104
      t508 = 0.188463850175263420266D-1 * t115 * t116 * t344
      t514 = t98 * t116
      t515 = t136 ** 2
      t516 = 0.1D1 / t515
      t521 = t122 / t123 / t112 * t126 * t129
      t522 = t131 * t63
      t528 = 0.1D1 / t128 / t93 * t131
      t533 = t121 ** 2
      t537 = 0.1D1 / t130 / t101
      t541 = 0.591977047048068965421D-4 * t533 * t124 * t126 * t129 * t5
     >37 * t63
      t544 = 0.355186228228841379252D-3 * t127 * t132 * t417
      t551 = 0.1D1 / t142
      t553 = t96 * (-0.188463850175263420266D-1 * t462 * t463 * t137 * t
     >457 - t476 - t480 + 0.942319250876317101329D-2 * t100 * t105 * (-t
     >495 - t499 - t505 - t508) * t137 - 0.942319250876317101329D-2 * t5
     >14 * t463 * t516 * (-t495 - t499 - t505 - t508 - 0.177593114114420
     >689627D-3 * t521 * t522 * t492 - 0.355186228228841379252D-3 * t127
     > * t528 * t63 * t457 - t541 - t544)) * t551
      t555 = t356 + t381 - t385 - t399 + t412 + t423 + t442 + t445 + t44
     >7 - t449 + 0.920558458320164071749D0 * t459 + 0.306852819440054690
     >583D0 * t553
      t562 = 0.174D1 * t400 + 0.200D1 * t414 + 0.1356D2 * t60 * t49 * t1
     >52
      t563 = t562 * t174
      t567 = 0.1D1 / t173 / t171
      t568 = t155 * t567
      t569 = t568 * t126
      t570 = t55 * sigmaaa
      t573 = t51 * sigmabb
      t582 = -0.2D1 * t570 * t17 + 0.2D1 * t573 * t17 + 0.2D1 * t17 * t5
     >1 * sigmaab - 0.2D1 * t55 * t17 * sigmaab
      t587 = 0.1D1 / t52 / t158
      t591 = 0.1D1 / t56 / t156
      t594 = -0.133333333333333333333D1 * t587 * t17 + 0.133333333333333
     >333333D1 * t591 * t17
      t598 = 0.151426716069344971574D0 * t582 * t104 * t168 + 0.15142671
     >6069344971574D0 * t164 * t131 * t594
      t613 = 0.1D1 * rhoa * t104 * t261
      t615 = tanh(-t251)
      t616 = t615 ** 2
      t618 = 0.1D1 - 0.1D1 * t616
      t619 = 0.402436643158883263334D28 * t355
      t620 = 0.100000000000000000000D31 * t380
      t621 = 0.128016206138671616206D28 * t384
      t622 = 0.112499995668310776915D31 * t398
      t623 = 0.379955235370239451738D29 * t411
      t624 = 0.379955235370239451738D29 * t422
      t625 = 0.192366105093153631974D31 * t441
      t626 = 0.192366105093153631974D31 * t444
      t627 = 0.769464420372614527896D31 * t446
      t628 = 0.769464420372614527896D31 * t448
      t629 = 0.920558458320164071749D30 * t459
      t630 = 0.306852819440054690583D30 * t553
      t631 = 0.1D1 / t202
      t632 = t631 * t16
      t633 = t217 * t208
      t634 = t632 * t633
      t636 = t204 ** 2
      t638 = t196 / t636
      t639 = t197 ** 2
      t640 = t639 ** 2
      t658 = (-0.224298561906574129855D1 / t640 / t197 * t16 * t217 - 0.
     >187699471636595866065D1 * t632 * t217 - 0.145760735710958868637D1
     >/ t200 * t16 * t217 - 0.344044309698575627326D0 / t194 * t16 * t21
     >7) / t207
      t659 = t638 * t658
      t663 = 0.1D1 / t213 / t212
      t668 = 0.1D1 / t216 / rhoa
      t675 = t214 * t217
      t683 = (-0.126105037207067989169D-1 * t631 * pi * t633 - 0.6517782
     >70654185890920D1 * t638 * t658 * t87) * t220
      t686 = 0.149583857013095144140D-1 * t97 * t231 * sigmaaa * t675 * 
     >t683
      t691 = 0.498612856710317147135D-2 * t468 * t222 * sigmaaa * t663 *
     > t217
      t694 = 0.299167714026190288281D-1 * t223 * t224 * t668
      t700 = t242 ** 2
      t701 = 0.1D1 / t700
      t729 = 0.1D1 / t247
      t730 = t88 * (-0.498612856710317147135D-2 * t468 * sigmaaa * t663 
     >* t244 - 0.299167714026190288281D-1 * t215 * t668 * t228 * t243 +
     >0.149583857013095144140D-1 * t215 * t217 * (-t686 - t691 - t694) *
     > t243 - 0.149583857013095144140D-1 * t215 * t229 * t701 * (-t686 -
     > t691 - t694 - 0.447506605578281866741D-3 * t122 / t230 / t221 * t
     >233 * t235 * t238 * t683 - 0.149168868526093955580D-3 * t533 * t23
     >1 * t233 / t234 / t212 * t238 - 0.895013211156563733482D-3 * t232
     >* t236 / t237 / rhoa)) * t729
      t732 = -t619 - t620 + t621 + t622 - t623 - t624 - t625 - t626 - t6
     >27 + t628 - t629 - t630 + 0.193478431062909061652D28 * t634 + 0.10
     >0000000000000000000D31 * t659 + 0.153426409720027345292D30 * t730
      t733 = t618 * t732
      t748 = 0.1D1 * rhob * t104 * t333
      t750 = tanh(-t323)
      t751 = t750 ** 2
      t753 = 0.1D1 - 0.1D1 * t751
      t754 = -t619 - t620 + t621 + t622 - t623 - t624 - t625 - t626 - t6
     >27 + t628 - t629 - t630
      t755 = t753 * t754
      t767 = t555 * t182 + t146 * (0.625000000000000000000D-1 * t563 * t
     >179 - 0.250000000000000000000D0 * t569 * t186 * t598) - 0.62500000
     >0000000000000D-1 * (t563 - 0.4D1 * t568 * t598) * t126 * t336 - 0.
     >625000000000000000000D-1 * t185 * t186 * (t17 * t261 - t613 + t187
     > * (-0.500000000000000000000D0 * t733 * t146 + t254 * t555 + 0.500
     >000000000000000000D0 * t733 * t259 + t256 * (0.1934784310629090616
     >52D-2 * t634 + 0.100000000000000000000D1 * t659 + 0.15342640972002
     >7345292D0 * t730)) - t748 + t263 * (-0.500000000000000000000D0 * t
     >755 * t146 + t326 * t555 + 0.500000000000000000000D0 * t755 * t331
     >))
      t770 = 0.500000000000000000000D0 * rho * t767 * t350
      t771 = t104 * t339
      t775 = 0.218750000000000000000D-1 * t771 * t767 * t341 * t346
      t776 = -t405 - t401
      t779 = t17 + t400
      t782 = 0.133333333333333333333D1 * t52 * t776 + 0.1333333333333333
     >33333D1 * t56 * t779
      t784 = t47 * t782 * t66
      t785 = 0.379955235370239451738D-1 * t784
      t788 = t47 * t58 * (t415 + t419)
      t789 = 0.379955235370239451738D-1 * t788
      t791 = t83 * t782 * t64
      t792 = 0.192366105093153631974D1 * t791
      t797 = 0.333333333333333333333D0 * t451 * t776 + 0.333333333333333
     >333333D0 * t454 * t779
      t799 = t450 * t143 * t797
      t813 = (-0.325889135327092945460D1 * (t356 + t381 - t385 - t399 + 
     >t785 + t789 + t442 + t792 - t447 - t449) * t87 * t108 + 0.97766740
     >5981278836380D1 * t107 * t129 * t797) * t111
      t816 = 0.942319250876317101329D-2 * t483 * t105 * t813
      t820 = 0.188463850175263420266D-1 * t115 * t461 * t104 * t797
      t840 = t96 * (-0.188463850175263420266D-1 * t462 * t463 * t137 * t
     >797 - t476 - t480 + 0.942319250876317101329D-2 * t100 * t105 * (-t
     >816 - t820 - t505 - t508) * t137 - 0.942319250876317101329D-2 * t5
     >14 * t463 * t516 * (-t816 - t820 - t505 - t508 - 0.177593114114420
     >689627D-3 * t521 * t522 * t813 - 0.355186228228841379252D-3 * t127
     > * t528 * t63 * t797 - t541 - t544)) * t551
      t842 = t356 + t381 - t385 - t399 + t785 + t789 + t442 + t792 - t44
     >7 - t449 + 0.920558458320164071749D0 * t799 + 0.306852819440054690
     >583D0 * t840
      t845 = -t562 * t174
      t856 = -0.151426716069344971574D0 * t582 * t104 * t168 - 0.1514267
     >16069344971574D0 * t164 * t131 * t594
      t868 = 0.379955235370239451738D29 * t784
      t869 = 0.379955235370239451738D29 * t788
      t870 = 0.192366105093153631974D31 * t791
      t871 = 0.920558458320164071749D30 * t799
      t872 = 0.306852819440054690583D30 * t840
      t873 = -t619 - t620 + t621 + t622 - t868 - t869 - t625 - t870 + t6
     >27 + t628 - t871 - t872
      t874 = t618 * t873
      t883 = 0.1D1 / t274
      t884 = t883 * t16
      t885 = t289 * t280
      t886 = t884 * t885
      t888 = t276 ** 2
      t890 = t268 / t888
      t891 = t269 ** 2
      t892 = t891 ** 2
      t910 = (-0.224298561906574129855D1 / t892 / t269 * t16 * t289 - 0.
     >187699471636595866065D1 * t884 * t289 - 0.145760735710958868637D1
     >/ t272 * t16 * t289 - 0.344044309698575627326D0 / t266 * t16 * t28
     >9) / t279
      t911 = t890 * t910
      t915 = 0.1D1 / t285 / t284
      t920 = 0.1D1 / t288 / rhob
      t927 = t286 * t289
      t935 = (-0.126105037207067989169D-1 * t883 * pi * t885 - 0.6517782
     >70654185890920D1 * t890 * t910 * t87) * t292
      t938 = 0.149583857013095144140D-1 * t97 * t303 * sigmabb * t927 * 
     >t935
      t943 = 0.498612856710317147135D-2 * t468 * t294 * sigmabb * t915 *
     > t289
      t946 = 0.299167714026190288281D-1 * t295 * t296 * t920
      t952 = t314 ** 2
      t953 = 0.1D1 / t952
      t981 = 0.1D1 / t319
      t982 = t88 * (-0.498612856710317147135D-2 * t468 * sigmabb * t915 
     >* t316 - 0.299167714026190288281D-1 * t287 * t920 * t300 * t315 +
     >0.149583857013095144140D-1 * t287 * t289 * (-t938 - t943 - t946) *
     > t315 - 0.149583857013095144140D-1 * t287 * t301 * t953 * (-t938 -
     > t943 - t946 - 0.447506605578281866741D-3 * t122 / t302 / t293 * t
     >305 * t307 * t310 * t935 - 0.149168868526093955580D-3 * t533 * t30
     >3 * t305 / t306 / t284 * t310 - 0.895013211156563733482D-3 * t304
     >* t308 / t309 / rhob)) * t981
      t984 = -t619 - t620 + t621 + t622 - t868 - t869 - t625 - t870 + t6
     >27 + t628 - t871 - t872 + 0.193478431062909061652D28 * t886 + 0.10
     >0000000000000000000D31 * t911 + 0.153426409720027345292D30 * t982
      t985 = t753 * t984
      t1002 = t842 * t182 + t146 * (0.625000000000000000000D-1 * t845 * 
     >t179 - 0.250000000000000000000D0 * t569 * t186 * t856) - 0.6250000
     >00000000000000D-1 * (t845 - 0.4D1 * t568 * t856) * t126 * t336 - 0
     >.625000000000000000000D-1 * t185 * t186 * (-t613 + t187 * (-0.5000
     >00000000000000000D0 * t874 * t146 + t254 * t842 + 0.50000000000000
     >0000000D0 * t874 * t259) + t17 * t333 - t748 + t263 * (-0.50000000
     >0000000000000D0 * t985 * t146 + t326 * t842 + 0.500000000000000000
     >000D0 * t985 * t331 + t328 * (0.193478431062909061652D-2 * t886 +
     >0.100000000000000000000D1 * t911 + 0.153426409720027345292D0 * t98
     >2)))
      t1005 = 0.500000000000000000000D0 * rho * t1002 * t350
      t1009 = 0.218750000000000000000D-1 * t771 * t1002 * t341 * t346
      t1019 = (-0.174D1 * t59 * t344 - 0.200D1 * t418 - 0.1356D2 * t150 
     >/ t62 / t343) * t174
      t1054 = 0.151426716069344971574D0 * (0.2D1 * t570 * t400 - 0.2D1 *
     > t573 * t400 - 0.2D1 * t400 * t51 * sigmaab + 0.2D1 * t55 * t49 *
     >t104 * sigmaab) * t104 * t168 - 0.302853432138689943149D0 * t163 *
     > t344 * t168 - 0.100951144046229981050D0 * t164 * t537 * t167 * t8
     >7 + 0.151426716069344971574D0 * t164 * t131 * (0.13333333333333333
     >3333D1 * t587 * t49 * t104 - 0.133333333333333333333D1 * t591 * t4
     >9 * t104)
      t1074 = t146 * (0.625000000000000000000D-1 * t1019 * t179 - 0.2500
     >00000000000000000D0 * t569 * t186 * t1054 - 0.12500000000000000000
     >0D0 * t175 * t126 * t344 * t178) - 0.625000000000000000000D-1 * (t
     >1019 - 0.4D1 * t568 * t1054) * t126 * t336 + 0.1250000000000000000
     >00D0 * t185 * t344 * t178 * t335
      vrhoc(i) = vrhoc(i) + t770 + t775 + t1005 + t1009 + t339 * t350 + 
     >rho * t1074 * t350 + t340 * (0.437500000000000000000D-1 * t1074 *
     >t341 * t347 - 0.131250000000000000000D0 * t342 * t63 * t346)
      vrhoo(i) = vrhoo(i) + t770 + t775 - t1005 - t1009
      t1087 = t146 * t155
      t1088 = t567 * t126
      t1089 = t1087 * t1088
      t1090 = t63 * t178
      t1100 = t168 * t126 * t178 * t335
      t1103 = t185 * t104
      t1105 = t618 * t88
      t1127 = (0.149583857013095144140D-1 * t97 * t214 * t244 + 0.223753
     >302789140933370D-3 * t122 * sigmaaa * t235 * t238 * t222 * t243 -
     >0.149583857013095144140D-1 * t215 * t229 * t701 * (0.1495838570130
     >95144140D-1 * t223 * t675 + 0.447506605578281866741D-3 * t232 * si
     >gmaaa * t235 * t238)) * t729
      t1142 = -0.378566790173362428935D-1 * t1089 * t1090 * t156 * t131 
     >* t167 + 0.378566790173362428935D-1 * t568 * t156 * t63 * t1100 -
     >0.625000000000000000000D-1 * t1103 * t178 * rhoa * t17 * (-0.76713
     >2048600136726458D29 * t1105 * t1127 * t146 + 0.7671320486001367264
     >58D29 * t1105 * t1127 * t259 + 0.153426409720027345292D0 * t256 *
     >t88 * t1127)
      t1144 = rho * t1142 * t350
      t1145 = 0.25D0 * t1144
      t1148 = t771 * t1142 * t341 * t346
      t1149 = 0.109375000000000000000D-1 * t1148
      t1160 = t753 * t88
      t1182 = (0.149583857013095144140D-1 * t97 * t286 * t316 + 0.223753
     >302789140933370D-3 * t122 * sigmabb * t307 * t310 * t294 * t315 -
     >0.149583857013095144140D-1 * t287 * t301 * t953 * (0.1495838570130
     >95144140D-1 * t295 * t927 + 0.447506605578281866741D-3 * t304 * si
     >gmabb * t307 * t310)) * t981
      t1197 = -0.378566790173362428935D-1 * t1089 * t1090 * t158 * t131 
     >* t167 + 0.378566790173362428935D-1 * t568 * t158 * t63 * t1100 -
     >0.625000000000000000000D-1 * t1103 * t178 * rhob * t17 * (-0.76713
     >2048600136726458D29 * t1160 * t1182 * t146 + 0.7671320486001367264
     >58D29 * t1160 * t1182 * t331 + 0.153426409720027345292D0 * t328 *
     >t88 * t1182)
      t1199 = rho * t1197 * t350
      t1200 = 0.25D0 * t1199
      t1203 = t771 * t1197 * t341 * t346
      t1204 = 0.109375000000000000000D-1 * t1203
      t1217 = 0.757133580346724857871D-1 * t1087 * t1088 * t63 * t178 * 
     >t55 * t51 * t131 * t167 - 0.757133580346724857871D-1 * t568 * t160
     > * t63 * t1100
      t1220 = 0.25D0 * rho * t1217 * t350
      t1224 = 0.109375000000000000000D-1 * t771 * t1217 * t341 * t346
      t1246 = 0.942319250876317101329D-2 * t97 * t99 * t103 * t463 * t13
     >7 + 0.887965570572103448136D-4 * t122 * sigma * t129 * t522 * t113
     > * t137 - 0.942319250876317101329D-2 * t514 * t463 * t516 * (0.942
     >319250876317101329D-2 * t114 * t117 + 0.177593114114420689627D-3 *
     > t125 * sigma * t133)
      t1247 = t1246 * t551
      t1251 = t1087 * t174
      t1259 = t1105 * t95
      t1260 = t1247 * t146
      t1265 = t95 * t1246 * t551
      t1273 = t1160 * t95
      t1288 = 0.306852819440054690583D0 * t96 * t1247 * t182 + 0.1250000
     >00000000000000D0 * t1251 * sigma * t104 * t178 - 0.125000000000000
     >000000D0 * t184 * sigma * t336 - 0.625000000000000000000D-1 * t185
     > * t186 * (t187 * (0.153426409720027345292D30 * t1259 * t1260 + 0.
     >306852819440054690583D0 * t254 * t88 * t1265 - 0.15342640972002734
     >5292D30 * t1259 * t1247 * t259) + t263 * (0.153426409720027345292D
     >30 * t1273 * t1260 + 0.306852819440054690583D0 * t326 * t88 * t126
     >5 - 0.153426409720027345292D30 * t1273 * t1247 * t331))
      vsigmacc(i) = vsigmacc(i) + t1145 + t1149 + t1200 + t1204 + t1220 
     >+ t1224 + rho * t1288 * t350 + t340 * (0.437500000000000000000D-1
     >* t1288 * t341 * t347 + 0.131250000000000000000D0 * t339 * t126 *
     >t347)
      vsigmaco(i) = vsigmaco(i) + 0.5D0 * t1144 + 0.21875000000000000000
     >0D-1 * t1148 - 0.5D0 * t1199 - 0.218750000000000000000D-1 * t1203
      vsigmaoo(i) = vsigmaoo(i) + t1145 + t1149 + t1200 + t1204 - t1220 
     >- t1224
      t1313 = -0.125000000000000000000D0 * t1251 * t176 * t346 + 0.12500
     >0000000000000000D0 * t185 * t104 * t346 * t335
      t1319 = t177 ** 2
      vtauc(i) = vtauc(i) + rho * t1313 * t350 + t340 * (0.4375000000000
     >00000000D-1 * t1313 * t341 * t347 - 0.131250000000000000000D0 * t3
     >42 * t344 / t1319)
      vtauo(i) = vtauo(i)

               endif
             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      sigmabb = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) - 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmabb
      taub = max(0.0D0, 0.500000000000000000000D0 * tauc(i) - 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taub
      t15 = 0.1D1 / pi / rhob
      t16 = t15 ** (0.1D1 / 0.3D1)
      t19 = t15 ** (0.1D1 / 0.6D1)
      t22 = sqrt(t15)
      t24 = t16 ** 2
      t30 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t19 + 0.563098414909787598194D1 * t16 + 0.291521471421917
     >737271D1 * t22 + 0.516066464547863440989D0 * t24))
      t31 = (0.1D1 + 0.186690969707574028554D0 * t16) * t30
      t33 = pi ** 2
      t35 = t33 * pi
      t38 = (t33 * rhob) ** (0.1D1 / 0.3D1)
      t39 = 0.1D1 / t38
      t41 = rhob ** 2
      t42 = 0.1D1 / t41
      t45 = exp(0.202642426794280972788D0 * t31 * t33)
      t46 = t45 - 0.1D1
      t52 = 0.149583857013095144140D-1 * t35 / t46 * sigmabb * t39 * t42
      t55 = t33 ** 2
      t57 = t46 ** 2
      t60 = sigmabb ** 2
      t61 = t38 ** 2
      t64 = t41 ** 2
      t75 = log(0.1D1 + 0.149583857013095144140D-1 * t35 * sigmabb * t39
     > * t42 * (0.1D1 + t52) / (0.1D1 + t52 + 0.223753302789140933370D-3
     > * t55 * t33 / t57 * t60 / t61 / t64))
      t78 = -0.3109070D-1 * t31 + 0.153426409720027345292D0 / t33 * t75
      t79 = t60 * t42
      t80 = taub ** 2
      t81 = 0.1D1 / t80
      t89 = t78 * (0.1D1 + 0.260000000000000000000D0 * t79 * t81) - 0.32
     >2500000000000000000D0 * t79 * t81 * t78
      zk(i) = rhob * t89 * (0.1D1 + 0.437500000000000000000D-1 * t89 * t
     >60 * sigmabb / t41 / rhob / t80 / taub)

               elseif(rhob.lt.tol) then
             rho = rhoa
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) + 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmaaa
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i) + 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taua
      t15 = 0.1D1 / pi / rhoa
      t16 = t15 ** (0.1D1 / 0.3D1)
      t19 = t15 ** (0.1D1 / 0.6D1)
      t22 = sqrt(t15)
      t24 = t16 ** 2
      t30 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t19 + 0.563098414909787598194D1 * t16 + 0.291521471421917
     >737271D1 * t22 + 0.516066464547863440989D0 * t24))
      t31 = (0.1D1 + 0.186690969707574028554D0 * t16) * t30
      t33 = pi ** 2
      t35 = t33 * pi
      t38 = (t33 * rhoa) ** (0.1D1 / 0.3D1)
      t39 = 0.1D1 / t38
      t41 = rhoa ** 2
      t42 = 0.1D1 / t41
      t45 = exp(0.202642426794280972788D0 * t31 * t33)
      t46 = t45 - 0.1D1
      t52 = 0.149583857013095144140D-1 * t35 / t46 * sigmaaa * t39 * t42
      t55 = t33 ** 2
      t57 = t46 ** 2
      t60 = sigmaaa ** 2
      t61 = t38 ** 2
      t64 = t41 ** 2
      t75 = log(0.1D1 + 0.149583857013095144140D-1 * t35 * sigmaaa * t39
     > * t42 * (0.1D1 + t52) / (0.1D1 + t52 + 0.223753302789140933370D-3
     > * t55 * t33 / t57 * t60 / t61 / t64))
      t78 = -0.3109070D-1 * t31 + 0.153426409720027345292D0 / t33 * t75
      t79 = t60 * t42
      t80 = taua ** 2
      t81 = 0.1D1 / t80
      t89 = t78 * (0.1D1 + 0.260000000000000000000D0 * t79 * t81) - 0.32
     >2500000000000000000D0 * t79 * t81 * t78
      zk(i) = rhoa * t89 * (0.1D1 + 0.437500000000000000000D-1 * t89 * t
     >60 * sigmaaa / t41 / rhoa / t80 / taua)

               else
             rho = rhoa + rhob
      t2 = 0.250000000000000000000D0 * sigmacc(i)
      t4 = 0.250000000000000000000D0 * sigmaoo(i)
      t6 = 0.500000000000000000000D0 * sigmaco(i)
      sigmaaa = max(0.0D0, t2 + t4 + t6)
      sigmaab = t2 - t4
      sigmabb = max(0.0D0, t2 + t4 - t6)
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t11 = 0.500000000000000000000D0 * tauc(i)
      t13 = 0.500000000000000000000D0 * tauo(i)
      taua = max(0.0D0, t11 + t13)
      taub = max(0.0D0, t11 - t13)
      tau = taua + taub
      t16 = 0.1D1 / pi
      t17 = 0.1D1 / rho
      t18 = t16 * t17
      t19 = t18 ** (0.1D1 / 0.3D1)
      t22 = t18 ** (0.1D1 / 0.6D1)
      t25 = sqrt(t18)
      t27 = t19 ** 2
      t33 = log(0.1D1 + 0.160819794986925350668D2 / (0.72401019343168311
     >3327D1 * t22 + 0.325955091942229212011D1 * t19 + 0.141872281647966
     >739112D1 * t25 + 0.406913004517529319387D0 * t27))
      t34 = (0.1D1 + 0.194159335344114122552D0 * t19) * t33
      t35 = 0.621814D-1 * t34
      t46 = log(0.1D1 + 0.296087499777934375166D2 / (0.98721297225692720
     >9438D1 * t22 + 0.329180480994506259905D1 * t19 + 0.762327521935289
     >963194D0 * t25 + 0.410025070949612505036D0 * t27))
      t49 = rhoa - 0.1D1 * rhob
      t50 = t49 * t17
      t51 = 0.1D1 + t50
      t52 = t51 ** (0.1D1 / 0.3D1)
      t53 = t52 * t51
      t55 = 0.1D1 - 0.1D1 * t50
      t56 = t55 ** (0.1D1 / 0.3D1)
      t57 = t56 * t55
      t58 = t53 + t57 - 0.2D1
      t59 = t49 ** 2
      t60 = t59 ** 2
      t61 = rho ** 2
      t62 = t61 ** 2
      t63 = 0.1D1 / t62
      t64 = t60 * t63
      t68 = (0.1D1 + 0.101077332976287768525D0 * t19) * t46 * t58 * (0.1
     >D1 - 0.1D1 * t64)
      t69 = 0.379955235370239451738D-1 * t68
      t80 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t22 + 0.563098414909787598194D1 * t19 + 0.291521471421917
     >737271D1 * t25 + 0.516066464547863440989D0 * t27))
      t85 = (-0.3109070D-1 * (0.1D1 + 0.186690969707574028554D0 * t19) *
     > t80 + t35) * t58 * t64
      t86 = 0.192366105093153631974D1 * t85
      t87 = pi ** 2
      t88 = 0.1D1 / t87
      t89 = t52 ** 2
      t91 = t56 ** 2
      t93 = 0.500000000000000000000D0 * t89 + 0.500000000000000000000D0 
     >* t91
      t94 = t93 ** 2
      t95 = t94 * t93
      t97 = t87 * pi
      t99 = 0.1D1 / t94
      t102 = (t87 * rho) ** (0.1D1 / 0.3D1)
      t103 = 0.1D1 / t102
      t104 = 0.1D1 / t61
      t111 = exp(-0.325889135327092945460D1 * (-t35 + t69 + t86) * t87 /
     > t95)
      t112 = t111 - 0.1D1
      t119 = 0.942319250876317101329D-2 * t97 / t112 * sigma * t99 * t10
     >3 * t104
      t121 = t87 ** 2
      t122 = t121 * t87
      t123 = t112 ** 2
      t126 = sigma ** 2
      t128 = t94 ** 2
      t130 = t102 ** 2
      t131 = 0.1D1 / t130
      t143 = log(0.1D1 + 0.942319250876317101329D-2 * t97 * sigma * t99 
     >* t103 * t104 * (0.1D1 + t119) / (0.1D1 + t119 + 0.887965570572103
     >448132D-4 * t122 / t123 * t126 / t128 * t131 * t63))
      t144 = t88 * t95 * t143
      t146 = -t35 + t69 + t86 + 0.306852819440054690583D0 * t144
      t156 = t55 ** 2
      t158 = t51 ** 2
      t172 = (0.1D1 + 0.151426716069344971574D0 * (t156 * sigmaaa + t158
     > * sigmabb - 0.2D1 * t55 * t51 * sigmaab) * t104 * t131 * (0.1D1 /
     > t53 + 0.1D1 / t57)) ** 2
      t173 = t172 ** 2
      t175 = (0.53D0 + 0.87D0 * t59 * t104 + 0.50D0 * t64 + 0.226D1 * t6
     >0 * t59 / t62 / t61) / t173
      t177 = tau ** 2
      t178 = 0.1D1 / t177
      t188 = 0.621814D29 * t34
      t189 = 0.379955235370239451738D29 * t68
      t190 = 0.192366105093153631974D31 * t85
      t191 = 0.306852819440054690583D30 * t144
      t193 = t16 / rhoa
      t194 = t193 ** (0.1D1 / 0.3D1)
      t197 = t193 ** (0.1D1 / 0.6D1)
      t200 = sqrt(t193)
      t202 = t194 ** 2
      t208 = log(0.1D1 + 0.321639589973850701335D2 / (0.1345791371439444
     >77912D2 * t197 + 0.563098414909787598194D1 * t194 + 0.291521471421
     >917737271D1 * t200 + 0.516066464547863440989D0 * t202))
      t209 = (0.1D1 + 0.186690969707574028554D0 * t194) * t208
      t213 = (t87 * rhoa) ** (0.1D1 / 0.3D1)
      t214 = 0.1D1 / t213
      t216 = rhoa ** 2
      t217 = 0.1D1 / t216
      t220 = exp(0.202642426794280972788D0 * t209 * t87)
      t221 = t220 - 0.1D1
      t227 = 0.149583857013095144140D-1 * t97 / t221 * sigmaaa * t214 * 
     >t217
      t230 = t221 ** 2
      t233 = sigmaaa ** 2
      t234 = t213 ** 2
      t237 = t216 ** 2
      t248 = log(0.1D1 + 0.149583857013095144140D-1 * t97 * sigmaaa * t2
     >14 * t217 * (0.1D1 + t227) / (0.1D1 + t227 + 0.2237533027891409333
     >70D-3 * t122 / t230 * t233 / t234 / t237))
      t249 = t88 * t248
      t252 = tanh(t188 - t189 - t190 - t191 - 0.3109070D29 * t209 + 0.15
     >3426409720027345292D30 * t249)
      t253 = 0.500000000000000000000D0 * t252
      t265 = t16 / rhob
      t266 = t265 ** (0.1D1 / 0.3D1)
      t269 = t265 ** (0.1D1 / 0.6D1)
      t272 = sqrt(t265)
      t274 = t266 ** 2
      t280 = log(0.1D1 + 0.321639589973850701335D2 / (0.1345791371439444
     >77912D2 * t269 + 0.563098414909787598194D1 * t266 + 0.291521471421
     >917737271D1 * t272 + 0.516066464547863440989D0 * t274))
      t281 = (0.1D1 + 0.186690969707574028554D0 * t266) * t280
      t285 = (t87 * rhob) ** (0.1D1 / 0.3D1)
      t286 = 0.1D1 / t285
      t288 = rhob ** 2
      t289 = 0.1D1 / t288
      t292 = exp(0.202642426794280972788D0 * t281 * t87)
      t293 = t292 - 0.1D1
      t299 = 0.149583857013095144140D-1 * t97 / t293 * sigmabb * t286 * 
     >t289
      t302 = t293 ** 2
      t305 = sigmabb ** 2
      t306 = t285 ** 2
      t309 = t288 ** 2
      t320 = log(0.1D1 + 0.149583857013095144140D-1 * t97 * sigmabb * t2
     >86 * t289 * (0.1D1 + t299) / (0.1D1 + t299 + 0.2237533027891409333
     >70D-3 * t122 / t302 * t305 / t306 / t309))
      t321 = t88 * t320
      t324 = tanh(t188 - t189 - t190 - t191 - 0.3109070D29 * t281 + 0.15
     >3426409720027345292D30 * t321)
      t325 = 0.500000000000000000000D0 * t324
      t339 = t146 * (0.1D1 + 0.625000000000000000000D-1 * t175 * t126 * 
     >t104 * t178) - 0.625000000000000000000D-1 * (0.1D1 + t175) * t126
     >* t104 * t178 * (rhoa * t17 * ((0.500000000000000000000D0 - t253)
     >* t146 + (0.500000000000000000000D0 + t253) * (-0.3109070D-1 * t20
     >9 + 0.153426409720027345292D0 * t249)) + rhob * t17 * ((0.50000000
     >0000000000000D0 - t325) * t146 + (0.500000000000000000000D0 + t325
     >) * (-0.3109070D-1 * t281 + 0.153426409720027345292D0 * t321)))
      zk(i) = rho * t339 * (0.1D1 + 0.437500000000000000000D-1 * t339 * 
     >t126 * sigma / t61 / rho / t177 / tau)

               endif
             endif
           enddo
         endif
       else
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i))
      sigmaab = sigmaaa
      sigmabb = sigmaab
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i))
      taub = taua
      tau = taua + taub
      t8 = 0.1D1 / pi
      t9 = 0.1D1 / rho
      t10 = t8 * t9
      t11 = t10 ** (0.1D1 / 0.3D1)
      t13 = 0.1D1 + 0.194159335344114122552D0 * t11
      t14 = t10 ** (0.1D1 / 0.6D1)
      t17 = sqrt(t10)
      t19 = t11 ** 2
      t21 = 0.724010193431683113327D1 * t14 + 0.325955091942229212011D1 
     >* t11 + 0.141872281647966739112D1 * t17 + 0.406913004517529319387D
     >0 * t19
      t24 = 0.1D1 + 0.160819794986925350668D2 / t21
      t25 = log(t24)
      t26 = t13 * t25
      t27 = 0.621814D-1 * t26
      t29 = 0.1D1 + 0.101077332976287768525D0 * t11
      t34 = 0.987212972256927209438D1 * t14 + 0.329180480994506259905D1 
     >* t11 + 0.762327521935289963194D0 * t17 + 0.410025070949612505036D
     >0 * t19
      t37 = 0.1D1 + 0.296087499777934375166D2 / t34
      t38 = log(t37)
      t39 = t29 * t38
      t41 = rhoa - 0.1D1 * rhob
      t42 = t41 * t9
      t43 = 0.1D1 + t42
      t44 = t43 ** (0.1D1 / 0.3D1)
      t45 = t44 * t43
      t47 = 0.1D1 - 0.1D1 * t42
      t48 = t47 ** (0.1D1 / 0.3D1)
      t49 = t48 * t47
      t50 = t45 + t49 - 0.2D1
      t51 = t41 ** 2
      t52 = t51 ** 2
      t53 = rho ** 2
      t54 = t53 ** 2
      t55 = 0.1D1 / t54
      t56 = t52 * t55
      t58 = 0.1D1 - 0.1D1 * t56
      t60 = t39 * t50 * t58
      t61 = 0.379955235370239451738D-1 * t60
      t63 = 0.1D1 + 0.186690969707574028554D0 * t11
      t68 = 0.134579137143944477912D2 * t14 + 0.563098414909787598194D1 
     >* t11 + 0.291521471421917737271D1 * t17 + 0.516066464547863440989D
     >0 * t19
      t71 = 0.1D1 + 0.321639589973850701335D2 / t68
      t72 = log(t71)
      t75 = -0.3109070D-1 * t63 * t72 + t27
      t76 = t75 * t50
      t77 = t76 * t56
      t78 = 0.192366105093153631974D1 * t77
      t79 = pi ** 2
      t80 = 0.1D1 / t79
      t81 = t44 ** 2
      t83 = t48 ** 2
      t85 = 0.500000000000000000000D0 * t81 + 0.500000000000000000000D0 
     >* t83
      t86 = t85 ** 2
      t87 = t86 * t85
      t88 = t80 * t87
      t89 = t79 * pi
      t90 = t89 * sigma
      t91 = 0.1D1 / t86
      t92 = t90 * t91
      t93 = t79 * rho
      t94 = t93 ** (0.1D1 / 0.3D1)
      t95 = 0.1D1 / t94
      t96 = 0.1D1 / t53
      t97 = t95 * t96
      t99 = (-t27 + t61 + t78) * t79
      t100 = 0.1D1 / t87
      t103 = exp(-0.325889135327092945460D1 * t99 * t100)
      t104 = t103 - 0.1D1
      t105 = 0.1D1 / t104
      t106 = t89 * t105
      t107 = t106 * sigma
      t108 = t91 * t95
      t109 = t108 * t96
      t111 = 0.942319250876317101329D-2 * t107 * t109
      t112 = 0.1D1 + t111
      t113 = t79 ** 2
      t114 = t113 * t79
      t115 = t104 ** 2
      t116 = 0.1D1 / t115
      t117 = t114 * t116
      t118 = sigma ** 2
      t119 = t117 * t118
      t120 = t86 ** 2
      t121 = 0.1D1 / t120
      t122 = t94 ** 2
      t123 = 0.1D1 / t122
      t124 = t121 * t123
      t125 = t124 * t55
      t128 = 0.1D1 + t111 + 0.887965570572103448132D-4 * t119 * t125
      t129 = 0.1D1 / t128
      t130 = t112 * t129
      t134 = 0.1D1 + 0.942319250876317101329D-2 * t92 * t97 * t130
      t135 = log(t134)
      t136 = t88 * t135
      t138 = -t27 + t61 + t78 + 0.306852819440054690583D0 * t136
      t142 = t52 * t51
      t144 = 0.1D1 / t54 / t53
      t147 = 0.53D0 + 0.87D0 * t51 * t96 + 0.50D0 * t56 + 0.226D1 * t142
     > * t144
      t148 = t47 ** 2
      t150 = t43 ** 2
      t152 = t47 * t43
      t155 = t148 * sigmaaa + t150 * sigmabb - 0.2D1 * t152 * sigmaab
      t156 = t155 * t96
      t159 = 0.1D1 / t45 + 0.1D1 / t49
      t160 = t123 * t159
      t163 = 0.1D1 + 0.151426716069344971574D0 * t156 * t160
      t164 = t163 ** 2
      t165 = t164 ** 2
      t166 = 0.1D1 / t165
      t167 = t147 * t166
      t168 = t118 * t96
      t169 = tau ** 2
      t170 = 0.1D1 / t169
      t171 = t168 * t170
      t174 = 0.1D1 + 0.625000000000000000000D-1 * t167 * t171
      t176 = 0.1D1 + t167
      t177 = t176 * t118
      t178 = t96 * t170
      t179 = rhoa * t9
      t180 = 0.621814D29 * t26
      t181 = 0.379955235370239451738D29 * t60
      t182 = 0.192366105093153631974D31 * t77
      t183 = 0.306852819440054690583D30 * t136
      t185 = t8 / rhoa
      t186 = t185 ** (0.1D1 / 0.3D1)
      t188 = 0.1D1 + 0.186690969707574028554D0 * t186
      t189 = t185 ** (0.1D1 / 0.6D1)
      t192 = sqrt(t185)
      t194 = t186 ** 2
      t196 = 0.134579137143944477912D2 * t189 + 0.563098414909787598194D
     >1 * t186 + 0.291521471421917737271D1 * t192 + 0.516066464547863440
     >989D0 * t194
      t199 = 0.1D1 + 0.321639589973850701335D2 / t196
      t200 = log(t199)
      t201 = t188 * t200
      t204 = t79 * rhoa
      t205 = t204 ** (0.1D1 / 0.3D1)
      t206 = 0.1D1 / t205
      t207 = t89 * sigmaaa * t206
      t208 = rhoa ** 2
      t209 = 0.1D1 / t208
      t212 = exp(0.202642426794280972788D0 * t201 * t79)
      t213 = t212 - 0.1D1
      t214 = 0.1D1 / t213
      t215 = t89 * t214
      t216 = sigmaaa * t206
      t219 = 0.149583857013095144140D-1 * t215 * t216 * t209
      t220 = 0.1D1 + t219
      t221 = t209 * t220
      t222 = t213 ** 2
      t223 = 0.1D1 / t222
      t224 = t114 * t223
      t225 = sigmaaa ** 2
      t226 = t205 ** 2
      t227 = 0.1D1 / t226
      t228 = t225 * t227
      t229 = t208 ** 2
      t230 = 0.1D1 / t229
      t234 = 0.1D1 + t219 + 0.223753302789140933370D-3 * t224 * t228 * t
     >230
      t235 = 0.1D1 / t234
      t236 = t221 * t235
      t239 = 0.1D1 + 0.149583857013095144140D-1 * t207 * t236
      t240 = log(t239)
      t241 = t80 * t240
      t243 = t180 - t181 - t182 - t183 - 0.3109070D29 * t201 + 0.1534264
     >09720027345292D30 * t241
      t244 = tanh(t243)
      t245 = 0.500000000000000000000D0 * t244
      t246 = 0.500000000000000000000D0 - t245
      t248 = 0.500000000000000000000D0 + t245
      t251 = -0.3109070D-1 * t201 + 0.153426409720027345292D0 * t241
      t253 = t246 * t138 + t248 * t251
      t255 = rhob * t9
      t257 = t8 / rhob
      t258 = t257 ** (0.1D1 / 0.3D1)
      t260 = 0.1D1 + 0.186690969707574028554D0 * t258
      t261 = t257 ** (0.1D1 / 0.6D1)
      t264 = sqrt(t257)
      t266 = t258 ** 2
      t268 = 0.134579137143944477912D2 * t261 + 0.563098414909787598194D
     >1 * t258 + 0.291521471421917737271D1 * t264 + 0.516066464547863440
     >989D0 * t266
      t271 = 0.1D1 + 0.321639589973850701335D2 / t268
      t272 = log(t271)
      t273 = t260 * t272
      t276 = t79 * rhob
      t277 = t276 ** (0.1D1 / 0.3D1)
      t278 = 0.1D1 / t277
      t279 = t89 * sigmabb * t278
      t280 = rhob ** 2
      t281 = 0.1D1 / t280
      t284 = exp(0.202642426794280972788D0 * t273 * t79)
      t285 = t284 - 0.1D1
      t286 = 0.1D1 / t285
      t287 = t89 * t286
      t288 = sigmabb * t278
      t291 = 0.149583857013095144140D-1 * t287 * t288 * t281
      t292 = 0.1D1 + t291
      t293 = t281 * t292
      t294 = t285 ** 2
      t295 = 0.1D1 / t294
      t296 = t114 * t295
      t297 = sigmabb ** 2
      t298 = t277 ** 2
      t299 = 0.1D1 / t298
      t300 = t297 * t299
      t301 = t280 ** 2
      t302 = 0.1D1 / t301
      t306 = 0.1D1 + t291 + 0.223753302789140933370D-3 * t296 * t300 * t
     >302
      t307 = 0.1D1 / t306
      t308 = t293 * t307
      t311 = 0.1D1 + 0.149583857013095144140D-1 * t279 * t308
      t312 = log(t311)
      t313 = t80 * t312
      t315 = t180 - t181 - t182 - t183 - 0.3109070D29 * t273 + 0.1534264
     >09720027345292D30 * t313
      t316 = tanh(t315)
      t317 = 0.500000000000000000000D0 * t316
      t318 = 0.500000000000000000000D0 - t317
      t320 = 0.500000000000000000000D0 + t317
      t323 = -0.3109070D-1 * t273 + 0.153426409720027345292D0 * t313
      t325 = t318 * t138 + t320 * t323
      t327 = t179 * t253 + t255 * t325
      t328 = t178 * t327
      t331 = t138 * t174 - 0.625000000000000000000D-1 * t177 * t328
      t332 = rho * t331
      t333 = t118 * sigma
      t334 = t331 * t333
      t335 = t53 * rho
      t336 = 0.1D1 / t335
      t338 = 0.1D1 / t169 / tau
      t339 = t336 * t338
      t342 = 0.1D1 + 0.437500000000000000000D-1 * t334 * t339
      zk(i) = t332 * t342
      t345 = 0.1D1 / t19 * t8
      t347 = t345 * t96 * t25
      t348 = 0.402436643158883263334D-2 * t347
      t349 = t21 ** 2
      t352 = t14 ** 2
      t353 = t352 ** 2
      t357 = 0.1D1 / t353 / t14 * t8 * t96
      t359 = t345 * t96
      t363 = 0.1D1 / t17 * t8 * t96
      t367 = 0.1D1 / t11 * t8 * t96
      t372 = t13 / t349 * (-0.120668365571947185555D1 * t357 - 0.1086516
     >97314076404004D1 * t359 - 0.709361408239833695563D0 * t363 - 0.271
     >275336345019546258D0 * t367) / t24
      t373 = 0.100000000000000000000D1 * t372
      t376 = t359 * t38 * t50 * t58
      t377 = 0.128016206138671616206D-2 * t376
      t378 = t34 ** 2
      t390 = t29 / t378 * (-0.164535495376154534908D1 * t357 - 0.1097268
     >26998168753302D1 * t359 - 0.381163760967644981600D0 * t363 - 0.273
     >350047299741670024D0 * t367) / t37 * t50 * t58
      t391 = 0.112499995668310776915D1 * t390
      t392 = t41 * t96
      t393 = 0.1D1 * t392
      t394 = t9 - t393
      t397 = 0.1D1 * t9
      t398 = -t397 + t392
      t401 = 0.133333333333333333333D1 * t44 * t394 + 0.1333333333333333
     >33333D1 * t48 * t398
      t403 = t39 * t401 * t58
      t404 = 0.379955235370239451738D-1 * t403
      t406 = t51 * t41 * t55
      t407 = 0.4D1 * t406
      t409 = 0.1D1 / t54 / rho
      t410 = t52 * t409
      t411 = 0.4D1 * t410
      t414 = t39 * t50 * (-t407 + t411)
      t415 = 0.379955235370239451738D-1 * t414
      t419 = t68 ** 2
      t433 = (0.193478431062909061652D-2 * t345 * t96 * t72 + 0.10000000
     >0000000000000D1 * t63 / t419 * (-0.224298561906574129855D1 * t357
     >- 0.187699471636595866065D1 * t359 - 0.145760735710958868637D1 * t
     >363 - 0.344044309698575627326D0 * t367) / t71 - t348 - t373) * t50
     > * t56
      t434 = 0.192366105093153631974D1 * t433
      t436 = t75 * t401 * t56
      t437 = 0.192366105093153631974D1 * t436
      t438 = t76 * t406
      t439 = 0.769464420372614527896D1 * t438
      t440 = t76 * t410
      t441 = 0.769464420372614527896D1 * t440
      t442 = t80 * t86
      t443 = 0.1D1 / t44
      t446 = 0.1D1 / t48
      t449 = 0.333333333333333333333D0 * t443 * t394 + 0.333333333333333
     >333333D0 * t446 * t398
      t451 = t442 * t135 * t449
      t453 = t100 * t95
      t454 = t90 * t453
      t455 = t96 * t112
      t460 = t113 * pi
      t464 = 0.1D1 / t94 / t93
      t468 = 0.314106416958772367110D-2 * t460 * sigma * t91 * t464 * t9
     >6 * t130
      t472 = 0.188463850175263420266D-1 * t92 * t95 * t336 * t130
      t475 = t89 * t116 * sigma * t91
      t484 = (-0.325889135327092945460D1 * (t348 + t373 - t377 - t391 + 
     >t404 + t415 + t434 + t437 + t439 - t441) * t79 * t100 + 0.97766740
     >5981278836380D1 * t99 * t121 * t449) * t103
      t487 = 0.942319250876317101329D-2 * t475 * t97 * t484
      t491 = 0.188463850175263420266D-1 * t107 * t453 * t96 * t449
      t497 = 0.314106416958772367110D-2 * t460 * t105 * sigma * t91 * t4
     >64 * t96
      t500 = 0.188463850175263420266D-1 * t107 * t108 * t336
      t506 = t90 * t108
      t507 = t128 ** 2
      t508 = 0.1D1 / t507
      t513 = t114 / t115 / t104 * t118 * t121
      t514 = t123 * t55
      t520 = 0.1D1 / t120 / t85 * t123
      t525 = t113 ** 2
      t529 = 0.1D1 / t122 / t93
      t533 = 0.591977047048068965421D-4 * t525 * t116 * t118 * t121 * t5
     >29 * t55
      t536 = 0.355186228228841379252D-3 * t119 * t124 * t409
      t543 = 0.1D1 / t134
      t545 = t88 * (-0.188463850175263420266D-1 * t454 * t455 * t129 * t
     >449 - t468 - t472 + 0.942319250876317101329D-2 * t92 * t97 * (-t48
     >7 - t491 - t497 - t500) * t129 - 0.942319250876317101329D-2 * t506
     > * t455 * t508 * (-t487 - t491 - t497 - t500 - 0.17759311411442068
     >9627D-3 * t513 * t514 * t484 - 0.355186228228841379252D-3 * t119 *
     > t520 * t55 * t449 - t533 - t536)) * t543
      t547 = t348 + t373 - t377 - t391 + t404 + t415 + t434 + t437 + t43
     >9 - t441 + 0.920558458320164071749D0 * t451 + 0.306852819440054690
     >583D0 * t545
      t554 = 0.174D1 * t392 + 0.200D1 * t406 + 0.1356D2 * t52 * t41 * t1
     >44
      t555 = t554 * t166
      t559 = 0.1D1 / t165 / t163
      t560 = t147 * t559
      t561 = t560 * t118
      t562 = t47 * sigmaaa
      t565 = t43 * sigmabb
      t574 = -0.2D1 * t562 * t9 + 0.2D1 * t565 * t9 + 0.2D1 * t9 * t43 *
     > sigmaab - 0.2D1 * t47 * t9 * sigmaab
      t579 = 0.1D1 / t44 / t150
      t583 = 0.1D1 / t48 / t148
      t586 = -0.133333333333333333333D1 * t579 * t9 + 0.1333333333333333
     >33333D1 * t583 * t9
      t590 = 0.151426716069344971574D0 * t574 * t96 * t160 + 0.151426716
     >069344971574D0 * t156 * t123 * t586
      t605 = 0.1D1 * rhoa * t96 * t253
      t607 = tanh(-t243)
      t608 = t607 ** 2
      t610 = 0.1D1 - 0.1D1 * t608
      t611 = 0.402436643158883263334D28 * t347
      t612 = 0.100000000000000000000D31 * t372
      t613 = 0.128016206138671616206D28 * t376
      t614 = 0.112499995668310776915D31 * t390
      t615 = 0.379955235370239451738D29 * t403
      t616 = 0.379955235370239451738D29 * t414
      t617 = 0.192366105093153631974D31 * t433
      t618 = 0.192366105093153631974D31 * t436
      t619 = 0.769464420372614527896D31 * t438
      t620 = 0.769464420372614527896D31 * t440
      t621 = 0.920558458320164071749D30 * t451
      t622 = 0.306852819440054690583D30 * t545
      t623 = 0.1D1 / t194
      t624 = t623 * t8
      t625 = t209 * t200
      t626 = t624 * t625
      t628 = t196 ** 2
      t630 = t188 / t628
      t631 = t189 ** 2
      t632 = t631 ** 2
      t650 = (-0.224298561906574129855D1 / t632 / t189 * t8 * t209 - 0.1
     >87699471636595866065D1 * t624 * t209 - 0.145760735710958868637D1 /
     > t192 * t8 * t209 - 0.344044309698575627326D0 / t186 * t8 * t209)
     >/ t199
      t651 = t630 * t650
      t655 = 0.1D1 / t205 / t204
      t660 = 0.1D1 / t208 / rhoa
      t667 = t206 * t209
      t675 = (-0.126105037207067989169D-1 * t623 * pi * t625 - 0.6517782
     >70654185890920D1 * t630 * t650 * t79) * t212
      t678 = 0.149583857013095144140D-1 * t89 * t223 * sigmaaa * t667 * 
     >t675
      t683 = 0.498612856710317147135D-2 * t460 * t214 * sigmaaa * t655 *
     > t209
      t686 = 0.299167714026190288281D-1 * t215 * t216 * t660
      t692 = t234 ** 2
      t693 = 0.1D1 / t692
      t721 = 0.1D1 / t239
      t722 = t80 * (-0.498612856710317147135D-2 * t460 * sigmaaa * t655 
     >* t236 - 0.299167714026190288281D-1 * t207 * t660 * t220 * t235 +
     >0.149583857013095144140D-1 * t207 * t209 * (-t678 - t683 - t686) *
     > t235 - 0.149583857013095144140D-1 * t207 * t221 * t693 * (-t678 -
     > t683 - t686 - 0.447506605578281866741D-3 * t114 / t222 / t213 * t
     >225 * t227 * t230 * t675 - 0.149168868526093955580D-3 * t525 * t22
     >3 * t225 / t226 / t204 * t230 - 0.895013211156563733482D-3 * t224
     >* t228 / t229 / rhoa)) * t721
      t724 = -t611 - t612 + t613 + t614 - t615 - t616 - t617 - t618 - t6
     >19 + t620 - t621 - t622 + 0.193478431062909061652D28 * t626 + 0.10
     >0000000000000000000D31 * t651 + 0.153426409720027345292D30 * t722
      t725 = t610 * t724
      t740 = 0.1D1 * rhob * t96 * t325
      t742 = tanh(-t315)
      t743 = t742 ** 2
      t745 = 0.1D1 - 0.1D1 * t743
      t746 = -t611 - t612 + t613 + t614 - t615 - t616 - t617 - t618 - t6
     >19 + t620 - t621 - t622
      t747 = t745 * t746
      t759 = t547 * t174 + t138 * (0.625000000000000000000D-1 * t555 * t
     >171 - 0.250000000000000000000D0 * t561 * t178 * t590) - 0.62500000
     >0000000000000D-1 * (t555 - 0.4D1 * t560 * t590) * t118 * t328 - 0.
     >625000000000000000000D-1 * t177 * t178 * (t9 * t253 - t605 + t179
     >* (-0.500000000000000000000D0 * t725 * t138 + t246 * t547 + 0.5000
     >00000000000000000D0 * t725 * t251 + t248 * (0.19347843106290906165
     >2D-2 * t626 + 0.100000000000000000000D1 * t651 + 0.153426409720027
     >345292D0 * t722)) - t740 + t255 * (-0.500000000000000000000D0 * t7
     >47 * t138 + t318 * t547 + 0.500000000000000000000D0 * t747 * t323)
     >)
      t763 = t96 * t331
      t768 = -t397 - t393
      t771 = t9 + t392
      t774 = 0.133333333333333333333D1 * t44 * t768 + 0.1333333333333333
     >33333D1 * t48 * t771
      t776 = t39 * t774 * t58
      t777 = 0.379955235370239451738D-1 * t776
      t780 = t39 * t50 * (t407 + t411)
      t781 = 0.379955235370239451738D-1 * t780
      t783 = t75 * t774 * t56
      t784 = 0.192366105093153631974D1 * t783
      t789 = 0.333333333333333333333D0 * t443 * t768 + 0.333333333333333
     >333333D0 * t446 * t771
      t791 = t442 * t135 * t789
      t805 = (-0.325889135327092945460D1 * (t348 + t373 - t377 - t391 + 
     >t777 + t781 + t434 + t784 - t439 - t441) * t79 * t100 + 0.97766740
     >5981278836380D1 * t99 * t121 * t789) * t103
      t808 = 0.942319250876317101329D-2 * t475 * t97 * t805
      t812 = 0.188463850175263420266D-1 * t107 * t453 * t96 * t789
      t832 = t88 * (-0.188463850175263420266D-1 * t454 * t455 * t129 * t
     >789 - t468 - t472 + 0.942319250876317101329D-2 * t92 * t97 * (-t80
     >8 - t812 - t497 - t500) * t129 - 0.942319250876317101329D-2 * t506
     > * t455 * t508 * (-t808 - t812 - t497 - t500 - 0.17759311411442068
     >9627D-3 * t513 * t514 * t805 - 0.355186228228841379252D-3 * t119 *
     > t520 * t55 * t789 - t533 - t536)) * t543
      t834 = t348 + t373 - t377 - t391 + t777 + t781 + t434 + t784 - t43
     >9 - t441 + 0.920558458320164071749D0 * t791 + 0.306852819440054690
     >583D0 * t832
      t837 = -t554 * t166
      t848 = -0.151426716069344971574D0 * t574 * t96 * t160 - 0.15142671
     >6069344971574D0 * t156 * t123 * t586
      t860 = 0.379955235370239451738D29 * t776
      t861 = 0.379955235370239451738D29 * t780
      t862 = 0.192366105093153631974D31 * t783
      t863 = 0.920558458320164071749D30 * t791
      t864 = 0.306852819440054690583D30 * t832
      t865 = -t611 - t612 + t613 + t614 - t860 - t861 - t617 - t862 + t6
     >19 + t620 - t863 - t864
      t866 = t610 * t865
      t875 = 0.1D1 / t266
      t876 = t875 * t8
      t877 = t281 * t272
      t878 = t876 * t877
      t880 = t268 ** 2
      t882 = t260 / t880
      t883 = t261 ** 2
      t884 = t883 ** 2
      t902 = (-0.224298561906574129855D1 / t884 / t261 * t8 * t281 - 0.1
     >87699471636595866065D1 * t876 * t281 - 0.145760735710958868637D1 /
     > t264 * t8 * t281 - 0.344044309698575627326D0 / t258 * t8 * t281)
     >/ t271
      t903 = t882 * t902
      t907 = 0.1D1 / t277 / t276
      t912 = 0.1D1 / t280 / rhob
      t919 = t278 * t281
      t927 = (-0.126105037207067989169D-1 * t875 * pi * t877 - 0.6517782
     >70654185890920D1 * t882 * t902 * t79) * t284
      t930 = 0.149583857013095144140D-1 * t89 * t295 * sigmabb * t919 * 
     >t927
      t935 = 0.498612856710317147135D-2 * t460 * t286 * sigmabb * t907 *
     > t281
      t938 = 0.299167714026190288281D-1 * t287 * t288 * t912
      t944 = t306 ** 2
      t945 = 0.1D1 / t944
      t973 = 0.1D1 / t311
      t974 = t80 * (-0.498612856710317147135D-2 * t460 * sigmabb * t907 
     >* t308 - 0.299167714026190288281D-1 * t279 * t912 * t292 * t307 +
     >0.149583857013095144140D-1 * t279 * t281 * (-t930 - t935 - t938) *
     > t307 - 0.149583857013095144140D-1 * t279 * t293 * t945 * (-t930 -
     > t935 - t938 - 0.447506605578281866741D-3 * t114 / t294 / t285 * t
     >297 * t299 * t302 * t927 - 0.149168868526093955580D-3 * t525 * t29
     >5 * t297 / t298 / t276 * t302 - 0.895013211156563733482D-3 * t296
     >* t300 / t301 / rhob)) * t973
      t976 = -t611 - t612 + t613 + t614 - t860 - t861 - t617 - t862 + t6
     >19 + t620 - t863 - t864 + 0.193478431062909061652D28 * t878 + 0.10
     >0000000000000000000D31 * t903 + 0.153426409720027345292D30 * t974
      t977 = t745 * t976
      t994 = t834 * t174 + t138 * (0.625000000000000000000D-1 * t837 * t
     >171 - 0.250000000000000000000D0 * t561 * t178 * t848) - 0.62500000
     >0000000000000D-1 * (t837 - 0.4D1 * t560 * t848) * t118 * t328 - 0.
     >625000000000000000000D-1 * t177 * t178 * (-t605 + t179 * (-0.50000
     >0000000000000000D0 * t866 * t138 + t246 * t834 + 0.500000000000000
     >000000D0 * t866 * t251) + t9 * t325 - t740 + t255 * (-0.5000000000
     >00000000000D0 * t977 * t138 + t318 * t834 + 0.50000000000000000000
     >0D0 * t977 * t323 + t320 * (0.193478431062909061652D-2 * t878 + 0.
     >100000000000000000000D1 * t903 + 0.153426409720027345292D0 * t974)
     >))
      t1011 = (-0.174D1 * t51 * t336 - 0.200D1 * t410 - 0.1356D2 * t142 
     >/ t54 / t335) * t166
      t1046 = 0.151426716069344971574D0 * (0.2D1 * t562 * t392 - 0.2D1 *
     > t565 * t392 - 0.2D1 * t392 * t43 * sigmaab + 0.2D1 * t47 * t41 *
     >t96 * sigmaab) * t96 * t160 - 0.302853432138689943149D0 * t155 * t
     >336 * t160 - 0.100951144046229981050D0 * t156 * t529 * t159 * t79
     >+ 0.151426716069344971574D0 * t156 * t123 * (0.1333333333333333333
     >33D1 * t579 * t41 * t96 - 0.133333333333333333333D1 * t583 * t41 *
     > t96)
      t1066 = t138 * (0.625000000000000000000D-1 * t1011 * t171 - 0.2500
     >00000000000000000D0 * t561 * t178 * t1046 - 0.12500000000000000000
     >0D0 * t167 * t118 * t336 * t170) - 0.625000000000000000000D-1 * (t
     >1011 - 0.4D1 * t560 * t1046) * t118 * t328 + 0.1250000000000000000
     >00D0 * t177 * t336 * t170 * t327
      vrhoc(i) = vrhoc(i) + 0.5D0 * rho * t759 * t342 + 0.21875000000000
     >0000000D-1 * t763 * t759 * t333 * t338 + 0.5D0 * rho * t994 * t342
     > + 0.218750000000000000000D-1 * t763 * t994 * t333 * t338 + t331 *
     > t342 + rho * t1066 * t342 + t332 * (0.437500000000000000000D-1 *
     >t1066 * t333 * t339 - 0.131250000000000000000D0 * t334 * t55 * t33
     >8)
      t1078 = t138 * t147
      t1079 = t559 * t118
      t1080 = t1078 * t1079
      t1081 = t55 * t170
      t1091 = t160 * t118 * t170 * t327
      t1094 = t177 * t96
      t1096 = t610 * t80
      t1118 = (0.149583857013095144140D-1 * t89 * t206 * t236 + 0.223753
     >302789140933370D-3 * t114 * sigmaaa * t227 * t230 * t214 * t235 -
     >0.149583857013095144140D-1 * t207 * t221 * t693 * (0.1495838570130
     >95144140D-1 * t215 * t667 + 0.447506605578281866741D-3 * t224 * si
     >gmaaa * t227 * t230)) * t721
      t1133 = -0.378566790173362428935D-1 * t1080 * t1081 * t148 * t123 
     >* t159 + 0.378566790173362428935D-1 * t560 * t148 * t55 * t1091 -
     >0.625000000000000000000D-1 * t1094 * t170 * rhoa * t9 * (-0.767132
     >048600136726458D29 * t1096 * t1118 * t138 + 0.76713204860013672645
     >8D29 * t1096 * t1118 * t251 + 0.153426409720027345292D0 * t248 * t
     >80 * t1118)
      t1151 = t745 * t80
      t1173 = (0.149583857013095144140D-1 * t89 * t278 * t308 + 0.223753
     >302789140933370D-3 * t114 * sigmabb * t299 * t302 * t286 * t307 -
     >0.149583857013095144140D-1 * t279 * t293 * t945 * (0.1495838570130
     >95144140D-1 * t287 * t919 + 0.447506605578281866741D-3 * t296 * si
     >gmabb * t299 * t302)) * t973
      t1188 = -0.378566790173362428935D-1 * t1080 * t1081 * t150 * t123 
     >* t159 + 0.378566790173362428935D-1 * t560 * t150 * t55 * t1091 -
     >0.625000000000000000000D-1 * t1094 * t170 * rhob * t9 * (-0.767132
     >048600136726458D29 * t1151 * t1173 * t138 + 0.76713204860013672645
     >8D29 * t1151 * t1173 * t323 + 0.153426409720027345292D0 * t320 * t
     >80 * t1173)
      t1208 = 0.757133580346724857871D-1 * t1078 * t1079 * t55 * t170 * 
     >t47 * t43 * t123 * t159 - 0.757133580346724857871D-1 * t560 * t152
     > * t55 * t1091
      t1237 = 0.942319250876317101329D-2 * t89 * t91 * t95 * t455 * t129
     > + 0.887965570572103448136D-4 * t114 * sigma * t121 * t514 * t105
     >* t129 - 0.942319250876317101329D-2 * t506 * t455 * t508 * (0.9423
     >19250876317101329D-2 * t106 * t109 + 0.177593114114420689627D-3 *
     >t117 * sigma * t125)
      t1238 = t1237 * t543
      t1242 = t1078 * t166
      t1250 = t1096 * t87
      t1251 = t1238 * t138
      t1256 = t87 * t1237 * t543
      t1264 = t1151 * t87
      t1279 = 0.306852819440054690583D0 * t88 * t1238 * t174 + 0.1250000
     >00000000000000D0 * t1242 * sigma * t96 * t170 - 0.1250000000000000
     >00000D0 * t176 * sigma * t328 - 0.625000000000000000000D-1 * t177
     >* t178 * (t179 * (0.153426409720027345292D30 * t1250 * t1251 + 0.3
     >06852819440054690583D0 * t246 * t80 * t1256 - 0.153426409720027345
     >292D30 * t1250 * t1238 * t251) + t255 * (0.153426409720027345292D3
     >0 * t1264 * t1251 + 0.306852819440054690583D0 * t318 * t80 * t1256
     > - 0.153426409720027345292D30 * t1264 * t1238 * t323))
      vsigmacc(i) = vsigmacc(i) + 0.25D0 * rho * t1133 * t342 + 0.109375
     >000000000000000D-1 * t763 * t1133 * t333 * t338 + 0.25D0 * rho * t
     >1188 * t342 + 0.109375000000000000000D-1 * t763 * t1188 * t333 * t
     >338 + 0.25D0 * rho * t1208 * t342 + 0.109375000000000000000D-1 * t
     >763 * t1208 * t333 * t338 + rho * t1279 * t342 + t332 * (0.4375000
     >00000000000000D-1 * t1279 * t333 * t339 + 0.131250000000000000000D
     >0 * t331 * t118 * t339)
      t1298 = -0.125000000000000000000D0 * t1242 * t168 * t338 + 0.12500
     >0000000000000000D0 * t177 * t96 * t338 * t327
      t1304 = t169 ** 2
      vtauc(i) = vtauc(i) + rho * t1298 * t342 + t332 * (0.4375000000000
     >00000000D-1 * t1298 * t333 * t339 - 0.131250000000000000000D0 * t3
     >34 * t336 / t1304)

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i))
      sigmaab = sigmaaa
      sigmabb = sigmaab
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i))
      taub = taua
      tau = taua + taub
      t8 = 0.1D1 / pi
      t9 = 0.1D1 / rho
      t10 = t8 * t9
      t11 = t10 ** (0.1D1 / 0.3D1)
      t14 = t10 ** (0.1D1 / 0.6D1)
      t17 = sqrt(t10)
      t19 = t11 ** 2
      t25 = log(0.1D1 + 0.160819794986925350668D2 / (0.72401019343168311
     >3327D1 * t14 + 0.325955091942229212011D1 * t11 + 0.141872281647966
     >739112D1 * t17 + 0.406913004517529319387D0 * t19))
      t26 = (0.1D1 + 0.194159335344114122552D0 * t11) * t25
      t27 = 0.621814D-1 * t26
      t38 = log(0.1D1 + 0.296087499777934375166D2 / (0.98721297225692720
     >9438D1 * t14 + 0.329180480994506259905D1 * t11 + 0.762327521935289
     >963194D0 * t17 + 0.410025070949612505036D0 * t19))
      t41 = rhoa - 0.1D1 * rhob
      t42 = t41 * t9
      t43 = 0.1D1 + t42
      t44 = t43 ** (0.1D1 / 0.3D1)
      t45 = t44 * t43
      t47 = 0.1D1 - 0.1D1 * t42
      t48 = t47 ** (0.1D1 / 0.3D1)
      t49 = t48 * t47
      t50 = t45 + t49 - 0.2D1
      t51 = t41 ** 2
      t52 = t51 ** 2
      t53 = rho ** 2
      t54 = t53 ** 2
      t55 = 0.1D1 / t54
      t56 = t52 * t55
      t60 = (0.1D1 + 0.101077332976287768525D0 * t11) * t38 * t50 * (0.1
     >D1 - 0.1D1 * t56)
      t61 = 0.379955235370239451738D-1 * t60
      t72 = log(0.1D1 + 0.321639589973850701335D2 / (0.13457913714394447
     >7912D2 * t14 + 0.563098414909787598194D1 * t11 + 0.291521471421917
     >737271D1 * t17 + 0.516066464547863440989D0 * t19))
      t77 = (-0.3109070D-1 * (0.1D1 + 0.186690969707574028554D0 * t11) *
     > t72 + t27) * t50 * t56
      t78 = 0.192366105093153631974D1 * t77
      t79 = pi ** 2
      t80 = 0.1D1 / t79
      t81 = t44 ** 2
      t83 = t48 ** 2
      t85 = 0.500000000000000000000D0 * t81 + 0.500000000000000000000D0 
     >* t83
      t86 = t85 ** 2
      t87 = t86 * t85
      t89 = t79 * pi
      t91 = 0.1D1 / t86
      t94 = (t79 * rho) ** (0.1D1 / 0.3D1)
      t95 = 0.1D1 / t94
      t96 = 0.1D1 / t53
      t103 = exp(-0.325889135327092945460D1 * (-t27 + t61 + t78) * t79 /
     > t87)
      t104 = t103 - 0.1D1
      t111 = 0.942319250876317101329D-2 * t89 / t104 * sigma * t91 * t95
     > * t96
      t113 = t79 ** 2
      t114 = t113 * t79
      t115 = t104 ** 2
      t118 = sigma ** 2
      t120 = t86 ** 2
      t122 = t94 ** 2
      t123 = 0.1D1 / t122
      t135 = log(0.1D1 + 0.942319250876317101329D-2 * t89 * sigma * t91 
     >* t95 * t96 * (0.1D1 + t111) / (0.1D1 + t111 + 0.88796557057210344
     >8132D-4 * t114 / t115 * t118 / t120 * t123 * t55))
      t136 = t80 * t87 * t135
      t138 = -t27 + t61 + t78 + 0.306852819440054690583D0 * t136
      t148 = t47 ** 2
      t150 = t43 ** 2
      t164 = (0.1D1 + 0.151426716069344971574D0 * (t148 * sigmaaa + t150
     > * sigmabb - 0.2D1 * t47 * t43 * sigmaab) * t96 * t123 * (0.1D1 /
     >t45 + 0.1D1 / t49)) ** 2
      t165 = t164 ** 2
      t167 = (0.53D0 + 0.87D0 * t51 * t96 + 0.50D0 * t56 + 0.226D1 * t52
     > * t51 / t54 / t53) / t165
      t169 = tau ** 2
      t170 = 0.1D1 / t169
      t180 = 0.621814D29 * t26
      t181 = 0.379955235370239451738D29 * t60
      t182 = 0.192366105093153631974D31 * t77
      t183 = 0.306852819440054690583D30 * t136
      t185 = t8 / rhoa
      t186 = t185 ** (0.1D1 / 0.3D1)
      t189 = t185 ** (0.1D1 / 0.6D1)
      t192 = sqrt(t185)
      t194 = t186 ** 2
      t200 = log(0.1D1 + 0.321639589973850701335D2 / (0.1345791371439444
     >77912D2 * t189 + 0.563098414909787598194D1 * t186 + 0.291521471421
     >917737271D1 * t192 + 0.516066464547863440989D0 * t194))
      t201 = (0.1D1 + 0.186690969707574028554D0 * t186) * t200
      t205 = (t79 * rhoa) ** (0.1D1 / 0.3D1)
      t206 = 0.1D1 / t205
      t208 = rhoa ** 2
      t209 = 0.1D1 / t208
      t212 = exp(0.202642426794280972788D0 * t201 * t79)
      t213 = t212 - 0.1D1
      t219 = 0.149583857013095144140D-1 * t89 / t213 * sigmaaa * t206 * 
     >t209
      t222 = t213 ** 2
      t225 = sigmaaa ** 2
      t226 = t205 ** 2
      t229 = t208 ** 2
      t240 = log(0.1D1 + 0.149583857013095144140D-1 * t89 * sigmaaa * t2
     >06 * t209 * (0.1D1 + t219) / (0.1D1 + t219 + 0.2237533027891409333
     >70D-3 * t114 / t222 * t225 / t226 / t229))
      t241 = t80 * t240
      t244 = tanh(t180 - t181 - t182 - t183 - 0.3109070D29 * t201 + 0.15
     >3426409720027345292D30 * t241)
      t245 = 0.500000000000000000000D0 * t244
      t257 = t8 / rhob
      t258 = t257 ** (0.1D1 / 0.3D1)
      t261 = t257 ** (0.1D1 / 0.6D1)
      t264 = sqrt(t257)
      t266 = t258 ** 2
      t272 = log(0.1D1 + 0.321639589973850701335D2 / (0.1345791371439444
     >77912D2 * t261 + 0.563098414909787598194D1 * t258 + 0.291521471421
     >917737271D1 * t264 + 0.516066464547863440989D0 * t266))
      t273 = (0.1D1 + 0.186690969707574028554D0 * t258) * t272
      t277 = (t79 * rhob) ** (0.1D1 / 0.3D1)
      t278 = 0.1D1 / t277
      t280 = rhob ** 2
      t281 = 0.1D1 / t280
      t284 = exp(0.202642426794280972788D0 * t273 * t79)
      t285 = t284 - 0.1D1
      t291 = 0.149583857013095144140D-1 * t89 / t285 * sigmabb * t278 * 
     >t281
      t294 = t285 ** 2
      t297 = sigmabb ** 2
      t298 = t277 ** 2
      t301 = t280 ** 2
      t312 = log(0.1D1 + 0.149583857013095144140D-1 * t89 * sigmabb * t2
     >78 * t281 * (0.1D1 + t291) / (0.1D1 + t291 + 0.2237533027891409333
     >70D-3 * t114 / t294 * t297 / t298 / t301))
      t313 = t80 * t312
      t316 = tanh(t180 - t181 - t182 - t183 - 0.3109070D29 * t273 + 0.15
     >3426409720027345292D30 * t313)
      t317 = 0.500000000000000000000D0 * t316
      t331 = t138 * (0.1D1 + 0.625000000000000000000D-1 * t167 * t118 * 
     >t96 * t170) - 0.625000000000000000000D-1 * (0.1D1 + t167) * t118 *
     > t96 * t170 * (rhoa * t9 * ((0.500000000000000000000D0 - t245) * t
     >138 + (0.500000000000000000000D0 + t245) * (-0.3109070D-1 * t201 +
     > 0.153426409720027345292D0 * t241)) + rhob * t9 * ((0.500000000000
     >000000000D0 - t317) * t138 + (0.500000000000000000000D0 + t317) *
     >(-0.3109070D-1 * t273 + 0.153426409720027345292D0 * t313)))
      zk(i) = rho * t331 * (0.1D1 + 0.437500000000000000000D-1 * t331 * 
     >t118 * sigma / t53 / rho / t169 / tau)

             endif
           enddo
         endif
       endif

      return
      end

c:TPSSCsubrend
c:TPSSXsubrstart

c    Generated: Mon Jun 20 16:05:44 CEST 2016

      subroutine dftacg_tpssx
     > (name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   tauc,tauo,upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,
     >                   vtauc,vtauo,vupsilonc,vupsilono)
      implicit double precision (a-h,o-z)
      logical fderiv,open
      integer igrad,npt
      character*(*) name
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)
      double precision tauc(*),tauo(*)
      double precision upsilonc(*),upsilono(*)
      double precision zk(*),vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)
      double precision vtauc(*),vtauo(*)
      double precision vupsilonc(*),vupsilono(*)
      include "common/cdft"
      include "common/tapes"
      parameter(tol=1d-12)
      pi=acos(-1d0)
      name='Automatically generated TPSSX'
      igrad=2
       if(open) then
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      sigmabb = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) - 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmabb
      taub = max(0.0D0, 0.500000000000000000000D0 * tauc(i) - 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taub
      t13 = pi ** 2
      t14 = t13 * rhob
      t15 = t14 ** (0.1D1 / 0.3D1)
      t16 = rhob * t15
      t17 = 0.1D1 / pi
      t18 = sigmabb ** 2
      t19 = rhob ** 2
      t20 = 0.1D1 / t19
      t21 = t18 * t20
      t22 = taub ** 2
      t23 = 0.1D1 / t22
      t24 = t21 * t23
      t26 = 0.1D1 + 0.625000000000000000000D-1 * t24
      t27 = t26 ** 2
      t28 = 0.1D1 / t27
      t29 = t23 * t28
      t32 = 0.123456790123456790123D0 + 0.994350000000000000000D-1 * t21
     > * t29
      t33 = t32 * sigmabb
      t34 = t15 ** 2
      t35 = 0.1D1 / t34
      t36 = t35 * t20
      t39 = sigmabb * t35
      t40 = 0.1D1 / sigmabb
      t44 = 0.4D1 * t40 * rhob * taub - 0.1D1
      t45 = t20 * t44
      t48 = 0.126188930057787476312D0 * t39 * t45 - 0.1D1
      t52 = 0.1D1 + 0.504755720231149905248D-1 * t39 * t45 * t48
      t53 = sqrt(t52)
      t54 = 0.1D1 / t53
      t57 = t39 * t20
      t59 = 0.450000000000000000000D0 * t48 * t54 + 0.504755720231149905
     >248D-1 * t57
      t60 = t59 ** 2
      t64 = 0.1D1 / t15 / t14
      t65 = t18 * t64
      t66 = t19 ** 2
      t67 = 0.1D1 / t66
      t68 = t65 * t67
      t71 = sqrt(0.648D3 * t24 + 0.165096362444731334194D3 * t68)
      t76 = 0.1D1 / t13
      t77 = t18 * sigmabb
      t78 = t76 * t77
      t79 = t66 ** 2
      t80 = 0.1D1 / t79
      t83 = 0.757133580346724857871D-1 * t33 * t36 + 0.72098765432098765
     >4321D-1 * t60 - 0.751028806584362139918D-3 * t59 * t71 + 0.1086723
     >17896997724749D-3 * t68 + 0.688754467171997887761D-2 * t24 + 0.148
     >374312789351851852D-4 * t78 * t80
      t85 = 0.1D1 + 0.938662444277523956252D-1 * t57
      t86 = t85 ** 2
      t87 = 0.1D1 / t86
      t90 = 0.1D1 + 0.124378109452736318408D1 * t83 * t87
      t93 = 0.1804D1 - 0.804D0 / t90
      zk(i) = -0.136284044462410474417D1 * t16 * t17 * t93
      t100 = 0.681420222312052372084D0 * t15 * t17 * t93
      t104 = 0.227140074104017457361D0 * rhob * t35 * pi * t93
      t105 = t90 ** 2
      t106 = 0.1D1 / t105
      t107 = t17 * t106
      t109 = 0.1D1 / t19 / rhob
      t110 = t18 * t109
      t113 = t18 ** 2
      t115 = 0.1D1 / t66 / rhob
      t117 = t22 ** 2
      t120 = 0.1D1 / t27 / t26
      t121 = 0.1D1 / t117 * t120
      t129 = 0.1D1 / t34 / t14
      t137 = sigmabb * t129
      t141 = t109 * t44
      t146 = -0.841259533718583175412D-1 * t137 * t45 * t13 - 0.25237786
     >0115574952624D0 * t39 * t141 + 0.504755720231149905248D0 * t36 * t
     >aub
      t151 = t48 / t53 / t52
      t153 = t44 * t48
      t170 = t137 * t20 * t13
      t172 = t39 * t109
      t174 = 0.450000000000000000000D0 * t146 * t54 - 0.2250000000000000
     >00000D0 * t151 * (-0.336503813487433270164D-1 * t137 * t20 * t153
     >* t13 - 0.100951144046229981050D0 * t39 * t141 * t48 + 0.201902288
     >092459962099D0 * t36 * taub * t48 + 0.504755720231149905248D-1 * t
     >39 * t45 * t146) - 0.336503813487433270164D-1 * t170 - 0.100951144
     >046229981050D0 * t172
      t180 = t59 / t71
      t181 = t110 * t23
      t183 = t13 ** 2
      t189 = t18 / t15 / t183 / t19 * t67 * t13
      t191 = t65 * t115
      t208 = t83 / t86 / t85
      t217 = 0.547861858738890107155D0 * t16 * t107 * (0.124378109452736
     >318408D1 * (0.757133580346724857871D-1 * (-0.198870000000000000000
     >D0 * t110 * t29 + 0.248587500000000000000D-1 * t113 * t115 * t121)
     > * sigmabb * t36 - 0.504755720231149905248D-1 * t33 * t129 * t20 *
     > t13 - 0.151426716069344971574D0 * t33 * t35 * t109 + 0.1441975308
     >64197530864D0 * t59 * t174 - 0.751028806584362139918D-3 * t174 * t
     >71 - 0.375514403292181069959D-3 * t180 * (-0.1296D4 * t181 - 0.220
     >128483259641778925D3 * t189 - 0.660385449778925336774D3 * t191) -
     >0.144896423862663632999D-3 * t189 - 0.434689271587990898997D-3 * t
     >191 - 0.137750893434399577552D-1 * t181 - 0.118699450231481481482D
     >-3 * t78 / t79 / rhob) * t87 - 0.248756218905472636816D1 * t208 *
     >(-0.625774962851682637502D-1 * t170 - 0.187732488855504791250D0 *
     >t172))
      vrhoc(i) = vrhoc(i) - t100 - t104 - t217
      vrhoo(i) = vrhoo(i) + t100 + t104 + t217
      t220 = sigmabb * t20
      t235 = t40 * t35
      t236 = 0.1D1 / rhob
      t237 = t236 * taub
      t240 = 0.126188930057787476312D0 * t36 * t44 - 0.50475572023114990
     >5248D0 * t235 * t237
      t255 = 0.450000000000000000000D0 * t240 * t54 - 0.2250000000000000
     >00000D0 * t151 * (0.504755720231149905248D-1 * t36 * t153 - 0.2019
     >02288092459962099D0 * t235 * t237 * t48 + 0.504755720231149905248D
     >-1 * t39 * t45 * t240) + 0.504755720231149905248D-1 * t36
      t260 = t220 * t23
      t262 = sigmabb * t64
      t263 = t262 * t67
      t280 = t16 * t107 * (0.124378109452736318408D1 * (0.75713358034672
     >4857871D-1 * (0.198870000000000000000D0 * t220 * t29 - 0.248587500
     >000000000000D-1 * t77 * t67 * t121) * sigmabb * t36 + 0.7571335803
     >46724857871D-1 * t32 * t35 * t20 + 0.144197530864197530864D0 * t59
     > * t255 - 0.751028806584362139918D-3 * t255 * t71 - 0.375514403292
     >181069959D-3 * t180 * (0.1296D4 * t260 + 0.330192724889462668387D3
     > * t263) + 0.217344635793995449498D-3 * t263 + 0.13775089343439957
     >7552D-1 * t260 + 0.445122938368055555556D-4 * t76 * t18 * t80) * t
     >87 - 0.233498120467045760261D0 * t208 * t36)
      t281 = 0.273930929369445053578D0 * t280
      vsigmacc(i) = vsigmacc(i) - t281
      vsigmaco(i) = vsigmaco(i) + 0.547861858738890107155D0 * t280
      vsigmaoo(i) = vsigmaoo(i) - t281
      t288 = 0.1D1 / t22 / taub
      t302 = t35 * t236
      t312 = 0.227140074104017457361D0 * t302 * t54 - 0.2250000000000000
     >00000D0 * t151 * (0.201902288092459962099D0 * t302 * t48 + 0.25477
     >8337106066873756D-1 * t262 * t141)
      t317 = t21 * t288
      t325 = 0.681420222312052372084D0 * t16 * t17 * t106 * (0.757133580
     >346724857871D-1 * (-0.198870000000000000000D0 * t21 * t288 * t28 +
     > 0.248587500000000000000D-1 * t113 * t67 / t117 / taub * t120) * s
     >igmabb * t36 + 0.144197530864197530864D0 * t59 * t312 - 0.75102880
     >6584362139918D-3 * t312 * t71 + 0.486666666666666666667D0 * t180 *
     > t317 - 0.137750893434399577552D-1 * t317) * t87
      vtauc(i) = vtauc(i) - t325
      vtauo(i) = vtauo(i) + t325

               elseif(rhob.lt.tol) then
             rho = rhoa
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) + 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmaaa
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i) + 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taua
      t13 = pi ** 2
      t14 = t13 * rhoa
      t15 = t14 ** (0.1D1 / 0.3D1)
      t16 = rhoa * t15
      t17 = 0.1D1 / pi
      t18 = sigmaaa ** 2
      t19 = rhoa ** 2
      t20 = 0.1D1 / t19
      t21 = t18 * t20
      t22 = taua ** 2
      t23 = 0.1D1 / t22
      t24 = t21 * t23
      t26 = 0.1D1 + 0.625000000000000000000D-1 * t24
      t27 = t26 ** 2
      t28 = 0.1D1 / t27
      t29 = t23 * t28
      t32 = 0.123456790123456790123D0 + 0.994350000000000000000D-1 * t21
     > * t29
      t33 = t32 * sigmaaa
      t34 = t15 ** 2
      t35 = 0.1D1 / t34
      t36 = t35 * t20
      t39 = sigmaaa * t35
      t40 = 0.1D1 / sigmaaa
      t44 = 0.4D1 * t40 * rhoa * taua - 0.1D1
      t45 = t20 * t44
      t48 = 0.126188930057787476312D0 * t39 * t45 - 0.1D1
      t52 = 0.1D1 + 0.504755720231149905248D-1 * t39 * t45 * t48
      t53 = sqrt(t52)
      t54 = 0.1D1 / t53
      t57 = t39 * t20
      t59 = 0.450000000000000000000D0 * t48 * t54 + 0.504755720231149905
     >248D-1 * t57
      t60 = t59 ** 2
      t64 = 0.1D1 / t15 / t14
      t65 = t18 * t64
      t66 = t19 ** 2
      t67 = 0.1D1 / t66
      t68 = t65 * t67
      t71 = sqrt(0.648D3 * t24 + 0.165096362444731334194D3 * t68)
      t76 = 0.1D1 / t13
      t77 = t18 * sigmaaa
      t78 = t76 * t77
      t79 = t66 ** 2
      t80 = 0.1D1 / t79
      t83 = 0.757133580346724857871D-1 * t33 * t36 + 0.72098765432098765
     >4321D-1 * t60 - 0.751028806584362139918D-3 * t59 * t71 + 0.1086723
     >17896997724749D-3 * t68 + 0.688754467171997887761D-2 * t24 + 0.148
     >374312789351851852D-4 * t78 * t80
      t85 = 0.1D1 + 0.938662444277523956252D-1 * t57
      t86 = t85 ** 2
      t87 = 0.1D1 / t86
      t90 = 0.1D1 + 0.124378109452736318408D1 * t83 * t87
      t93 = 0.1804D1 - 0.804D0 / t90
      zk(i) = -0.136284044462410474417D1 * t16 * t17 * t93
      t100 = 0.681420222312052372084D0 * t15 * t17 * t93
      t104 = 0.227140074104017457361D0 * rhoa * t35 * pi * t93
      t105 = t90 ** 2
      t106 = 0.1D1 / t105
      t107 = t17 * t106
      t109 = 0.1D1 / t19 / rhoa
      t110 = t18 * t109
      t113 = t18 ** 2
      t115 = 0.1D1 / t66 / rhoa
      t117 = t22 ** 2
      t120 = 0.1D1 / t27 / t26
      t121 = 0.1D1 / t117 * t120
      t129 = 0.1D1 / t34 / t14
      t137 = sigmaaa * t129
      t141 = t109 * t44
      t146 = -0.841259533718583175412D-1 * t137 * t45 * t13 - 0.25237786
     >0115574952624D0 * t39 * t141 + 0.504755720231149905248D0 * t36 * t
     >aua
      t151 = t48 / t53 / t52
      t153 = t44 * t48
      t170 = t137 * t20 * t13
      t172 = t39 * t109
      t174 = 0.450000000000000000000D0 * t146 * t54 - 0.2250000000000000
     >00000D0 * t151 * (-0.336503813487433270164D-1 * t137 * t20 * t153
     >* t13 - 0.100951144046229981050D0 * t39 * t141 * t48 + 0.201902288
     >092459962099D0 * t36 * taua * t48 + 0.504755720231149905248D-1 * t
     >39 * t45 * t146) - 0.336503813487433270164D-1 * t170 - 0.100951144
     >046229981050D0 * t172
      t180 = t59 / t71
      t181 = t110 * t23
      t183 = t13 ** 2
      t189 = t18 / t15 / t183 / t19 * t67 * t13
      t191 = t65 * t115
      t208 = t83 / t86 / t85
      t217 = 0.547861858738890107155D0 * t16 * t107 * (0.124378109452736
     >318408D1 * (0.757133580346724857871D-1 * (-0.198870000000000000000
     >D0 * t110 * t29 + 0.248587500000000000000D-1 * t113 * t115 * t121)
     > * sigmaaa * t36 - 0.504755720231149905248D-1 * t33 * t129 * t20 *
     > t13 - 0.151426716069344971574D0 * t33 * t35 * t109 + 0.1441975308
     >64197530864D0 * t59 * t174 - 0.751028806584362139918D-3 * t174 * t
     >71 - 0.375514403292181069959D-3 * t180 * (-0.1296D4 * t181 - 0.220
     >128483259641778925D3 * t189 - 0.660385449778925336774D3 * t191) -
     >0.144896423862663632999D-3 * t189 - 0.434689271587990898997D-3 * t
     >191 - 0.137750893434399577552D-1 * t181 - 0.118699450231481481482D
     >-3 * t78 / t79 / rhoa) * t87 - 0.248756218905472636816D1 * t208 *
     >(-0.625774962851682637502D-1 * t170 - 0.187732488855504791250D0 *
     >t172))
      vrhoc(i) = vrhoc(i) - t100 - t104 - t217
      vrhoo(i) = vrhoo(i) - t100 - t104 - t217
      t220 = sigmaaa * t20
      t235 = t40 * t35
      t236 = 0.1D1 / rhoa
      t237 = t236 * taua
      t240 = 0.126188930057787476312D0 * t36 * t44 - 0.50475572023114990
     >5248D0 * t235 * t237
      t255 = 0.450000000000000000000D0 * t240 * t54 - 0.2250000000000000
     >00000D0 * t151 * (0.504755720231149905248D-1 * t36 * t153 - 0.2019
     >02288092459962099D0 * t235 * t237 * t48 + 0.504755720231149905248D
     >-1 * t39 * t45 * t240) + 0.504755720231149905248D-1 * t36
      t260 = t220 * t23
      t262 = sigmaaa * t64
      t263 = t262 * t67
      t280 = t16 * t107 * (0.124378109452736318408D1 * (0.75713358034672
     >4857871D-1 * (0.198870000000000000000D0 * t220 * t29 - 0.248587500
     >000000000000D-1 * t77 * t67 * t121) * sigmaaa * t36 + 0.7571335803
     >46724857871D-1 * t32 * t35 * t20 + 0.144197530864197530864D0 * t59
     > * t255 - 0.751028806584362139918D-3 * t255 * t71 - 0.375514403292
     >181069959D-3 * t180 * (0.1296D4 * t260 + 0.330192724889462668387D3
     > * t263) + 0.217344635793995449498D-3 * t263 + 0.13775089343439957
     >7552D-1 * t260 + 0.445122938368055555556D-4 * t76 * t18 * t80) * t
     >87 - 0.233498120467045760261D0 * t208 * t36)
      t281 = 0.273930929369445053578D0 * t280
      vsigmacc(i) = vsigmacc(i) - t281
      vsigmaco(i) = vsigmaco(i) - 0.547861858738890107155D0 * t280
      vsigmaoo(i) = vsigmaoo(i) - t281
      t288 = 0.1D1 / t22 / taua
      t302 = t35 * t236
      t312 = 0.227140074104017457361D0 * t302 * t54 - 0.2250000000000000
     >00000D0 * t151 * (0.201902288092459962099D0 * t302 * t48 + 0.25477
     >8337106066873756D-1 * t262 * t141)
      t317 = t21 * t288
      t325 = 0.681420222312052372084D0 * t16 * t17 * t106 * (0.757133580
     >346724857871D-1 * (-0.198870000000000000000D0 * t21 * t288 * t28 +
     > 0.248587500000000000000D-1 * t113 * t67 / t117 / taua * t120) * s
     >igmaaa * t36 + 0.144197530864197530864D0 * t59 * t312 - 0.75102880
     >6584362139918D-3 * t312 * t71 + 0.486666666666666666667D0 * t180 *
     > t317 - 0.137750893434399577552D-1 * t317) * t87
      vtauc(i) = vtauc(i) - t325
      vtauo(i) = vtauo(i) - t325

               else
             rho = rhoa + rhob
      t2 = 0.250000000000000000000D0 * sigmacc(i)
      t4 = 0.250000000000000000000D0 * sigmaoo(i)
      t6 = 0.500000000000000000000D0 * sigmaco(i)
      sigmaaa = max(0.0D0, t2 + t4 + t6)
      sigmaab = t2 - t4
      sigmabb = max(0.0D0, t2 + t4 - t6)
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t11 = 0.500000000000000000000D0 * tauc(i)
      t13 = 0.500000000000000000000D0 * tauo(i)
      taua = max(0.0D0, t11 + t13)
      taub = max(0.0D0, t11 - t13)
      tau = taua + taub
      t16 = pi ** 2
      t17 = t16 * rhoa
      t18 = t17 ** (0.1D1 / 0.3D1)
      t19 = rhoa * t18
      t20 = 0.1D1 / pi
      t21 = sigmaaa ** 2
      t22 = rhoa ** 2
      t23 = 0.1D1 / t22
      t24 = t21 * t23
      t25 = taua ** 2
      t26 = 0.1D1 / t25
      t27 = t24 * t26
      t29 = 0.1D1 + 0.625000000000000000000D-1 * t27
      t30 = t29 ** 2
      t31 = 0.1D1 / t30
      t32 = t26 * t31
      t35 = 0.123456790123456790123D0 + 0.994350000000000000000D-1 * t24
     > * t32
      t36 = sigmaaa * t35
      t37 = t18 ** 2
      t38 = 0.1D1 / t37
      t39 = t38 * t23
      t42 = sigmaaa * t38
      t43 = 0.1D1 / sigmaaa
      t47 = 0.4D1 * t43 * rhoa * taua - 0.1D1
      t48 = t23 * t47
      t51 = 0.126188930057787476312D0 * t42 * t48 - 0.1D1
      t55 = 0.1D1 + 0.504755720231149905248D-1 * t42 * t48 * t51
      t56 = sqrt(t55)
      t57 = 0.1D1 / t56
      t60 = t42 * t23
      t62 = 0.450000000000000000000D0 * t51 * t57 + 0.504755720231149905
     >248D-1 * t60
      t63 = t62 ** 2
      t67 = 0.1D1 / t18 / t17
      t68 = t21 * t67
      t69 = t22 ** 2
      t70 = 0.1D1 / t69
      t71 = t68 * t70
      t74 = sqrt(0.648D3 * t27 + 0.165096362444731334194D3 * t71)
      t79 = 0.1D1 / t16
      t80 = t21 * sigmaaa
      t81 = t79 * t80
      t82 = t69 ** 2
      t83 = 0.1D1 / t82
      t86 = 0.757133580346724857871D-1 * t36 * t39 + 0.72098765432098765
     >4321D-1 * t63 - 0.751028806584362139918D-3 * t62 * t74 + 0.1086723
     >17896997724749D-3 * t71 + 0.688754467171997887761D-2 * t27 + 0.148
     >374312789351851852D-4 * t81 * t83
      t88 = 0.1D1 + 0.938662444277523956252D-1 * t60
      t89 = t88 ** 2
      t90 = 0.1D1 / t89
      t93 = 0.1D1 + 0.124378109452736318408D1 * t86 * t90
      t96 = 0.1804D1 - 0.804D0 / t93
      t100 = t16 * rhob
      t101 = t100 ** (0.1D1 / 0.3D1)
      t102 = rhob * t101
      t103 = sigmabb ** 2
      t104 = rhob ** 2
      t105 = 0.1D1 / t104
      t106 = t103 * t105
      t107 = taub ** 2
      t108 = 0.1D1 / t107
      t109 = t106 * t108
      t111 = 0.1D1 + 0.625000000000000000000D-1 * t109
      t112 = t111 ** 2
      t113 = 0.1D1 / t112
      t114 = t108 * t113
      t117 = 0.123456790123456790123D0 + 0.994350000000000000000D-1 * t1
     >06 * t114
      t118 = t117 * sigmabb
      t119 = t101 ** 2
      t120 = 0.1D1 / t119
      t121 = t120 * t105
      t124 = sigmabb * t120
      t125 = 0.1D1 / sigmabb
      t129 = 0.4D1 * t125 * rhob * taub - 0.1D1
      t130 = t105 * t129
      t133 = 0.126188930057787476312D0 * t124 * t130 - 0.1D1
      t137 = 0.1D1 + 0.504755720231149905248D-1 * t124 * t130 * t133
      t138 = sqrt(t137)
      t139 = 0.1D1 / t138
      t142 = t124 * t105
      t144 = 0.450000000000000000000D0 * t133 * t139 + 0.504755720231149
     >905248D-1 * t142
      t145 = t144 ** 2
      t149 = 0.1D1 / t101 / t100
      t150 = t103 * t149
      t151 = t104 ** 2
      t152 = 0.1D1 / t151
      t153 = t150 * t152
      t156 = sqrt(0.648D3 * t109 + 0.165096362444731334194D3 * t153)
      t161 = t103 * sigmabb
      t162 = t79 * t161
      t163 = t151 ** 2
      t164 = 0.1D1 / t163
      t167 = 0.757133580346724857871D-1 * t118 * t121 + 0.72098765432098
     >7654321D-1 * t145 - 0.751028806584362139918D-3 * t144 * t156 + 0.1
     >08672317896997724749D-3 * t153 + 0.688754467171997887761D-2 * t109
     > + 0.148374312789351851852D-4 * t162 * t164
      t169 = 0.1D1 + 0.938662444277523956252D-1 * t142
      t170 = t169 ** 2
      t171 = 0.1D1 / t170
      t174 = 0.1D1 + 0.124378109452736318408D1 * t167 * t171
      t177 = 0.1804D1 - 0.804D0 / t174
      zk(i) = -0.136284044462410474417D1 * t19 * t20 * t96 - 0.136284044
     >462410474417D1 * t102 * t20 * t177
      t184 = 0.681420222312052372084D0 * t18 * t20 * t96
      t188 = 0.227140074104017457361D0 * rhoa * t38 * pi * t96
      t189 = t93 ** 2
      t190 = 0.1D1 / t189
      t191 = t20 * t190
      t193 = 0.1D1 / t22 / rhoa
      t194 = t21 * t193
      t197 = t21 ** 2
      t199 = 0.1D1 / t69 / rhoa
      t201 = t25 ** 2
      t204 = 0.1D1 / t30 / t29
      t205 = 0.1D1 / t201 * t204
      t213 = 0.1D1 / t37 / t17
      t221 = sigmaaa * t213
      t225 = t193 * t47
      t230 = -0.841259533718583175412D-1 * t221 * t48 * t16 - 0.25237786
     >0115574952624D0 * t42 * t225 + 0.504755720231149905248D0 * t39 * t
     >aua
      t235 = t51 / t56 / t55
      t237 = t47 * t51
      t254 = t221 * t23 * t16
      t256 = t42 * t193
      t258 = 0.450000000000000000000D0 * t230 * t57 - 0.2250000000000000
     >00000D0 * t235 * (-0.336503813487433270164D-1 * t221 * t23 * t237
     >* t16 - 0.100951144046229981050D0 * t42 * t225 * t51 + 0.201902288
     >092459962099D0 * t39 * taua * t51 + 0.504755720231149905248D-1 * t
     >42 * t48 * t230) - 0.336503813487433270164D-1 * t254 - 0.100951144
     >046229981050D0 * t256
      t264 = t62 / t74
      t265 = t194 * t26
      t267 = t16 ** 2
      t273 = t21 / t18 / t267 / t22 * t70 * t16
      t275 = t68 * t199
      t292 = t86 / t89 / t88
      t301 = 0.547861858738890107155D0 * t19 * t191 * (0.124378109452736
     >318408D1 * (0.757133580346724857871D-1 * (-0.198870000000000000000
     >D0 * t194 * t32 + 0.248587500000000000000D-1 * t197 * t199 * t205)
     > * sigmaaa * t39 - 0.504755720231149905248D-1 * t36 * t213 * t23 *
     > t16 - 0.151426716069344971574D0 * t36 * t38 * t193 + 0.1441975308
     >64197530864D0 * t62 * t258 - 0.751028806584362139918D-3 * t258 * t
     >74 - 0.375514403292181069959D-3 * t264 * (-0.1296D4 * t265 - 0.220
     >128483259641778925D3 * t273 - 0.660385449778925336774D3 * t275) -
     >0.144896423862663632999D-3 * t273 - 0.434689271587990898997D-3 * t
     >275 - 0.137750893434399577552D-1 * t265 - 0.118699450231481481482D
     >-3 * t81 / t82 / rhoa) * t90 - 0.248756218905472636816D1 * t292 *
     >(-0.625774962851682637502D-1 * t254 - 0.187732488855504791250D0 *
     >t256))
      t304 = 0.681420222312052372084D0 * t101 * t20 * t177
      t308 = 0.227140074104017457361D0 * rhob * t120 * pi * t177
      t309 = t174 ** 2
      t310 = 0.1D1 / t309
      t311 = t20 * t310
      t313 = 0.1D1 / t104 / rhob
      t314 = t103 * t313
      t317 = t103 ** 2
      t319 = 0.1D1 / t151 / rhob
      t321 = t107 ** 2
      t324 = 0.1D1 / t112 / t111
      t325 = 0.1D1 / t321 * t324
      t333 = 0.1D1 / t119 / t100
      t341 = sigmabb * t333
      t345 = t313 * t129
      t350 = -0.841259533718583175412D-1 * t341 * t130 * t16 - 0.2523778
     >60115574952624D0 * t124 * t345 + 0.504755720231149905248D0 * t121
     >* taub
      t355 = t133 / t138 / t137
      t357 = t129 * t133
      t374 = t341 * t105 * t16
      t376 = t124 * t313
      t378 = 0.450000000000000000000D0 * t350 * t139 - 0.225000000000000
     >000000D0 * t355 * (-0.336503813487433270164D-1 * t341 * t105 * t35
     >7 * t16 - 0.100951144046229981050D0 * t124 * t345 * t133 + 0.20190
     >2288092459962099D0 * t121 * taub * t133 + 0.504755720231149905248D
     >-1 * t124 * t130 * t350) - 0.336503813487433270164D-1 * t374 - 0.1
     >00951144046229981050D0 * t376
      t384 = t144 / t156
      t385 = t314 * t108
      t392 = t103 / t101 / t267 / t104 * t152 * t16
      t394 = t150 * t319
      t411 = t167 / t170 / t169
      t420 = 0.547861858738890107155D0 * t102 * t311 * (0.12437810945273
     >6318408D1 * (0.757133580346724857871D-1 * (-0.19887000000000000000
     >0D0 * t314 * t114 + 0.248587500000000000000D-1 * t317 * t319 * t32
     >5) * sigmabb * t121 - 0.504755720231149905248D-1 * t118 * t333 * t
     >105 * t16 - 0.151426716069344971574D0 * t118 * t120 * t313 + 0.144
     >197530864197530864D0 * t144 * t378 - 0.751028806584362139918D-3 *
     >t378 * t156 - 0.375514403292181069959D-3 * t384 * (-0.1296D4 * t38
     >5 - 0.220128483259641778925D3 * t392 - 0.660385449778925336774D3 *
     > t394) - 0.144896423862663632999D-3 * t392 - 0.4346892715879908989
     >97D-3 * t394 - 0.137750893434399577552D-1 * t385 - 0.1186994502314
     >81481482D-3 * t162 / t163 / rhob) * t171 - 0.248756218905472636816
     >D1 * t411 * (-0.625774962851682637502D-1 * t374 - 0.18773248885550
     >4791250D0 * t376))
      vrhoc(i) = vrhoc(i) - t184 - t188 - t301 - t304 - t308 - t420
      vrhoo(i) = vrhoo(i) - t184 - t188 - t301 + t304 + t308 + t420
      t423 = sigmaaa * t23
      t438 = t43 * t38
      t439 = 0.1D1 / rhoa
      t440 = t439 * taua
      t443 = 0.126188930057787476312D0 * t39 * t47 - 0.50475572023114990
     >5248D0 * t438 * t440
      t458 = 0.450000000000000000000D0 * t443 * t57 - 0.2250000000000000
     >00000D0 * t235 * (0.504755720231149905248D-1 * t39 * t237 - 0.2019
     >02288092459962099D0 * t438 * t440 * t51 + 0.504755720231149905248D
     >-1 * t42 * t48 * t443) + 0.504755720231149905248D-1 * t39
      t463 = t423 * t26
      t465 = sigmaaa * t67
      t466 = t465 * t70
      t483 = t19 * t191 * (0.124378109452736318408D1 * (0.75713358034672
     >4857871D-1 * (0.198870000000000000000D0 * t423 * t32 - 0.248587500
     >000000000000D-1 * t80 * t70 * t205) * sigmaaa * t39 + 0.7571335803
     >46724857871D-1 * t35 * t38 * t23 + 0.144197530864197530864D0 * t62
     > * t458 - 0.751028806584362139918D-3 * t458 * t74 - 0.375514403292
     >181069959D-3 * t264 * (0.1296D4 * t463 + 0.330192724889462668387D3
     > * t466) + 0.217344635793995449498D-3 * t466 + 0.13775089343439957
     >7552D-1 * t463 + 0.445122938368055555556D-4 * t79 * t21 * t83) * t
     >90 - 0.233498120467045760261D0 * t292 * t39)
      t484 = 0.273930929369445053578D0 * t483
      t485 = sigmabb * t105
      t500 = t125 * t120
      t501 = 0.1D1 / rhob
      t502 = t501 * taub
      t505 = 0.126188930057787476312D0 * t121 * t129 - 0.504755720231149
     >905248D0 * t500 * t502
      t520 = 0.450000000000000000000D0 * t505 * t139 - 0.225000000000000
     >000000D0 * t355 * (0.504755720231149905248D-1 * t121 * t357 - 0.20
     >1902288092459962099D0 * t500 * t502 * t133 + 0.5047557202311499052
     >48D-1 * t124 * t130 * t505) + 0.504755720231149905248D-1 * t121
      t525 = t485 * t108
      t527 = sigmabb * t149
      t528 = t527 * t152
      t545 = t102 * t311 * (0.124378109452736318408D1 * (0.7571335803467
     >24857871D-1 * (0.198870000000000000000D0 * t485 * t114 - 0.2485875
     >00000000000000D-1 * t161 * t152 * t325) * sigmabb * t121 + 0.75713
     >3580346724857871D-1 * t117 * t120 * t105 + 0.144197530864197530864
     >D0 * t144 * t520 - 0.751028806584362139918D-3 * t520 * t156 - 0.37
     >5514403292181069959D-3 * t384 * (0.1296D4 * t525 + 0.3301927248894
     >62668387D3 * t528) + 0.217344635793995449498D-3 * t528 + 0.1377508
     >93434399577552D-1 * t525 + 0.445122938368055555556D-4 * t79 * t103
     > * t164) * t171 - 0.233498120467045760261D0 * t411 * t121)
      t546 = 0.273930929369445053578D0 * t545
      vsigmacc(i) = vsigmacc(i) - t484 - t546
      vsigmaco(i) = vsigmaco(i) - 0.547861858738890107155D0 * t483 + 0.5
     >47861858738890107155D0 * t545
      vsigmaoo(i) = vsigmaoo(i) - t484 - t546
      t554 = 0.1D1 / t25 / taua
      t568 = t38 * t439
      t578 = 0.227140074104017457361D0 * t568 * t57 - 0.2250000000000000
     >00000D0 * t235 * (0.201902288092459962099D0 * t568 * t51 + 0.25477
     >8337106066873756D-1 * t465 * t225)
      t583 = t24 * t554
      t591 = 0.681420222312052372084D0 * t19 * t20 * t190 * (0.757133580
     >346724857871D-1 * (-0.198870000000000000000D0 * t24 * t554 * t31 +
     > 0.248587500000000000000D-1 * t197 * t70 / t201 / taua * t204) * s
     >igmaaa * t39 + 0.144197530864197530864D0 * t62 * t578 - 0.75102880
     >6584362139918D-3 * t578 * t74 + 0.486666666666666666667D0 * t264 *
     > t583 - 0.137750893434399577552D-1 * t583) * t90
      t594 = 0.1D1 / t107 / taub
      t608 = t120 * t501
      t618 = 0.227140074104017457361D0 * t608 * t139 - 0.225000000000000
     >000000D0 * t355 * (0.201902288092459962099D0 * t608 * t133 + 0.254
     >778337106066873756D-1 * t527 * t345)
      t623 = t106 * t594
      t631 = 0.681420222312052372084D0 * t102 * t20 * t310 * (0.75713358
     >0346724857871D-1 * (-0.198870000000000000000D0 * t106 * t594 * t11
     >3 + 0.248587500000000000000D-1 * t317 * t152 / t321 / taub * t324)
     > * sigmabb * t121 + 0.144197530864197530864D0 * t144 * t618 - 0.75
     >1028806584362139918D-3 * t618 * t156 + 0.486666666666666666667D0 *
     > t384 * t623 - 0.137750893434399577552D-1 * t623) * t171
      vtauc(i) = vtauc(i) - t591 - t631
      vtauo(i) = vtauo(i) - t591 + t631

               endif
             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

               if(rhoa.lt.tol) then
             rho = rhob
      sigmabb = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) - 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmabb
      taub = max(0.0D0, 0.500000000000000000000D0 * tauc(i) - 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taub
      t13 = pi ** 2
      t14 = t13 * rhob
      t15 = t14 ** (0.1D1 / 0.3D1)
      t18 = sigmabb ** 2
      t19 = rhob ** 2
      t20 = 0.1D1 / t19
      t21 = t18 * t20
      t22 = taub ** 2
      t23 = 0.1D1 / t22
      t24 = t21 * t23
      t27 = (0.1D1 + 0.625000000000000000000D-1 * t24) ** 2
      t34 = t15 ** 2
      t35 = 0.1D1 / t34
      t39 = sigmabb * t35
      t45 = t20 * (0.4D1 / sigmabb * rhob * taub - 0.1D1)
      t48 = 0.126188930057787476312D0 * t39 * t45 - 0.1D1
      t53 = sqrt(0.1D1 + 0.504755720231149905248D-1 * t39 * t45 * t48)
      t57 = t39 * t20
      t59 = 0.450000000000000000000D0 * t48 / t53 + 0.504755720231149905
     >248D-1 * t57
      t60 = t59 ** 2
      t66 = t19 ** 2
      t68 = t18 / t15 / t14 / t66
      t71 = sqrt(0.648D3 * t24 + 0.165096362444731334194D3 * t68)
      t79 = t66 ** 2
      t86 = (0.1D1 + 0.938662444277523956252D-1 * t57) ** 2
      zk(i) = -0.136284044462410474417D1 * rhob * t15 / pi * (0.1804D1 -
     > 0.804D0 / (0.1D1 + 0.124378109452736318408D1 * (0.757133580346724
     >857871D-1 * (0.123456790123456790123D0 + 0.994350000000000000000D-
     >1 * t21 * t23 / t27) * sigmabb * t35 * t20 + 0.7209876543209876543
     >21D-1 * t60 - 0.751028806584362139918D-3 * t59 * t71 + 0.108672317
     >896997724749D-3 * t68 + 0.688754467171997887761D-2 * t24 + 0.14837
     >4312789351851852D-4 / t13 * t18 * sigmabb / t79) / t86))

               elseif(rhob.lt.tol) then
             rho = rhoa
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i) + 0.25
     >0000000000000000000D0 * sigmaoo(i) + 0.500000000000000000000D0 * s
     >igmaco(i))
      sigma = sigmaaa
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i) + 0.50000000
     >0000000000000D0 * tauo(i))
      tau = taua
      t13 = pi ** 2
      t14 = t13 * rhoa
      t15 = t14 ** (0.1D1 / 0.3D1)
      t18 = sigmaaa ** 2
      t19 = rhoa ** 2
      t20 = 0.1D1 / t19
      t21 = t18 * t20
      t22 = taua ** 2
      t23 = 0.1D1 / t22
      t24 = t21 * t23
      t27 = (0.1D1 + 0.625000000000000000000D-1 * t24) ** 2
      t34 = t15 ** 2
      t35 = 0.1D1 / t34
      t39 = sigmaaa * t35
      t45 = t20 * (0.4D1 / sigmaaa * rhoa * taua - 0.1D1)
      t48 = 0.126188930057787476312D0 * t39 * t45 - 0.1D1
      t53 = sqrt(0.1D1 + 0.504755720231149905248D-1 * t39 * t45 * t48)
      t57 = t39 * t20
      t59 = 0.450000000000000000000D0 * t48 / t53 + 0.504755720231149905
     >248D-1 * t57
      t60 = t59 ** 2
      t66 = t19 ** 2
      t68 = t18 / t15 / t14 / t66
      t71 = sqrt(0.648D3 * t24 + 0.165096362444731334194D3 * t68)
      t79 = t66 ** 2
      t86 = (0.1D1 + 0.938662444277523956252D-1 * t57) ** 2
      zk(i) = -0.136284044462410474417D1 * rhoa * t15 / pi * (0.1804D1 -
     > 0.804D0 / (0.1D1 + 0.124378109452736318408D1 * (0.757133580346724
     >857871D-1 * (0.123456790123456790123D0 + 0.994350000000000000000D-
     >1 * t21 * t23 / t27) * sigmaaa * t35 * t20 + 0.7209876543209876543
     >21D-1 * t60 - 0.751028806584362139918D-3 * t59 * t71 + 0.108672317
     >896997724749D-3 * t68 + 0.688754467171997887761D-2 * t24 + 0.14837
     >4312789351851852D-4 / t13 * t18 * sigmaaa / t79) / t86))

               else
             rho = rhoa + rhob
      t2 = 0.250000000000000000000D0 * sigmacc(i)
      t4 = 0.250000000000000000000D0 * sigmaoo(i)
      t6 = 0.500000000000000000000D0 * sigmaco(i)
      sigmaaa = max(0.0D0, t2 + t4 + t6)
      sigmaab = t2 - t4
      sigmabb = max(0.0D0, t2 + t4 - t6)
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      t11 = 0.500000000000000000000D0 * tauc(i)
      t13 = 0.500000000000000000000D0 * tauo(i)
      taua = max(0.0D0, t11 + t13)
      taub = max(0.0D0, t11 - t13)
      tau = taua + taub
      t16 = pi ** 2
      t17 = t16 * rhoa
      t18 = t17 ** (0.1D1 / 0.3D1)
      t20 = 0.1D1 / pi
      t21 = sigmaaa ** 2
      t22 = rhoa ** 2
      t23 = 0.1D1 / t22
      t24 = t21 * t23
      t25 = taua ** 2
      t26 = 0.1D1 / t25
      t27 = t24 * t26
      t30 = (0.1D1 + 0.625000000000000000000D-1 * t27) ** 2
      t37 = t18 ** 2
      t38 = 0.1D1 / t37
      t42 = sigmaaa * t38
      t48 = t23 * (0.4D1 / sigmaaa * rhoa * taua - 0.1D1)
      t51 = 0.126188930057787476312D0 * t42 * t48 - 0.1D1
      t56 = sqrt(0.1D1 + 0.504755720231149905248D-1 * t42 * t48 * t51)
      t60 = t42 * t23
      t62 = 0.450000000000000000000D0 * t51 / t56 + 0.504755720231149905
     >248D-1 * t60
      t63 = t62 ** 2
      t69 = t22 ** 2
      t71 = t21 / t18 / t17 / t69
      t74 = sqrt(0.648D3 * t27 + 0.165096362444731334194D3 * t71)
      t79 = 0.1D1 / t16
      t82 = t69 ** 2
      t89 = (0.1D1 + 0.938662444277523956252D-1 * t60) ** 2
      t100 = t16 * rhob
      t101 = t100 ** (0.1D1 / 0.3D1)
      t103 = sigmabb ** 2
      t104 = rhob ** 2
      t105 = 0.1D1 / t104
      t106 = t103 * t105
      t107 = taub ** 2
      t108 = 0.1D1 / t107
      t109 = t106 * t108
      t112 = (0.1D1 + 0.625000000000000000000D-1 * t109) ** 2
      t119 = t101 ** 2
      t120 = 0.1D1 / t119
      t124 = sigmabb * t120
      t130 = t105 * (0.4D1 / sigmabb * rhob * taub - 0.1D1)
      t133 = 0.126188930057787476312D0 * t124 * t130 - 0.1D1
      t138 = sqrt(0.1D1 + 0.504755720231149905248D-1 * t124 * t130 * t13
     >3)
      t142 = t124 * t105
      t144 = 0.450000000000000000000D0 * t133 / t138 + 0.504755720231149
     >905248D-1 * t142
      t145 = t144 ** 2
      t151 = t104 ** 2
      t153 = t103 / t101 / t100 / t151
      t156 = sqrt(0.648D3 * t109 + 0.165096362444731334194D3 * t153)
      t163 = t151 ** 2
      t170 = (0.1D1 + 0.938662444277523956252D-1 * t142) ** 2
      zk(i) = -0.136284044462410474417D1 * rhoa * t18 * t20 * (0.1804D1 
     >- 0.804D0 / (0.1D1 + 0.124378109452736318408D1 * (0.75713358034672
     >4857871D-1 * (0.123456790123456790123D0 + 0.994350000000000000000D
     >-1 * t24 * t26 / t30) * sigmaaa * t38 * t23 + 0.720987654320987654
     >321D-1 * t63 - 0.751028806584362139918D-3 * t62 * t74 + 0.10867231
     >7896997724749D-3 * t71 + 0.688754467171997887761D-2 * t27 + 0.1483
     >74312789351851852D-4 * t79 * t21 * sigmaaa / t82) / t89)) - 0.1362
     >84044462410474417D1 * rhob * t101 * t20 * (0.1804D1 - 0.804D0 / (0
     >.1D1 + 0.124378109452736318408D1 * (0.757133580346724857871D-1 * (
     >0.123456790123456790123D0 + 0.994350000000000000000D-1 * t106 * t1
     >08 / t112) * sigmabb * t120 * t105 + 0.720987654320987654321D-1 *
     >t145 - 0.751028806584362139918D-3 * t144 * t156 + 0.10867231789699
     >7724749D-3 * t153 + 0.688754467171997887761D-2 * t109 + 0.14837431
     >2789351851852D-4 * t79 * t103 * sigmabb / t163) / t170))

               endif
             endif
           enddo
         endif
       else
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i))
      sigmaab = sigmaaa
      sigmabb = sigmaab
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i))
      taub = taua
      tau = taua + taub
      t8 = pi ** 2
      t9 = t8 * rhoa
      t10 = t9 ** (0.1D1 / 0.3D1)
      t11 = rhoa * t10
      t12 = 0.1D1 / pi
      t13 = sigmaaa ** 2
      t14 = rhoa ** 2
      t15 = 0.1D1 / t14
      t16 = t13 * t15
      t17 = taua ** 2
      t18 = 0.1D1 / t17
      t19 = t16 * t18
      t21 = 0.1D1 + 0.625000000000000000000D-1 * t19
      t22 = t21 ** 2
      t23 = 0.1D1 / t22
      t24 = t18 * t23
      t27 = 0.123456790123456790123D0 + 0.994350000000000000000D-1 * t16
     > * t24
      t28 = t27 * sigmaaa
      t29 = t10 ** 2
      t30 = 0.1D1 / t29
      t31 = t30 * t15
      t34 = sigmaaa * t30
      t35 = 0.1D1 / sigmaaa
      t39 = 0.4D1 * t35 * rhoa * taua - 0.1D1
      t40 = t15 * t39
      t43 = 0.126188930057787476312D0 * t34 * t40 - 0.1D1
      t47 = 0.1D1 + 0.504755720231149905248D-1 * t34 * t40 * t43
      t48 = sqrt(t47)
      t49 = 0.1D1 / t48
      t52 = t34 * t15
      t54 = 0.450000000000000000000D0 * t43 * t49 + 0.504755720231149905
     >248D-1 * t52
      t55 = t54 ** 2
      t59 = 0.1D1 / t10 / t9
      t60 = t13 * t59
      t61 = t14 ** 2
      t62 = 0.1D1 / t61
      t63 = t60 * t62
      t66 = sqrt(0.648D3 * t19 + 0.165096362444731334194D3 * t63)
      t71 = 0.1D1 / t8
      t72 = t13 * sigmaaa
      t73 = t71 * t72
      t74 = t61 ** 2
      t75 = 0.1D1 / t74
      t78 = 0.757133580346724857871D-1 * t28 * t31 + 0.72098765432098765
     >4321D-1 * t55 - 0.751028806584362139918D-3 * t54 * t66 + 0.1086723
     >17896997724749D-3 * t63 + 0.688754467171997887761D-2 * t19 + 0.148
     >374312789351851852D-4 * t73 * t75
      t80 = 0.1D1 + 0.938662444277523956252D-1 * t52
      t81 = t80 ** 2
      t82 = 0.1D1 / t81
      t85 = 0.1D1 + 0.124378109452736318408D1 * t78 * t82
      t88 = 0.1804D1 - 0.804D0 / t85
      t92 = t8 * rhob
      t93 = t92 ** (0.1D1 / 0.3D1)
      t94 = rhob * t93
      t95 = sigmabb ** 2
      t96 = rhob ** 2
      t97 = 0.1D1 / t96
      t98 = t95 * t97
      t99 = taub ** 2
      t100 = 0.1D1 / t99
      t101 = t98 * t100
      t103 = 0.1D1 + 0.625000000000000000000D-1 * t101
      t104 = t103 ** 2
      t105 = 0.1D1 / t104
      t106 = t100 * t105
      t109 = 0.123456790123456790123D0 + 0.994350000000000000000D-1 * t9
     >8 * t106
      t110 = t109 * sigmabb
      t111 = t93 ** 2
      t112 = 0.1D1 / t111
      t113 = t112 * t97
      t116 = sigmabb * t112
      t117 = 0.1D1 / sigmabb
      t121 = 0.4D1 * t117 * rhob * taub - 0.1D1
      t122 = t97 * t121
      t125 = 0.126188930057787476312D0 * t116 * t122 - 0.1D1
      t129 = 0.1D1 + 0.504755720231149905248D-1 * t116 * t122 * t125
      t130 = sqrt(t129)
      t131 = 0.1D1 / t130
      t134 = t116 * t97
      t136 = 0.450000000000000000000D0 * t125 * t131 + 0.504755720231149
     >905248D-1 * t134
      t137 = t136 ** 2
      t141 = 0.1D1 / t93 / t92
      t142 = t95 * t141
      t143 = t96 ** 2
      t144 = 0.1D1 / t143
      t145 = t142 * t144
      t148 = sqrt(0.648D3 * t101 + 0.165096362444731334194D3 * t145)
      t153 = t95 * sigmabb
      t154 = t71 * t153
      t155 = t143 ** 2
      t156 = 0.1D1 / t155
      t159 = 0.757133580346724857871D-1 * t110 * t113 + 0.72098765432098
     >7654321D-1 * t137 - 0.751028806584362139918D-3 * t136 * t148 + 0.1
     >08672317896997724749D-3 * t145 + 0.688754467171997887761D-2 * t101
     > + 0.148374312789351851852D-4 * t154 * t156
      t161 = 0.1D1 + 0.938662444277523956252D-1 * t134
      t162 = t161 ** 2
      t163 = 0.1D1 / t162
      t166 = 0.1D1 + 0.124378109452736318408D1 * t159 * t163
      t169 = 0.1804D1 - 0.804D0 / t166
      zk(i) = -0.136284044462410474417D1 * t11 * t12 * t88 - 0.136284044
     >462410474417D1 * t94 * t12 * t169
      t181 = t85 ** 2
      t182 = 0.1D1 / t181
      t183 = t12 * t182
      t185 = 0.1D1 / t14 / rhoa
      t186 = t13 * t185
      t189 = t13 ** 2
      t191 = 0.1D1 / t61 / rhoa
      t193 = t17 ** 2
      t196 = 0.1D1 / t22 / t21
      t197 = 0.1D1 / t193 * t196
      t205 = 0.1D1 / t29 / t9
      t213 = sigmaaa * t205
      t217 = t185 * t39
      t222 = -0.841259533718583175412D-1 * t213 * t40 * t8 - 0.252377860
     >115574952624D0 * t34 * t217 + 0.504755720231149905248D0 * t31 * ta
     >ua
      t227 = t43 / t48 / t47
      t229 = t39 * t43
      t246 = t213 * t15 * t8
      t248 = t34 * t185
      t250 = 0.450000000000000000000D0 * t222 * t49 - 0.2250000000000000
     >00000D0 * t227 * (-0.336503813487433270164D-1 * t213 * t15 * t229
     >* t8 - 0.100951144046229981050D0 * t34 * t217 * t43 + 0.2019022880
     >92459962099D0 * t31 * taua * t43 + 0.504755720231149905248D-1 * t3
     >4 * t40 * t222) - 0.336503813487433270164D-1 * t246 - 0.1009511440
     >46229981050D0 * t248
      t256 = t54 / t66
      t257 = t186 * t18
      t259 = t8 ** 2
      t265 = t13 / t10 / t259 / t14 * t62 * t8
      t267 = t60 * t191
      t284 = t78 / t81 / t80
      t301 = t166 ** 2
      t302 = 0.1D1 / t301
      t303 = t12 * t302
      t305 = 0.1D1 / t96 / rhob
      t306 = t95 * t305
      t309 = t95 ** 2
      t311 = 0.1D1 / t143 / rhob
      t313 = t99 ** 2
      t316 = 0.1D1 / t104 / t103
      t317 = 0.1D1 / t313 * t316
      t325 = 0.1D1 / t111 / t92
      t333 = sigmabb * t325
      t337 = t305 * t121
      t342 = -0.841259533718583175412D-1 * t333 * t122 * t8 - 0.25237786
     >0115574952624D0 * t116 * t337 + 0.504755720231149905248D0 * t113 *
     > taub
      t347 = t125 / t130 / t129
      t349 = t121 * t125
      t366 = t333 * t97 * t8
      t368 = t116 * t305
      t370 = 0.450000000000000000000D0 * t342 * t131 - 0.225000000000000
     >000000D0 * t347 * (-0.336503813487433270164D-1 * t333 * t97 * t349
     > * t8 - 0.100951144046229981050D0 * t116 * t337 * t125 + 0.2019022
     >88092459962099D0 * t113 * taub * t125 + 0.504755720231149905248D-1
     > * t116 * t122 * t342) - 0.336503813487433270164D-1 * t366 - 0.100
     >951144046229981050D0 * t368
      t376 = t136 / t148
      t377 = t306 * t100
      t384 = t95 / t93 / t259 / t96 * t144 * t8
      t386 = t142 * t311
      t403 = t159 / t162 / t161
      vrhoc(i) = vrhoc(i) - 0.681420222312052372084D0 * t10 * t12 * t88 
     >- 0.227140074104017457361D0 * rhoa * t30 * pi * t88 - 0.5478618587
     >38890107155D0 * t11 * t183 * (0.124378109452736318408D1 * (0.75713
     >3580346724857871D-1 * (-0.198870000000000000000D0 * t186 * t24 + 0
     >.248587500000000000000D-1 * t189 * t191 * t197) * sigmaaa * t31 -
     >0.504755720231149905248D-1 * t28 * t205 * t15 * t8 - 0.15142671606
     >9344971574D0 * t28 * t30 * t185 + 0.144197530864197530864D0 * t54
     >* t250 - 0.751028806584362139918D-3 * t250 * t66 - 0.3755144032921
     >81069959D-3 * t256 * (-0.1296D4 * t257 - 0.220128483259641778925D3
     > * t265 - 0.660385449778925336774D3 * t267) - 0.144896423862663632
     >999D-3 * t265 - 0.434689271587990898997D-3 * t267 - 0.137750893434
     >399577552D-1 * t257 - 0.118699450231481481482D-3 * t73 / t74 / rho
     >a) * t82 - 0.248756218905472636816D1 * t284 * (-0.6257749628516826
     >37502D-1 * t246 - 0.187732488855504791250D0 * t248)) - 0.681420222
     >312052372084D0 * t93 * t12 * t169 - 0.227140074104017457361D0 * rh
     >ob * t112 * pi * t169 - 0.547861858738890107155D0 * t94 * t303 * (
     >0.124378109452736318408D1 * (0.757133580346724857871D-1 * (-0.1988
     >70000000000000000D0 * t306 * t106 + 0.248587500000000000000D-1 * t
     >309 * t311 * t317) * sigmabb * t113 - 0.504755720231149905248D-1 *
     > t110 * t325 * t97 * t8 - 0.151426716069344971574D0 * t110 * t112
     >* t305 + 0.144197530864197530864D0 * t136 * t370 - 0.7510288065843
     >62139918D-3 * t370 * t148 - 0.375514403292181069959D-3 * t376 * (-
     >0.1296D4 * t377 - 0.220128483259641778925D3 * t384 - 0.66038544977
     >8925336774D3 * t386) - 0.144896423862663632999D-3 * t384 - 0.43468
     >9271587990898997D-3 * t386 - 0.137750893434399577552D-1 * t377 - 0
     >.118699450231481481482D-3 * t154 / t155 / rhob) * t163 - 0.2487562
     >18905472636816D1 * t403 * (-0.625774962851682637502D-1 * t366 - 0.
     >187732488855504791250D0 * t368))
      t414 = sigmaaa * t15
      t429 = t35 * t30
      t430 = 0.1D1 / rhoa
      t431 = t430 * taua
      t434 = 0.126188930057787476312D0 * t31 * t39 - 0.50475572023114990
     >5248D0 * t429 * t431
      t449 = 0.450000000000000000000D0 * t434 * t49 - 0.2250000000000000
     >00000D0 * t227 * (0.504755720231149905248D-1 * t31 * t229 - 0.2019
     >02288092459962099D0 * t429 * t431 * t43 + 0.504755720231149905248D
     >-1 * t34 * t40 * t434) + 0.504755720231149905248D-1 * t31
      t454 = t414 * t18
      t456 = sigmaaa * t59
      t457 = t456 * t62
      t476 = sigmabb * t97
      t491 = t117 * t112
      t492 = 0.1D1 / rhob
      t493 = t492 * taub
      t496 = 0.126188930057787476312D0 * t113 * t121 - 0.504755720231149
     >905248D0 * t491 * t493
      t511 = 0.450000000000000000000D0 * t496 * t131 - 0.225000000000000
     >000000D0 * t347 * (0.504755720231149905248D-1 * t113 * t349 - 0.20
     >1902288092459962099D0 * t491 * t493 * t125 + 0.5047557202311499052
     >48D-1 * t116 * t122 * t496) + 0.504755720231149905248D-1 * t113
      t516 = t476 * t100
      t518 = sigmabb * t141
      t519 = t518 * t144
      vsigmacc(i) = vsigmacc(i) - 0.273930929369445053578D0 * t11 * t183
     > * (0.124378109452736318408D1 * (0.757133580346724857871D-1 * (0.1
     >98870000000000000000D0 * t414 * t24 - 0.248587500000000000000D-1 *
     > t72 * t62 * t197) * sigmaaa * t31 + 0.757133580346724857871D-1 *
     >t27 * t30 * t15 + 0.144197530864197530864D0 * t54 * t449 - 0.75102
     >8806584362139918D-3 * t449 * t66 - 0.375514403292181069959D-3 * t2
     >56 * (0.1296D4 * t454 + 0.330192724889462668387D3 * t457) + 0.2173
     >44635793995449498D-3 * t457 + 0.137750893434399577552D-1 * t454 +
     >0.445122938368055555556D-4 * t71 * t13 * t75) * t82 - 0.2334981204
     >67045760261D0 * t284 * t31) - 0.273930929369445053578D0 * t94 * t3
     >03 * (0.124378109452736318408D1 * (0.757133580346724857871D-1 * (0
     >.198870000000000000000D0 * t476 * t106 - 0.248587500000000000000D-
     >1 * t153 * t144 * t317) * sigmabb * t113 + 0.757133580346724857871
     >D-1 * t109 * t112 * t97 + 0.144197530864197530864D0 * t136 * t511
     >- 0.751028806584362139918D-3 * t511 * t148 - 0.3755144032921810699
     >59D-3 * t376 * (0.1296D4 * t516 + 0.330192724889462668387D3 * t519
     >) + 0.217344635793995449498D-3 * t519 + 0.137750893434399577552D-1
     > * t516 + 0.445122938368055555556D-4 * t71 * t95 * t156) * t163 -
     >0.233498120467045760261D0 * t403 * t113)
      t541 = 0.1D1 / t17 / taua
      t555 = t30 * t430
      t565 = 0.227140074104017457361D0 * t555 * t49 - 0.2250000000000000
     >00000D0 * t227 * (0.201902288092459962099D0 * t555 * t43 + 0.25477
     >8337106066873756D-1 * t456 * t217)
      t570 = t16 * t541
      t581 = 0.1D1 / t99 / taub
      t595 = t112 * t492
      t605 = 0.227140074104017457361D0 * t595 * t131 - 0.225000000000000
     >000000D0 * t347 * (0.201902288092459962099D0 * t595 * t125 + 0.254
     >778337106066873756D-1 * t518 * t337)
      t610 = t98 * t581
      vtauc(i) = vtauc(i) - 0.681420222312052372084D0 * t11 * t12 * t182
     > * (0.757133580346724857871D-1 * (-0.198870000000000000000D0 * t16
     > * t541 * t23 + 0.248587500000000000000D-1 * t189 * t62 / t193 / t
     >aua * t196) * sigmaaa * t31 + 0.144197530864197530864D0 * t54 * t5
     >65 - 0.751028806584362139918D-3 * t565 * t66 + 0.48666666666666666
     >6667D0 * t256 * t570 - 0.137750893434399577552D-1 * t570) * t82 -
     >0.681420222312052372084D0 * t94 * t12 * t302 * (0.7571335803467248
     >57871D-1 * (-0.198870000000000000000D0 * t98 * t581 * t105 + 0.248
     >587500000000000000D-1 * t309 * t144 / t313 / taub * t316) * sigmab
     >b * t113 + 0.144197530864197530864D0 * t136 * t605 - 0.75102880658
     >4362139918D-3 * t605 * t148 + 0.486666666666666666667D0 * t376 * t
     >610 - 0.137750893434399577552D-1 * t610) * t163

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      sigmaaa = max(0.0D0, 0.250000000000000000000D0 * sigmacc(i))
      sigmaab = sigmaaa
      sigmabb = sigmaab
      sigma = sigmaaa + sigmabb + 0.2D1 * sigmaab
      taua = max(0.0D0, 0.500000000000000000000D0 * tauc(i))
      taub = taua
      tau = taua + taub
      t8 = pi ** 2
      t9 = t8 * rhoa
      t10 = t9 ** (0.1D1 / 0.3D1)
      t12 = 0.1D1 / pi
      t13 = sigmaaa ** 2
      t14 = rhoa ** 2
      t15 = 0.1D1 / t14
      t16 = t13 * t15
      t17 = taua ** 2
      t18 = 0.1D1 / t17
      t19 = t16 * t18
      t22 = (0.1D1 + 0.625000000000000000000D-1 * t19) ** 2
      t29 = t10 ** 2
      t30 = 0.1D1 / t29
      t34 = sigmaaa * t30
      t40 = t15 * (0.4D1 / sigmaaa * rhoa * taua - 0.1D1)
      t43 = 0.126188930057787476312D0 * t34 * t40 - 0.1D1
      t48 = sqrt(0.1D1 + 0.504755720231149905248D-1 * t34 * t40 * t43)
      t52 = t34 * t15
      t54 = 0.450000000000000000000D0 * t43 / t48 + 0.504755720231149905
     >248D-1 * t52
      t55 = t54 ** 2
      t61 = t14 ** 2
      t63 = t13 / t10 / t9 / t61
      t66 = sqrt(0.648D3 * t19 + 0.165096362444731334194D3 * t63)
      t71 = 0.1D1 / t8
      t74 = t61 ** 2
      t81 = (0.1D1 + 0.938662444277523956252D-1 * t52) ** 2
      t92 = t8 * rhob
      t93 = t92 ** (0.1D1 / 0.3D1)
      t95 = sigmabb ** 2
      t96 = rhob ** 2
      t97 = 0.1D1 / t96
      t98 = t95 * t97
      t99 = taub ** 2
      t100 = 0.1D1 / t99
      t101 = t98 * t100
      t104 = (0.1D1 + 0.625000000000000000000D-1 * t101) ** 2
      t111 = t93 ** 2
      t112 = 0.1D1 / t111
      t116 = sigmabb * t112
      t122 = t97 * (0.4D1 / sigmabb * rhob * taub - 0.1D1)
      t125 = 0.126188930057787476312D0 * t116 * t122 - 0.1D1
      t130 = sqrt(0.1D1 + 0.504755720231149905248D-1 * t116 * t122 * t12
     >5)
      t134 = t116 * t97
      t136 = 0.450000000000000000000D0 * t125 / t130 + 0.504755720231149
     >905248D-1 * t134
      t137 = t136 ** 2
      t143 = t96 ** 2
      t145 = t95 / t93 / t92 / t143
      t148 = sqrt(0.648D3 * t101 + 0.165096362444731334194D3 * t145)
      t155 = t143 ** 2
      t162 = (0.1D1 + 0.938662444277523956252D-1 * t134) ** 2
      zk(i) = -0.136284044462410474417D1 * rhoa * t10 * t12 * (0.1804D1 
     >- 0.804D0 / (0.1D1 + 0.124378109452736318408D1 * (0.75713358034672
     >4857871D-1 * (0.123456790123456790123D0 + 0.994350000000000000000D
     >-1 * t16 * t18 / t22) * sigmaaa * t30 * t15 + 0.720987654320987654
     >321D-1 * t55 - 0.751028806584362139918D-3 * t54 * t66 + 0.10867231
     >7896997724749D-3 * t63 + 0.688754467171997887761D-2 * t19 + 0.1483
     >74312789351851852D-4 * t71 * t13 * sigmaaa / t74) / t81)) - 0.1362
     >84044462410474417D1 * rhob * t93 * t12 * (0.1804D1 - 0.804D0 / (0.
     >1D1 + 0.124378109452736318408D1 * (0.757133580346724857871D-1 * (0
     >.123456790123456790123D0 + 0.994350000000000000000D-1 * t98 * t100
     > / t104) * sigmabb * t112 * t97 + 0.720987654320987654321D-1 * t13
     >7 - 0.751028806584362139918D-3 * t136 * t148 + 0.10867231789699772
     >4749D-3 * t145 + 0.688754467171997887761D-2 * t101 + 0.14837431278
     >9351851852D-4 * t71 * t95 * sigmabb / t155) / t162))

             endif
           enddo
         endif
       endif

      return
      end

c:TPSSXsubrend
