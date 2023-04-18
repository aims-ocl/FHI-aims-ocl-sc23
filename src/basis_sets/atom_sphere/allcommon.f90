!  This file is part of ATOM_SPHERE
!
!  ATOM_SPHERE was originally written by Stefan Goedecker
!  while he was at Cornell University, the Max Planck Institut Stuttgart
!  the CEA Grenoble and the university of Basel.
!  Santanu Saha at Basel university has interfaced it to LIBXC and
!  improved a few other things

! ATOM_SPHERE is free software: you can redistribute it and/or modify
! it under the terms of the version 3 of the license of the
! GNU Lesser General Public License as published by the Free
! Software Foundation.
!
! ATOM_SPHERE is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with ATOM_SPHERE. If not, see 
! https://www.gnu.org/licenses/lgpl-3.0.en.html
!
! ATOM_SPHERE reflects a substantial effort on the part of the
! developers, and we ask you to respect the spirit of the
! license that we chose and in particular to  keep any derivatives
! of ATOM_SPHERE under the same license that we chose for
! the original distribution, the GNU Lesser General Public License.

! We use atom_sphere in FHI-aims with permission from Stefan Goedecker 
! and his group. Please respect the license - thanks! VB

    subroutine contract_clebsch(lmax,lmaxp,ccleb)
    use localorb_io, only: localorb_info
    implicit real*8 (a-h,o-z)
    dimension ccleb(lmaxp+1,lmax+1,lmax+1)
    real*8, allocatable :: ap(:,:),aw(:),ylm(:,:),wylm(:,:),cleb(:,:,:),func(:),funcnew(:)
    character*300 :: info_str
! nang=12,32,72,110
    nang=72
    nlm=(lmax+1)**2
    nlmp=(lmaxp+1)**2

       if (lmax.ge.2 .and. nang.lt.32) stop 'nang at least 32'
       if (lmax.ge.3 .and. nang.lt.72) stop 'nang at least 72'
       if (lmax.ge.4 .and. nang.lt.110) stop 'nang at least 110'
       if (lmax.ge.5) stop 'no sufficient nang'

    allocate(ap(3,nang),aw(nang),ylm(nang,nlmp),wylm(nang,nlmp), & 
              cleb(nlmp,nlm,nlm),func(nang),funcnew(nang))

! angular grid
       call anggrid(nang,ap,aw)
! spherical harmonics
       call sphhar(nang,ap,aw,lmaxp,nlmp,ylm,wylm)
       call hartest(nang,nlmp,ylm,wylm)

       ilm1=0
       do l1=0,lmax
       do m1=-l1,l1
       ilm1=ilm1+1

         ilm2=0
         do l2=0,lmax
         do m2=-l2,l2
         ilm2=ilm2+1

           do ilmp=1,nlmp
           cleb(ilmp,ilm1,ilm2)=0.d0
           enddo
           
           do iang=1,nang
           func(iang)=ylm(iang,ilm1)*ylm(iang,ilm2)
           do ilmp=1,nlmp
           cleb(ilmp,ilm1,ilm2)=cleb(ilmp,ilm1,ilm2)+func(iang)*wylm(iang,ilmp)
           enddo
           enddo

           do ilmp=1,nlmp
           if (abs(cleb(ilmp,ilm1,ilm2)).lt.1.d-14) cleb(ilmp,ilm1,ilm2)=0.d0
           enddo

     !      write(6,*) 'product of l,m= ',l1,m1,' and  ',l2,m2
     !      ilmp=0
     !      do lp=0,lmaxp
     !      do mp=-lp,lp
     !      ilmp=ilmp+1
     !      write(6,*) lp,mp,cleb(ilmp,ilm1,ilm2)
     !      enddo
     !      enddo

!reconstruct function to check whether original function is regained
           do iang=1,nang
           funcnew(iang)=0.d0
           do ilmp=1,nlmp
           funcnew(iang)=funcnew(iang)+cleb(ilmp,ilm1,ilm2)*ylm(iang,ilmp)
           enddo
           if (abs(func(iang)-funcnew(iang)).gt. 1.d-14) then
           write(6,*) 'ERROR ',iang,func(iang),funcnew(iang)
           endif
           enddo

         enddo
         enddo
       enddo
       enddo

! contracted clebsch coefficients
         do lp=0,lmaxp
         do l1=0,lmax
         do l2=0,lmax
         ccleb(lp+1,l1+1,l2+1)=0.d0
         enddo
         enddo
         enddo

         ilmp=0
         do lp=0,lmaxp
         do mp=-lp,lp
         ilmp=ilmp+1
           ilm1=0
           do l1=0,lmax
           do m1=-l1,l1
           ilm1=ilm1+1
             ilm2=0
             do l2=0,lmax
             do m2=-l2,l2
             ilm2=ilm2+1
               ccleb(lp+1,l1+1,l2+1)=ccleb(lp+1,l1+1,l2+1)+cleb(ilmp,ilm1,ilm2)**2
               !write(6,'(6(i3),2(1x,e10.3))') lp,l1,l2,mp,m1,m2,cleb(ilmp,ilm1,ilm2) !,cleb(ilmp,ilm1,ilm2)**2
             enddo
             enddo
           enddo
           enddo
         enddo
         enddo

         do l1=0,lmax
         do l2=0,lmax
           write(info_str,'(2X,I2,I2,14(1X,E12.5))') l1,l2,(ccleb(lp+1,l1+1,l2+1),lp=0,lmaxp)
           call localorb_info( info_str )
         enddo
         enddo

         deallocate(ap,aw,ylm,wylm,cleb,func,funcnew)

       end

        subroutine sphhar(nang,ap,aw,lmaxp,nlmp,ylm,wylm)
        implicit real*8 (a-h,o-z)
        dimension ap(3,nang),aw(nang),ylm(nang,nlmp),wylm(nang,nlmp)
!
! REAL spherical harmonics up to l=5
!
       do 120 iang=1,nang
       x=ap(1,iang)
       y=ap(2,iang)
       z=ap(3,iang)

        ylm(iang,1)=0.2820947917738781d0
! l=1
        if (lmaxp.ge.1) then
        ylm(iang,2)=0.4886025119029199d0*x
        ylm(iang,3)=0.4886025119029199d0*y
        ylm(iang,4)=0.4886025119029199d0*z
        endif
! l=2
        if (lmaxp.ge.2) then
        ylm(iang,5)=0.5462742152960396d0*(x**2 - y**2)
        ylm(iang,6)=1.092548430592079d0*x*y
        ylm(iang,7)=1.092548430592079d0*x*z
        ylm(iang,8)=1.092548430592079d0*y*z
        ylm(iang,9)=0.3153915652525201d0*(-1.d0 + 3.d0*z**2)
        endif
! l=3
        if (lmaxp.ge.3) then
        ylm(iang,10)=0.5900435899266436d0*(x**3 - 3.d0*x*y**2)
        ylm(iang,11)=0.5900435899266436d0*(-3.d0*x**2*y + y**3)
        ylm(iang,12)=1.445305721320277d0*(x**2 - y**2)*z
        ylm(iang,13)=2.890611442640554d0*x*y*z
        ylm(iang,14)=0.4570457994644658d0*x*(-1.d0 + 5.d0*z**2)
        ylm(iang,15)=0.4570457994644658d0*y*(-1.d0 + 5.d0*z**2)
        ylm(iang,16)=0.3731763325901154d0*(-3.d0*z + 5.d0*z**3)
        endif
! l=4
        if (lmaxp.ge.4) then
        ylm(iang,17)=0.6258357354491761d0*(x**4 - 6.d0*x**2*y**2 + y**4)
        ylm(iang,18)=2.503342941796704d0*(x**3*y - x*y**3)
        ylm(iang,19)=1.77013076977993d0*(x**3 - 3.d0*x*y**2)*z
        ylm(iang,20)=1.77013076977993d0*(3.d0*x**2*y - y**3)*z
        ylm(iang,21)=0.4730873478787801d0*(x**2 - y**2)*(7.d0*z**2-1.d0)
        ylm(iang,22)=0.94617469575756d0*x*y*(-1.d0 + 7.d0*z**2)
        ylm(iang,23)=0.6690465435572892d0*x*(-3.d0*z + 7.d0*z**3)
        ylm(iang,24)=0.6690465435572892d0*y*(-3.d0*z + 7.d0*z**3)
        ylm(iang,25)=0.1057855469152043d0*(3.d0-30.d0*z**2+35.d0*z**4)
        endif
! l=5
        if (lmaxp.ge.5) then
        ylm(iang,26)=0.65638205684017d0*y* & 
                     (5.d0*x**4 - 10.d0*x**2*y**2 + y**4)
        ylm(iang,27)=8.30264925952417d0*x*y* & 
                     (x**2 - y**2)*z
        ylm(iang,28)=0.4892382994352503d0*y* & 
                     (-3.d0*x**2 + y**2)*(x**2 + y**2 - 8.d0*z**2)
        ylm(iang,29)=4.793536784973324d0*x*y*z* & 
                     (x**2 + y**2 - 2.d0*z**2)
        ylm(iang,30)=0.4529466511956969d0*y* & 
                     (x**4 + y**4 - 12.d0*y**2*z**2 + 8.d0*z**4 +  & 
                                       2.d0*x**2*(y**2 - 6.d0*z**2))
        ylm(iang,31)=0.1169503224534235d0*z* & 
                     (15.d0*x**4 + 15.d0*y**4 - 40.d0*y**2*z**2 +  & 
                      8.d0*z**4 + 10.d0*x**2*(3.d0*y**2 - 4.d0*z**2))
        ylm(iang,32)=0.4529466511956969d0*x* & 
                     (x**4 + y**4 - 12.d0*y**2*z**2 + 8.d0*z**4 +  & 
                                        2.d0*x**2*(y**2 - 6.d0*z**2))
        ylm(iang,33)=-2.396768392486662d0* & 
                     (x**2 - y**2)*z*(x**2 + y**2 - 2.d0*z**2)
        ylm(iang,34)=-0.4892382994352503d0*x* & 
                     (x**2 - 3.d0*y**2)*(x**2 + y**2 - 8.d0*z**2)
        ylm(iang,35)=2.075662314881041d0* & 
                     (x**4 - 6.d0*x**2*y**2 + y**4)*z
        ylm(iang,36)=0.65638205684017d0* & 
                     (x**5 - 10.d0*x**3*y**2 + 5.d0*x*y**4)
        endif

        if (lmaxp.ge.6) then
        ylm(iang,37)=1.366368210383828d0*x*y*(3.d0*x**4 -  & 
                     10.d0*x**2*y**2 + 3.d0*y**4)
        ylm(iang,38)=2.366619162231752d0*y*(5.d0*x**4 -  & 
                     10.d0*x**2*y**2 + y**4)*z
        ylm(iang,39)=-2.018259602914896d0*x*y*(x**2 - y**2)* & 
                     (x**2 + y**2 - 10.d0*z**2)
        ylm(iang,40)=0.921205259514924d0*y*(-3.d0*x**2 + y**2)*z* & 
                     (3.d0*x**2 + 3.d0*y**2 - 8.d0*z**2)
        ylm(iang,41)=0.921205259514924d0*x*y* & 
                     (x**4 + y**4 - 16.d0*y**2*z**2 + 16.d0*z**4 +  &
                     2.d0*x**2*(y**2 - 8.d0*z**2))
        ylm(iang,42)=0.5826213625187314d0*y*z* &
                     (5.d0*x**4 + 5.d0*y**4 - 20.d0*y**2*z**2 +  &
                     8.d0*z**4 + 10.d0*x**2*(y**2 - 2.d0*z**2))
        ylm(iang,43)=0.06356920226762841d0* &
                     (231.d0*z**6 - 315.d0*z**4*(x**2 + y**2 + z**2) +  &
                      105.d0*z**2*(x**2 + y**2 + z**2)**2 -  &
                      5.d0*(x**2 + y**2 + z**2)**3)
        ylm(iang,44)=0.5826213625187314d0*x*z* &
                     (5.d0*x**4 + 5.d0*y**4 - 20.d0*y**2*z**2 +  &
                     8.d0*z**4 + 10.d0*x**2*(y**2 - 2.d0*z**2))
        ylm(iang,45)=0.4606026297574618d0*(x**2 - y**2)* &
                     (x**4 + y**4 - 16.d0*y**2*z**2 + 16.d0*z**4 +  &
                     2.d0*x**2*(y**2 - 8.d0*z**2))
        ylm(iang,46)=-0.921205259514924d0*x*(x**2 - 3.d0*y**2)*z* &
                     (3.d0*x**2 + 3.d0*y**2 - 8.d0*z**2)
        ylm(iang,47)=-0.5045649007287241d0* &
                     (x**4 - 6.d0*x**2*y**2 + y**4)* &
                     (x**2 + y**2 - 10.d0*z**2)
        ylm(iang,48)=2.366619162231752d0*(x**5 -  &
                     10.d0*x**3*y**2 + 5.d0*x*y**4)*z
        ylm(iang,49)=0.6831841051919143d0* &
                     (x**6 - 15.d0*x**4*y**2 + 15.d0*x**2*y**4 - y**6)
        endif

        if (lmaxp.ge.7) then
        ylm(iang,50)=-0.7071627325245961d0*y*(-7.d0*x**6 +  &
                     35.d0*x**4*y**2 - 21.d0*x**2*y**4 + y**6)
        ylm(iang,51)=5.2919213236038d0*x*y*(3.d0*x**4 -  &
                     10.d0*x**2*y**2 + 3.d0*y**4)*z
        ylm(iang,52)=-0.5189155787202603d0*y*(5.d0*x**4 -  &
                     10.d0*x**2*y**2 + y**4)*(x**2 + y**2 - 12.d0*z**2)
        ylm(iang,53)=4.151324629762083d0*x*y*(-x**2 + y**2)*z* &
                     (3.d0*x**2 + 3.d0*y**2 - 10.d0*z**2)
        ylm(iang,54)=-0.156458933862294d0*y*(-3.d0*x**2 + y**2)* &
                     (3.d0*x**4 + 3.d0*y**4 - 60.d0*y**2*z**2 +  &
                     80.d0*z**4 + 6.d0*x**2*(y**2 - 10.d0*z**2))
        ylm(iang,55)=0.4425326924449826d0*x*y*z*(15.d0*x**4 +  &
                     15.d0*y**4 - 80.d0*y**2*z**2 + 48.d0*z**4 +  &
                     10.d0*x**2*(3.d0*y**2 - 8.d0*z**2))
        ylm(iang,56)=-0.0903316075825173d0*y*(5.d0*x**6 + 5.d0*y**6 -  &
                     120.d0*y**4*z**2 + 240.d0*y**2*z**4 - 64.d0*z**6 +  &
                     15.d0*x**4*(y**2 - 8.d0*z**2) + 15.d0*x**2*(y**4 -  &
                     16.d0*y**2*z**2 + 16.d0*z**4))
        ylm(iang,57)=0.06828427691200494d0*z*(-35.d0*x**6 - 35.d0*y**6 +  &
                     210.d0*y**4*z**2 - 168.d0*y**2*z**4 + 16.d0*z**6 -  &
                     105.d0*x**4*(y**2 - 2.d0*z**2) - 21.d0*x**2* &
                     (5.d0*y**4 - 20.d0*y**2*z**2 + 8.d0*z**4))
        ylm(iang,58)=-0.0903316075825173d0*x*(5.d0*x**6 + 5.d0*y**6 -  &
                     120.d0*y**4*z**2 + 240.d0*y**2*z**4 - 64.d0*z**6 +  &
                     15.d0*x**4*(y**2 - 8.d0*z**2) + 15.d0*x**2*(y**4 -  &
                     16.d0*y**2*z**2 + 16.d0*z**4))
        ylm(iang,59)=0.2212663462224913d0*(x**2 - y**2)*z*(15.d0*x**4 +  &
                     15.d0*y**4 - 80.d0*y**2*z**2 + 48.d0*z**4 +  &
                     10.d0*x**2*(3.d0*y**2 - 8.d0*z**2))
        ylm(iang,60)=0.156458933862294d0*x*(x**2 - 3.d0*y**2)* &
                     (3.d0*x**4 + 3.d0*y**4 - 60.d0*y**2*z**2 +  &
                     80.d0*z**4 + 6.d0*x**2*(y**2 - 10.d0*z**2))
        ylm(iang,61)=-1.03783115744052d0*(x**4 - 6.d0*x**2*y**2 + y**4)* &
                     z*(3.d0*x**2 + 3.d0*y**2 - 10.d0*z**2)
        ylm(iang,62)=-0.5189155787202603d0*x*(x**4 - 10.d0*x**2*y**2 +  &
                     5.d0*y**4)*(x**2 + y**2 - 12.d0*z**2)
        ylm(iang,63)=2.6459606618019d0*(x**6 - 15.d0*x**4*y**2 +  &
                     15.d0*x**2*y**4 - y**6)*z
        ylm(iang,64)=0.7071627325245961d0*x*(x**6 - 21.d0*x**4*y**2 +  &
                     35.d0*x**2*y**4 - 7.d0*y**6)
        endif

        if (lmaxp.ge.8) then
        ylm(iang,65)=5.831413281398639d0*x*y*(x**6 - 7.d0*x**4*y**2 +  &
                     7.d0*x**2*y**4 - y**6)
        ylm(iang,66)=-2.915706640699319d0*y*(-7.d0*x**6 +  &
                     35.d0*x**4*y**2 - 21.d0*x**2*y**4 + y**6)*z
        ylm(iang,67)=-1.064665532119085d0*x*y*(3.d0*x**4 -  &
                     10.d0*x**2*y**2 + 3.d0*y**4)* &
                     (x**2 + y**2 - 14.d0*z**2)
        ylm(iang,68)=-3.449910622098108d0*y*(5.d0*x**4 -  &
                     10.d0*x**2*y**2 + y**4)*z* &
                     (x**2 + y**2 - 4.d0*z**2)
        ylm(iang,69)=1.913666099037322d0*x*y*(x**2 - y**2)*(x**4 +  &
                     y**4 - 24.d0*y**2*z**2 + 40.d0*z**4 + 2.d0*x**2* &
                     (y**2 - 12.d0*z**2))
        ylm(iang,70)=-1.235266155295544d0*y*(-3.d0*x**2 + y**2)*z* &
                     (3.d0*x**4 + 3.d0*y**4 - 20.d0*y**2*z**2 +  &
                     16.d0*z**4 + x**2*(6.d0*y**2 - 20.d0*z**2))
        ylm(iang,71)=-0.912304516869819d0*x*y*(x**6 + y**6 -  &
                     30.d0*y**4*z**2 + 80.d0*y**2*z**4 - 32.d0*z**6 +  &
                     3.d0*x**4*(y**2 - 10.d0*z**2) + x**2*(3.d0*y**4 -  &
                     60.d0*y**2*z**2 + 80.d0*z**4))
        ylm(iang,72)=-0.1090412458987799d0*y*z*(35.d0*x**6 +  &
                     35.d0*y**6 - 280.d0*y**4*z**2 + 336.d0*y**2*z**4 -  &
                     64.d0*z**6 + 35.d0*x**4*(3.d0*y**2 - 8.d0*z**2) +  &
                     7.d0*x**2*(15.d0*y**4 -  &
                     80.d0*y**2*z**2 + 48.d0*z**4))
        ylm(iang,73)=0.009086770491565d0*(35.d0*x**8 + 35.d0*y**8 -  &
                     1120.d0*y**6*z**2 + 3360.d0*y**4*z**4 -  &
                     1792.d0*y**2*z**6 + 128.d0*z**8 + 140.d0*x**6* &
                     (y**2 - 8.d0*z**2) + 210.d0*x**4*(y**4 -  &
                     16.d0*y**2*z**2 + 16.d0*z**4) + 28.d0*x**2* &
                     (5.d0*y**6 - 120.d0*y**4*z**2 +  &
                     240.d0*y**2*z**4 - 64.d0*z**6))
        ylm(iang,74)=-0.1090412458987799d0*x*z*(35.d0*x**6 +  &
                     35.d0*y**6 - 280.d0*y**4*z**2 +  &
                     336.d0*y**2*z**4 - 64.d0*z**6 + 35.d0*x**4* &
                     (3.d0*y**2 - 8.d0*z**2) + 7.d0*x**2*(15.d0*y**4 -  &
                     80.d0*y**2*z**2 + 48.d0*z**4))
        ylm(iang,75)=-0.4561522584349094d0*(x**2 - y**2)*(x**6 + y**6 -  &
                     30.d0*y**4*z**2 + 80.d0*y**2*z**4 - 32.d0*z**6 +  &
                     3.d0*x**4*(y**2 - 10.d0*z**2) + x**2*(3.d0*y**4 -  &
                     60.d0*y**2*z**2 + 80.d0*z**4))
        ylm(iang,76)=1.235266155295544d0*x*(x**2 - 3.d0*y**2)*z* &
                     (3.d0*x**4 + 3.d0*y**4 - 20.d0*y**2*z**2 +  &
                     16.d0*z**4 + x**2*(6.d0*y**2 - 20.d0*z**2))
        ylm(iang,77)=0.4784165247593306d0*(x**4 - 6.d0*x**2*y**2 + y**4) &
                     *(x**4 + y**4 - 24.d0*y**2*z**2 + 40.d0*z**4 +  &
                     2.d0*x**2*(y**2 - 12.d0*z**2))
        ylm(iang,78)=-3.449910622098108d0*x*(x**4 - 10.d0*x**2*y**2 +  &
                     5.d0*y**4)*z*(x**2 + y**2 - 4.d0*z**2)
        ylm(iang,79)=0.5323327660595426d0*(-x**6 + 15.d0*x**4*y**2 -  &
                     15.d0*x**2*y**4 + y**6)*(x**2 + y**2 - 14.d0*z**2) 
        ylm(iang,80)=2.915706640699319d0*x*(x**6 - 21.d0*x**4*y**2 +  & 
                     35.d0*x**2*y**4 - 7.d0*y**6)*z 
        ylm(iang,81)=0.7289266601748298d0*(x**8 - 28.d0*x**6*y**2 +  &
                     70.d0*x**4*y**4 - 28.d0*x**2*y**6 + y**8)
        endif

        if (lmaxp.gt.8) stop 'sphhar'
! wylm
        ilm=0
        do 15,l=0,lmaxp
        do 15,m=-l,l
        ilm=ilm+1
        wylm(iang,ilm)=ylm(iang,ilm)*aw(iang)
15      continue

120    continue

       return
       end



        subroutine anggrid(nang,ap,aw)
        use localorb_io, only: localorb_info
        implicit real*8 (a-h,o-z)
        dimension ap(3,nang),aw(nang)
        character*300 :: info_str
       if (nang.eq.12) then
! A. H. Stroud, Prentice Hall 1971, page 296
! integrates up to order 5 (.eq. up to x^5 or x^1*y^2*z*2)
       write(info_str,'(2X,A)')'      '
       call localorb_info( info_str )
       write(info_str,'(2X,A)') '5-th order angular integration grid'
       call localorb_info( info_str )
!     
       r2=(5.d0+dsqrt(5.d0))/10.d0
       rr=dsqrt(r2)
       s2=(5.d0-dsqrt(5.d0))/10.d0
       ss=dsqrt(s2)
!
!  angular integration points
!
       ap(1,1)=rr
       ap(2,1)=ss
       ap(3,1)=0.d0

       ap(1,2)=-rr
       ap(2,2)=ss
       ap(3,2)=0.d0

       ap(1,3)=rr
       ap(2,3)=-ss
       ap(3,3)=0.d0

       ap(1,4)=-rr
       ap(2,4)=-ss
       ap(3,4)=0.d0

       ap(1,5)=0.d0
       ap(2,5)=rr
       ap(3,5)=ss

       ap(1,6)=0.d0
       ap(2,6)=-rr
       ap(3,6)=ss

       ap(1,7)=0.d0
       ap(2,7)=rr
       ap(3,7)=-ss

       ap(1,8)=0.d0
       ap(2,8)=-rr
       ap(3,8)=-ss

       ap(1,9)=ss
       ap(2,9)=0.d0
       ap(3,9)=rr

       ap(1,10)=-ss
       ap(2,10)=0.d0
       ap(3,10)=rr

       ap(1,11)=ss
       ap(2,11)=0.d0
       ap(3,11)=-rr

       ap(1,12)=-ss
       ap(2,12)=0.d0
       ap(3,12)=-rr

! angular integration weights

       fpi=4.d0*dacos(-1.d0)
       do iang=1,nang
       aw(iang)=(1.d0/12.d0)*fpi
       enddo


       else if (nang.eq.32) then
! A. H. Stroud, Prentice Hall 1971, page 299
! integrates up to order 9 (.eq. up to x^9 or x^3*y^3*z*3
       write(info_str,'(2X,A)') '      '
       call localorb_info( info_str )
       write(info_str,'(2X,A)') '9-th order angular integration grid'
       call localorb_info( info_str )
     
       r2=(5.d0+dsqrt(5.d0))/10.d0
       rr=dsqrt(r2)
       s2=(5.d0-dsqrt(5.d0))/10.d0
       ss=dsqrt(s2)
       u2=(3.d0-dsqrt(5.d0))/6.d0
       uu=dsqrt(u2)
       v2=(3.d0+dsqrt(5.d0))/6.d0
       vv=dsqrt(v2)
       t2=1.d0/3.d0
       tt=dsqrt(t2)
!  angular integration points

       ap(1,1)=rr
       ap(2,1)=ss
       ap(3,1)=0.d0

       ap(1,2)=-rr
       ap(2,2)=ss
       ap(3,2)=0.d0

       ap(1,3)=rr
       ap(2,3)=-ss
       ap(3,3)=0.d0

       ap(1,4)=-rr
       ap(2,4)=-ss
       ap(3,4)=0.d0

       ap(1,5)=0.d0
       ap(2,5)=rr
       ap(3,5)=ss

       ap(1,6)=0.d0
       ap(2,6)=-rr
       ap(3,6)=ss

       ap(1,7)=0.d0
       ap(2,7)=rr
       ap(3,7)=-ss

       ap(1,8)=0.d0
       ap(2,8)=-rr
       ap(3,8)=-ss

       ap(1,9)=ss
       ap(2,9)=0.d0
       ap(3,9)=rr

       ap(1,10)=-ss
       ap(2,10)=0.d0
       ap(3,10)=rr

       ap(1,11)=ss
       ap(2,11)=0.d0
       ap(3,11)=-rr

       ap(1,12)=-ss
       ap(2,12)=0.d0
       ap(3,12)=-rr

       ap(1,13)=uu
       ap(2,13)=vv
       ap(3,13)=0.d0

       ap(1,14)=-uu
       ap(2,14)=vv
       ap(3,14)=0.d0

       ap(1,15)=uu
       ap(2,15)=-vv
       ap(3,15)=0.d0

       ap(1,16)=-uu
       ap(2,16)=-vv
       ap(3,16)=0.d0

       ap(1,17)=0.d0
       ap(2,17)=uu
       ap(3,17)=vv

       ap(1,18)=0.d0
       ap(2,18)=-uu
       ap(3,18)=vv

       ap(1,19)=0.d0
       ap(2,19)=uu
       ap(3,19)=-vv

       ap(1,20)=0.d0
       ap(2,20)=-uu
       ap(3,20)=-vv

       ap(1,21)=vv
       ap(2,21)=0.d0
       ap(3,21)=uu

       ap(1,22)=-vv
       ap(2,22)=0.d0
       ap(3,22)=uu

       ap(1,23)=vv
       ap(2,23)=0.d0
       ap(3,23)=-uu

       ap(1,24)=-vv
       ap(2,24)=0.d0
       ap(3,24)=-uu

       ap(1,25)=tt
       ap(2,25)=tt
       ap(3,25)=tt

       ap(1,26)=-tt
       ap(2,26)=tt
       ap(3,26)=tt

       ap(1,27)=tt
       ap(2,27)=-tt
       ap(3,27)=tt

       ap(1,28)=tt
       ap(2,28)=tt
       ap(3,28)=-tt

       ap(1,29)=tt
       ap(2,29)=-tt
       ap(3,29)=-tt

       ap(1,30)=-tt
       ap(2,30)=tt
       ap(3,30)=-tt

       ap(1,31)=-tt
       ap(2,31)=-tt
       ap(3,31)=tt

       ap(1,32)=-tt
       ap(2,32)=-tt
       ap(3,32)=-tt


! angular integration weights

       fpi=4.d0*dacos(-1.d0)
       do iang=1,nang
       if(iang.le.12)then
       aw(iang)=(25.d0/840.d0)*fpi
       else
       aw(iang)=(27.d0/840.d0)*fpi
       end if
       enddo

       else if (nang.eq.72) then
! A. H. Stroud, Prentice Hall 1971, page 302
! integrates up to order 14 (.eg. up to x^14 or x^6*y^4*z*4)
       write(info_str,'(2X,A)')'      '
       call localorb_info( info_str )
       write(info_str,'(2X,A)') '14-th order angular integration grid'
       call localorb_info( info_str )


       r2=(5.d0-dsqrt(5.d0))/10.d0
       r=dsqrt(r2)
       s2=(5.d0+dsqrt(5.d0))/10.d0
       s=dsqrt(s2)

       z6=0.52613220610787372459d-1
       z5=0.211497863190394288d0
       z4=0.3847053176951983194d0
       z3=0.641786069679381079d0
       z2=0.748834163668201058d0
       z1=0.9120637902629030996d0

       u1=(-z3+z4)/(2*s)
       v1=( z5+z6)/(2*s)
       w1=( z1+z2)/(2*s)

       u2=(-z5+z2)/(2*s)
       v2=( z6+z4)/(2*s)
       w2=( z1+z3)/(2*s)

       u3=(-z2+z6)/(2*s)
       v3=( z3+z5)/(2*s)
       w3=( z1+z4)/(2*s)

       u4=(-z6+z3)/(2*s)
       v4=( z4+z2)/(2*s)
       w4=( z1+z5)/(2*s)

       u5=(-z4+z5)/(2*s)
       v5=( z2+z3)/(2*s)
       w5=( z1+z6)/(2*s)

      ic=0
      ic=ic+1
       ap(1,ic)=r
       ap(2,ic)=s
       ap(3,ic)=0
      ic=ic+1
       ap(1,ic)=r
       ap(2,ic)=-s
       ap(3,ic)=0
      ic=ic+1
       ap(1,ic)=-r
       ap(2,ic)=s
       ap(3,ic)=0
      ic=ic+1
       ap(1,ic)=-r
       ap(2,ic)=-s
       ap(3,ic)=0
      ic=ic+1
       ap(1,ic)=0
       ap(2,ic)=r
       ap(3,ic)=s
      ic=ic+1
       ap(1,ic)=0
       ap(2,ic)=r
       ap(3,ic)=-s
      ic=ic+1
       ap(1,ic)=0
       ap(2,ic)=-r
       ap(3,ic)=s
      ic=ic+1
       ap(1,ic)=0
       ap(2,ic)=-r
       ap(3,ic)=-s
      ic=ic+1
       ap(1,ic)=s
       ap(2,ic)=0
       ap(3,ic)=r
      ic=ic+1
       ap(1,ic)=-s
       ap(2,ic)=0
       ap(3,ic)=r
      ic=ic+1
       ap(1,ic)=s
       ap(2,ic)=0
       ap(3,ic)=-r
      ic=ic+1
       ap(1,ic)=-s
       ap(2,ic)=0
       ap(3,ic)=-r
      ic=ic+1
       ap(1,ic)=u1
       ap(2,ic)=v1
       ap(3,ic)=w1
      ic=ic+1
       ap(1,ic)=u2
       ap(2,ic)=v2
       ap(3,ic)=w2
      ic=ic+1
       ap(1,ic)=u3
       ap(2,ic)=v3
       ap(3,ic)=w3
      ic=ic+1
       ap(1,ic)=u4
       ap(2,ic)=v4
       ap(3,ic)=w4
      ic=ic+1
       ap(1,ic)=u5
       ap(2,ic)=v5
       ap(3,ic)=w5
      ic=ic+1
       ap(1,ic)=u1
       ap(2,ic)=-v1
       ap(3,ic)=-w1
      ic=ic+1
       ap(1,ic)=u2
       ap(2,ic)=-v2
       ap(3,ic)=-w2
      ic=ic+1
       ap(1,ic)=u3
       ap(2,ic)=-v3
       ap(3,ic)=-w3
      ic=ic+1
       ap(1,ic)=u4
       ap(2,ic)=-v4
       ap(3,ic)=-w4
      ic=ic+1
       ap(1,ic)=u5
       ap(2,ic)=-v5
       ap(3,ic)=-w5
      ic=ic+1
       ap(1,ic)=-u1
       ap(2,ic)=-v1
       ap(3,ic)=w1
      ic=ic+1
       ap(1,ic)=-u2
       ap(2,ic)=-v2
       ap(3,ic)=w2
      ic=ic+1
       ap(1,ic)=-u3
       ap(2,ic)=-v3
       ap(3,ic)=w3
      ic=ic+1
       ap(1,ic)=-u4
       ap(2,ic)=-v4
       ap(3,ic)=w4
      ic=ic+1
       ap(1,ic)=-u5
       ap(2,ic)=-v5
       ap(3,ic)=w5
      ic=ic+1
       ap(1,ic)=-u1
       ap(2,ic)=v1
       ap(3,ic)=-w1
      ic=ic+1
       ap(1,ic)=-u2
       ap(2,ic)=v2
       ap(3,ic)=-w2
      ic=ic+1
       ap(1,ic)=-u3
       ap(2,ic)=v3
       ap(3,ic)=-w3
      ic=ic+1
       ap(1,ic)=-u4
       ap(2,ic)=v4
       ap(3,ic)=-w4
      ic=ic+1
       ap(1,ic)=-u5
       ap(2,ic)=v5
       ap(3,ic)=-w5
      ic=ic+1
       ap(1,ic)=v1
       ap(2,ic)=w1
       ap(3,ic)=u1
      ic=ic+1
       ap(1,ic)=v2
       ap(2,ic)=w2
       ap(3,ic)=u2
      ic=ic+1
       ap(1,ic)=v3
       ap(2,ic)=w3
       ap(3,ic)=u3
      ic=ic+1
       ap(1,ic)=v4
       ap(2,ic)=w4
       ap(3,ic)=u4
      ic=ic+1
       ap(1,ic)=v5
       ap(2,ic)=w5
       ap(3,ic)=u5
      ic=ic+1
       ap(1,ic)=v1
       ap(2,ic)=-w1
       ap(3,ic)=-u1
      ic=ic+1
       ap(1,ic)=v2
       ap(2,ic)=-w2
       ap(3,ic)=-u2
      ic=ic+1
       ap(1,ic)=v3
       ap(2,ic)=-w3
       ap(3,ic)=-u3
      ic=ic+1
       ap(1,ic)=v4
       ap(2,ic)=-w4
       ap(3,ic)=-u4
      ic=ic+1
       ap(1,ic)=v5
       ap(2,ic)=-w5
       ap(3,ic)=-u5
      ic=ic+1
       ap(1,ic)=-v1
       ap(2,ic)=-w1
       ap(3,ic)=u1
      ic=ic+1
       ap(1,ic)=-v2
       ap(2,ic)=-w2
       ap(3,ic)=u2
      ic=ic+1
       ap(1,ic)=-v3
       ap(2,ic)=-w3
       ap(3,ic)=u3
      ic=ic+1
       ap(1,ic)=-v4
       ap(2,ic)=-w4
       ap(3,ic)=u4
      ic=ic+1
       ap(1,ic)=-v5
       ap(2,ic)=-w5
       ap(3,ic)=u5
      ic=ic+1
       ap(1,ic)=-v1
       ap(2,ic)=w1
       ap(3,ic)=-u1
      ic=ic+1
       ap(1,ic)=-v2
       ap(2,ic)=w2
       ap(3,ic)=-u2
      ic=ic+1
       ap(1,ic)=-v3
       ap(2,ic)=w3
       ap(3,ic)=-u3
      ic=ic+1
       ap(1,ic)=-v4
       ap(2,ic)=w4
       ap(3,ic)=-u4
      ic=ic+1
       ap(1,ic)=-v5
       ap(2,ic)=w5
       ap(3,ic)=-u5
      ic=ic+1
       ap(1,ic)=w1
       ap(2,ic)=u1
       ap(3,ic)=v1
      ic=ic+1
       ap(1,ic)=w2
       ap(2,ic)=u2
       ap(3,ic)=v2
      ic=ic+1
       ap(1,ic)=w3
       ap(2,ic)=u3
       ap(3,ic)=v3
      ic=ic+1
       ap(1,ic)=w4
       ap(2,ic)=u4
       ap(3,ic)=v4
      ic=ic+1
       ap(1,ic)=w5
       ap(2,ic)=u5
       ap(3,ic)=v5
      ic=ic+1
       ap(1,ic)=w1
       ap(2,ic)=-u1
       ap(3,ic)=-v1
      ic=ic+1
       ap(1,ic)=w2
       ap(2,ic)=-u2
       ap(3,ic)=-v2
      ic=ic+1
       ap(1,ic)=w3
       ap(2,ic)=-u3
       ap(3,ic)=-v3
      ic=ic+1
       ap(1,ic)=w4
       ap(2,ic)=-u4
       ap(3,ic)=-v4
      ic=ic+1
       ap(1,ic)=w5
       ap(2,ic)=-u5
       ap(3,ic)=-v5
      ic=ic+1
       ap(1,ic)=-w1
       ap(2,ic)=-u1
       ap(3,ic)=v1
      ic=ic+1
       ap(1,ic)=-w2
       ap(2,ic)=-u2
       ap(3,ic)=v2
      ic=ic+1
       ap(1,ic)=-w3
       ap(2,ic)=-u3
       ap(3,ic)=v3
      ic=ic+1
       ap(1,ic)=-w4
       ap(2,ic)=-u4
       ap(3,ic)=v4
      ic=ic+1
       ap(1,ic)=-w5
       ap(2,ic)=-u5
       ap(3,ic)=v5
      ic=ic+1
       ap(1,ic)=-w1
       ap(2,ic)=u1
       ap(3,ic)=-v1
      ic=ic+1
       ap(1,ic)=-w2
       ap(2,ic)=u2
       ap(3,ic)=-v2
      ic=ic+1
       ap(1,ic)=-w3
       ap(2,ic)=u3
       ap(3,ic)=-v3
      ic=ic+1
       ap(1,ic)=-w4
       ap(2,ic)=u4
       ap(3,ic)=-v4
      ic=ic+1
       ap(1,ic)=-w5
       ap(2,ic)=u5
       ap(3,ic)=-v5
       if (ic.ne.72) stop 'stroud 72'

! angular integration weights

       fpi=4.d0*dacos(-1.d0)
       do iang=1,nang
       if(iang.le.12)then
       aw(iang)=(125.d0/10080.d0)*fpi
       else
       aw(iang)=(143.d0/10080.d0)*fpi
       end if
       enddo

       else if (nang.eq.110) then
       open(unit=50,file='lebedev110.dat',status='unknown')
       do iang=1,nang
       read(50,*) ap(1,iang),ap(2,iang),ap(3,iang),aw(iang)
       enddo
       close(50)

       else 
       write(6,*) 'NANG not permitted'
       stop
       endif

       return
       end


       subroutine hartest(nang,nlmp,ylm,wylm)
! test spherical harmonics and radial grid
        implicit real*8 (a-h,o-z)
       parameter(eps=1.d-8)
        dimension ylm(nang,nlmp),wylm(nang,nlmp)

       do 150,ilm=1,nlmp
       do 150,jlm=1,nlmp
       tt=0.d0
       do 120, iang=1,nang
       tt=tt+ylm(iang,ilm)*wylm(iang,jlm)
120     continue
!       write(6,*) ilm,jlm,tt
       if (ilm.eq.jlm) then
        if (abs(tt-1.d0).gt.eps) write(6,*) 'ERROR',ilm,jlm,tt
       else
        if (abs(tt).gt.eps) write(6,*) 'ERROR',ilm,jlm,tt
       endif
150     continue

       return
       end


!***THE KINETIC ENERGY IS CALCULATED ON 4*RR RADIAL GRID INSTED OF************!
!***RR RADIAL GRID FOR BETTER ACCURACY. THIS SUBROUTINE SUM UP THE************!
!***CONTRIBUTION OF THE ADDITIONAL GRIDS IN 4*RR AND PROJECT IT TO************!
!***RR RADIAL GRID. **********************************************************!
          subroutine mult_t(nrad,w,y)
          implicit real*8 (a-h,o-z)
          dimension w(4,nrad),y(nrad)
          dimension c41( 0:8),c42( 0:8),c43( 0:8),c44( 0:8)
          dimension c31(-1:7),c32(-1:7),c33(-1:7),c34(-1:7)
          dimension c21(-2:6),c22(-2:6),c23(-2:6),c24(-2:6)
          dimension c11(-3:5),c12(-3:5),c13(-3:5),c14(-3:5)
          dimension c01(-4:4),c02(-4:4),c03(-4:4),c04(-4:4)


         c41(0)=-2.717857142857143d0
         c41(1)=8.d0
         c41(2)=-14.d0
         c41(3)=18.66666666666667d0
         c41(4)=-17.5d0
         c41(5)=11.2d0
         c41(6)=-4.666666666666667d0
         c41(7)=1.142857142857143d0
         c41(8)=-0.125d0

         c42(0)=-1.533237675258092d0
         c42(1)=2.732823108491443d0
         c42(2)=-2.637493896484375d0
         c42(3)=2.849429321289063d0
         c42(4)=-2.394930521647135d0
         c42(5)=1.433224487304688d0
         c42(6)=-0.57060546875d0
         c42(7)=0.1352159772600446d0
         c42(8)=-0.01442533220563616d0

         c43(0)=-0.7940848214285714d0
         c43(1)=-0.06849888392857143d0
         c43(2)=2.523763020833333d0
         c43(3)=-3.6150390625d0
         c43(4)=3.4521484375d0
         c43(5)=-2.2255859375d0
         c43(6)=0.9306640625d0
         c43(7)=-0.2283761160714286d0
         c43(8)=0.0250093005952381d0

         c44(0)=-0.359611329578218d0
         c44(1)=-1.296153477260045d0
         c44(2)=3.98780517578125d0
         c44(3)=-4.811203002929687d0
         c44(4)=4.290122985839844d0
         c44(5)=-2.665542602539062d0
         c44(6)=1.089182535807292d0
         c44(7)=-0.2630802699497768d0
         c44(8)=0.02847998482840402d0



         c31(-1)=-0.125d0
         c31(0)=-1.592857142857143d0
         c31(1)=3.5d0
         c31(2)=-3.5d0
         c31(3)=2.916666666666667d0
         c31(4)=-1.75d0
         c31(5)=0.7d0
         c31(6)=-0.1666666666666667d0
         c31(7)=0.01785714285714286d0

         c32(-1)=-0.01442533220563616d0
         c32(0)=-1.403409685407366d0
         c32(1)=2.213511149088542d0
         c32(2)=-1.425765991210937d0
         c32(3)=1.031837463378906d0
         c32(4)=-0.5773386637369792d0
         c32(5)=0.22149658203125d0
         c32(6)=-0.05129350934709821d0
         c32(7)=0.005387987409319196d0

         c33(-1)=0.0250093005952381d0
         c33(0)=-1.019168526785714d0
         c33(1)=0.8318359375d0
         c33(2)=0.4229817708333333d0
         c33(3)=-0.4638671875d0
         c33(4)=0.3009765625d0
         c33(5)=-0.1248046875d0
         c33(6)=0.03032924107142857d0
         c33(7)=-0.003292410714285714d0

         c34(-1)=0.02847998482840402d0
         c34(0)=-0.6159311930338542d0
         c34(1)=-0.2708740234375d0
         c34(2)=1.595486450195312d0
         c34(3)=-1.222724914550781d0
         c34(4)=0.7016448974609375d0
         c34(5)=-0.273223876953125d0
         c34(6)=0.06390308198474702d0
         c34(7)=-0.006760406494140625d0

         c21(-2)=0.01785714285714286d0
         c21(-1)=-0.2857142857142857d0
         c21(0)=-0.95d0
         c21(1)=2.d0
         c21(2)=-1.25d0
         c21(3)=0.6666666666666667d0
         c21(4)=-0.25d0
         c21(5)=0.05714285714285714d0
         c21(6)=-0.005952380952380952d0

         c22(-2)=0.005387987409319196d0
         c22(-1)=-0.06291721888950893d0
         c22(0)=-1.209442138671875d0
         c22(1)=1.760920206705729d0
         c22(2)=-0.7468795776367188d0
         c22(3)=0.3529510498046875d0
         c22(4)=-0.1247477213541667d0
         c22(5)=0.02752903529575893d0
         c22(6)=-0.002801622663225446d0

         c23(-2)=-0.003292410714285714d0
         c23(-1)=0.05464099702380952d0
         c23(0)=-1.1376953125d0
         c23(1)=1.1083984375d0
         c23(2)=0.008138020833333333d0
         c23(3)=-0.0490234375d0
         c23(4)=0.0244140625d0
         c23(5)=-0.006277901785714286d0
         c23(6)=0.0006975446428571429d0

         c24(-2)=-0.006760406494140625d0
         c24(-1)=0.08932364327566964d0
         c24(0)=-0.8593058268229167d0
         c24(1)=0.2970001220703125d0
         c24(2)=0.7436752319335938d0
         c24(3)=-0.3709136962890625d0
         c24(4)=0.133770751953125d0
         c24(5)=-0.0298492431640625d0
         c24(6)=0.003059423537481399d0



         c11(-3)=-0.005952380952380952d0
         c11(-2)=0.07142857142857143d0
         c11(-1)=-0.5d0
         c11(0)=-0.45d0
         c11(1)=1.25d0
         c11(2)=-0.5d0
         c11(3)=0.1666666666666667d0
         c11(4)=-0.03571428571428571d0
         c11(5)=0.003571428571428571d0

         c12(-3)=-0.002801622663225446d0
         c12(-2)=0.03060259137834821d0
         c12(-1)=-0.163775634765625d0
         c12(0)=-0.9741058349609375d0
         c12(1)=1.407915751139323d0
         c12(2)=-0.3938751220703125d0
         c12(3)=0.11761474609375d0
         c12(4)=-0.0238893054780506d0
         c12(5)=0.002314431326729911d0

         c13(-3)=0.0006975446428571429d0
         c13(-2)=-0.0095703125d0
         c13(-1)=0.07975260416666667d0
         c13(0)=-1.1962890625d0
         c13(1)=1.1962890625d0
         c13(2)=-0.07975260416666667d0
         c13(3)=0.0095703125d0
         c13(4)=-0.0006975446428571429d0
         c13(5)=0.d0

         c14(-3)=0.003059423537481399d0
         c14(-2)=-0.03429521833147321d0
         c14(-1)=0.199462890625d0
         c14(0)=-1.116297403971354d0
         c14(1)=0.6824874877929688d0
         c14(2)=0.3581878662109375d0
         c14(3)=-0.113922119140625d0
         c14(4)=0.02363150460379464d0
         c14(5)=-0.002314431326729911d0



         c01(-4)=0.003571428571428571d0
         c01(-3)=-0.0380952380952381d0
         c01(-2)=0.2d0
         c01(-1)=-0.8d0
         c01(0)=0.d0
         c01(1)=0.8d0
         c01(2)=-0.2d0
         c01(3)=0.0380952380952381d0
         c01(4)=-0.003571428571428571d0

         c02(-4)=0.002314431326729911d0
         c02(-3)=-0.02363150460379464d0
         c02(-2)=0.113922119140625d0
         c02(-1)=-0.3581878662109375d0
         c02(0)=-0.6824874877929688d0
         c02(1)=1.116297403971354d0
         c02(2)=-0.199462890625d0
         c02(3)=0.03429521833147321d0
         c02(4)=-0.003059423537481399d0

         c03(-4)=0.d0
         c03(-3)=0.0006975446428571429d0
         c03(-2)=-0.0095703125d0
         c03(-1)=0.07975260416666667d0
         c03(0)=-1.1962890625d0
         c03(1)=1.1962890625d0
         c03(2)=-0.07975260416666667d0
         c03(3)=0.0095703125d0
         c03(4)=-0.0006975446428571429d0

         c04(-4)=-0.002314431326729911d0
         c04(-3)=0.0238893054780506d0
         c04(-2)=-0.11761474609375d0
         c04(-1)=0.3938751220703125d0
         c04(0)=-1.407915751139323d0
         c04(1)=0.9741058349609375d0
         c04(2)=0.163775634765625d0
         c04(3)=-0.03060259137834821d0
         c04(4)=0.002801622663225446d0

         !do j=1,nrad
         !y(j)=0.d0
         !enddo

         j=1
         do i=-0,8
         y(j+i)=y(j+i)+c41(i)*w(1,j)+c42(i)*w(2,j)+c43(i)*w(3,j)+c44(i)*w(4,j)
         enddo

         j=2
         do i=-1,7
         y(j+i)=y(j+i)+c31(i)*w(1,j)+c32(i)*w(2,j)+c33(i)*w(3,j)+c34(i)*w(4,j)
         enddo

         j=3
         do i=-2,6
         y(j+i)=y(j+i)+c21(i)*w(1,j)+c22(i)*w(2,j)+c23(i)*w(3,j)+c24(i)*w(4,j)
         enddo

         j=4
         do i=-3,5
         y(j+i)=y(j+i)+c11(i)*w(1,j)+c12(i)*w(2,j)+c13(i)*w(3,j)+c14(i)*w(4,j)
         enddo

         do j=5,nrad-4
         do i=-4,4
         y(j+i)=y(j+i)+c01(i)*w(1,j)+c02(i)*w(2,j)+c03(i)*w(3,j)+c04(i)*w(4,j)
         enddo
         enddo

         j=nrad-3
         do i=-4,3
         y(j+i)=y(j+i)+c01(i)*w(1,j)+c02(i)*w(2,j)+c03(i)*w(3,j)+c04(i)*w(4,j)
         enddo

         j=nrad-2
         do i=-4,2
         y(j+i)=y(j+i)+c01(i)*w(1,j)+c02(i)*w(2,j)+c03(i)*w(3,j)+c04(i)*w(4,j)
         enddo

         j=nrad-1
         do i=-4,1
         y(j+i)=y(j+i)+c01(i)*w(1,j)+c02(i)*w(2,j)+c03(i)*w(3,j)+c04(i)*w(4,j)
         enddo

         j=nrad
         do i=-4,0
         y(j+i)=y(j+i)+c01(i)*w(1,j)+c02(i)*w(2,j)+c03(i)*w(3,j)+c04(i)*w(4,j)
         enddo


         end subroutine


!***THE KINETIC ENERGY IS CALCULATED ON 4*RR RADIAL GRID INSTED OF************!
!***RR RADIAL GRID FOR BETTER ACCURACY. THIS SUBROUTINE CALCULATES************!
!***THE WEIGHT FACTOR FOR THE KINETIC ENERGY OPERATOR IN 4*RR ****************!
!***RADIAL GRID. *************************************************************!
          subroutine mult(nrad,x,w)
          implicit real*8 (a-h,o-z)
          dimension x(nrad)!,y(nrad)
          dimension w(4,nrad)
          dimension c41( 0:8),c42( 0:8),c43( 0:8),c44( 0:8)
          dimension c31(-1:7),c32(-1:7),c33(-1:7),c34(-1:7)
          dimension c21(-2:6),c22(-2:6),c23(-2:6),c24(-2:6)
          dimension c11(-3:5),c12(-3:5),c13(-3:5),c14(-3:5)
          dimension c01(-4:4),c02(-4:4),c03(-4:4),c04(-4:4)


         c41(0)=-2.717857142857143d0
         c41(1)=8.d0
         c41(2)=-14.d0
         c41(3)=18.66666666666667d0
         c41(4)=-17.5d0
         c41(5)=11.2d0
         c41(6)=-4.666666666666667d0
         c41(7)=1.142857142857143d0
         c41(8)=-0.125d0

         c42(0)=-1.533237675258092d0
         c42(1)=2.732823108491443d0
         c42(2)=-2.637493896484375d0
         c42(3)=2.849429321289063d0
         c42(4)=-2.394930521647135d0
         c42(5)=1.433224487304688d0
         c42(6)=-0.57060546875d0
         c42(7)=0.1352159772600446d0
         c42(8)=-0.01442533220563616d0


         c43(0)=-0.7940848214285714d0
         c43(1)=-0.06849888392857143d0
         c43(2)=2.523763020833333d0
         c43(3)=-3.6150390625d0
         c43(4)=3.4521484375d0
         c43(5)=-2.2255859375d0
         c43(6)=0.9306640625d0
         c43(7)=-0.2283761160714286d0
         c43(8)=0.0250093005952381d0

         c44(0)=-0.359611329578218d0
         c44(1)=-1.296153477260045d0
         c44(2)=3.98780517578125d0
         c44(3)=-4.811203002929687d0
         c44(4)=4.290122985839844d0
         c44(5)=-2.665542602539062d0
         c44(6)=1.089182535807292d0
         c44(7)=-0.2630802699497768d0
         c44(8)=0.02847998482840402d0



         c31(-1)=-0.125d0
         c31(0)=-1.592857142857143d0
         c31(1)=3.5d0
         c31(2)=-3.5d0
         c31(3)=2.916666666666667d0
         c31(4)=-1.75d0
         c31(5)=0.7d0
         c31(6)=-0.1666666666666667d0
         c31(7)=0.01785714285714286d0


         c32(-1)=-0.01442533220563616d0
         c32(0)=-1.403409685407366d0
         c32(1)=2.213511149088542d0
         c32(2)=-1.425765991210937d0
         c32(3)=1.031837463378906d0
         c32(4)=-0.5773386637369792d0
         c32(5)=0.22149658203125d0
         c32(6)=-0.05129350934709821d0
         c32(7)=0.005387987409319196d0

         c33(-1)=0.0250093005952381d0
         c33(0)=-1.019168526785714d0
         c33(1)=0.8318359375d0
         c33(2)=0.4229817708333333d0
         c33(3)=-0.4638671875d0
         c33(4)=0.3009765625d0
         c33(5)=-0.1248046875d0
         c33(6)=0.03032924107142857d0
         c33(7)=-0.003292410714285714d0

         c34(-1)=0.02847998482840402d0
         c34(0)=-0.6159311930338542d0
         c34(1)=-0.2708740234375d0
         c34(2)=1.595486450195312d0
         c34(3)=-1.222724914550781d0
         c34(4)=0.7016448974609375d0
         c34(5)=-0.273223876953125d0
         c34(6)=0.06390308198474702d0
         c34(7)=-0.006760406494140625d0

         c21(-2)=0.01785714285714286d0
         c21(-1)=-0.2857142857142857d0
         c21(0)=-0.95d0
         c21(1)=2.d0
         c21(2)=-1.25d0
         c21(3)=0.6666666666666667d0
         c21(4)=-0.25d0
         c21(5)=0.05714285714285714d0
         c21(6)=-0.005952380952380952d0

         c22(-2)=0.005387987409319196d0
         c22(-1)=-0.06291721888950893d0
         c22(0)=-1.209442138671875d0
         c22(1)=1.760920206705729d0
         c22(2)=-0.7468795776367188d0
         c22(3)=0.3529510498046875d0
         c22(4)=-0.1247477213541667d0
         c22(5)=0.02752903529575893d0
         c22(6)=-0.002801622663225446d0

         c23(-2)=-0.003292410714285714d0
         c23(-1)=0.05464099702380952d0
         c23(0)=-1.1376953125d0
         c23(1)=1.1083984375d0
         c23(2)=0.008138020833333333d0
         c23(3)=-0.0490234375d0
         c23(4)=0.0244140625d0
         c23(5)=-0.006277901785714286d0
         c23(6)=0.0006975446428571429d0

         c24(-2)=-0.006760406494140625d0
         c24(-1)=0.08932364327566964d0
         c24(0)=-0.8593058268229167d0
         c24(1)=0.2970001220703125d0
         c24(2)=0.7436752319335938d0
         c24(3)=-0.3709136962890625d0
         c24(4)=0.133770751953125d0
         c24(5)=-0.0298492431640625d0
         c24(6)=0.003059423537481399d0



         c11(-3)=-0.005952380952380952d0
         c11(-2)=0.07142857142857143d0
         c11(-1)=-0.5d0
         c11(0)=-0.45d0
         c11(1)=1.25d0
         c11(2)=-0.5d0
         c11(3)=0.1666666666666667d0
         c11(4)=-0.03571428571428571d0
         c11(5)=0.003571428571428571d0

         c12(-3)=-0.002801622663225446d0
         c12(-2)=0.03060259137834821d0
         c12(-1)=-0.163775634765625d0
         c12(0)=-0.9741058349609375d0
         c12(1)=1.407915751139323d0
         c12(2)=-0.3938751220703125d0
         c12(3)=0.11761474609375d0
         c12(4)=-0.0238893054780506d0
         c12(5)=0.002314431326729911d0

         c13(-3)=0.0006975446428571429d0
         c13(-2)=-0.0095703125d0
         c13(-1)=0.07975260416666667d0
         c13(0)=-1.1962890625d0
         c13(1)=1.1962890625d0
         c13(2)=-0.07975260416666667d0
         c13(3)=0.0095703125d0
         c13(4)=-0.0006975446428571429d0
         c13(5)=0.d0

         c14(-3)=0.003059423537481399d0
         c14(-2)=-0.03429521833147321d0
         c14(-1)=0.199462890625d0
         c14(0)=-1.116297403971354d0
         c14(1)=0.6824874877929688d0
         c14(2)=0.3581878662109375d0
         c14(3)=-0.113922119140625d0
         c14(4)=0.02363150460379464d0
         c14(5)=-0.002314431326729911d0



         c01(-4)=0.003571428571428571d0
         c01(-3)=-0.0380952380952381d0
         c01(-2)=0.2d0
         c01(-1)=-0.8d0
         c01(0)=0.d0
         c01(1)=0.8d0
         c01(2)=-0.2d0
         c01(3)=0.0380952380952381d0
         c01(4)=-0.003571428571428571d0

         c02(-4)=0.002314431326729911d0
         c02(-3)=-0.02363150460379464d0
         c02(-2)=0.113922119140625d0
         c02(-1)=-0.3581878662109375d0
         c02(0)=-0.6824874877929688d0
         c02(1)=1.116297403971354d0
         c02(2)=-0.199462890625d0
         c02(3)=0.03429521833147321d0
         c02(4)=-0.003059423537481399d0

         c03(-4)=0.d0
         c03(-3)=0.0006975446428571429d0
         c03(-2)=-0.0095703125d0
         c03(-1)=0.07975260416666667d0
         c03(0)=-1.1962890625d0
         c03(1)=1.1962890625d0
         c03(2)=-0.07975260416666667d0
         c03(3)=0.0095703125d0
         c03(4)=-0.0006975446428571429d0

         c04(-4)=-0.002314431326729911d0
         c04(-3)=0.0238893054780506d0
         c04(-2)=-0.11761474609375d0
         c04(-1)=0.3938751220703125d0
         c04(0)=-1.407915751139323d0
         c04(1)=0.9741058349609375d0
         c04(2)=0.163775634765625d0
         c04(3)=-0.03060259137834821d0
         c04(4)=0.002801622663225446d0

          j=1
          w(1,j)=0.d0
          w(2,j)=0.d0
          w(3,j)=0.d0
          w(4,j)=0.d0
          do l=0,8
          w(1,j)=w(1,j)+c41(l)*x(j+l)
          w(2,j)=w(2,j)+c42(l)*x(j+l)
          w(3,j)=w(3,j)+c43(l)*x(j+l)
          w(4,j)=w(4,j)+c44(l)*x(j+l)
          enddo

          j=2
          w(1,j)=0.d0
          w(2,j)=0.d0
          w(3,j)=0.d0
          w(4,j)=0.d0
          do l=-1,7
          w(1,j)=w(1,j)+c31(l)*x(j+l)
          w(2,j)=w(2,j)+c32(l)*x(j+l)
          w(3,j)=w(3,j)+c33(l)*x(j+l)
          w(4,j)=w(4,j)+c34(l)*x(j+l)
          enddo

          j=3
          w(1,j)=0.d0
          w(2,j)=0.d0
          w(3,j)=0.d0
          w(4,j)=0.d0
          do l=-2,6
          w(1,j)=w(1,j)+c21(l)*x(j+l)
          w(2,j)=w(2,j)+c22(l)*x(j+l)
          w(3,j)=w(3,j)+c23(l)*x(j+l)
          w(4,j)=w(4,j)+c24(l)*x(j+l)
          enddo

          j=4
          w(1,j)=0.d0
          w(2,j)=0.d0
          w(3,j)=0.d0
          w(4,j)=0.d0
          do l=-3,5
          w(1,j)=w(1,j)+c11(l)*x(j+l)
          w(2,j)=w(2,j)+c12(l)*x(j+l)
          w(3,j)=w(3,j)+c13(l)*x(j+l)
          w(4,j)=w(4,j)+c14(l)*x(j+l)
          enddo



       do j=5,nrad-4
       w(1,j)=c01(-4)*x(j-4)+c01(-3)*x(j-3)+c01(-2)*x(j-2)+c01(-1)*x(j-1)+c01(0)*x(j)+c01(1)*x(j+1)+c01(2)*x(j+2)+c01(3)*x(j+3)+c01(4)*x(j+4)
       w(2,j)=c02(-4)*x(j-4)+c02(-3)*x(j-3)+c02(-2)*x(j-2)+c02(-1)*x(j-1)+c02(0)*x(j)+c02(1)*x(j+1)+c02(2)*x(j+2)+c02(3)*x(j+3)+c02(4)*x(j+4)
       w(3,j)=c03(-4)*x(j-4)+c03(-3)*x(j-3)+c03(-2)*x(j-2)+c03(-1)*x(j-1)+c03(0)*x(j)+c03(1)*x(j+1)+c03(2)*x(j+2)+c03(3)*x(j+3)+c03(4)*x(j+4)
       w(4,j)=c04(-4)*x(j-4)+c04(-3)*x(j-3)+c04(-2)*x(j-2)+c04(-1)*x(j-1)+c04(0)*x(j)+c04(1)*x(j+1)+c04(2)*x(j+2)+c04(3)*x(j+3)+c04(4)*x(j+4)
       enddo

       j=nrad-3
       w(1,j)=c01(-4)*x(j-4)+c01(-3)*x(j-3)+c01(-2)*x(j-2)+c01(-1)*x(j-1)+c01(0)*x(j)+c01(1)*x(j+1)+c01(2)*x(j+2)+c01(3)*x(j+3)
       w(2,j)=c02(-4)*x(j-4)+c02(-3)*x(j-3)+c02(-2)*x(j-2)+c02(-1)*x(j-1)+c02(0)*x(j)+c02(1)*x(j+1)+c02(2)*x(j+2)+c02(3)*x(j+3)
       w(3,j)=c03(-4)*x(j-4)+c03(-3)*x(j-3)+c03(-2)*x(j-2)+c03(-1)*x(j-1)+c03(0)*x(j)+c03(1)*x(j+1)+c03(2)*x(j+2)+c03(3)*x(j+3)
       w(4,j)=c04(-4)*x(j-4)+c04(-3)*x(j-3)+c04(-2)*x(j-2)+c04(-1)*x(j-1)+c04(0)*x(j)+c04(1)*x(j+1)+c04(2)*x(j+2)+c04(3)*x(j+3)

       j=nrad-2
       w(1,j)=c01(-4)*x(j-4)+c01(-3)*x(j-3)+c01(-2)*x(j-2)+c01(-1)*x(j-1)+c01(0)*x(j)+c01(1)*x(j+1)+c01(2)*x(j+2)
       w(2,j)=c02(-4)*x(j-4)+c02(-3)*x(j-3)+c02(-2)*x(j-2)+c02(-1)*x(j-1)+c02(0)*x(j)+c02(1)*x(j+1)+c02(2)*x(j+2)
       w(3,j)=c03(-4)*x(j-4)+c03(-3)*x(j-3)+c03(-2)*x(j-2)+c03(-1)*x(j-1)+c03(0)*x(j)+c03(1)*x(j+1)+c03(2)*x(j+2)
       w(4,j)=c04(-4)*x(j-4)+c04(-3)*x(j-3)+c04(-2)*x(j-2)+c04(-1)*x(j-1)+c04(0)*x(j)+c04(1)*x(j+1)+c04(2)*x(j+2)

       j=nrad-1
       w(1,j)=c01(-4)*x(j-4)+c01(-3)*x(j-3)+c01(-2)*x(j-2)+c01(-1)*x(j-1)+c01(0)*x(j)+c01(1)*x(j+1)
       w(2,j)=c02(-4)*x(j-4)+c02(-3)*x(j-3)+c02(-2)*x(j-2)+c02(-1)*x(j-1)+c02(0)*x(j)+c02(1)*x(j+1)
       w(3,j)=c03(-4)*x(j-4)+c03(-3)*x(j-3)+c03(-2)*x(j-2)+c03(-1)*x(j-1)+c03(0)*x(j)+c03(1)*x(j+1)
       w(4,j)=c04(-4)*x(j-4)+c04(-3)*x(j-3)+c04(-2)*x(j-2)+c04(-1)*x(j-1)+c04(0)*x(j)+c04(1)*x(j+1)

       j=nrad
       w(1,j)=c01(-4)*x(j-4)+c01(-3)*x(j-3)+c01(-2)*x(j-2)+c01(-1)*x(j-1)+c01(0)*x(j)
       w(2,j)=c02(-4)*x(j-4)+c02(-3)*x(j-3)+c02(-2)*x(j-2)+c02(-1)*x(j-1)+c02(0)*x(j)
       w(3,j)=c03(-4)*x(j-4)+c03(-3)*x(j-3)+c03(-2)*x(j-2)+c03(-1)*x(j-1)+c03(0)*x(j)
       w(4,j)=c04(-4)*x(j-4)+c04(-3)*x(j-3)+c04(-2)*x(j-2)+c04(-1)*x(j-1)+c04(0)*x(j)


!        do i=1,nrad
!        do l=1,4
!         write(11,*) -10*w(l,i)
!        enddo
!        enddo

       end  subroutine mult


        subroutine gsortho(nrad,nspol,nprinx,nprin,lmax,rw,wrky,psi,psiout)
!      Gram Schmidt orthogonalization
        implicit real*8 (a-h,o-z)
        dimension psi(nrad,nspol,nprinx,lmax+1),rw(nrad),wrky(nrad,nprinx), &
                  nprin(lmax+1,nspol),&
               psiout(nrad,nspol,nprinx,lmax+1)

      do 111,l=0,lmax
      do 111,isp=1,nspol

        do 100,iprin=1,nprin(l+1,isp)
        do j=1,nrad
        wrky(j,iprin)=psi(j,isp,iprin,l+1)*rw(j)
        enddo
        do 200,jprin=1,iprin-1
        tt=DDOT(nrad,psi(1,isp,jprin,l+1),1,wrky(1,iprin),1)
        do j=1,nrad
        psi(j,isp,iprin,l+1)=psi(j,isp,iprin,l+1)-tt*psi(j,isp,jprin,l+1)
        enddo
200     continue
        tt1=0.d0
        tt2=0.d0
        tt3=0.d0
        tt4=0.d0
        do j=1,nrad
        tt1=tt1+rw(j)*psi(j,isp,iprin,l+1)**2
        enddo
        tt=1.d0/sqrt(tt1+tt2+tt3+tt4)
        do j=1,nrad
        psiout(j,isp,iprin,l+1)=psi(j,isp,iprin,l+1)*tt
        enddo
100     continue

111     continue

        return
        end

!****** SUBROUTINE TO SET UP THE RADIAL GRID AND THE WEIGHT FACTOR ***********! 
!****** FOR INTEGRATION BASED ON INPUT PARAMETERS ****************************!
        subroutine point(xj,rj,rw1j,rw2j,c,b,bi,dr1,rist)
              !call point(1.d0,rj,rw1j,rw2j,c,b,bi,0.d0,rist)
        implicit real*8 (a-h,o-z)

        aj=xj-(rist+1.d0)
        if (aj*bi.lt.-4.d0) then
        arg=exp(aj*bi)
!       grid point
        rj=b*c*arg*(1.d0+arg*(-.5d0+arg*(.3333333333333336d0+  &
        arg*(-.25d0+arg*(.2d0+arg*(-.166666666666667d0+  &
        arg*(.1428571428571428d0+arg*(-.125d0+arg*(.1111111111111111d0+ &
        arg*(-.1d0+arg*(.9090909090909091d-1+arg))))))))))) - dr1
!       weight
        rw1j=c*arg/(1.d0+arg)
        rw2j=rj*rj*rw1j
        else
        arg=exp(-aj*bi)
!       grid point
        rj=c*(aj + b*log(1.d0+arg))-dr1
!       weight
        rw1j=c/(1.d0+arg)
        rw2j=rj*rj*rw1j
        endif
100     continue

        return
        end

!****TRIDIAGONAL OVERLAP MATRIX SSP FOR LINEAR FINITE ELEMENTS****************!
        subroutine crtssp(nrad,rr,ssp)
        implicit real*8 (a-h,o-z)
        dimension rr(nrad),ssp(2,nrad)
!  Calculates tridiagonal overlap matrix ssp for linear finite elements

        const=1.d0/30.d0
        r2=rr(1)
        ssold=0.d0
        do 7092,i=1,nrad-1
        r1=r2
        r2=rr(i+1)
        ssp(1,i)=const*((r2-r1)*( 6.d0*r1**2+3.d0*r1*r2+      r2**2) &
        +ssold)
        ssp(2,i)=const*(r2-r1)* (1.5d0*r1**2+2.d0*r1*r2+1.5d0*r2**2)
        ssold=(r2-r1)*(r1**2 + 3.d0*r1*r2 + 6.d0*r2**2)
7092    continue
        ssp(1,nrad)=const*ssold
        ssp(2,nrad)=0.d0
        return
        end
        subroutine zero(n,x)
        implicit real*8 (a-h,o-z)
        dimension x(n)
        do 10,i=1,n-1,4
        x(i)=0.d0
        x(i+1)=0.d0
        x(i+2)=0.d0
        x(i+3)=0.d0
10      continue
        istart=i
        do 11,i=istart,n
        x(i)=0.d0
11      continue
        return
        end

!*****TRIDIAGONAL MATRIX OF THE HAMILTONIAN MATRIX FOR NON SPIN CASE**********!
      SUBROUTINE tridag(n,hhp,r,u)
      INTEGER n,NMAX
      REAL*8 hhp(2,n),r(n),u(n)
!      PARAMETER (NMAX=532)
      INTEGER j
      REAL*8 bet

      real*8, allocatable :: gam(:)

      NMAX = n
      allocate( gam(NMAX) )

      if (n.gt.nmax) stop 'tridag'
      if(hhp(1,1).eq.0.d0) stop  'tridag: rewrite equations'
      bet=hhp(1,1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=hhp(2,j-1)/bet
        bet=hhp(1,j)-hhp(2,j-1)*gam(j)
        if(bet.eq.0.d0) stop  'tridag failed'
        u(j)=(r(j)-hhp(2,j-1)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue

      return
      END

!*****TRIDIAGONAL MATRIX OF THE HAMILTONIAN MATRIX FOR SPIN CASE**************!
      SUBROUTINE tridag_two(n,hhp,r,u)
      INTEGER n,NMAX
      REAL*8 hhp(2,n),r(n,2),u(n,2)
!      PARAMETER (NMAX=532)
      INTEGER j
      REAL*8 bet1,bet2

      real*8, allocatable :: gam(:,:)

      NMAX = n
      allocate( gam(NMAX,2) )

      if (n.gt.nmax) stop 'tridag'
!      if(hhp(1,1).eq.0.d0) stop  'tridag: rewrite equations'
      bet1=hhp(1,1)
      bet2=hhp(1,1)
      u(1,1)=r(1,1)/bet1
      u(1,2)=r(1,2)/bet2
      do 11 j=2,n
        gam(j,1)=hhp(2,j-1)/bet1
        gam(j,2)=hhp(2,j-1)/bet2
        bet1=hhp(1,j)-hhp(2,j-1)*gam(j,1)
        bet2=hhp(1,j)-hhp(2,j-1)*gam(j,2)
!        if(bet.eq.0.d0) stop  'tridag failed'
        u(j,1)=(r(j,1)-hhp(2,j-1)*u(j-1,1))/bet1
        u(j,2)=(r(j,2)-hhp(2,j-1)*u(j-1,2))/bet2
11    continue
      do 12 j=n-1,1,-1
        u(j,1)=u(j,1)-gam(j+1,1)*u(j+1,1)
        u(j,2)=u(j,2)-gam(j+1,2)*u(j+1,2)
12    continue

      return
      END

!*** SUBROUTINE TO SET RADIAL GRID. THE INPUT PARAMETER IS HARD CODED HERE ***!
!*** BASED ON PARAMETERS SUBTOUTINE POINT SET UP THE RADIAL GRID *************!
!*** AND OTHER FACTORS TO CALCULATE RW ***************************************!
        subroutine radgrid(nrad,rr,rw)
        implicit real*8 (a-h,o-z)
        dimension rr(nrad),rw(nrad,3)
!   the infinite intervall x={-Infinity,Infinity} is mapped onto the 
!   semi-infinite intervall 
!   y=t(x)=c*(x+(1/Log[fac])*Log[1+(1/fac)**x]={0,Infinity}.
!   The integral \int_0^Infinity f(y) dy = 
!   \int_{-Infinity}^{Infinity} f(t(x)) (dy/dx) dx 
!   is then approximated by \sum_i f(t(x_i)) tp(x_i)
!        if (nrad.ne.192) stop 'radgrid'
!        write(6,*) 'radgrid called'

! begin parameterset for radial grid. 
        rad=10.d0
!       geometric progression factor for radial grid has to be same as in RHOTOSHO
        fac=2.d0**(1.d0/6.d0)
        rist=351
        fac=sqrt(fac)
!       ratio of geometric progression versus origin
        bi=log(fac)
        b=1.d0/bi
 !       grid spacing far away from origin: rad/10
        c=rad*.05d0
        call point(1.d0,r1,rw1j,rw2j,c,b,bi,0.d0,rist)
        dr1=r1
! end parameterset for radial grid

        do 100,j=1,nrad
        xj=j
        call point(xj,rj,rw1j,rw2j,c,b,bi,dr1,rist)
        rr(j)=rj
        rw(j,1)=rw1j
        rw(j,2)=rw2j
        rw(j,3)=1.d0/rw1j
!        write(*,*) rj,rw1j,rw2j
100        continue

! modify weights at en point
        rw(1,2)=rw(1,2)*17.d0/48.d0
        rw(2,2)=rw(2,2)*59.d0/48.d0
        rw(3,2)=rw(3,2)*43.d0/48.d0
        rw(4,2)=rw(4,2)*49.d0/48.d0

        rw(nrad-0,2)=rw(nrad-0,2)*17.d0/48.d0
        rw(nrad-1,2)=rw(nrad-1,2)*59.d0/48.d0
        rw(nrad-2,2)=rw(nrad-2,2)*43.d0/48.d0
        rw(nrad-3,2)=rw(nrad-3,2)*49.d0/48.d0
        
        return
        end

!*** SUBROUTINE TO SET 4*RR RADIAL GRID FOR THE KINETIC ENERGY CALCULATION ***!
!*** FOLLOWS SAME PROTOCOL AS SUBROUTINE RADGRID *****************************!
        subroutine radgrid4(nrad,rr4,rw4)
        implicit real*8 (a-h,o-z)
        dimension rr4(4,nrad),rw4(4,nrad,3)
!   the infinite intervall x={-Infinity,Infinity} is mapped onto the 
!   semi-infinite intervall 
!   y=t(x)=c*(x+(1/Log[fac])*Log[1+(1/fac)**x]={0,Infinity}.
!   The integral \int_0^Infinity f(y) dy = 
!   \int_{-Infinity}^{Infinity} f(t(x)) (dy/dx) dx 
!   is then approximated by \sum_i f(t(x_i)) tp(x_i)
!        if (nrad.ne.192) stop 'radgrid4'
        !!write(6,*) 'radgrid called'

! begin parameterset for radial grid. 
        rad=10.d0
!       geometric progression factor for radial grid 
        fac=2.d0**(1.d0/6.d0)
        rist=351
        fac=sqrt(fac)
!       ratio of geometric progression versus origin
        bi=log(fac)
        b=1.d0/bi
!       grid spacing far away from origin: rad/10
        c=rad*.05d0
        call point(1.d0,r1,rw1j,rw2j,c,b,bi,0.d0,rist)
        dr1=r1
! end parameterset for radial grid

        do 100,j=1,nrad
        do 100,i=1,4
        xj=j+.25d0*(i-1)
        call point(xj,rj,rw1j,rw2j,c,b,bi,dr1,rist)
        rr4(i,j)=rj
        rw4(i,j,1)=rw1j*.25d0
        rw4(i,j,2)=rw2j*.25d0
        rw4(i,j,3)=1.d0/rw1j
!        write(6,*) j,rj,rw1j,rw2j
 100        continue

! modify weights at en point
        rw4(1,1,2)=rw4(1,1,2)*17.d0/48.d0
        rw4(2,1,2)=rw4(2,1,2)*59.d0/48.d0
        rw4(3,1,2)=rw4(3,1,2)*43.d0/48.d0
        rw4(4,1,2)=rw4(4,1,2)*49.d0/48.d0

        rw4(4,nrad,2)=rw4(1,nrad,2)*17.d0/48.d0
        rw4(3,nrad,2)=rw4(2,nrad,2)*59.d0/48.d0
        rw4(2,nrad,2)=rw4(3,nrad,2)*43.d0/48.d0
        rw4(1,nrad,2)=rw4(4,nrad,2)*49.d0/48.d0
        
        return
        end


        subroutine loggrid(nrad,rrmin,f,rr,rw)
        use localorb_io, only: localorb_info
        implicit real*8 (a-h,o-z)
        dimension rr(nrad),rw(nrad,3)
        character*300 :: info_str
        write(info_str,'(2X,A)') 'loggrid called'
        call localorb_info( info_str )

!         rrmin=1.d-3  ! grid point closest to origin
!         f=1.05d0  ! geometric progression factor
        do 100,j=1,nrad
        xj=j-1
        rj=rrmin*exp(xj*log(f))
        rr(j)=rj-rrmin
        rw1j=rj*log(f)
        rw(j,1)=rw1j
        rw(j,2)=rw1j*rr(j)**2
        rw(j,3)=1.d0/rw1j
!        write(*,*) rj,rw1j,rw2j
100        continue

! modify weights at en point
        rw(1,2)=rw(1,2)*17.d0/48.d0
        rw(2,2)=rw(2,2)*59.d0/48.d0
        rw(3,2)=rw(3,2)*43.d0/48.d0
        rw(4,2)=rw(4,2)*49.d0/48.d0

        rw(nrad-0,2)=rw(nrad-0,2)*17.d0/48.d0
        rw(nrad-1,2)=rw(nrad-1,2)*59.d0/48.d0
        rw(nrad-2,2)=rw(nrad-2,2)*43.d0/48.d0
        rw(nrad-3,2)=rw(nrad-3,2)*49.d0/48.d0

        return
        end

        subroutine loggrid4(nrad,rrmin,f,rr4,rw4)
        use localorb_io, only: localorb_info
        implicit real*8 (a-h,o-z)
        dimension rr4(4,nrad),rw4(4,nrad,3)
        character*300 :: info_str
        write(info_str,'(2X,A)') 'loggrid4 called'
        call localorb_info( info_str )
        
!         rrmin=1.d-3  ! grid point closest to origin
!         f=1.05d0  ! geometric progression factor

        do 100,j=1,nrad
        do 100,i=1,4
        xj=(j-1)+.25d0*(i-1)
        rj=rrmin*exp(xj*log(f))
        rr4(i,j)=rj-rrmin
        rw1j=rj*log(f)*.25d0
        rw4(i,j,1)=rw1j
        rw4(i,j,2)=rw1j*rr4(i,j)**2
        rw4(i,j,3)=0.25d0/rw1j
!        write(*,*) j,rj,rw1j,rw2j
100        continue

! modify weights at en point
        rw4(1,1,2)=rw4(1,1,2)*17.d0/48.d0
        rw4(2,1,2)=rw4(2,1,2)*59.d0/48.d0
        rw4(3,1,2)=rw4(3,1,2)*43.d0/48.d0
        rw4(4,1,2)=rw4(4,1,2)*49.d0/48.d0

        rw4(4,nrad,2)=rw4(1,nrad,2)*17.d0/48.d0
        rw4(3,nrad,2)=rw4(2,nrad,2)*59.d0/48.d0
        rw4(2,nrad,2)=rw4(3,nrad,2)*43.d0/48.d0
        rw4(1,nrad,2)=rw4(4,nrad,2)*49.d0/48.d0

        return
        end



!**** POISSON SOLVER FOR A GIVEN CHARGE DENSITY ******************************!
!**** WITH ASSUMPTION OF SPHERICAL SYMMETRY **********************************!
        subroutine poissoll(nrad,lmaxp,lmax,ccleb,ppois,rho,rw,pot)
!       Solves the Poisson equation for the radial charge distribution
!       rho
!       multiplied by all spherical harmonics up to lmax and multiplies
!       the
!       result with the radial weight rw(*,2)
!        parameter(nradx=532)
        implicit real*8 (a-h,o-z)
        dimension ppois(nrad,-(lmaxp+1):lmaxp,2),  &
                  rho(nrad),rw(nrad),pot(nrad,lmaxp+1),ccleb(lmaxp+1)

        integer :: nradx
        real*8, allocatable :: a(:)
        real*8, allocatable :: b(:)

        nradx = nrad
        allocate( a(nradx) )
        allocate( b(nradx) )
!        if (nrad.gt.nradx) stop 'nrad>nradx'

   do l=0,lmax
   if (ccleb(l+1).ne.0.d0) then

!  outward integration
        j=1
        sum=0.d0
        a(j)=sum

        j=2
        ym1=ppois(j-1,l,2)*rho(j-1)
        y0 =ppois(j  ,l,2)*rho(j  )
        yp1=ppois(j+1,l,2)*rho(j+1)
        yp2=ppois(j+2,l,2)*rho(j+2)
        yp3=ppois(j+3,l,2)*rho(j+3)
        yp4=ppois(j+4,l,2)*rho(j+4)
        yp5=ppois(j+5,l,2)*rho(j+5)
        yp6=ppois(j+6,l,2)*rho(j+6)
        sum=sum+  &
        0.3042245370370371d0*ym1 - 1.006919642857143d0*yp1 +   &
         1.017964616402116d0*yp2 - 0.7320353835978836d0*yp3 +   &
         0.3430803571428571d0*yp4 - 0.0938409391534391d0*yp5 +   &
         0.01136739417989418d0*yp6 + 1.156159060846561d0*y0
        a(j)=sum

        j=3
        ym2=ppois(j-2,l,2)*rho(j-2)
        ym1=ppois(j-1,l,2)*rho(j-1)
        y0 =ppois(j  ,l,2)*rho(j  )
        yp1=ppois(j+1,l,2)*rho(j+1)
        yp2=ppois(j+2,l,2)*rho(j+2)
        yp3=ppois(j+3,l,2)*rho(j+3)
        yp4=ppois(j+4,l,2)*rho(j+4)
        yp5=ppois(j+5,l,2)*rho(j+5)
        sum=sum+  &
        0.3951636904761905d0*ym1 - 0.01136739417989418d0*ym2 -   &
         0.3703455687830688d0*yp1 + 0.2222470238095238d0*yp2 -   &
         0.0954613095238095d0*yp3 + 0.02479332010582011d0*yp4 -   &
         0.002901785714285714d0*yp5 + 0.837872023809524d0*y0
        a(j)=sum

        j=4
        ym3=ppois(j-3,l,2)*rho(j-3)
        ym2=ppois(j-2,l,2)*rho(j-2)
        ym1=ppois(j-1,l,2)*rho(j-1)
        y0 =ppois(j  ,l,2)*rho(j  )
        yp1=ppois(j+1,l,2)*rho(j+1)
        yp2=ppois(j+2,l,2)*rho(j+2)
        yp3=ppois(j+3,l,2)*rho(j+3)
        yp4=ppois(j+4,l,2)*rho(j+4)
        sum=sum+  &
        0.4764136904761905d0*ym1 - 0.0345816798941799d0*ym2 +   &
         0.002901785714285714d0*ym3 - 0.1672205687830688d0*yp1 +   &
         0.05974702380952382d0*yp2 - 0.01421130952380952d0*yp3 +   &
         0.001579034391534391d0*yp4 + 0.6753720238095238d0*y0
        a(j)=sum


      do j=5,nrad-3
        ym4=ppois(j-4,l,2)*rho(j-4)
        ym3=ppois(j-3,l,2)*rho(j-3)
        ym2=ppois(j-2,l,2)*rho(j-2)
        ym1=ppois(j-1,l,2)*rho(j-1)
        y0 =ppois(j  ,l,2)*rho(j  )
        yp1=ppois(j+1,l,2)*rho(j+1)
        yp2=ppois(j+2,l,2)*rho(j+2)
        yp3=ppois(j+3,l,2)*rho(j+3)
        sum=sum+  &
        0.5648396164021165d0*ym1 - 0.07879464285714287d0*ym2 +   &
         0.01553406084656085d0*ym3 - 0.001579034391534391d0*ym4 -   &
         0.07879464285714287d0*yp1 + 0.01553406084656085d0*yp2 -   &
         0.001579034391534391d0*yp3 + 0.5648396164021165d0*y0
        a(j)=sum
      enddo

        j=nrad-2
        ym4=ppois(j-4,l,2)*rho(j-4)
        ym3=ppois(j-3,l,2)*rho(j-3)
        ym2=ppois(j-2,l,2)*rho(j-2)
        ym1=ppois(j-1,l,2)*rho(j-1)
        y0 =ppois(j  ,l,2)*rho(j  )
        yp1=ppois(j+1,l,2)*rho(j+1)
        yp2=ppois(j+2,l,2)*rho(j+2)
        yp3=0.d0
        sum=sum+  &
        0.5648396164021165d0*ym1 - 0.07879464285714287d0*ym2 +   &
         0.01553406084656085d0*ym3 - 0.001579034391534391d0*ym4 -   &
         0.07879464285714287d0*yp1 + 0.01553406084656085d0*yp2 -   &
         0.001579034391534391d0*yp3 + 0.5648396164021165d0*y0
        a(j)=sum

        j=nrad-1
        ym4=ppois(j-4,l,2)*rho(j-4)
        ym3=ppois(j-3,l,2)*rho(j-3)
        ym2=ppois(j-2,l,2)*rho(j-2)
        ym1=ppois(j-1,l,2)*rho(j-1)
        y0 =ppois(j  ,l,2)*rho(j  )
        yp1=ppois(j+1,l,2)*rho(j+1)
        yp2=0.d0
        yp3=0.d0
        sum=sum+  &
        0.5648396164021165d0*ym1 - 0.07879464285714287d0*ym2 +   &
         0.01553406084656085d0*ym3 - 0.001579034391534391d0*ym4 -   &
         0.07879464285714287d0*yp1 + 0.01553406084656085d0*yp2 -   &
         0.001579034391534391d0*yp3 + 0.5648396164021165d0*y0
        a(j)=sum

        j=nrad
        ym4=ppois(j-4,l,2)*rho(j-4)
        ym3=ppois(j-3,l,2)*rho(j-3)
        ym2=ppois(j-2,l,2)*rho(j-2)
        ym1=ppois(j-1,l,2)*rho(j-1)
        y0 =ppois(j  ,l,2)*rho(j  )
        yp1=0.d0
        yp2=0.d0
        yp3=0.d0
        sum=sum+  &
        0.5648396164021165d0*ym1 - 0.07879464285714287d0*ym2 +   &
         0.01553406084656085d0*ym3 - 0.001579034391534391d0*ym4 -   &
         0.07879464285714287d0*yp1 + 0.01553406084656085d0*yp2 -   &
         0.001579034391534391d0*yp3 + 0.5648396164021165d0*y0
        a(j)=sum

! inward integration

        j=nrad
        sum=0.d0
        b(j)=sum

        j=nrad-1
        ym3=ppois(j-3,-(l+1),2)*rho(j-3)
        ym2=ppois(j-2,-(l+1),2)*rho(j-2)
        ym1=ppois(j-1,-(l+1),2)*rho(j-1)
        y0 =ppois(j  ,-(l+1),2)*rho(j  )
        yp1=ppois(j+1,-(l+1),2)*rho(j+1)
        yp2=0.d0
        yp3=0.d0
        yp4=0.d0
        sum=sum  &
        -0.07879464285714287d0*ym1 + 0.01553406084656085d0*ym2 -   &
         0.001579034391534391d0*ym3 + 0.5648396164021165d0*yp1 -   &
         0.07879464285714287d0*yp2 + 0.01553406084656085d0*yp3 -   &
         0.001579034391534391d0*yp4 + 0.5648396164021165d0*y0
        b(j)=sum

        j=nrad-2
        ym3=ppois(j-3,-(l+1),2)*rho(j-3)
        ym2=ppois(j-2,-(l+1),2)*rho(j-2)
        ym1=ppois(j-1,-(l+1),2)*rho(j-1)
        y0 =ppois(j  ,-(l+1),2)*rho(j  )
        yp1=ppois(j+1,-(l+1),2)*rho(j+1)
        yp2=ppois(j+2,-(l+1),2)*rho(j+2)
        yp3=0.d0
        yp4=0.d0
        sum=sum  &
        -0.07879464285714287d0*ym1 + 0.01553406084656085d0*ym2 -   &
         0.001579034391534391d0*ym3 + 0.5648396164021165d0*yp1 -   &
         0.07879464285714287d0*yp2 + 0.01553406084656085d0*yp3 -   &
         0.001579034391534391d0*yp4 + 0.5648396164021165d0*y0
        b(j)=sum

        j=nrad-3
        ym3=ppois(j-3,-(l+1),2)*rho(j-3)
        ym2=ppois(j-2,-(l+1),2)*rho(j-2)
        ym1=ppois(j-1,-(l+1),2)*rho(j-1)
        y0 =ppois(j  ,-(l+1),2)*rho(j  )
        yp1=ppois(j+1,-(l+1),2)*rho(j+1)
        yp2=ppois(j+2,-(l+1),2)*rho(j+2)
        yp3=ppois(j+3,-(l+1),2)*rho(j+3)
        yp4=0.d0
        sum=sum  &
        -0.07879464285714287d0*ym1 + 0.01553406084656085d0*ym2 -   &
         0.001579034391534391d0*ym3 + 0.5648396164021165d0*yp1 -   &
         0.07879464285714287d0*yp2 + 0.01553406084656085d0*yp3 -   &
         0.001579034391534391d0*yp4 + 0.5648396164021165d0*y0
        b(j)=sum


        do 9321,j=nrad-4,4,-1
        ym3=ppois(j-3,-(l+1),2)*rho(j-3)
        ym2=ppois(j-2,-(l+1),2)*rho(j-2)
        ym1=ppois(j-1,-(l+1),2)*rho(j-1)
        y0 =ppois(j  ,-(l+1),2)*rho(j  )
        yp1=ppois(j+1,-(l+1),2)*rho(j+1)
        yp2=ppois(j+2,-(l+1),2)*rho(j+2)
        yp3=ppois(j+3,-(l+1),2)*rho(j+3)
        yp4=ppois(j+4,-(l+1),2)*rho(j+4)
        sum=sum  &
        -0.07879464285714287d0*ym1 + 0.01553406084656085d0*ym2 -   &
         0.001579034391534391d0*ym3 + 0.5648396164021165d0*yp1 -   &
         0.07879464285714287d0*yp2 + 0.01553406084656085d0*yp3 -   &
         0.001579034391534391d0*yp4 + 0.5648396164021165d0*y0
        b(j)=sum
9321        continue

        j=3
        ym2=ppois(j-2,-(l+1),2)*rho(j-2)
        ym1=ppois(j-1,-(l+1),2)*rho(j-1)
        y0 =ppois(j  ,-(l+1),2)*rho(j  )
        yp1=ppois(j+1,-(l+1),2)*rho(j+1)
        yp2=ppois(j+2,-(l+1),2)*rho(j+2)
        yp3=ppois(j+3,-(l+1),2)*rho(j+3)
        yp4=ppois(j+4,-(l+1),2)*rho(j+4)
        yp5=ppois(j+5,-(l+1),2)*rho(j+5)
        sum=sum  &
        -0.0345816798941799d0*ym1 + 0.002901785714285714d0*ym2 +   &
         0.6753720238095238d0*yp1 - 0.1672205687830688d0*yp2 +   &
         0.05974702380952382d0*yp3 - 0.01421130952380952d0*yp4 +   &
         0.001579034391534391d0*yp5 + 0.4764136904761905d0*y0
        b(j)=sum

        j=2
        ym1=ppois(j-1,-(l+1),2)*rho(j-1)
        y0 =ppois(j  ,-(l+1),2)*rho(j  )
        yp1=ppois(j+1,-(l+1),2)*rho(j+1)
        yp2=ppois(j+2,-(l+1),2)*rho(j+2)
        yp3=ppois(j+3,-(l+1),2)*rho(j+3)
        yp4=ppois(j+4,-(l+1),2)*rho(j+4)
        yp5=ppois(j+5,-(l+1),2)*rho(j+5)
        yp6=ppois(j+6,-(l+1),2)*rho(j+6)
        sum=sum  &
        -0.01136739417989418d0*ym1 + 0.837872023809524d0*yp1 -   &
         0.3703455687830688d0*yp2 + 0.2222470238095238d0*yp3 -   &
         0.0954613095238095d0*yp4 + 0.02479332010582011d0*yp5 -   &
         0.002901785714285714d0*yp6 + 0.3951636904761905d0*y0
        b(j)=sum

        j=1
        y0 =ppois(j  ,-(l+1),2)*rho(j  )
        yp1=ppois(j+1,-(l+1),2)*rho(j+1)
        yp2=ppois(j+2,-(l+1),2)*rho(j+2)
        yp3=ppois(j+3,-(l+1),2)*rho(j+3)
        yp4=ppois(j+4,-(l+1),2)*rho(j+4)
        yp5=ppois(j+5,-(l+1),2)*rho(j+5)
        yp6=ppois(j+6,-(l+1),2)*rho(j+6)
        yp7=ppois(j+7,-(l+1),2)*rho(j+7)
        sum=sum+  &
        1.156159060846561d0*yp1 - 1.006919642857143d0*yp2 +   &
         1.017964616402116d0*yp3 - 0.7320353835978836d0*yp4 +   &
         0.3430803571428571d0*yp5 - 0.0938409391534391d0*yp6 +   &
         0.01136739417989418d0*yp7 + 0.3042245370370371d0*y0
        b(j)=sum

        do j=1,nrad
        pot(j,l+1)=(ppois(j,-(l+1),1)*a(j)+ppois(j,l,1)*b(j))*rw(j)
        enddo
   else

        do j=1,nrad
        pot(j,l+1)=0.d0
        enddo

   endif
   enddo

        return
        end


!**** DETERMINE THE CUT OFF POINT FOR PLOTTING A FUNCTION ********************!
!> Determine NPSP
subroutine detnp(nn,r,rad0,npsp)
   implicit none
   !Arguments
   integer, intent(in) :: nn
   real(kind=8), dimension(nn), intent(in) :: r
   real(kind=8), intent(in) :: rad0
   integer, intent(out) :: npsp
   !Local variables
   real(kind=8) :: rmin
   integer :: i
   rmin=1.d10
   do i=1,nn
      if (abs(r(i)-rad0) < rmin) then
      rmin=abs(r(i)-rad0)
      npsp=i
      end if
   end do
end subroutine detnp


!**** ESTIMATE THE COEFFICIENTS FOR THE POSSION SOLVER TO BE USED ************!
!**** BY SUBROUTINE POISSOL **************************************************!
       subroutine crtpois(nrad,lmaxp,rr,rw,ppois)
           implicit real*8 (a-h,o-z)
           dimension rr(nrad),rw(nrad),ppois(nrad,-(lmaxp+1):lmaxp,2)

        do 134,l=0,lmaxp
        const=16.d0*atan(1.d0)/(2*l+1)
        j=1

        if (l.eq.0) then
        ppois(j,l,1)=const
        else
        ppois(j,l,1)=0.d0
        endif

        ppois(j,l,2)=0.d0

        ppois(j,-(l+1),1)=0.d0

        if (-l+1.eq.0) then
        ppois(j,-(l+1),2)=rw(j)
        else
        ppois(j,-(l+1),2)=0.d0
        endif
134        continue

        do 234,l=0,lmaxp
        const=16.d0*atan(1.d0)/(2*l+1)
        do 234,j=2,nrad
        ppois(j,l,1)=const*rr(j)**l
        ppois(j,l,2)=rw(j)*rr(j)**(l+2)
        ppois(j,-(l+1),1)=const*rr(j)**(-l-1)
        ppois(j,-(l+1),2)=rw(j)*rr(j)**(-l+1)
234        continue

        return
        end


!*** CANONICAL TRANSFORMATION OF THE ORBITALS ********************************!
        subroutine  canonical(nrad,nspol,nprinx,nprin,lmax,epslag,work,psi,grad,grads,eval)
!   transform to canonical orbitals (i.e. eigenorbitals)
        implicit real*8 (a-h,o-z)
        dimension psi(nrad,nspol,nprinx,lmax+1), grad(nrad,nspol,nprinx,lmax+1),  & 
                  grads(nrad,nspol,nprinx,lmax+1),work(nrad,nprinx),nprin(lmax+1,nspol), & 
                  epslag(nprinx,nprinx,nspol,lmax+1),worklp(200),eval(nprinx,nspol,lmax+1)

              if (nprinx.gt.10) stop 'eval too small'
      do 111,l=0,lmax
      do 111,isp=1,nspol


       !do iprin=1,nprin(l+1,isp)
       !write(*,'(10(1x,e16.9))') (epslag(iprin,jprin,isp,l+1),jprin=1,nprin(l+1,isp))
       !enddo


        call DSYEV('V','L',nprin(l+1,isp),epslag(1,1,isp,l+1),nprinx,eval(1,isp,l+1),worklp,200,info)
        if (info.ne.0) stop 'Can DSYEV'
        !if(nprin(l+1,isp).gt.0) write(*,'(a,10(1x,e12.5))') 'can. EVs',(eval(i,isp,l+1),i=1,nprin(l+1,isp))

 ! transform wavefunctions
        do iprin=1,nprin(l+1,isp)
        do j=1,nrad
        work(j,iprin)=0.d0
        enddo
        enddo

        do iprin=1,nprin(l+1,isp)
        do jprin=1,nprin(l+1,isp)
        do j=1,nrad
        work(j,iprin)=work(j,iprin)+psi(j+0,isp,jprin,l+1)*epslag(jprin,iprin,isp,l+1)
        enddo
        enddo
        enddo

        do iprin=1,nprin(l+1,isp)
        do j=1,nrad
        psi(j,isp,iprin,l+1)=work(j,iprin)
        enddo
        enddo

 ! transform grad
        do iprin=1,nprin(l+1,isp)
        do j=1,nrad
        work(j,iprin)=0.d0
        enddo
        enddo

        do iprin=1,nprin(l+1,isp)
        do jprin=1,nprin(l+1,isp)
        do j=1,nrad
        work(j,iprin)=work(j,iprin)+grad(j,isp,jprin,l+1)*epslag(jprin,iprin,isp,l+1)
        enddo
        enddo
        enddo

        do iprin=1,nprin(l+1,isp)
        do j=1,nrad
        grad(j,isp,iprin,l+1)=work(j,iprin)
        enddo
        enddo

 ! transform grads
        do iprin=1,nprin(l+1,isp)
        do j=1,nrad
        work(j,iprin)=0.d0
        enddo
        enddo

        do iprin=1,nprin(l+1,isp)
        do jprin=1,nprin(l+1,isp)
        do j=1,nrad
        work(j,iprin)=work(j,iprin)+grads(j,isp,jprin,l+1)*epslag(jprin,iprin,isp,l+1)
        enddo
        enddo
        enddo

        do iprin=1,nprin(l+1,isp)
        do j=1,nrad
        grads(j,isp,iprin,l+1)=work(j,iprin)
        enddo
        enddo

111     continue

        return
        end


!** TRIDIAGONAL SYMMETRIC PRECONDITIONING HAMILTONIAN MATRIX *****************!
        subroutine crthhp(nspol,nprinx,nprin,nrad,lmax,rr,znuc,pot0,hhp)
!  Calculates a tridiagonal symmetric preconditioning hamiltonian matrix
!  hhp in
!  a basis of linear finite elements
        implicit real*8 (a-h,o-z)
        dimension  rr(nrad),hhp(2,nrad,nprinx,lmax+1),pot0(nrad,nprinx,lmax+1),nprin(lmax+1,nspol)

   do l=0,lmax
        const=.5d0*l*(l+1)*.333333333333333d0
   do iprin=1,max(nprin(l+1,1),nprin(l+1,2))

! tt kinetic energy
        ttold=0.d0
! ss SC potential energy
        ssold=0.d0
! ee potential energy from z/r
        eeold=0.d0
! potential from angular term
        aold=0.d0
        r2=rr(1)
        q=pot0(1,iprin,l+1)
        do 5092,i=1,nrad-1
        aa=const*(rr(i+1)-rr(i))
        r1=r2
        r2=rr(i+1)
        dr=r2-r1
        r1q=r1*r1
        r2q=r2*r2
        r1r2=r1*r2
        p=q
        q=pot0(i+1,iprin,l+1)
        tt=(r1q + r1r2 + r2q)/dr*.16666666666666d0
        ss=dr*(10.d0*p*r1q + 2.d0*q*r1q + 4.d0*p*r1r2 +   &
                     2.d0*q*r1r2 + p*r2q + q*r2q)*.16666666666666d-1
        ssof=(r2q*r2*(p + 2.d0*q) - r1q*r1*(2.d0*p + q)   &
              - q*r1q*r2 + p*r1*r2q )*.16666666666666d-1
        ee=dr*(3*r1 + r2)*.83333333333333d-1
        eeof=dr*(r1 + r2)*.83333333333333d-1

        hhp(1,i,iprin,l+1)=(tt+ttold - znuc*(ee+eeold)) + ss+ssold + (aa+aold)
        hhp(2,i,iprin,l+1)=(-tt - znuc*eeof) + ssof + .5d0*aa
        aold=aa

        ttold=tt
        ssold =dr*  &
             (    p*r1q +      q*r1q +  2.d0*p*r1r2 +   &
             4.d0*q*r1r2 + 2.d0*p*r2q + 10.d0*q*r2q)*.16666666666666d-1
        eeold=dr*(r1 + 3*r2)*.83333333333333d-1
5092    continue
        hhp(1,nrad,iprin,l+1)=ttold + ssold + aold
        hhp(2,nrad,iprin,l+1)=0.d0

   enddo
   enddo
        return
        end


        subroutine checkorthogonalization(nrad,nspol,nprinx,nprin,lmax,rw,wrky,psi)
        use localorb_io, only: localorb_info
        implicit real*8 (a-h,o-z)
        dimension  psi(nrad,nspol,nprinx,lmax+1),rw(nrad),wrky(nrad,nprinx),nprin(lmax+1,nspol)
        character*300 :: info_str
      do 111,l=0,lmax
      do 111,isp=1,nspol

        do 100,iprin=1,nprin(l+1,isp)
        do j=1,nrad !-3,4
        wrky(j,iprin)=psi(j,isp,iprin,l+1)*rw(j)
        !wrky(j+1,iprin)=psi(j+1,isp,iprin,l+1)*rw(j+1)
        !wrky(j+2,iprin)=psi(j+2,isp,iprin,l+1)*rw(j+2)
        !wrky(j+3,iprin)=psi(j+3,isp,iprin,l+1)*rw(j+3)
        enddo
        do jprin=1,nprin(l+1,isp)
        tt=DDOT(nrad,psi(1,isp,jprin,l+1),1,wrky(1,iprin),1)
        write(info_str,'(2X,A,4(I3),1X,E12.5)') 'ORTHO: l,isp,iprin,jprin',l,isp,iprin,jprin,tt
        call localorb_info( info_str )
        enddo
100     continue

111     continue

        return
        end


!*** CONVOLUTION OF THE SR/LR KERNEL FOR THE HYBRID FUNCTIONAL ***************!
        subroutine convolut_kernel(nrad,rr,rw,width,convkernel)
        implicit real*8 (a-h,o-z)
        dimension rr(nrad),rw(nrad),convkernel(nrad,nrad)

        twopi=8.d0*atan(1.d0)
        fact=twopi/(width*sqrt(twopi))**3
        !fact=1.d0*fact*width**2
        do i=1,nrad
        ri=rr(i)/width
        do j=1,nrad
        rj=rr(j)/width
        if (ri*rj .eq. 0.d0) then
        convkernel(i,j)=2.d0*fact*exp(-.5d0*ri**2)*exp(-.5d0*rj**2)*rw(j)
        else if (ri*rj .lt. 600.d0) then
        convkernel(i,j)=2.d0*fact*exp(-.5d0*ri**2)*(sinh(ri*rj)/(ri*rj))*exp(-.5d0*rj**2)*rw(j)
        else
        convkernel(i,j)=fact*(exp(-.5d0*(ri-rj)**2)-exp(-.5d0*(ri+rj)**2))*rw(j)/(ri*rj)
        endif
        enddo
        enddo

        return
        end

!**** CREATION OF THE SCREENING KERNEL******************************!
      subroutine create_screenexkernel(nrad,rr,rw,lmax,lmaxp,omega,screenexkernel)
          implicit real*8 (a-h,o-z)
          !!parameter(nradxx=232)
!          parameter(nradxx=532)

          dimension clebone(lmaxp+1), screenexkernel(nrad,nrad,lmaxp+1)

          real*8, dimension(nrad)   :: rr
          real*8, dimension(nrad,3) :: rw

          integer :: nradxx
!          real*8, allocatable :: rr(:)
!          real*8, allocatable :: rw(:,:)
          real*8, allocatable :: rhoij(:)
          real*8, allocatable :: smoothrhoij(:)
          real*8, allocatable :: pothf(:,:)
          real*8, allocatable :: convkernel(:,:) ! Kernel for the convolution with aGaussian
          real*8, allocatable :: ppois(:,:,:) !auxiliary array for solving Poisson's equation

          nradxx = nrad
!          allocate( rr(nradxx) )
!          allocate( rw(nradxx,3) )
          allocate( rhoij(nradxx) )
          allocate( smoothrhoij(nradxx) )
          allocate( pothf(nradxx,lmaxp+1) )
          allocate( convkernel(nradxx,nradxx) )
          allocate( ppois(nradxx,-(lmaxp+1):lmaxp,2) )

!            if (mod(nradxx,4).ne.0) stop 'inconsistent loop unrolling'
            do l=0,lmaxp
            clebone(l+1)=1.d0
            enddo

!             call  radgrid(nradxx,rr,rw)

             width=sqrt(.5d0)/omega
             call convolut_kernel(nradxx,rr,rw(1,2),width,convkernel)
             call crtpois(nradxx,lmaxp,rr,rw(1,1),ppois)

           if (nrad.gt.nradxx) stop 'nrad> nradxx in create_screenexkernel'

          do i=1,nrad
             do j=1,nradxx
             rhoij(j)=0.d0
             enddo
             rhoij(i)=1.d0

             call DGEMV('N', nradxx, nradxx, 1.d0, convkernel, nradxx, rhoij, 1,0.d0, smoothrhoij, 1)
             do j=1,nradxx
             rhoij(j)=rhoij(j) -smoothrhoij(j)
             enddo
             call poissoll(nradxx,lmaxp,lmaxp,clebone,ppois,rhoij,rw(1,2),pothf)
             do l=0,lmaxp
             do j=1,nrad
             screenexkernel(j,i,l+1)=pothf(j,l+1)
             enddo
             enddo

           enddo
           ! symmetrize matrix
             do l=0,lmaxp
             do j=1,nrad
             do i=j+1,nrad
             t1=screenexkernel(i,j,l+1)
             t2=screenexkernel(j,i,l+1)
             screenexkernel(i,j,l+1)=.5d0*(t1+t2)
             screenexkernel(j,i,l+1)=.5d0*(t1+t2)
             enddo
             enddo
             enddo

           end subroutine


!*** CALCULATE THE RESULT OF KINETIC ENERGY OPERATOR *************************!
!*** AND LOCAL POTENTIAL *****************************************************!
        subroutine applykinpot(nrad,rr,rw,rw4,kinonly,zion,ll,potloc,x,y)
        implicit real*8 (a-h,o-z)
!        parameter(nradx=532)
        logical kinonly
        dimension x(nrad),y(nrad),rr(nrad),potloc(nrad)
        dimension rw4(4*nrad,3),rw(nrad,3)

        integer :: nradx
        real*8, allocatable :: w(:)

        nradx = nrad
        allocate( w(4*nradx) )
        
!        if (nrad.gt.nradx) stop 'nrad>nradx'


         call mult(nrad,x,w)
         do i=1,4*nrad
         w(i)=.5d0*rw4(i,2)*w(i)*rw4(i,3)**2
         enddo
         if (kinonly) then
         do i=1,nrad
         y(i)=(.5d0*ll*(ll+1))*rw(i,1)*x(i)
         enddo
         else
         do i=1,nrad
         y(i)=( ( (.5d0*ll*(ll+1)) - zion*rr(i) )*rw(i,1) + potloc(i)*rw(i,2))*x(i)
         enddo
         endif
         call mult_t(nrad,w,y)

        end subroutine


        subroutine applyzorakinpot(nrad,rr,rw,rw4,kinonly,zion,ll,potloc,damprel,damprel4, x,y)
          implicit real*8 (a-h,o-z)
          !parameter(nradx=532) ! Set largest value of nradx
          !parameter(nradx=1236) ! Largest value from For loggrid
          parameter(c=137.0359895d0)
          logical kinonly
          dimension x(nrad),y(nrad),rr(nrad),potloc(nrad), damprel(nrad),damprel4(4*nrad)
          dimension w(4*nrad),rw4(4*nrad,3),rw(nrad,3)
         !if (nrad.gt.nradx) stop 'nrad>nradx'
          
          call mult(nrad,x,w)
        
          do i=1,4*nrad
             w(i)=damprel4(i)*rw4(i,2)*w(i)*rw4(i,3)**2
             !w(i)=.5d0*rw4(i,2)*w(i)*rw4(i,3)**2
          enddo
        
          if (kinonly) then
            do i=1,nrad
               y(i)=damprel(i)*ll*(ll+1)*rw(i,1)*x(i)
              !y(i)=(.5d0*ll*(ll+1))*rw(i,1)*x(i)
            enddo
          else
            do i=1,nrad
               y(i)=( ( damprel(i)*ll*(ll+1) - zion*rr(i))*rw(i,1) + potloc(i)*rw(i,2))*x(i)
              !y(i)=( ( (.5d0*ll*(ll+1)) - zion*rr(i) )*rw(i,1) +
              !potloc(i)*rw(i,2))*x(i)
            enddo
          endif

          call mult_t(nrad,w,y)
        end subroutine applyzorakinpot


        subroutine intrapolate_nrad(pmatrix, pmatrix4, nrad)
          !subroutine intrapolate_nrad(pmatrix, pmatrix4, nrad, pmatrixname)
          implicit none
        
          integer, intent(in) :: nrad
          real(8), dimension(nrad), intent(in) :: pmatrix
        
          !character(len=7) :: pmatrixname, filenumber
          !character(len=27) :: filename
          integer :: i
          !integer, save :: j
          !data j /0/
        
          real(8), dimension(4*nrad), intent(out) :: pmatrix4
        
          !j = j+1
        
          pmatrix4(1) = pmatrix(1)
          pmatrix4(2) = (77.d0*pmatrix(1) + 77.d0*pmatrix(2) - 33.d0*pmatrix(3) + 7.d0*pmatrix(4)) / 128.d0
          pmatrix4(3) = (5.d0*pmatrix(1) + 15.d0*pmatrix(2) - 5.d0*pmatrix(3) + pmatrix(4)) / 16.d0
          pmatrix4(4) = (15.d0*pmatrix(1) + 135.d0*pmatrix(2) - 27.d0*pmatrix(3) + 5.d0*pmatrix(4)) / 128.d0
        
          do i = 1, nrad -3
            pmatrix4(4*i+1) = pmatrix(i+1)
            pmatrix4(4*i+2) = (-7.d0*pmatrix(i) + 105.d0*pmatrix(i+1) + 35.d0*pmatrix(i+2) &
                              - 5.d0*pmatrix(i+3)) / 128.d0
            pmatrix4(4*i+3) = (-pmatrix(i) + 9.d0*pmatrix(i+1) + 9.d0*pmatrix(i+2) - pmatrix(i+3)) &
                              / 16.d0
            pmatrix4(4*i+4) = (-5.d0*pmatrix(i) + 35.d0*pmatrix(i+1) + 105.d0*pmatrix(i+2) &
                              - 7.d0*pmatrix(i+3)) / 128.d0
          end do
        
          pmatrix4(4*nrad-7) = pmatrix(nrad-1)
          pmatrix4(4*nrad-6) = (15.d0*pmatrix(nrad-3) + 135.d0*pmatrix(nrad-2) - 27.d0*pmatrix(nrad-1) + 5.d0*pmatrix(nrad)) & 
                               / 128.d0
          pmatrix4(4*nrad-5) = (5.d0*pmatrix(nrad-3) + 15.d0*pmatrix(nrad-2) - 5.d0*pmatrix(nrad-1) + pmatrix(nrad)) &
                               / 16.d0
          pmatrix4(4*nrad-4) = (77.d0*pmatrix(nrad-3) + 77.d0*pmatrix(nrad-2) - 33.d0*pmatrix(nrad-1) + 7.d0*pmatrix(nrad)) &
                               / 128.d0
          pmatrix4(4*nrad-3) = pmatrix(nrad)
          pmatrix4(4*nrad-2) = pmatrix(nrad)
          pmatrix4(4*nrad-1) = pmatrix(nrad)
          pmatrix4(4*nrad)   = pmatrix(nrad)
        
        end subroutine intrapolate_nrad


        subroutine zerotail(lmax,nspol,nprinx,nprin,nrad,nprb,psi)
          implicit real*8 (a-h,o-z)
          dimension nprin(lmax+1,nspol),psi(nrad,nspol,nprinx,lmax+1)
          do l=0,lmax
          do isp=1,nspol
          do iprin=1,nprin(l+1,isp)
          do j=nprb,nrad
          psi(j,isp,iprin,l+1)=0.d0
          enddo
          enddo
          enddo
          enddo
       end subroutine zerotail
