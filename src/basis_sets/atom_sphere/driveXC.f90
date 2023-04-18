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



        subroutine initlibXC(method,nspol)
!       Initalizes libXC once.

!       Method 1 is backwards compatible: Teter 93 LDA 
!       Methods 0 and 2 are without libXC, so nothing is done here.
!       Method 3,4 and 5 are abbreviations for PBE
!                        and polarized Teter93 LSDA or PBE.

!       General ixc and nspol may be accessed by simply
!       concatenating the two integers.
!       Example: Spin polarized VWN LDA:
!       ixc = - 001 007 nspol=2  -> method = -0010072

        use libxcModule
        integer:: method,ixc,nspol,family
        real*8 omega,hfmix
        nspol=1

        write(*,'(a)')' ************************************'
!       general case
        if(method<0)then
           ixc=method/10
           if(abs(method-10*ixc).eq.2) nspol=2
           !!write(*,*)"method-10*ixc",method-10*ixc
           write(*,'(a)')' *** method variable is negative  ***'
           write(*,'(a)')' *** method   ->     ixc   nspol  ***' 
           write(*,'(a,i8.7,a,i7.6,6x,i2.1,a)')&
                         ' *** ',method,&
                                      ' -> ',ixc,  nspol,'  ***'
           family=-ixc/1000
           if ((family.lt.100).and.(nspol.eq.1)) then                               
           write(*,*)'***      LDA (closed shell)      ***'
           else if ((family.lt.100).and.(nspol.eq.2)) then                               
           write(*,*)'***      LDA (spin polariz)      ***'
           else if ((family.ge.100).and.(family.lt.200).and.  &
              (nspol.eq.1)) then                               
           write(*,*)'***      PBE (closed shell)      ***'
           else if ((family.ge.100).and.(family.lt.200).and.   &
               (nspol.eq.2)) then                               
           write(*,*)'***      PBE (spin polariz)      ***'
           else if ((family.ge.200).and.(family.lt.300).and.   &
               (nspol.eq.1)) then   
           write(*,*)'***     MGGA (closed shell)      ***'                           
           else if ((family.ge.200).and.(family.lt.300).and.   &
               (nspol.eq.2)) then                               
           write(*,*)'***     MGGA (spin polariz)      ***'
           else if ((family.ge.300).and.(family.lt.400).and.   &
               (nspol.eq.1)) then                               
           write(*,*)'***      LCA (closed shell)      ***'
           else if ((family.ge.300).and.(family.lt.400).and.   &
               (nspol.eq.2)) then                               
           write(*,*)'***      LCA (spin polariz)      ***'
           else if ((family.gt.400).and.(nspol.eq.1)) then                               
           write(*,*)'***    HYB_GGA (closed shell)    ***'
           else if ((family.gt.400).and.(nspol.eq.2)) then                               
           write(*,*)'***    HYB_GGA (spin polariz)    ***'
           end if
           

           call libxc_functionals_init(ixc,nspol,omega,hfmix)
        end if

!       predefined cases
        select case(method)
           case(0)
           write(*,*)'***   TEST case (method=0)       ***'
           case(1)
           write(*,*)'*** LDA (Teter 93  closed shell) ***'
           call libxc_functionals_init(-20,1,omega,hfmix)
           case(2)
           write(*,*)'***       GHF method             ***'
           case(3)
           write(*,*)'***     PBE (in closed shell)    ***'
           call libxc_functionals_init(-101130,1,omega,hfmix)
           case(4)
           write(*,*)'*** LSDA (Teter 93 spin polariz) ***'
!          This requires some modifications in the atomic program
           call libxc_functionals_init(-20,2,omega,hfmix)
           nspol=2
           case(5)
           write(*,*)'***    PBE (spin  polarized)     ***'
!          This requires some modifications in the atomic program
           call libxc_functionals_init(-101130,2,omega,hfmix)
           nspol=2
        end select
        write(*,'(a)')' ************************************'

!       below convention may change
        !if(method/=2 .and. method /= 0) method = 1

        end subroutine initlibXC

!        subroutine driveXC(nrad,nang,nspol,rw,aw,rhor,potential,eexcu)
!!       same as mikeu in atom.f, but calls libXC routines.

!!       The radial part of the integration over the angular
!!       grid is done in below routines that support libXCs
!!       spin polarized GGAs. This routine here just does
!!       the outer loop over the angular part ignoring
!!       angular derivatives for now.

!        use libxcModule
!        implicit real*8 (a-h,o-z)


!!       arguments from calling ldagga (former ldapot) routine 
!        dimension  rhor(nrad,nang,3),potential(nrad,nang,nspol),rw(nrad,3),aw(nang)

!!       arguments to radial integration routine
!        real(8):: rd(nrad), enexc, pot(nrad,nspol)

!!       local variables: None for now

!!       set up the weights for integrating the derivatives
!!       according to the conventions of the radial routines:
!        do irad=1,nrad
!           rd(irad) = 1d0 / rw(irad,1)
!          rw(irad) =       rw(irad,2)
!        end do

!!       initialize exc to zero, vxc is added to pot from psolver
!        eexcu = 0d0

!!       in the polarized case, we will pass rhor(:,:, 1:2),
!!       otherwise we need the total charge  rhor(:,:, 3:3),
!!       let this be                         rhor(:,:,i1:i2), then
!        i1=5-2*nspol
!        i2=4-nspol

!        do iang=1,nang
!              call driveXCradial(nspol,nrad,rw(:,2),rd,rhor(:,iang,i1:i2),enexc,pot)

!!             outer integration over angular component
!              eexcu=eexcu+enexc*aw(iang)

!!             the radial vxc term pot is added to the potential 
!              potential(:,iang,:)=potential(:,iang,:)+pot(:,:)
!!             White Bird terms for angular derivatives?!?
!!             This will require a scheme to correct the gradients
!!             before passing them to libXC, too!
!        end do

!        return
!        end subroutine driveXC



subroutine driveXCradial(nspol,nrad,rw,rd,rho,enexc,pot)
   ! ****************************************************************
   ! This is essentially the same wrapper as used by pseudos gatom.f
   ! and the spherically symmetric implementation of atom.f
   ! ****************************************************************
   use libxcModule
   implicit real*8 (a-h,o-z)
   !Dummy arguments
   real(kind=8) :: rho(nrad,nspol),rw(nrad),rd(nrad),pot(nrad,nspol)

   if (libxc_functionals_isgga()) then
      select case(nspol)
         case (1)
         call driveGGAsimple(nspol,nrad,rw,rd,rho,enexc,pot)
         !write(*,*)"calling GGAsimple"
         case (2)
         call driveGGApolarized(nspol,nrad,rw,rd,rho,enexc,pot)
         !write(*,*)"calling GGApolarized"
      end select
   else
      select case(nspol)
         case (1)
         call driveLDAsimple(nspol,nrad,rw,rd,rho,enexc,pot)
         !write(*,*)"calling LDAsimple"
         case (2)
         call driveLDApolarized(nspol,nrad,rw,rd,rho,enexc,pot)
         !write(6,*)"calling LDApolarized"
      end select
   end if

   return
end subroutine driveXCradial


subroutine driveLDAsimple(nspol,nrad,rw,rd,rho,enexc,pot)
   ! ****************************************************************
   ! greatly simplified version of ggaenergy_15 for LDA XC
   ! ****************************************************************
   use libxcModule
   implicit real*8 (a-h,o-z)
   real(kind=8) :: rho(nrad,nspol),rw(nrad),rd(nrad),pot(nrad,nspol) ! dummy arguments
   real(kind=8) :: Exc,Vxc(nspol),dEdg(nspol), grad(nspol)  ! arguments to XCfunction

   enexc=0.d0
   pot=0d0
   grad=0d0

   do j=1,nrad
      call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
      enexc=enexc+Exc*rw(j)*rho(j,1)
      pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   enddo


!   do j=2,nrad
!      pot(j,1)=pot(j,1)/rw(j)
!   enddo
!   pot(1,1)=pot(2,1)

   return
end subroutine


subroutine driveLDApolarized(nspol,nrad,rw,rd,rho,enexc,pot)
   ! ****************************************************************
   !   the simple LDA XC driver generalized for two spin channels
   ! ****************************************************************
   use libxcModule
   implicit real*8 (a-h,o-z)
   real(kind=8) :: rho(nrad,nspol),rw(nrad),rd(nrad),pot(nrad,nspol) ! dummy arguments
   real(kind=8) :: Exc,Vxc(nspol),dEdg(nspol), grad(nspol),rhoupdown(2)      ! arguments to XCfunction

   enexc=0.d0
   pot=0d0
   grad=0d0

   do j=1,nrad
      rhoupdown(1)=rho(j,1)
      rhoupdown(2)=rho(j,2)
      call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
     ! call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
      enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
      pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
      pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   enddo

!   do j=2,nrad
!      pot(j,1)=pot(j,1)/rw(j)
!      pot(j,2)=pot(j,2)/rw(j)
!   enddo
!   pot(1,:)=pot(2,:)

   return
end subroutine



subroutine driveGGAsimple(nspol,nrad,rw,rd,rho,enexc,pot)
   ! ****************************************************************
   ! this one is very close to the original version of ggaenergy_15 
   ! ****************************************************************
   use libxcModule
   implicit real*8 (a-h,o-z)
   real(8):: rho(nrad,nspol),rw(nrad),rd(nrad),pot(nrad,nspol),& ! dummy arguments
             Exc,Vxc(nspol),dEdg(nspol), grad(nspol),&      ! arguments to XCfunction
             c(-8:8), sign(nspol) !                           ! local varialbes for derivatives 

   enexc=0.d0
   pot=0d0
   grad=0d0

   j=1
    c(0)=-2.717857142857143d0
    c(1)=8.d0
    c(2)=-14.d0
    c(3)=18.66666666666667d0
    c(4)=-17.5d0
    c(5)=11.2d0
    c(6)=-4.666666666666667d0
    c(7)=1.142857142857143d0
    c(8)=-0.125d0
   grad=0.d0
   do i=-0,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-0,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=2
    c(-1)=-0.1111111111111111d0
    c(0)=-1.717857142857143d0
    c(1)=4.d0
    c(2)=-4.666666666666667d0
    c(3)=4.666666666666667d0
    c(4)=-3.5d0
    c(5)=1.866666666666667d0
    c(6)=-0.6666666666666666d0
    c(7)=0.1428571428571428d0
    c(8)=-0.01388888888888889d0
   grad=0.d0
   do i=-1,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-1,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=3
    c(-2)=0.01111111111111111d0
    c(-1)=-0.2222222222222222d0
    c(0)=-1.217857142857143d0
    c(1)=2.666666666666666d0
    c(2)=-2.333333333333333d0
    c(3)=1.866666666666667d0
    c(4)=-1.166666666666667d0
    c(5)=0.5333333333333333d0
    c(6)=-0.1666666666666666d0
    c(7)=0.03174603174603174d0
    c(8)=-0.2777777777777778d-2
   grad=0.d0
   do i=-2,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-2,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=4
    c(-3)=-0.202020202020202d-2
    c(-2)=0.03333333333333333d0
    c(-1)=-0.3333333333333333d0
    c(0)=-0.88452380952381d0
    c(1)=2.d0
    c(2)=-1.4d0
    c(3)=0.933333333333333d0
    c(4)=-0.5d0
    c(5)=0.2d0
    c(6)=-0.05555555555555556d0
    c(7)=0.952380952380952d-2
    c(8)=-0.7575757575757577d-3
   grad=0.d0
   do i=-3,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-3,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=5
    c(-4)=0.5050505050505051d-3
    c(-3)=-0.808080808080808d-2
    c(-2)=0.06666666666666666d0
    c(-1)=-0.4444444444444445d0
    c(0)=-0.6345238095238095d0
    c(1)=1.6d0
    c(2)=-0.933333333333333d0
    c(3)=0.5333333333333333d0
    c(4)=-0.25d0
    c(5)=0.0888888888888889d0
    c(6)=-0.02222222222222222d0
    c(7)=0.3463203463203463d-2
    c(8)=-0.2525252525252525d-3
   grad=0.d0
   do i=-4,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-4,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=6
    c(-5)=-0.1554001554001554d-3
    c(-4)=0.2525252525252525d-2
    c(-3)=-0.0202020202020202d0
    c(-2)=0.1111111111111111d0
    c(-1)=-0.5555555555555556d0
    c(0)=-0.4345238095238095d0
    c(1)=1.333333333333333d0
    c(2)=-0.6666666666666666d0
    c(3)=0.3333333333333333d0
    c(4)=-0.1388888888888889d0
    c(5)=0.04444444444444445d0
    c(6)=-0.0101010101010101d0
    c(7)=0.1443001443001443d-2
    c(8)=-0.971250971250971d-4
   grad=0.d0
   do i=-5,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-5,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=7
    c(-6)=0.555000555000555d-4
    c(-5)=-0.932400932400932d-3
    c(-4)=0.7575757575757577d-2
    c(-3)=-0.04040404040404041d0
    c(-2)=0.1666666666666666d0
    c(-1)=-0.6666666666666666d0
    c(0)=-0.2678571428571428d0
    c(1)=1.142857142857143d0
    c(2)=-0.5d0
    c(3)=0.2222222222222222d0
    c(4)=-0.0833333333333333d0
    c(5)=0.02424242424242424d0
    c(6)=-0.5050505050505051d-2
    c(7)=0.6660006660006659d-3
    c(8)=-0.4162504162504162d-4
   grad=0.d0
   do i=-6,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-6,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=8
    c(-7)=-0.222000222000222d-4
    c(-6)=0.3885003885003884d-3
    c(-5)=-0.3263403263403263d-2
    c(-4)=0.01767676767676768d0
    c(-3)=-0.07070707070707071d0
    c(-2)=0.2333333333333333d0
    c(-1)=-0.7777777777777778d0
    c(0)=-0.125d0
    c(1)=1.d0
    c(2)=-0.3888888888888889d0
    c(3)=0.1555555555555556d0
    c(4)=-0.05303030303030303d0
    c(5)=0.01414141414141414d0
    c(6)=-0.2719502719502719d-2
    c(7)=0.3330003330003329d-3
    c(8)=-0.1942501942501942d-4
   grad=0.d0
   do i=-7,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-7,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo


    c(-8)=9.71250971250971d-6
    c(-7)=-0.1776001776001776d-3
    c(-6)=0.1554001554001554d-2
    c(-5)=-0.87024087024087d-2
    c(-4)=0.3535353535353535d-1
    c(-3)=-0.1131313131313131d0
    c(-2)=0.3111111111111111d0
    c(-1)=-0.888888888888889d0
    c(0)=0.d0
    c(1)=0.888888888888889d0
    c(2)=-0.3111111111111111d0
    c(3)=0.1131313131313131d0
    c(4)=-0.3535353535353535d-1
    c(5)=0.87024087024087d-2
    c(6)=-0.1554001554001554d-2
    c(7)=0.1776001776001776d-3
    c(8)=-9.71250971250971d-6
   do 100,j=9,nrad-8
   grad=0.d0
   do i=-8,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo
100   continue

   j=nrad-7
   grad=0.d0
   do i=-8,7
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,7
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-6
   grad=0.d0
   do i=-8,6
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,6
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-5
   grad=0.d0
   do i=-8,5
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,5
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-4
   grad=0.d0
   do i=-8,4
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,4
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-3
   grad=0.d0
   do i=-8,3
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,3
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-2
   grad=0.d0
   do i=-8,2
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,2
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-1
   grad=0.d0
   do i=-8,1
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,1
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-0
   grad=0.d0
   do i=-8,0
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,0
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

!   do j=2,nrad
!      pot(j,1)=pot(j,1)/rw(j)
!   enddo
!   pot(1,1)=pot(2,1)

   return
end subroutine




subroutine driveGGApolarized(nspol,nrad,rw,rd,rho,enexc,pot)
   ! ****************************************************************
   !    ggaenergy_15 generalized to colinear spin polarization
   ! ****************************************************************
   use libxcModule
   implicit real*8 (a-h,o-z)
   real(8):: rho(nrad,nspol),rw(nrad),rd(nrad),pot(nrad,nspol),& ! dummy arguments
             Exc,Vxc(nspol),dEdg(nspol), grad(nspol),&      ! arguments to XCfunction
             c(-8:16), sign(nspol),rhoupdown(2) !                           ! local varialbes for derivatives 

   enexc=0.d0
   pot=0d0
   grad=0d0

   j=1
    !c(0)=-2.717857142857143d0
    !c(1)=8.d0
    !c(2)=-14.d0
    !c(3)=18.66666666666667d0
    !c(4)=-17.5d0
    !c(5)=11.2d0
    !c(6)=-4.666666666666667d0
    !c(7)=1.142857142857143d0
    !c(8)=-0.125d0
         c(0)=-3.380728993228993d0
         c(1)=16.d0
         c(2)=-60.d0
         c(3)=186.6666666666667d0
         c(4)=-455.d0
         c(5)=873.6d0
         c(6)=-1334.666666666667d0
         c(7)=1634.285714285714d0
         c(8)=-1608.75d0
         c(9)=1271.111111111111d0
         c(10)=-800.8d0
         c(11)=397.0909090909091d0
         c(12)=-151.6666666666667d0
         c(13)=43.07692307692308d0
         c(14)=-8.571428571428571d0
         c(15)=1.066666666666667d0
         c(16)=-0.0625d0
   grad=0.d0
   do i=-0,16
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-0,16
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=2
    !c(-1)=-0.1111111111111111d0
    !c(0)=-1.717857142857143d0
    !c(1)=4.d0
    !c(2)=-4.666666666666667d0
    !c(3)=4.666666666666667d0
    !c(4)=-3.5d0
    !c(5)=1.866666666666667d0
    !c(6)=-0.6666666666666666d0
    !c(7)=0.1428571428571428d0
    !c(8)=-0.01388888888888889d0
         c(-1)=-0.0625d0
         c(0)=-2.318228993228993d0
         c(1)=7.5d0
         c(2)=-17.5d0
         c(3)=37.91666666666667d0
         c(4)=-68.25d0
         c(5)=100.1d0
         c(6)=-119.1666666666667d0
         c(7)=114.9107142857143d0
         c(8)=-89.375d0
         c(9)=55.61111111111111d0
         c(10)=-27.3d0
         c(11)=10.34090909090909d0
         c(12)=-2.916666666666667d0
         c(13)=0.5769230769230769d0
         c(14)=-0.07142857142857143d0
         c(15)=0.4166666666666667d-2
   grad=0.d0
   do i=-1,15
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-1,15
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=3
    !c(-2)=0.01111111111111111d0
    !c(-1)=-0.2222222222222222d0
    !c(0)=-1.217857142857143d0
    !c(1)=2.666666666666666d0
    !c(2)=-2.333333333333333d0
    !c(3)=1.866666666666667d0
    !c(4)=-1.166666666666667d0
    !c(5)=0.5333333333333333d0
    !c(6)=-0.1666666666666666d0
    !c(7)=0.03174603174603174d0
    !c(8)=-0.2777777777777778d-2
         c(-2)=0.4166666666666667d-2
         c(-1)=-0.1333333333333333d0
         c(0)=-1.751562326562327d0
         c(1)=4.666666666666667d0
         c(2)=-7.583333333333333d0
         c(3)=12.13333333333333d0
         c(4)=-16.68333333333333d0
         c(5)=19.06666666666667d0
         c(6)=-17.875d0
         c(7)=13.61904761904762d0
         c(8)=-8.341666666666667d0
         c(9)=4.044444444444444d0
         c(10)=-1.516666666666667d0
         c(11)=0.4242424242424242d0
         c(12)=-0.08333333333333333d0
         c(13)=0.01025641025641026d0
         c(14)=-0.5952380952380952d-3
   grad=0.d0
   do i=-2,14
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-2,14
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=4
    !c(-3)=-0.202020202020202d-2
    !c(-2)=0.03333333333333333d0
    !c(-1)=-0.3333333333333333d0
    !c(0)=-0.88452380952381d0
    !c(1)=2.d0
    !c(2)=-1.4d0
    !c(3)=0.933333333333333d0
    !c(4)=-0.5d0
    !c(5)=0.2d0
    !c(6)=-0.05555555555555556d0
    !c(7)=0.952380952380952d-2
    !c(8)=-0.7575757575757577d-3
         c(-3)=-0.5952380952380952d-3
         c(-2)=0.01428571428571429d0
         c(-1)=-0.2142857142857143d0
         c(0)=-1.346800421800422d0
         c(1)=3.25d0
         c(2)=-3.9d0
         c(3)=4.766666666666667d0
         c(4)=-5.107142857142857d0
         c(5)=4.596428571428571d0
         c(6)=-3.404761904761905d0
         c(7)=2.042857142857143d0
         c(8)=-0.975d0
         c(9)=0.3611111111111111d0
         c(10)=-0.1d0
         c(11)=0.01948051948051948d0
         c(12)=-0.2380952380952381d-2
         c(13)=0.1373626373626374d-3
   grad=0.d0
   do i=-3,13
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-3,13
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=5
    !c(-4)=0.5050505050505051d-3
    !c(-3)=-0.808080808080808d-2
    !c(-2)=0.06666666666666666d0
    !c(-1)=-0.4444444444444445d0
    !c(0)=-0.6345238095238095d0
    !c(1)=1.6d0
    !c(2)=-0.933333333333333d0
    !c(3)=0.5333333333333333d0
    !c(4)=-0.25d0
    !c(5)=0.0888888888888889d0
    !c(6)=-0.02222222222222222d0
    !c(7)=0.3463203463203463d-2
    !c(8)=-0.2525252525252525d-3
         c(-4)=0.1373626373626374d-3
         c(-3)=-0.293040293040293d-2
         c(-2)=0.03296703296703297d0
         c(-1)=-0.3076923076923077d0
         c(0)=-1.019877344877345d0
         c(1)=2.4d0
         c(2)=-2.2d0
         c(3)=2.095238095238095d0
         c(4)=-1.767857142857143d0
         c(5)=1.257142857142857d0
         c(6)=-0.7333333333333333d0
         c(7)=0.3428571428571429d0
         c(8)=-0.125d0
         c(9)=0.03418803418803419d0
         c(10)=-0.6593406593406593d-2
         c(11)=0.7992007992007992d-3
         c(12)=-0.4578754578754579d-4
   grad=0.d0
   do i=-4,12
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-4,12
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=6
    !c(-5)=-0.1554001554001554d-3
    !c(-4)=0.2525252525252525d-2
    !c(-3)=-0.0202020202020202d0
    !c(-2)=0.1111111111111111d0
    !c(-1)=-0.5555555555555556d0
    !c(0)=-0.4345238095238095d0
    !c(1)=1.333333333333333d0
    !c(2)=-0.6666666666666666d0
    !c(3)=0.3333333333333333d0
    !c(4)=-0.1388888888888889d0
    !c(5)=0.04444444444444445d0
    !c(6)=-0.0101010101010101d0
    !c(7)=0.1443001443001443d-2
    !c(8)=-0.971250971250971d-4
         c(-5)=-0.4578754578754579d-4
         c(-4)=0.9157509157509158d-3
         c(-3)=-0.9157509157509158d-2
         c(-2)=0.0641025641025641d0
         c(-1)=-0.4166666666666667d0
         c(0)=-0.7365440115440115d0
         c(1)=1.833333333333333d0
         c(2)=-1.30952380952381d0
         c(3)=0.9821428571428571d0
         c(4)=-0.6547619047619048d0
         c(5)=0.3666666666666667d0
         c(6)=-0.1666666666666667d0
         c(7)=0.05952380952380952d0
         c(8)=-0.01602564102564103d0
         c(9)=0.3052503052503053d-2
         c(10)=-0.3663003663003663d-3
         c(11)=0.2081252081252081d-4
   grad=0.d0
   do i=-5,11
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-5,11
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=7
    !c(-6)=0.555000555000555d-4
    !c(-5)=-0.932400932400932d-3
    !c(-4)=0.7575757575757577d-2
    !c(-3)=-0.04040404040404041d0
    !c(-2)=0.1666666666666666d0
    !c(-1)=-0.6666666666666666d0
    !c(0)=-0.2678571428571428d0
    !c(1)=1.142857142857143d0
    !c(2)=-0.5d0
    !c(3)=0.2222222222222222d0
    !c(4)=-0.0833333333333333d0
    !c(5)=0.02424242424242424d0
    !c(6)=-0.5050505050505051d-2
    !c(7)=0.6660006660006659d-3
    !c(8)=-0.4162504162504162d-4
         c(-6)=0.2081252081252081d-4
         c(-5)=-0.3996003996003996d-3
         c(-4)=0.3746253746253746d-2
         c(-3)=-0.02331002331002331d0
         c(-2)=0.1136363636363636d0
         c(-1)=-0.5454545454545455d0
         c(0)=-0.478968253968254d0
         c(1)=1.428571428571429d0
         c(2)=-0.8035714285714286d0
         c(3)=0.4761904761904762d0
         c(4)=-0.25d0
         c(5)=0.1090909090909091d0
         c(6)=-0.03787878787878788d0
         c(7)=0.999000999000999d-2
         c(8)=-0.1873126873126873d-2
         c(9)=0.222000222000222d-3
         c(10)=-0.1248751248751249d-4
   grad=0.d0
   do i=-6,10
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-6,10
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=8
    !c(-7)=-0.222000222000222d-4
    !c(-6)=0.3885003885003884d-3
    !c(-5)=-0.3263403263403263d-2
    !c(-4)=0.01767676767676768d0
    !c(-3)=-0.07070707070707071d0
    !c(-2)=0.2333333333333333d0
    !c(-1)=-0.7777777777777778d0
    !c(0)=-0.125d0
    !c(1)=1.d0
    !c(2)=-0.3888888888888889d0
    !c(3)=0.1555555555555556d0
    !c(4)=-0.05303030303030303d0
    !c(5)=0.01414141414141414d0
    !c(6)=-0.2719502719502719d-2
    !c(7)=0.3330003330003329d-3
    !c(8)=-0.1942501942501942d-4
         c(-7)=-0.1248751248751249d-4
         c(-6)=0.2331002331002331d-3
         c(-5)=-0.2097902097902098d-2
         c(-4)=0.01223776223776224d0
         c(-3)=-0.05303030303030303d0
         c(-2)=0.1909090909090909d0
         c(-1)=-0.7d0
         c(0)=-0.2361111111111111d0
         c(1)=1.125d0
         c(2)=-0.5d0
         c(3)=0.2333333333333333d0
         c(4)=-0.09545454545454545d0
         c(5)=0.03181818181818182d0
         c(6)=-0.8158508158508159d-2
         c(7)=0.1498501498501499d-2
         c(8)=-0.1748251748251748d-3
         c(9)=9.712509712509713d-6
   grad=0.d0
   do i=-7,9
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-7,9
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo


    c(-8)=9.71250971250971d-6
    c(-7)=-0.1776001776001776d-3
    c(-6)=0.1554001554001554d-2
    c(-5)=-0.87024087024087d-2
    c(-4)=0.3535353535353535d-1
    c(-3)=-0.1131313131313131d0
    c(-2)=0.3111111111111111d0
    c(-1)=-0.888888888888889d0
    c(0)=0.d0
    c(1)=0.888888888888889d0
    c(2)=-0.3111111111111111d0
    c(3)=0.1131313131313131d0
    c(4)=-0.3535353535353535d-1
    c(5)=0.87024087024087d-2
    c(6)=-0.1554001554001554d-2
    c(7)=0.1776001776001776d-3
    c(8)=-9.71250971250971d-6
   do 100,j=9,nrad-8
   grad=0.d0
   do i=-8,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo
100   continue

   j=nrad-7
   grad=0.d0
   do i=-8,7
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,7
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-6
   grad=0.d0
   do i=-8,6
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,6
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-5
   grad=0.d0
   do i=-8,5
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,5
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-4
   grad=0.d0
   do i=-8,4
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,4
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-3
   grad=0.d0
   do i=-8,3
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,3
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-2
   grad=0.d0
   do i=-8,2
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,2
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-1
   grad=0.d0
   do i=-8,1
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,1
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-0
   grad=0.d0
   do i=-8,0
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   !write(44,*) j,grad(1),grad(2)
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   !write(45,*) j,grad(1),grad(2)
   rhoupdown(1)=rho(j,1)
   rhoupdown(2)=rho(j,2)
   if (rhoupdown(1).lt.1.d-16)grad(1)=0.d0
   if (rhoupdown(2).lt.1.d-16)grad(2)=0.d0
   call xcfunction(nspol,rhoupdown,grad,Exc,Vxc,dEdg)
!   call xcfunction(nspol,rho(j,:),grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,0
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

!   do j=2,nrad
!      pot(j,1)=pot(j,1)/rw(j)
!      pot(j,2)=pot(j,2)/rw(j)
!   enddo
!   pot(1,:)=pot(2,:)

!   do j=1,nrad
!      write(22,*) j,rho(j,1),rho(j,2),rw(j)
!      write(55,*) j,pot(j,1),pot(j,2)
!   enddo

   return
end subroutine




