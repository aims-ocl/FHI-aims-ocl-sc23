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


! libXC interfaces for pseudo.

! TEST version derived from BigDFTs libABINIT/src/56_xc/m_libxc_functionals.F90 


module libxcModule

  use xc_f90_types_m
  use libxc_funcs_m
  use xc_f90_lib_m
  use localorb_io, only: localorb_info
  use mpi_tasks, only: aims_stop 

  implicit none

  type libxc_functional
    private
    integer         :: family ! LDA, GGA, etc.
    integer         :: id     ! identifier

    type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
    type(xc_f90_pointer_t) :: info ! information about the functional
  end type libxc_functional

  type(libxc_functional) :: funcs(2)

  private
  public :: libxc_functionals_init, &
&      libxc_functionals_isgga, &
&      libxc_functionals_ismgga, &
&      libxc_functionals_exctXfac, &
&      libxc_functionals_end, &
!      libxc_functionals_getvxc, &
!      is replaced by a simple routine
&      xcfunction

contains
!!*** 

!!****f* libxc_functionals/libxc_functionals_init
!! NAME
!!  libxc_functionals_init
!!
!! FUNCTION
!!  Initialize the desired XC functional, from LibXC.
!!  * Call the LibXC initializer
!!  * Fill preliminary fields in module structures.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      xc_f90_gga_exc_vxc,xc_f90_gga_vxc,xc_f90_lda_exc_vxc,xc_f90_lda_vxc
!!      xc_f90_mgga_exc_vxc,xc_f90_mgga_vxc
!!
!! SOURCE

  subroutine libxc_functionals_init(ixc,nspden,omega,hfmix)



    implicit none

!Arguments ------------------------------------
!scalars

    integer, intent(in) :: nspden
    integer, intent(in) :: ixc
    real(8) :: omega,hfmix,alpha,beta

!Local variables-------------------------------
!scalars

    character*300 :: info_str
    integer :: i
    ! LibXC v4.0.2
    ! real(8), dimension (1:3) :: ext_params
! *************************************************************************

    if (ixc < 0) then
       funcs(1)%id = -ixc/1000
       funcs(2)%id = -ixc - funcs(1)%id*1000
    else
       ! DIFFERENCE to ABINIT: Only use LibXC, ignore the sign!
       funcs(1)%id = ixc/1000
       funcs(2)%id = ixc - funcs(1)%id*1000
    end if

    do i = 1, 2
      if (funcs(i)%id == 0) then
        funcs(i)%family = 0
        cycle
      end if

      ! Get XC functional family
      funcs(i)%family = xc_f90_family_from_id(funcs(i)%id)
      select case (funcs(i)%family)
      case (XC_FAMILY_LDA, XC_FAMILY_GGA,XC_FAMILY_MGGA, XC_FAMILY_HYB_GGA)
         call xc_f90_func_init(funcs(i)%conf,funcs(i)%info,funcs(i)%id,nspden)
         if (funcs(i)%family==XC_FAMILY_HYB_GGA) then
            if (funcs(i)%id == 406 ) then  !! PBE0
               omega=0.d0
               call xc_f90_hyb_exx_coef(funcs(i)%conf,hfmix)
            else if (funcs(i)%id == 402 ) then !! B3LYP
               omega=0.d0
               call xc_f90_hyb_exx_coef(funcs(i)%conf,hfmix)
               write(info_str,'(2X,A)')"B3LYP HYB FUNC"
               call localorb_info( info_str )
            else if (funcs(i)%id == 428 ) then
               if ((hfmix.lt.0.d0).and.(omega.lt.0.d0)) then
                  call xc_f90_hyb_cam_coef(funcs(i)%conf, omega, alpha, beta)
                  hfmix=beta
               else
                  write(info_str,'(2X,A)')"This is a stub in external/atom_sphere/xcfunctions.f90 and you shouldn't see me"
                  call localorb_info( info_str )
                  write(info_str,'(2X,A)')"Please contact the developers to explain your situation and we can investigate further"
                  call localorb_info( info_str )
                  call aims_stop("ERROR: Reached unexpected stub in xcfunction.f90")
               end if
               write(info_str,'(2X,A)')"HSE06 HYB FUNC" !!,funcs(i)%id,omega,alpha,beta
               call localorb_info( info_str )
            end if
         else
            omega=0.d0
            hfmix=0.d0
         end if

      case default

        write(info_str,'(2X,A)') ' libxc_functionals_init : ERROR -'
        call localorb_info( info_str )
        write(info_str,'(2X,A,I9)') '  Invalid IXC = ',ixc
        call localorb_info( info_str )
        write(info_str,'(A,I9,A)')'  The LibXC functional family ',funcs(i)%family,&
             &    '  is currently not supported by pseudo.'
        call localorb_info( info_str )
        write(info_str,'(2X,A)')'  (-1 means the family is unknown to the LibXC itself)'
        call localorb_info( info_str )
        write(info_str,'(2X,A)')'  Please consult the LibXC documentation'
        call localorb_info( info_str )
      end select

!!!!      ! Dump functional information
!!!!      call xc_f90_info_name(funcs(i)%info,message)
!!!!      call wrtout(std_out,message,'COLL')
!!!!      ii = 0
!!!!      call xc_f90_info_refs(funcs(i)%info,ii,str,message)
!!!!      do while (ii >= 0)
!!!!        call wrtout(std_out,message,'COLL')
!!!!        call xc_f90_info_refs(funcs(i)%info,ii,str,message)
!!!!      end do
    end do 
          
  end subroutine libxc_functionals_init
!!***

!!****f* libxc_functionals/libxc_functionals_end
!! NAME
!!  libxc_functionals_end
!!
!! FUNCTION
!!  End usage of LibXC functional. Call LibXC end function,
!!  and deallocate module contents.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      xc_f90_gga_exc_vxc,xc_f90_gga_vxc,xc_f90_lda_exc_vxc,xc_f90_lda_vxc
!!      xc_f90_mgga_exc_vxc,xc_f90_mgga_vxc
!!
!! SOURCE
  subroutine libxc_functionals_end()


    implicit none

    integer :: i

    do i = 1, 2
      if (funcs(i)%id == 0) cycle
      call xc_f90_func_end(funcs(i)%conf)
    end do
  end subroutine libxc_functionals_end
!!*** 

!!****f* libxc_functionals/libxc_functionals_isgga
!! NAME
!!  libxc_functionals_isgga
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is a GGA or not
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  function libxc_functionals_isgga()


    implicit none

!Arguments ------------------------------------

!Local variables-------------------------------

    logical :: libxc_functionals_isgga

! *************************************************************************

    if (any(funcs%family == XC_FAMILY_GGA) .or. any(funcs%family == XC_FAMILY_HYB_GGA)) then
      libxc_functionals_isgga = .true.
    else
      libxc_functionals_isgga = .false.
    end if
  end function libxc_functionals_isgga
!!*** 

!!****f* libxc_functionals/libxc_functionals_exctXfac
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 real(kind=8) function libxc_functionals_exctXfac()


    implicit none

!Arguments ------------------------------------

!Local variables-------------------------------
    character*300 :: info_str


! *************************************************************************

    if (any(funcs%family == XC_FAMILY_HYB_GGA)) then
       !factors for the exact exchange contribution of different hybrid functionals
       if (any(funcs%id == XC_HYB_GGA_XC_PBEH)) then
          libxc_functionals_exctXfac = 0.25d0 
       end if
    else
      libxc_functionals_exctXfac = 0.d0
    end if
    write(info_str,'(2X,A,F15.8)')"fact",libxc_functionals_exctXfac
    call localorb_info( info_str )
  end function libxc_functionals_exctXfac
!!*** 

!!****f* libxc_functionals/libxc_functionals_ismgga
!! NAME
!!  libxc_functionals_ismgga
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is a Meta-GGA or not
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  function libxc_functionals_ismgga()


    implicit none

!Arguments ------------------------------------

!Local variables-------------------------------

    logical :: libxc_functionals_ismgga

! *************************************************************************

    if (any(funcs%family == XC_FAMILY_MGGA)) then
      libxc_functionals_ismgga = .true.
    else
      libxc_functionals_ismgga = .false.
    end if

  end function libxc_functionals_ismgga
!!*** 

!!****f* libxc_functionals/libxc_functionals_getvxc
!! NAME
!!  libxc_functionals_getvxc
!!
!! FUNCTION
!!  Return XC potential and energy, from input density (event gradient etc...)
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!      xc_f90_gga_exc_vxc,xc_f90_gga_vxc,xc_f90_lda_exc_vxc,xc_f90_lda_vxc
!!      xc_f90_mgga_exc_vxc,xc_f90_mgga_vxc
!!
!! SOURCE

!end module 
!!***


SUBROUTINE XCFUNCTION(nspol,rho,grad,EXC,VXC,dEdg)
! this version calls BigDFTs XC drivers, which access LibXC as part of the ABINIT XC functions.
! I should really try to put this apart from BigDFT and ABINIT, and directly call libXC. 

!     use libxcModule

      implicit none
      integer :: nspol,i
      real(8) :: alpha,exchange,correlation
      real(8) :: EXC,rho(nspol),VXC(nspol),dEdg(nspol),grad(nspol)  ! dummy ARGUMENTS
      real(8) :: EXCi,VXCi(nspol), sigma(3),vsigma(3)!  ! summands and libxc arg

!     These output quantities may be summed over two functionals (X and C)
      exchange=0.d0
      correlation=0.d0
      EXC  =0.0d0
      VXC  =0.0d0
      dEdg =0.0d0


       
!     convert the gradient to sigma if needed 
      if (libxc_functionals_isgga()) then
         sigma(1)=grad(1)*grad(1)
         if(nspol==2)then
            sigma(2)=grad(1)*grad(2)
            sigma(3)=grad(2)*grad(2)
         end if
      end if

!     libXC can use rather independent parts for exchange and correlation
!     the outer loop goes over both functionals from the two 3 digit codes
!     Santanu
     !do i=1,2
     !   if (funcs(i)%id == 0) cycle
     !   select case (funcs(i)%family)
     !   case (XC_FAMILY_HYB_GGA)
     !     write(6,*)"Hybrid Accepted"
     !     call xc_f90_hyb_exx_coef(funcs(i)%conf,alpha)
     !     write(6,*)"Alpha",alpha
     !     !call xc_f90_hyb_gga_xc_pbeh_init(funcs(i)%conf)
     !     !!call xc_f90_hyb_gga_pbeh_init(funcs(i)%conf)
     !   end select
     !end do
      do i = 1,2
          !write(6,*)"xcfunction func_id",i,funcs(i)%id
          if (funcs(i)%id == 0) cycle
          select case (funcs(i)%family)
!-------------------------------------------------------------------------------
!         LDA XC 
          case (XC_FAMILY_LDA)
             !write(6,*)"LDA FAMILY SELECTED"
             call xc_f90_lda_exc_vxc(funcs(i)%conf,1,rho(1),EXCi,VXCi(1))
             EXC=EXC+EXCi
             VXC=VXC+VXCi

!-------------------------------------------------------------------------------
!         GGA XC
          case (XC_FAMILY_GGA)
             !write(6,*)"GGA FAMILY SELECTED",i,funcs(i)%id
             call xc_f90_gga_exc_vxc(funcs(i)%conf,1,rho(1),sigma(1),&
                                     EXCi,VXCi(1),vsigma(1))
             EXC=EXC+EXCi
             VXC=VXC+VXCi

!            vsigma  are derivatives with respect to some products of gradients,
!            we want the derivatives with respect to the up and down gradients.
             if(nspol==1)then
               dEdg(1) = dEdg(1) +  vsigma(1)*grad(1)*2d0

             elseif(nspol==2)then
!              dE/dgup  =          dE/d(gup**2) *2*gup  + dE/d(gup*gdn) * gdn 
               dEdg(1)=dEdg(1) + vsigma(1) *2d0*grad(1) + vsigma(2)*grad(2)
!              dE/dgdn  =          dE/d(gdn**2) *2*gd   + dE/d(gup*gdn) * gup 
               dEdg(2)=dEdg(2) + vsigma(3) *2d0*grad(2) + vsigma(2)*grad(1)
             end if
!-------------------------------------------------------------------------------
!         HYB_GGA XC
          case (XC_FAMILY_HYB_GGA)
             !write(6,*)"HYB_GGA FAMILY SELECTED",i,funcs(i)%id
             call xc_f90_gga_exc_vxc(funcs(i)%conf,1,rho(1),sigma(1),&
                                     EXCi,VXCi(1),vsigma(1))
               EXC=EXC+EXCi
               VXC=VXC+VXCi

!            vsigma  are derivatives with respect to some products of gradients,
!            we want the derivatives with respect to the up and down gradients.
             if(nspol==1)then
               dEdg(1) = dEdg(1) +  vsigma(1)*grad(1)*2d0

             elseif(nspol==2)then
!              dE/dgup  =          dE/d(gup**2) *2*gup  + dE/d(gup*gdn) * gdn 
               dEdg(1)=dEdg(1) + vsigma(1) *2d0*grad(1) + vsigma(2)*grad(2)
!              dE/dgdn  =          dE/d(gdn**2) *2*gd   + dE/d(gup*gdn) * gup 
               dEdg(2)=dEdg(2) + vsigma(3) *2d0*grad(2) + vsigma(2)*grad(1)
             end if
          end select
      end do



!     this part is to plot the Exc(rho) function for debugging purposes
!      do j=1,nspol
!       write(17,'(i3,4f20.12)')j,rho(j),grad(j),EXC,VXC(j)
!      end do

END SUBROUTINE XCFUNCTION

end module libxcModule
