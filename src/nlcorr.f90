subroutine nlcorr(coord_current,rho, rho_gradient, eps1,eps2,eps3)

!!$  Calculates nonlocal correlation terms[1]: 
!!$  eps1 :     nonlocal correlation energy density e(x)
!!$  eps2 :     partial derivative d(rho*e)/d(rho)
!!$  eps3 :     partial derivative d(rho*e)/d[|grad(rho)|^2]
!!$  rho and rho_gradient solved from multipoles or octree interpolation. 
!!$  [1] A.Gulans,M.J.Puska and R.M.Nieminen, Phys. Rev. B 79 201105(R) (2009).
!!$  30.11.2010
!!$  Simiam Ghan (SAG), Department of Applied Physics, COMP, Aalto University, Espoo, Finland.
!!$  sag@fyslab.hut.fi

  use octree_routines
  use nlcorr_routines
  use vdwkernel
  use dimensions
  use grids
  use constants
  use pbc_lists, only: centers_basis_integrals 
  
  implicit none
  
 
  real*8  rho
  real*8, dimension(3) :: rho_gradient
  real*8 ::  coord_current(3)
  integer rho_flag, grad_flag
!   real*8 xx,rm,w,E, E2, E3
  real*8 rm,E,E2,E3
  
  !Integration params:  
  real*8 nradfactor
  integer GridSize , i_lebedev
  integer nradmax, nradmin, nrad, Nquad
  integer i,j, i_coord 
  parameter (nradfactor=4d0,nradmin=10,nradmax=40)
  
  real*8 DensityThreshold
  parameter (DensityThreshold=1d-5)
  real*8, dimension(3) :: rho_gradient_mp
  real*8, dimension(3) :: coord_current1,rho_gradient_mp2,xr
  real*8, dimension(4) :: temp1, temp2
  real*8, dimension(3) :: cc
  real*8  xrabs, rho_mp2
  real*8 eps1, eps2, eps3
  real*8 q01,q02,D,delta
  real*8 sum_ang,x,y,z,xp, sum2, sum3
  real*8 phi,DPhiD1,DPhiD2
  real*8 rho2, aux_squared_grad_rho, rho_grad_abs1, rho_grad_abs2, np1
  real*8  dqdn, dqdg
  logical  within_box

  real*8 :: dist_tab_sq(n_centers_basis_integrals)
  real*8 :: dir_tab(3,n_centers_basis_integrals)
  real*8 :: dist_tab(n_centers_basis_integrals)
  integer :: i_center, current_center
  real*8 :: coord_unmapped(3)

  real*8, dimension(:,:), allocatable, save :: ang, deriv
  real*8, dimension(:), allocatable, save :: w_ang, xx, r, w, temp
  real*8, external :: ddot

  !start work
  
  
  rho_flag = treeflag   !1 for mp, 2 for linear interpolation, 3 for hermite interpolation
  grad_flag = treeflag  !1 for mp, 2 for linear interpolation, 3 for hermite interpolation
  
  !linear interpolation available for testing (faster than hermite):
  !rho_flag = 2
  !grad_flag = 2
  
  
  
  !set integration params
  
  if(flag_nrad_nlcorr)then
     Nquad = nrad_nlcorr           
  else
     nrad=nradmax+dint(nradfactor*log(rho))  !adaptive radial grid
     if (nrad.le.nradmin) then
        nrad=nradmin
     elseif (nrad.gt.nradmax) then
        nrad=nradmax
     endif
     Nquad = nrad
     !Nquad = 10   !testing ...
  endif
  
  
  if(flag_i_leb_nlcorr)then
     i_lebedev = i_leb_nlcorr    !optional input
  else
     i_lebedev = 7  !7 for productions. No adaptive radial grids implemented yet. 
     !i_lebedev = 1 !testing
  endif
  
  
  if((i_lebedev.eq.6).or.(i_lebedev.eq.7))then  !patch these grids. ( Not available ).
     i_lebedev = 7
  elseif((i_lebedev.eq.12).or.(i_lebedev.eq.13).or.(i_lebedev.eq.14)) then  
     i_lebedev=14
  endif
 
  GridSize = n_ang_lebedev(i_lebedev)
  !write(21,*)"flag_i_leb_nlcorr, flag_nrad_nlcorr",  flag_i_leb_nlcorr, flag_nrad_nlcorr
  !write(21,*)"i_lebedev, gridsize", i_lebedev, GridSize
  !write(21,*)"Nquad", Nquad
  !write(21,*)" "
  if (.not.allocated(ang)) allocate(ang(GridSize,Nquad-1))
  if (.not.allocated(deriv)) allocate(deriv(GridSize,Nquad-1))
  if (.not.allocated(w_ang)) allocate(w_ang(GridSize))
  if (.not.allocated(xx)) allocate(xx(Nquad-1))
  if (.not.allocated(r)) allocate(r(Nquad-1))
  if (.not.allocated(w)) allocate(w(Nquad-1))
  if (.not.allocated(temp)) allocate(temp(Nquad-1))

  w_ang(:) = w_ang_lebedev(1:GridSize,i_lebedev)

  rho_grad_abs1 = sqrt(rho_gradient(1)**2.d0+rho_gradient(2)**2.d0+rho_gradient(3)**2.d0   )
  xr(:) = coord_current(:)
  xrabs = sqrt(xr(1)**2+xr(2)**2+xr(3)**2)

 
  if(  (xr(1).le.coord_master_full(1,2)).and.(xr(1).ge.coord_master_full(1,1)).and.&
       (xr(2).le.coord_master_full(2,2)).and.(xr(2).ge.coord_master_full(2,1)).and.&
       (xr(3).le.coord_master_full(3,2)).and.(xr(3).ge.coord_master_full(3,1))  ) then
     within_box = .true.
  else
     within_box = .false.
  endif
  
    
  if (rho.gt.DensityThreshold.and.within_box) then 
     
     dqdn = 0d0
     dqdg = 0d0
     call Calc_q0_full(rho,rho_gradient(:),q01,dqdn,dqdg)

     
     E=0d0
     E2 = 0d0
     eps1 = 0d0
     eps2 = 0d0
     eps3 = 0d0
     
     rm=1d0/q01


     do i=1,Nquad-1
        xx(i)=cos(i*pi/Nquad)
        r(i)=rm*sqrt((1+xx(i))/(1-xx(i)))
        w(i)=pi/dble(Nquad)*(sin(pi*dble(i)/Nquad))**2
        sum_ang=0d0
        sum2 = 0d0

        do j=1,GridSize
           x=xr(1)+r(i)*r_ang_lebedev(1,j,i_lebedev)
           y=xr(2)+r(i)*r_ang_lebedev(2,j,i_lebedev)
           z=xr(3)+r(i)*r_ang_lebedev(3,j,i_lebedev)
           ! xp=sqrt(x**2+y**2+z**2)
                      
           coord_current1(:) = 0d0
           coord_current1(1) = x
           coord_current1(2) = y
           coord_current1(3) = z

           coord_unmapped = coord_current1

           if(n_periodic > 0)then
              call map_to_center_cell(coord_current1(1:3) )
           end if

           cc(:) = coord_current1(:)
           
           if(  (cc(1).le.coord_master_full(1,2)).and.(cc(1).ge.coord_master_full(1,1)).and.&
                (cc(2).le.coord_master_full(2,2)).and.(cc(2).ge.coord_master_full(2,1)).and.&
                (cc(3).le.coord_master_full(3,2)).and.(cc(3).ge.coord_master_full(3,1))  ) then
              within_box = .true.
           else
              within_box = .false.
           endif
           
           
           if(within_box)then    
              
              if( (rho_flag.eq.1).or.(grad_flag.eq.1) )then
                 call get_rho_mp(coord_unmapped,rho_mp2,rho_gradient_mp2)   !from multipoles. 
              endif
              
              if( (rho_flag.eq.2).or.(grad_flag.eq.2) )then
                 temp1 = tree_interp_v2( coord_current1,root_full, coord_master_full)   !linear interp
              endif
              
              if( (rho_flag.eq.3).or.(grad_flag.eq.3) )then
                 temp2 = tree_interp_v3( coord_current1,root_full, coord_master_full)   !hermite interpolation
              endif

              if( (rho_flag.eq.4).or.(grad_flag.eq.4) )then
                 call tab_atom_centered_coords_p0 &
                      ( coord_unmapped, &
                      dist_tab_sq, &
                      dir_tab, &
                      n_centers_basis_integrals, centers_basis_integrals )
                 dist_tab = sqrt(dist_tab_sq)
                 current_center = -1
                 do i_center = 1, n_centers_basis_integrals
                    if (dist_tab(i_center) < min_atom_atom_tab(i_center)) then
                       current_center = i_center
                       exit
                    end if
                 end do
                 if (current_center > 0) then
                    !print *, 'Myid: ', myid, ' Current center: ', current_center
                    call get_rho_mp_single_center(coord_unmapped,temp2(1),temp2(2:4),current_center)   !from multipoles. 
                    !call get_rho_mp(coord_unmapped,temp2(1),temp2(2:4))   !from multipoles. 
                 else
                    temp2 = tree_interp_v3( coord_current1,root_full, coord_master_full)   !hermite interpolation
                 end if
              endif
              

              if(rho_flag.eq.1)then
                 rho2 = rho_mp2
              elseif(rho_flag.eq.2)then
                 rho2 = temp1(1)
              else
                 rho2 = temp2(1)
              endif

              if(grad_flag.eq.1)then
                 np1 = sqrt( rho_gradient_mp2(1)**2 + rho_gradient_mp2(2)**2 + rho_gradient_mp2(3)**2)
              elseif(grad_flag.eq.2)then
                 np1 = sqrt((temp1(2))**2 +(temp1(3))**2 +(temp1(4))**2  )   
              else
                 np1 = sqrt((temp2(2))**2 +(temp2(3))**2 +(temp2(4))**2  )   
              endif

              
              if (rho2.gt.DensityThreshold) then  
                 
                 
                 q02 = qq00(rho2,np1)
                 
                 
                 D=0.5d0*(q01+q02)*r(i)
                 delta=(q01-q02)/(q01+q02)
                 
                 call GetKernelValue(D,delta,phi,DPhiD1,DPhiD2)
                 
                 ang(j,i)=phi*rho2
                 
                 ! sum_ang=sum_ang+phi*rho2*w_ang_lebedev(j,i_lebedev)
                 deriv(j,i)=(DPhiD1+r(i)*DPhiD2*q02/D**2)*rho2
                 ! sum2=sum2+(DPhiD1+r*DPhiD2*q02/D**2)*rho2*w_ang_lebedev(j,i_lebedev) 

              else ! not large enough rho2
                 ang(j,i) = 0.0d0
                 deriv(j,i) = 0.0d0
              endif
              
           else ! not within the box
              ang(j,i) = 0.0d0
              deriv(j,i) = 0.0d0
           endif
        enddo
        
!!$        sum_ang = ddot(GridSize,ang,1,w_ang,1)
!!$        E = E + sum_ang*w*(1-xx)**(-3)*rm**3
!!$        !potential:
!!$        sum2 = ddot(GridSize,deriv,1,w_ang,1)
!!$        E2 = E2 + r*sum2*w*(1-xx)**(-3)*rm**3 

        
     enddo

     call dgemv('T',GridSize,Nquad-1,1.0d0,ang,GridSize,w_ang,1,0.0d0,temp,1)
     temp = temp*w*(1-xx)**(-3)
     temp = temp*rm**3
     E = SUM(temp)

     call dgemv('T',GridSize,Nquad-1,1.0d0,deriv,GridSize,w_ang,1,0.0d0,temp,1)
     temp = r*temp*w*(1-xx)**(-3)
     temp = temp*rm**3
     E2 = SUM(temp)

     eps1 = E*(4d0/pi)    !TRUE ENERGY DENSITY (atomic units)
     
     !Potential terms:
     eps2 = ( 2d0*E + dqdn*E2) * (4d0/pi) !density_deriv.  
     eps3 = ( dqdg * E2 ) * (4d0/pi)                  !grad_deriv.
    
  else   
     eps1=0d0
     eps2 = 0d0
     eps3 = 0d0
  endif
  
  deallocate(ang)
  deallocate(deriv)
  deallocate(w_ang)
  deallocate(xx)
  deallocate(r)
  deallocate(w)
  deallocate(temp)
  
end subroutine nlcorr



      
