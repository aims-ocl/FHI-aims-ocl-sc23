module octree_routines
  
!!$ Contains routines needed to generate octrees for all-electron density interpolation. 
!!$ 30.11.2010
!!$ Simiam Ghan (SAG), Department of Applied Physics, COMP, Aalto University, Espoo, Finland.
!!$ sag@fyslab.hut.fi
  
  use nlcorr_routines  
  use localorb_io
  use dimensions
  implicit none
  
  
  type treenode
     type(treenode), pointer,dimension(:) :: branch
     real*8, dimension(4,8) :: val   
  end type treenode
   
  type(treenode), pointer :: root_full
  real*8, dimension(3,2) :: coord_master_full  
  real*8 r_lim
  logical cube_ready
  real*8 epsfinal
  integer treeflag
  
  integer nrad_nlcorr
  integer i_leb_nlcorr
  logical flag_nrad_nlcorr
  logical flag_i_leb_nlcorr
  
  real*8, dimension(3,8) :: unit_cell_vertices
  real*8, dimension(:, :), allocatable :: atom_atom_tab
  real*8, dimension(:),    allocatable :: min_atom_atom_tab
  
contains 
  
  
  
  subroutine gen_master_cube()  
    
!!$ Define master cubic area where vdw interactions (nlcorr.f90) will be calculated.
!!$ Build an octree if required.
    
    
    use geometry
    use dimensions
    use runtime_choices
    use pbc_lists
    use mpi_tasks, only: check_allocation, myid
    implicit none
    
    integer count,count2,i_point , count_ended,count_started, acc_not_achieved
    integer acc_achieved,divide_overly_large, to_be_divided 
    integer exact_evals, deall, reclevel, deepest
    integer interpol_evals
    real*8 R,n,np 
    integer ii, count_f11
    real*8, dimension(3,8):: cc, dummy
    real*8, dimension(3) :: x_m,x, dx,dy,dz, origin, cmax, cmin
    real*8, dimension(2) :: temp
    real*8  x1,x2,dxx,dd,abso, ddx, ddy,ddz
    real*8 diff, xd, yd, zd, maxdx
    logical freeonly, hermite, debug
    integer i,j,k, num, i_coord, i_atom, i_atom_2 
    real*8 x_val_exact1,x_val_exact2, x_val_int, offset,offset_origins
    real*8 diff1, diff2, diff3, cubicness, epsabss, epss

    character*100 :: info_str

!begin work

    ! First, let's tabulate the interatomic distances. When we are close enough to one of them
    ! the density will simply follow from that one due to the partition tab and no octree is needed.
    if (.not.allocated(atom_atom_tab)) then
       if ( (partition_type.eq.5) .or. (partition_type.eq.7) .or. (partition_type.eq.8) .or. (partition_type.eq.9) ) then
          ! for stratmann partitioning scheme, make available only if necessary
          ! tabulate the interatomic distances for all relevant atoms
          allocate(atom_atom_tab(n_centers_basis_integrals,n_centers_basis_integrals),stat=i_point)
          call check_allocation(i_point, 'atom_atom_tab                 ')
       else
          ! Not needed, allocate dummy
          allocate(atom_atom_tab(1,1))
       endif
    end if
       
    if (.not.allocated(min_atom_atom_tab)) then
       allocate(min_atom_atom_tab(n_centers_basis_integrals),stat=i_point)
       call check_allocation(i_point, 'min_atom_atom_tab             ')
    end if

    if ( (partition_type.eq.5) .or. (partition_type.eq.7) .or. (partition_type.eq.8) .or. (partition_type.eq.9) ) then
       ! if (myid==0) print *, 'HERE'
       call tab_interatomic_distances(n_centers_basis_integrals, centers_basis_integrals, atom_atom_tab)
       min_atom_atom_tab = 1.0d100
       do i_atom = 1, n_centers_basis_integrals
          do i_atom_2 = 1, n_centers_basis_integrals
             if (i_atom /= i_atom_2) then
                min_atom_atom_tab(i_atom) = MIN(min_atom_atom_tab(i_atom),atom_atom_tab(i_atom,i_atom_2))
             end if
          end do
       end do
!!$       if (myid==0) print *, min_atom_atom_tab
!!$       if (myid==0) print *, atom_atom_tab
       ! Scaled for direct comparison in the future.
       min_atom_atom_tab(:) = (1.0d0-stratmann_a)*min_atom_atom_tab(:)/2.0d0
    end if
!!$    if (myid == 0) print *, min_atom_atom_tab
!!$    if (myid == 0) print *, atom_atom_tab

    
!Create master cube based on system coords.
    
    debug = .false.

    if (n_periodic > 0) then
       unit_cell_vertices(:,1) = &
            0.5*(lattice_vector(:,1) + lattice_vector(:,2) + lattice_vector(:,3))
       unit_cell_vertices(:,2) = &
            0.5*(-lattice_vector(:,1) + lattice_vector(:,2) + lattice_vector(:,3))
       unit_cell_vertices(:,3) = &
            0.5*(lattice_vector(:,1) - lattice_vector(:,2) + lattice_vector(:,3))
       unit_cell_vertices(:,4) = &
            0.5*(lattice_vector(:,1) + lattice_vector(:,2) - lattice_vector(:,3))
       unit_cell_vertices(:,5) = &
            0.5*(- lattice_vector(:,1) - lattice_vector(:,2) + lattice_vector(:,3))
       unit_cell_vertices(:,6) = &
            0.5*(- lattice_vector(:,1) + lattice_vector(:,2) - lattice_vector(:,3))
       unit_cell_vertices(:,7) = &
            0.5*(lattice_vector(:,1) - lattice_vector(:,2) - lattice_vector(:,3))
       unit_cell_vertices(:,8) = &
            0.5*(-lattice_vector(:,1) - lattice_vector(:,2) - lattice_vector(:,3))
       if ((debug).and.(myid==0)) write(use_unit,*) unit_cell_vertices
       ! We need to have a master cube that contains the entire unit cell.
       ! Still, the fact that master cube is a cube is problematic for very non-cubic
       ! unit cells
 
       coord_master_full(:,1) = MINVAL(unit_cell_vertices(:,:))
       coord_master_full(:,2) = MAXVAL(unit_cell_vertices(:,:))
    else ! non-periodic case
       cmax(:) = coords(:,1)
       cmin(:) = coords(:,1)
    
       do i_atom = 1, n_atoms
          do i_coord = 1,3
             if(coords(i_coord,i_atom).gt.cmax(i_coord)) cmax(i_coord) = coords(i_coord,i_atom)
             if(coords(i_coord,i_atom).lt.cmin(i_coord)) cmin(i_coord) = coords(i_coord,i_atom)
          enddo
       enddo
    
       R = 6d0   !minimum distance between atom and edge of vdw region.  
       coord_master_full(:,2) = cmax(:) + R
       coord_master_full(:,1) = cmin(:) - R

       !that's a rectangle...
    
       xd = abs( coord_master_full(1,2) - coord_master_full(1,1) ) 
       yd = abs( coord_master_full(2,2) - coord_master_full(2,1) )
       zd = abs( coord_master_full(3,2) - coord_master_full(3,1) )
       
    
       maxdx = xd
       if((xd.gt.yd).and.(xd.gt.zd)) maxdx = xd
       if((yd.gt.xd).and.(yd.gt.zd)) maxdx = yd
       if((zd.gt.yd).and.(zd.gt.xd)) maxdx = zd
    
       if(xd.lt.maxdx)then
          diff = maxdx - xd
          coord_master_full(1,2) =  coord_master_full(1,2) + 0.5d0*diff
          coord_master_full(1,1) =  coord_master_full(1,1) - 0.5d0*diff
       endif
       
       if(yd.lt.maxdx)then
          diff = maxdx - yd
          coord_master_full(2,2) =  coord_master_full(2,2) + 0.5d0*diff
          coord_master_full(2,1) =  coord_master_full(2,1) - 0.5d0*diff
       endif
       
       if(zd.lt.maxdx)then
          diff = maxdx - zd
          coord_master_full(3,2) =  coord_master_full(3,2) + 0.5d0*diff
          coord_master_full(3,1) =  coord_master_full(3,1) - 0.5d0*diff
       endif
  
    end if ! periodic vs. non-periodic case
    
    xd = abs( coord_master_full(1,2) - coord_master_full(1,1) ) 
    yd = abs( coord_master_full(2,2) - coord_master_full(2,1) )
    zd = abs( coord_master_full(3,2) - coord_master_full(3,1) )
    
    diff1 = abs(xd-yd)
    diff2 = abs(xd-zd)
    diff3 = abs(yd-zd)
    
    cubicness = 0.0001d0
    
    if(diff1.gt.cubicness.or.diff2.gt.cubicness.or.diff3.gt.cubicness)then
       write(use_unit,*)"Warning: master octree cube not cubic! Consequences for interpolation accuracy!"
    else
       !write(use_unit,*)"Master cube is indeed cubic."
    endif
       
    count_f11 = 0
   
   
    ddx =abs( coord_master_full(1,2) -  coord_master_full(1,1) )
    ddy =abs( coord_master_full(2,2) -  coord_master_full(2,1) )
    ddz =abs( coord_master_full(3,2) -  coord_master_full(3,1) )
    
    dx = (/ddx,0d0,0d0/)
    dy = (/0d0,ddy,0d0/)
    dz = (/0d0,0d0,ddz/)
    
    origin(:) = coord_master_full(:,1)
    
    num=0
    do k = 1,2
       do j = 1,2
          do i = 1,2
             num = num+1
             cc(:,num) = origin(:) + (i-1)*dx(:) + (j-1)*dy(:) + (k-1)*dz(:)
          enddo
       enddo
    enddo
        
    do i = 1,8
       root_full%val(:,i) = fmp (cc(:,i),count_f11)  
    enddo
    
    if(debug.and.myid==0)then
       write(use_unit,*)"coord_master_full(:,:)"
       write(use_unit,*)coord_master_full(:,:)
    endif
    
    cube_ready = .true.   !nlcorr can now be used since the cube is initialized.


! Now create octree if required. 
        
    count_ended=0
    count_started= 0 
    acc_not_achieved = 0
    acc_achieved = 0 
    divide_overly_large  = 0
    to_be_divided = 0
    exact_evals = 0
    interpol_evals = 0
    reclevel = 0
    deepest = 0
    count = 0
    count2 = 0
    dummy(:,:) = 0d0
    

  
    if(treeflag.eq.1)then
       write(use_unit,*)" Nonlocal correlation using multipoles."
    endif
      
   
    if((treeflag==3).or.(treeflag==4))then    
            
       write(info_str,'(2X,A)') "Building charge-based octree for nonlocal correlation:"
       call localorb_info(info_str, use_unit,'(A)', OL_norm)
       write(info_str,'(2X,A,F11.8)') "epsilon: ",epsfinal
       call localorb_info(info_str, use_unit,'(A)', OL_norm)
             
       call add_child(epsfinal,coord_master_full, root_full,&
            reclevel,deepest,count_f11 )    
    endif
    
    
   

    if(.false.)then  !An inferior interpolation-based tree. Not used - remains for development purposes.  
       write(use_unit,*)"Calling interpolation-based octree."
       hermite = .true.  
       epss=0.0005d0  
       epsabss = 0.0001d0
       call add_child_interpolation(hermite,epsabss,epss,coord_master_full, root_full,count_ended,&
            count_started, acc_not_achieved, &
            acc_achieved,divide_overly_large,to_be_divided,exact_evals,reclevel,deepest, &
            interpol_evals,count_f11 )  
    endif
    
    
    call count_tree_nodes(root_full,count,count2)
        
    if(((treeflag==3).or.(treeflag==4)).and.debug.and.myid==0)then
       write(use_unit, *) "number of calls to get_rho_mp:", count_f11
       write(use_unit, *) "tree nodes total:             ", count
       write(use_unit, *) "lowest level nodes:           ", count2
       write(use_unit, *) "deepest level:                ", deepest
       write(use_unit,*) " "
    endif


  end subroutine gen_master_cube











 function fmp(coord,count_f11) result(val)
    
    use dimensions
    use constants
    use physics
    
    implicit none
    real*8 R
    integer i_coord
    real*8, dimension(3):: coord, rho_gradient_mp
    real*8, dimension(4) :: val
    real*8 rho_mp, aux_squared_grad_rho, rho_grad_abs
    integer,optional ::count_f11
    logical freeonly
    integer i
   
    if(Present(count_f11))    count_f11 = count_f11 +1

    if(.true.)then    
              
       call get_rho_mp(coord,rho_mp,rho_gradient_mp)
       
       val(1) = rho_mp
       do i= 1,3
          val(i+1) = rho_gradient_mp(i)
       enddo
              
    else !testing gaussian
       
       val(1) = exp(-1d0*(coord(1)**2+coord(2)**2+coord(3)**2))
       val(2) = -2d0*coord(1)*exp(-1d0*(coord(1)**2+coord(2)**2+coord(3)**2))
       val(3) = -2d0*coord(2)*exp(-1d0*(coord(1)**2+coord(2)**2+coord(3)**2))
       val(4) = -2d0*coord(3)*exp(-1d0*(coord(1)**2+coord(2)**2+coord(3)**2))       
              
    endif
    
  end function fmp




  recursive subroutine count_tree_nodes(root,count,count2)
    implicit none
    type(treenode) :: root

    integer count,count2,i,i_child
    real*8, dimension(3,8,8) :: c_coord
    real*8, dimension(3,8) :: c_origin
    
    count = count + 1
    if(associated(root%branch)) then
 
       
       do i_child = 1,8
          call count_tree_nodes(root%branch(i_child),count,count2)
       enddo
    else
       count2 = count2+1
    endif
    
  end subroutine count_tree_nodes
  







  
  recursive subroutine deallocate_tree(root, count)
    implicit none
    
    type(treenode) :: root
    integer count,count2,i,i_child, no_associated
    
    
    if(associated(root%branch))then
       
       do i_child = 1,8
          call deallocate_tree(root%branch(i_child),count)
       enddo
       
       count = count + 1  !number of deallocations.
       deallocate(root%branch)  
       
    endif
    
  end subroutine deallocate_tree
  





  
  subroutine test_interpolation(root,coord_m,R )
    !only for testing multipoles, interpolation.
    
    implicit none
    real*8, dimension(3) :: x_m,x
    type(treenode) :: root
    
    integer i ,ii
    double precision n,np
    real*8 x_val_exact1,x_val_exact2, x_val_int, offset,offset_origins,abso,dd, np_temp3, np_temp1, npexact
    real*8 R,q_val_exact1,q_val_exact2, q_val_int,testq, q0exact, q0v2, q0v3, q0v3_e
    real*8, dimension(4) :: x_exact, x_int, temp1,temp2, x_int2, temp3
    real*8, dimension(3,2) :: coord_m
    integer count_f11  
    
    real*8  q0_exrho_v2grad 
    real*8  q0_exrho_v3grad 
    real*8  q0_v2rho_exgrad 
    real*8  q0_v3rho_exgrad , step
    
    
    count_f11 = 0
    
    step = 0.02
    
    do i= 1,1000
       
       x_m(:) = 0.d0
       x_m(1) = -6d0 
       
       x_m(1) = x_m(1)+ dble(i)*step
       temp1(:) = tree_interp_v2(x_m,root,coord_m)
       temp3(:) = tree_interp_v3(x_m,root,coord_m)
       temp2(:) = fmp(x_m, count_f11) 
       write(10,*) x_m(1),abs(temp2(1)-temp1(1)), abs(temp2(1)-temp3(1))
       
       
       
    enddo
    
    
  end subroutine test_interpolation
  
     
     
     
     
     
  recursive  function tree_interp_v2(x,root,coord_v2) result(vals)  
    !linear interpolation of density from octree.  
    
    
    !  use pointertypes_3D
    implicit none
    real*8, dimension(3,8,8) :: c_coord
    real*8, dimension(3,8):: coord,c_origin
    real*8, dimension(3,2) :: coord_v2         
    real*8, dimension(3,2,8) :: cc_
    
    real*8, dimension(3) :: x
    real*8, dimension(4) :: vals
    real*8, dimension(4) :: y12,y34,y56,y78,y1234,y5678,yfinal
    real*8 dx,dy,dz,u0,v0,w0,u1,v1,w1, invdx
    real*8, dimension(3) :: ddx, ddy, ddz
    type(treenode) :: root

    integer i, i_child, j,k, num
    logical writ
    writ = .false.
    
    if(associated(root%branch)) then
       
       if(writ)then
          write(use_unit,*)"branches are associated"
          write(use_unit,*)"coord(:,:) of current branch"
          write(use_unit,*) coord(:,:)
          write(use_unit,*)"coord_v2"
          write(use_unit,*)coord_v2(:,:)
       endif
      
     
       !new child vertice coordinates
       
       ddx = abs(coord_v2(1,2)-coord_v2(1,1))*0.5d0*(/1d0,0d0,0d0/)
       ddy = abs(coord_v2(2,2)-coord_v2(2,1))*0.5d0*(/0d0,1d0,0d0/)
       ddz = abs(coord_v2(3,2)-coord_v2(3,1))*0.5d0*(/0d0,0d0,1d0/)
       
       num = 0
       do k = 1,2
          do j = 1,2
             do i = 1,2
                num = num + 1
                cc_(:,1,num) = coord_v2(:,1) + (i-1)*ddx + (j-1)*ddy + (k-1)*ddz           
             enddo
          enddo
       enddo
       num = 0
       do k = 2,3
          do j = 2,3
             do i = 2,3
                num = num + 1
                cc_(:,2,num) = coord_v2(:,1) + (i-1)*ddx + (j-1)*ddy + (k-1)*ddz           
             enddo
          enddo
       enddo
      
     

       !find which child has x

       if( (x(1).ge.cc_(1,1,1)).and.(x(1).le.cc_(1,2,1)))then
          
          if( (x(2).ge.cc_(2,1,1)).and.(x(2).le.cc_(2,2,1)) )then
          
             if( (x(3).ge.cc_(3,1,1)).and.(x(3).le.cc_(3,2,1)))then

                !search 1
                if(writ) write(use_unit,*)"going into branch no: 1"
                vals=tree_interp_v2(x,root%branch(1), cc_(:,:,1))!
                
             else
                ! search 5
                if(writ) write(use_unit,*)"going into branch no: 5"
                vals=tree_interp_v2(x,root%branch(5), cc_(:,:,5))!)
             
             endif
          else
             !search children 3,7
             if( (x(3).ge.cc_(3,1,3)).and.(x(3).le.cc_(3,2,3)))then

                !search 3
                if(writ)write(use_unit,*)"going into branch no: 3"
                vals=tree_interp_v2(x,root%branch(3), cc_(:,:,3))!)

             else
                
                ! search 7
                if(writ)write(use_unit,*)"going into branch no:7"
                vals=tree_interp_v2(x,root%branch(7), cc_(:,:,7))!)

             endif
          endif
          
          
       else
          !search children 2,4,6,8
          if( (x(2).ge.cc_(2,1,2)).and.(x(2).le.cc_(2,2,2)) )then

             !search children 2,6
             if( (x(3).ge.cc_(3,1,2)).and.(x(3).le.cc_(3,2,2)))then

                !search 2
                if(writ)write(use_unit,*)"going into branch no: 2"
                vals=tree_interp_v2(x,root%branch(2), cc_(:,:,2))!)

             else
                ! search 6
                if(writ)write(use_unit,*)"going into branch no: 6"
                vals=tree_interp_v2(x,root%branch(6), cc_(:,:,6))!)

             endif
             
          else
             !search 4,8
             if( (x(3).ge.cc_(3,1,4)).and.(x(3).le.cc_(3,2,4)))then
 
                !search 4
                if(writ)write(use_unit,*)"going into branch no: 4"
                vals=tree_interp_v2(x,root%branch(4), cc_(:,:,4))!)

             else
                ! search 8
                if(writ)write(use_unit,*)"going into branch no: 8"
                vals=tree_interp_v2(x,root%branch(8), cc_(:,:,8))!)

             endif
          endif
          
       endif
       
    else
       
       if(writ)then
          write(use_unit,*)"branches not associated"
          write(use_unit,*)"coord(:,:)"
          write(use_unit,*) coord(:,:)
          write(use_unit,*)"root%yval(:)"
          write(use_unit,*)root%val(:,:) 
       endif
       
       
    
      dx = abs(coord_v2(1,2) - coord_v2(1,1))
       invdx = 1d0/dx   !assume dx=dy=dz (cube).


       u0 = (coord_v2(1,2)-x(1)) * invdx 
       u1 = 1d0 - u0
       v0 = (coord_v2(2,2)-x(2)) * invdx 
       v1 =  1d0 - v0
       w0 = (coord_v2(3,2)-x(3)) * invdx 
       w1 =  1d0 - w0

       vals(:) = root%val(:,1)*u0*v0*w0 + root%val(:,2)*u1*v0*w0 +&
            root%val(:,3)*u0*v1*w0 + root%val(:,4)*u1*v1*w0 + &
            root%val(:,5)*u0*v0*w1 + root%val(:,6)*u1*v0*w1 + &
            root%val(:,7)*u0*v1*w1 + root%val(:,8)*u1*v1*w1 
       
       if(writ)   write(use_unit,*)"tree_interp: yfinal is",yfinal   
       
    endif

    if(writ)   write(use_unit,*)"tree_interp recursion ended"
    
  end function tree_interp_v2





!!$
!!$

  
  recursive  function tree_interp_v3(x,root,coord_v2) result(vals)  
    ! 3D Cubic Hermite spline interpolation of density and gradient from octree. 
    ! Follows Numerical Recipies, 3rd Ed. Chapter 3.3 "Cubic Spline Interpolation."
    ! Contains some code by Andris Gulans. 
    ! SAG
    
    
    
    
    implicit none
    real*8, dimension(3,8,8) :: c_coord
    real*8, dimension(3,8):: coord,c_origin
    real*8, dimension(3,2) :: coord_v2        
    real*8, dimension(3,2,8) :: cc_
    real*8, dimension(3) :: x, x1
    real*8, dimension(4) :: vals
    real*8, dimension(4,4) :: xedge
    real*8, dimension(4,2) :: yedge
    real*8, dimension(4) :: y12,y34,y56,y78,y1234,y5678,yfinal
    real*8 dx,dy,dz,u0,v0,w0,v1,w1, invdx, m0,m1,p1,p0,t, xx, xk, xk1, xk0, f0
    real*8, dimension(3) :: ddx, ddy, ddz
    type(treenode) :: root
    real*8 density(0:3,0:1,0:1,0:1)
    real*8, dimension(3) :: nabla
    real*8 answer, invdx3
   real*8 ffx( 0:1,0:1,0:1), ffy( 0:1,0:1,0:1), ffz( 0:1,0:1,0:1)
   real*8 f(4,0:1,0:1,0:1)
    real*8 a,b,c,d
    real*8 a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,f1,f2,f3,t1,t2,t3,u1,u2,u3
    real*8 a1a2, a1a3, a2a3, b1b2,b1b3,b2b3, a1b2, a1b3, a2b3, b1a2, b1a3, b2a3
    real*8 sixth  
    integer klo1,klo2,klo3
    parameter(sixth=0.16666666666666666666666666666667d0)
    
    integer i, i_child, j,k, num, count
    logical writ
    
    writ = .false.
    
    if(associated(root%branch)) then
       if(writ)then
          write(use_unit,*)"branches are associated"
          write(use_unit,*)"coord_v2 of current cube"
          write(use_unit,*)coord_v2(:,:)
       endif
       
       
       
       !new child vertice coordinates
       
       ddx = abs(coord_v2(1,2)-coord_v2(1,1))*0.5d0*(/1d0,0d0,0d0/)
       ddy = abs(coord_v2(2,2)-coord_v2(2,1))*0.5d0*(/0d0,1d0,0d0/)
       ddz = abs(coord_v2(3,2)-coord_v2(3,1))*0.5d0*(/0d0,0d0,1d0/)
       
       num = 0
       do k = 1,2
          do j = 1,2
             do i = 1,2
                num = num + 1
                cc_(:,1,num) = coord_v2(:,1) + (i-1)*ddx + (j-1)*ddy + (k-1)*ddz           
                
                
             enddo
          enddo
       enddo
       num = 0
       do k = 2,3
          do j = 2,3
             do i = 2,3
                num = num + 1
                cc_(:,2,num) = coord_v2(:,1) + (i-1)*ddx + (j-1)*ddy + (k-1)*ddz           
                
              
             enddo
          enddo
       enddo
       


       !Search tree for x

     
        if( (x(1).ge.cc_(1,1,1)).and.(x(1).le.cc_(1,2,1)))then
       !search children 1,3,5,7
                 
            if( (x(2).ge.cc_(2,1,1)).and.(x(2).le.cc_(2,2,1)) )then
               
               !search children 1,5
             
               if( (x(3).ge.cc_(3,1,1)).and.(x(3).le.cc_(3,2,1)))then

               !search 1
                if(writ) write(use_unit,*)"going into branch no: 1"
                vals=tree_interp_v3(x,root%branch(1), cc_(:,:,1))
                
             else
                ! search 5
                if(writ) write(use_unit,*)"going into branch no: 5"
                vals=tree_interp_v3(x,root%branch(5), cc_(:,:,5))

             endif
          else
             !search children 3,7
             if( (x(3).ge.cc_(3,1,3)).and.(x(3).le.cc_(3,2,3)))then

                !search 3
                if(writ)write(use_unit,*)"going into branch no: 3"
                vals=tree_interp_v3(x,root%branch(3), cc_(:,:,3))

             else
                ! search 7
                if(writ)write(use_unit,*)"going into branch no:7"
                vals=tree_interp_v3(x,root%branch(7), cc_(:,:,7))

             endif
          endif
          
          
       else
          !search children 2,4,6,8
          if( (x(2).ge.cc_(2,1,2)).and.(x(2).le.cc_(2,2,2)) )then

             !search children 2,6
             
             if( (x(3).ge.cc_(3,1,2)).and.(x(3).le.cc_(3,2,2)))then

                !search 2
                if(writ)write(use_unit,*)"going into branch no: 2"
                vals=tree_interp_v3(x,root%branch(2), cc_(:,:,2))

             else
                ! search 6
                if(writ)write(use_unit,*)"going into branch no: 6"
                vals=tree_interp_v3(x,root%branch(6), cc_(:,:,6))

             endif
             
          else
             !search 4,8
             if( (x(3).ge.cc_(3,1,4)).and.(x(3).le.cc_(3,2,4)))then
 
                !search 4
                if(writ)write(use_unit,*)"going into branch no: 4"
                vals=tree_interp_v3(x,root%branch(4), cc_(:,:,4))

             else
                ! search 8
                if(writ)write(use_unit,*)"going into branch no: 8"
                vals=tree_interp_v3(x,root%branch(8), cc_(:,:,8))

             endif
             
          endif
          
          
       endif
       
    else
       
       if(writ)then
          write(use_unit,*)"branches not associated. performing hermite interpolation."
          write(use_unit,*)"coord_v2(:,:) of current cube"
          write(use_unit,*) coord_v2(:,:)
          write(use_unit,*)"root%yval(:)"
          write(use_unit,*)root%val(:,:) 
       endif
       
       
       if(.true.)then  
          
          
          !function vals, gradients
          count = 0
          do k = 0,1
             do j = 0,1
                do i = 0,1
                   count = count +1
                   f(1,i,j,k) = root%val(1,count)
                enddo
             enddo
          enddo
          
          !approximate double-derivatives of rho along axis direction using a 1D cubic polynomial.
          !rho, gradient at node vertices serve as boundary conditions. 
          
          !x
          
          !first x edge
          
          !assume cube          
          
          xk = coord_v2(1,1)
          xk1 = coord_v2(1,2)
          f0 = root%val(1,1)
          f1 = root%val(1,2)
          m0 = root%val(2,1)
          m1 = root%val(2,2)
          invdx3 = (1 / (xk - xk1)**3)   !assume cube: use this for all directions. 
          A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3
          B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
          ffx(0,0,0) = 6*A*xk + 2*B
          ffx(1,0,0) = 6*A*xk1 + 2*B
          
          !second x edge
          
          xk = coord_v2(1,1)
          xk1 = coord_v2(1,2)
          f0 = root%val(1,3)
          f1 = root%val(1,4)
          m0 = root%val(2,3)
          m1 = root%val(2,4)
          A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3
          B =   invdx3* ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
          ffx(0,1,0) = 6*A*xk + 2*B
          ffx(1,1,0) = 6*A*xk1 + 2*B
          
          !third x edge
          
          xk = coord_v2(1,1)
          xk1 = coord_v2(1,2)
          f0 = root%val(1,5)
          f1 = root%val(1,6)
          m0 = root%val(2,5)
          m1 = root%val(2,6)
          A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3
          B =   invdx3* ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
          ffx(0,0,1) = 6*A*xk + 2*B
          ffx(1,0,1) = 6*A*xk1 + 2*B
          
          !fourth x edge
          
          xk = coord_v2(1,1)
          xk1 = coord_v2(1,2)
          f0 = root%val(1,7)
          f1 = root%val(1,8)
          m0 = root%val(2,7)
          m1 = root%val(2,8)
          A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3
          B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
          ffx(0,1,1) = 6*A*xk + 2*B
          ffx(1,1,1) = 6*A*xk1 + 2*B
          
          
          !y
          
          !first y edge
          
          xk = coord_v2(2,1)
          xk1 = coord_v2(2,2)
          f0 = root%val(1,1)
          f1 = root%val(1,3)
          m0 = root%val(3,1)
          m1 = root%val(3,3)
          A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3
          B =   invdx3* ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
          ffy(0,0,0) = 6*A*xk + 2*B
          ffy(0,1,0) = 6*A*xk1 + 2*B
          
          !second y edge
          
          
          f0 = root%val(1,2)
          f1 = root%val(1,4)
          m0 = root%val(3,2)
          m1 = root%val(3,4)
          A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3
          B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
          ffy(1,0,0) = 6*A*xk + 2*B
          ffy(1,1,0) = 6*A*xk1 + 2*B
          
          !third y edge
          
          
          f0 = root%val(1,5)
          f1 = root%val(1,7)
          m0 = root%val(3,5)
          m1 = root%val(3,7)
          A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3 !/(xk- xk1)**3
          B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
          ffy(0,0,1) = 6*A*xk + 2*B
          ffy(0,1,1) = 6*A*xk1 + 2*B
          
          !fourth y edge
          
          
          f0 = root%val(1,6)
          f1 = root%val(1,8)
          m0 = root%val(3,6)
          m1 = root%val(3,8)
          A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3 !/(xk- xk1)**3
          B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
          ffy(1,0,1) = 6*A*xk + 2*B
          ffy(1,1,1) = 6*A*xk1 + 2*B
          
          
          !write(use_unit,*)"ffy", ffy(0:1,0:1,0:1)
          !write(use_unit,*)" "
          !z
          
          !first z edge
          
          xk = coord_v2(3,1)
          xk1 = coord_v2(3,2)
          f0 = root%val(1,1)
          f1 = root%val(1,5)
          m0 = root%val(4,1)
          m1 = root%val(4,5)
          A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3 !/(xk- xk1)**3
          B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
          ffz(0,0,0) = 6*A*xk + 2*B
          ffz(0,0,1) = 6*A*xk1 + 2*B
          
          !second z edge
          
          f0 = root%val(1,2)
          f1 = root%val(1,6)
          m0 = root%val(4,2)
          m1 = root%val(4,6)
          A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3 !/(xk- xk1)**3
          B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
          ffz(1,0,0) = 6*A*xk  + 2*B
          ffz(1,0,1) = 6*A*xk1 + 2*B
          
          !third z edge
          
          f0 = root%val(1,3)
          f1 = root%val(1,7)
          m0 = root%val(4,3)
          m1 = root%val(4,7)
          A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3 !/(xk- xk1)**3
          B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
          ffz(0,1,0) = 6*A*xk + 2*B
          ffz(0,1,1) = 6*A*xk1 + 2*B
          
          !fourth z edge
          
          f0 = root%val(1,4)
          f1 = root%val(1,8)
          m0 = root%val(4,4)
          m1 = root%val(4,8)
          A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3 !/(xk- xk1)**3
          B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
          ffz(1,1,0) = 6*A*xk + 2*B
          ffz(1,1,1) = 6*A*xk1 + 2*B
  
          do k = 0,1
             do j = 0,1
                do i = 0,1
                   density(0,i,j,k) = f(1,i,j,k)
                   density(1,i,j,k) = ffx(i,j,k)
                   density(2,i,j,k) = ffy(i,j,k)
                   density(3,i,j,k) = ffz(i,j,k)
                enddo
             enddo
          enddo
          
          
          !now evaluate the cubic hermite polynomial for rho, and derivatives for rho_gradient.   
          !rho values and double derivatives now serve as boundary conditions.
          !Thanks to Andris Gulans for the following code:
          
          x1(:) = x(:)
          xk0 =  coord_v2(1,1)
          xk1 = coord_v2(1,2)
          dx = xk1 - xk0  !assume cube
          invdx = 1d0 / dx
          
          a1=(xk1-x1(1)) * invdx
          b1=(x1(1)-xk0) * invdx
          c1=(a1**3 - a1)*sixth*dx**2
          d1=(b1**3 - b1)*sixth*dx**2
          e1=-(3*a1**2 - 1d0)*sixth*(dx**2)!*(-invdx)
          f1=(3*b1**2 - 1d0)*sixth*(dx**2)!*(invdx)
          
          xk0 =  coord_v2(2,1)
          xk1 = coord_v2(2,2)
          a2=(xk1-x1(2)) * invdx
          b2=(x1(2)-xk0) * invdx
          c2=(a2**3 - a2)*sixth*dx**2
          d2=(b2**3 - b2)*sixth*dx**2
          e2=-(3*a2**2 - 1d0)*sixth*(dx**2)!*(-invdx)
          f2=(3*b2**2 - 1d0)*sixth*(dx**2)!*(invdx)
          
          xk0 =  coord_v2(3,1)
          xk1 = coord_v2(3,2)
          a3=(xk1-x1(3)) * invdx
          b3=(x1(3)-xk0) * invdx
          c3=(a3**3-a3)*sixth*dx**2
          d3=(b3**3-b3)*sixth*dx**2
          e3=-(3*a3**2 - 1d0)*sixth*(dx**2)!*(-invdx)
          f3=(3*b3**2 - 1d0)*sixth*(dx**2)!*(invdx)
          
          a1a2=a1*a2
          a1a3=a1*a3
          a2a3=a2*a3
          b1b2=b1*b2
          b1b3=b1*b3
          b2b3=b2*b3
          a1b2=a1*b2
          a1b3=a1*b3
          a2b3=a2*b3
          b1a2=b1*a2
          b1a3=b1*a3
          b2a3=b2*a3
          
          t1=(a1**2-1d0)*sixth*dx**2
          t2=(a2**2-1d0)*sixth*dx**2
          t3=(a3**2-1d0)*sixth*dx**2
          
          u1=(b1**2-1d0)*sixth*dx**2
          u2=(b2**2-1d0)*sixth*dx**2
          u3=(b3**2-1d0)*sixth*dx**2
          
          !write(use_unit,*)" "
          !write(use_unit,*)"a's", a1,a2,a3
          !write(use_unit,*)"b's", b1, b2,b3
          !write(use_unit,*)"c's", c1, c2,c3
          !write(use_unit,*)"d's", d1, d2,d3
          !write(use_unit,*)"e's", e1, e2,e3
          !write(use_unit,*)"f's", f1, f2,f3
          !
          !write(use_unit,*)" "
         
          
          klo1=0
          klo2=0
          klo3=0
          
          answer = 0d0
          nabla(:) = 0d0
          
         
         answer=       a1*a2*a3*density(0,klo1,klo2,klo3)+b1*a2*a3*density(0,klo1+1,klo2,klo3)
         answer=answer+c1*a2*a3*density(1,klo1,klo2,klo3)+d1*a2*a3*density(1,klo1+1,klo2,klo3)
         answer=answer+a1*c2*a3*density(2,klo1,klo2,klo3)+b1*c2*a3*density(2,klo1+1,klo2,klo3)
         answer=answer+a1*a2*c3*density(3,klo1,klo2,klo3)+b1*a2*c3*density(3,klo1+1,klo2,klo3)
         !   write(use_unit,*)"answer1", answer    
      
         nabla(1)=          - a2a3*density(0,klo1,klo2,klo3)+   a2a3*density(0,klo1+1,klo2,klo3)
         ! nabla(1)= e1*a2a3*density(1,klo1,klo2,klo3)+f1*a2a3*density(1,klo1+1,klo2,klo3) !Just for testing
         nabla(1)=nabla(1)+e1*a2a3*density(1,klo1,klo2,klo3)+f1*a2a3*density(1,klo1+1,klo2,klo3)
         nabla(1)=nabla(1)-t2*a2a3*density(2,klo1,klo2,klo3)+t2*a2a3*density(2,klo1+1,klo2,klo3)
         nabla(1)=nabla(1)-t3*a2a3*density(3,klo1,klo2,klo3)+t3*a2a3*density(3,klo1+1,klo2,klo3)

         nabla(2)=          -  a1a3*density(0,klo1,klo2,klo3)-   b1a3*density(0,klo1+1,klo2,klo3)
        ! nabla(2)= -t1*a1a3*density(1,klo1,klo2,klo3)-u1*b1a3*density(1,klo1+1,klo2,klo3) !just for testing
         nabla(2)=nabla(2)-t1*a1a3*density(1,klo1,klo2,klo3)-u1*b1a3*density(1,klo1+1,klo2,klo3)
         nabla(2)=nabla(2)+e2*a1a3*density(2,klo1,klo2,klo3)+e2*b1a3*density(2,klo1+1,klo2,klo3)
         nabla(2)=nabla(2)-t3*a1a3*density(3,klo1,klo2,klo3)-t3*b1a3*density(3,klo1+1,klo2,klo3)

         nabla(3)=        -   a1a2*density(0,klo1,klo2,klo3)-   b1a2*density(0,klo1+1,klo2,klo3)
       ! nabla(3)= -t1*a1a2*density(1,klo1,klo2,klo3)-u1*b1a2*density(1,klo1+1,klo2,klo3)!just for testing
         nabla(3)=nabla(3)-t1*a1a2*density(1,klo1,klo2,klo3)-u1*b1a2*density(1,klo1+1,klo2,klo3)
         nabla(3)=nabla(3)-t2*a1a2*density(2,klo1,klo2,klo3)-t2*b1a2*density(2,klo1+1,klo2,klo3)
         nabla(3)=nabla(3)+e3*a1a2*density(3,klo1,klo2,klo3)+e3*b1a2*density(3,klo1+1,klo2,klo3)


         answer=answer+a1*b2*a3*density(0,klo1,klo2+1,klo3)+b1*b2*a3*density(0,klo1+1,klo2+1,klo3)
         answer=answer+c1*b2*a3*density(1,klo1,klo2+1,klo3)+d1*b2*a3*density(1,klo1+1,klo2+1,klo3)
         answer=answer+a1*d2*a3*density(2,klo1,klo2+1,klo3)+b1*d2*a3*density(2,klo1+1,klo2+1,klo3)
         answer=answer+a1*b2*c3*density(3,klo1,klo2+1,klo3)+b1*b2*c3*density(3,klo1+1,klo2+1,klo3)
   !  write(use_unit,*)"answer2", answer    

         nabla(1)=nabla(1)-   b2a3*density(0,klo1,klo2+1,klo3)+   b2a3*density(0,klo1+1,klo2+1,klo3)
         nabla(1)=nabla(1)+e1*b2a3*density(1,klo1,klo2+1,klo3)+f1*b2a3*density(1,klo1+1,klo2+1,klo3)
         nabla(1)=nabla(1)-u2*b2a3*density(2,klo1,klo2+1,klo3)+u2*b2a3*density(2,klo1+1,klo2+1,klo3)
         nabla(1)=nabla(1)-t3*b2a3*density(3,klo1,klo2+1,klo3)+t3*b2a3*density(3,klo1+1,klo2+1,klo3)

         nabla(2)=nabla(2)+   a1a3*density(0,klo1,klo2+1,klo3)+   b1a3*density(0,klo1+1,klo2+1,klo3)
         nabla(2)=nabla(2)+t1*a1a3*density(1,klo1,klo2+1,klo3)+u1*b1a3*density(1,klo1+1,klo2+1,klo3)
         nabla(2)=nabla(2)+f2*a1a3*density(2,klo1,klo2+1,klo3)+f2*b1a3*density(2,klo1+1,klo2+1,klo3)
         nabla(2)=nabla(2)+t3*a1a3*density(3,klo1,klo2+1,klo3)+t3*b1a3*density(3,klo1+1,klo2+1,klo3)

         nabla(3)=nabla(3)-   a1b2*density(0,klo1,klo2+1,klo3)-   b1b2*density(0,klo1+1,klo2+1,klo3)
         nabla(3)=nabla(3)-t1*a1b2*density(1,klo1,klo2+1,klo3)-u1*b1b2*density(1,klo1+1,klo2+1,klo3)
         nabla(3)=nabla(3)-u2*a1b2*density(2,klo1,klo2+1,klo3)-u2*b1b2*density(2,klo1+1,klo2+1,klo3)
         nabla(3)=nabla(3)+e3*a1b2*density(3,klo1,klo2+1,klo3)+e3*b1b2*density(3,klo1+1,klo2+1,klo3)


         answer=answer+a1*a2*b3*density(0,klo1,klo2,klo3+1)+b1*a2*b3*density(0,klo1+1,klo2,klo3+1)
         answer=answer+c1*a2*b3*density(1,klo1,klo2,klo3+1)+d1*a2*b3*density(1,klo1+1,klo2,klo3+1)
         answer=answer+a1*c2*b3*density(2,klo1,klo2,klo3+1)+b1*c2*b3*density(2,klo1+1,klo2,klo3+1)
         answer=answer+a1*a2*d3*density(3,klo1,klo2,klo3+1)+b1*a2*d3*density(3,klo1+1,klo2,klo3+1)
  !   write(use_unit,*)a1*a2*b3, density(0,klo1,klo2,klo3+1) 
  !       write(use_unit,*)"answer3", answer    

         nabla(1)=nabla(1)-   a2b3*density(0,klo1,klo2,klo3+1)+   a2b3*density(0,klo1+1,klo2,klo3+1)
         nabla(1)=nabla(1)+e1*a2b3*density(1,klo1,klo2,klo3+1)+f1*a2b3*density(1,klo1+1,klo2,klo3+1)
         nabla(1)=nabla(1)-t2*a2b3*density(2,klo1,klo2,klo3+1)+t2*a2b3*density(2,klo1+1,klo2,klo3+1)
         nabla(1)=nabla(1)-u3*a2b3*density(3,klo1,klo2,klo3+1)+u3*a2b3*density(3,klo1+1,klo2,klo3+1)

         nabla(2)=nabla(2)-   a1b3*density(0,klo1,klo2,klo3+1)-   b1b3*density(0,klo1+1,klo2,klo3+1)
         nabla(2)=nabla(2)-t1*a1b3*density(1,klo1,klo2,klo3+1)-u1*b1b3*density(1,klo1+1,klo2,klo3+1)
         nabla(2)=nabla(2)+e2*a1b3*density(2,klo1,klo2,klo3+1)+e2*b1b3*density(2,klo1+1,klo2,klo3+1)
         nabla(2)=nabla(2)-u3*a1b3*density(3,klo1,klo2,klo3+1)-u3*b1b3*density(3,klo1+1,klo2,klo3+1)

         nabla(3)=nabla(3)+   a1a2*density(0,klo1,klo2,klo3+1)+   b1a2*density(0,klo1+1,klo2,klo3+1)
         nabla(3)=nabla(3)+t1*a1a2*density(1,klo1,klo2,klo3+1)+u1*b1a2*density(1,klo1+1,klo2,klo3+1)
         nabla(3)=nabla(3)+t2*a1a2*density(2,klo1,klo2,klo3+1)+t2*b1a2*density(2,klo1+1,klo2,klo3+1)
         nabla(3)=nabla(3)+f3*a1a2*density(3,klo1,klo2,klo3+1)+f3*b1a2*density(3,klo1+1,klo2,klo3+1)
         
        
         answer=answer+a1*b2*b3*density(0,klo1,klo2+1,klo3+1)+b1*b2*b3*density(0,klo1+1,klo2+1,klo3+1)
         answer=answer+c1*b2*b3*density(1,klo1,klo2+1,klo3+1)+d1*b2*b3*density(1,klo1+1,klo2+1,klo3+1)
         answer=answer+a1*d2*b3*density(2,klo1,klo2+1,klo3+1)+b1*d2*b3*density(2,klo1+1,klo2+1,klo3+1)
         answer=answer+a1*b2*d3*density(3,klo1,klo2+1,klo3+1)+b1*b2*d3*density(3,klo1+1,klo2+1,klo3+1)
    ! write(use_unit,*)"answer4", answer    

         nabla(1)=nabla(1)-   b2b3*density(0,klo1,klo2+1,klo3+1)+   b2b3*density(0,klo1+1,klo2+1,klo3+1)
         nabla(1)=nabla(1)+e1*b2b3*density(1,klo1,klo2+1,klo3+1)+f1*b2b3*density(1,klo1+1,klo2+1,klo3+1)
         nabla(1)=nabla(1)-u2*b2b3*density(2,klo1,klo2+1,klo3+1)+u2*b2b3*density(2,klo1+1,klo2+1,klo3+1)
         nabla(1)=nabla(1)-u3*b2b3*density(3,klo1,klo2+1,klo3+1)+u3*b2b3*density(3,klo1+1,klo2+1,klo3+1)

         nabla(2)=nabla(2)+   a1b3*density(0,klo1,klo2+1,klo3+1)+   b1b3*density(0,klo1+1,klo2+1,klo3+1)
         nabla(2)=nabla(2)+t1*a1b3*density(1,klo1,klo2+1,klo3+1)+u1*b1b3*density(1,klo1+1,klo2+1,klo3+1)
         nabla(2)=nabla(2)+f2*a1b3*density(2,klo1,klo2+1,klo3+1)+f2*b1b3*density(2,klo1+1,klo2+1,klo3+1)
         nabla(2)=nabla(2)+u3*a1b3*density(3,klo1,klo2+1,klo3+1)+u3*b1b3*density(3,klo1+1,klo2+1,klo3+1)

         nabla(3)=nabla(3)+   a1b2*density(0,klo1,klo2+1,klo3+1)+   b1b2*density(0,klo1+1,klo2+1,klo3+1)
         nabla(3)=nabla(3)+t1*a1b2*density(1,klo1,klo2+1,klo3+1)+u1*b1b2*density(1,klo1+1,klo2+1,klo3+1)
         nabla(3)=nabla(3)+u2*a1b2*density(2,klo1,klo2+1,klo3+1)+u2*b1b2*density(2,klo1+1,klo2+1,klo3+1)
         nabla(3)=nabla(3)+f3*a1b2*density(3,klo1,klo2+1,klo3+1)+f3*b1b2*density(3,klo1+1,klo2+1,klo3+1)

         nabla(:) = nabla(:) * invdx
         !write(use_unit,*)"answerfinal", answer
         !write(use_unit,*)"nabla(:)",nabla(:)
         
         vals(1) = answer
         vals(2:4) = nabla(:)
         

      endif


       if(writ)   write(use_unit,*)"tree_interp_v3: final vals",vals   
       if(writ)   write(use_unit,*)"x coordinate for interp:", x(:)    
       
       
    endif
    
    if(writ)   write(use_unit,*)"tree_interp recursion ended"
    
  end function tree_interp_v3
  








subroutine get_herm_charge( root, coord_v2, q  )
!integrate hermite polynomial to get charge in one tree node. 

  implicit none
  
  real*8, dimension(3,2) :: coord_v2 
  real*8 q  !integral of density for the whole cube volume.
  real*8, dimension(3,8,8) :: c_coord
    real*8, dimension(3,8):: coord,c_origin
    real*8, dimension(3,2,8) :: cc_

    real*8, dimension(3) :: x, x1
    real*8, dimension(4) :: vals
    real*8, dimension(4,4) :: xedge
    real*8, dimension(4,2) :: yedge

    real*8, dimension(4) :: y12,y34,y56,y78,y1234,y5678,yfinal
    real*8 dx,dy,dz,u0,v0,w0,v1,w1, invdx, m0,m1,p1,p0,t, xx, xk, xk1, xk0, f0
    real*8, dimension(3) :: ddx, ddy, ddz
    type(treenode) :: root

    real*8 density(0:3,0:1,0:1,0:1)
    real*8, dimension(3) :: nabla
    real*8 answer, invdx3
   real*8 ffx( 0:1,0:1,0:1), ffy( 0:1,0:1,0:1), ffz( 0:1,0:1,0:1)
   real*8 f(4,0:1,0:1,0:1)
    real*8 a,b,c,d
    real*8 a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,f1,f2,f3,t1,t2,t3,u1,u2,u3
    real*8 a1a2, a1a3, a2a3, b1b2,b1b3,b2b3, a1b2, a1b3, a2b3, b1a2, b1a3, b2a3
    real*8 sixth  , pref
    integer klo1,klo2,klo3
    parameter(sixth=0.16666666666666666666666666666667d0)
    
    integer i, i_child, j,k, num, count
    logical writ
    
    
    
    !function vals, gradients
    count = 0
    do k = 0,1
       do j = 0,1
          do i = 0,1
             count = count +1
             f(1,i,j,k) = root%val(1,count)
          enddo
       enddo
    enddo
    
    
    !approximate double-derivatives ( a 1D case ).  
    
    
    xk = coord_v2(1,1)
    xk1 = coord_v2(1,2)
    f0 = root%val(1,1)
    f1 = root%val(1,2)
    m0 = root%val(2,1)
    m1 = root%val(2,2)
    invdx3 = (1 / (xk - xk1)**3)   !assume cube: use this for all directions. 
    
    A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3
    B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
    ffx(0,0,0) = 6*A*xk + 2*B
    ffx(1,0,0) = 6*A*xk1 + 2*B
    
    !second x edge
    
    xk = coord_v2(1,1)
    xk1 = coord_v2(1,2)
    f0 = root%val(1,3)
    f1 = root%val(1,4)
    m0 = root%val(2,3)
    m1 = root%val(2,4)
    A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3
    B =   invdx3* ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
    ffx(0,1,0) = 6*A*xk + 2*B
    ffx(1,1,0) = 6*A*xk1 + 2*B
    
    !third x edge
    
    xk = coord_v2(1,1)
    xk1 = coord_v2(1,2)
    f0 = root%val(1,5)
    f1 = root%val(1,6)
    m0 = root%val(2,5)
    m1 = root%val(2,6)
    A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3
    B =   invdx3* ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
    
    
    ffx(0,0,1) = 6*A*xk + 2*B
    ffx(1,0,1) = 6*A*xk1 + 2*B
    
    !fourth x edge
    
    xk = coord_v2(1,1)
    xk1 = coord_v2(1,2)
    f0 = root%val(1,7)
    f1 = root%val(1,8)
    m0 = root%val(2,7)
    m1 = root%val(2,8)
    A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3
    B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
    ffx(0,1,1) = 6*A*xk + 2*B
    ffx(1,1,1) = 6*A*xk1 + 2*B
    
    
    !y
    
    !first y edge
    
    xk = coord_v2(2,1)
    xk1 = coord_v2(2,2)
    f0 = root%val(1,1)
    f1 = root%val(1,3)
    m0 = root%val(3,1)
    m1 = root%val(3,3)
    A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3
    B =   invdx3* ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
    ffy(0,0,0) = 6*A*xk + 2*B
    ffy(0,1,0) = 6*A*xk1 + 2*B
    
    !second y edge
    
    
    f0 = root%val(1,2)
    f1 = root%val(1,4)
    m0 = root%val(3,2)
    m1 = root%val(3,4)
    A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3
    B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
    ffy(1,0,0) = 6*A*xk + 2*B
    ffy(1,1,0) = 6*A*xk1 + 2*B
    
    !third y edge
    
    
    f0 = root%val(1,5)
    f1 = root%val(1,7)
    m0 = root%val(3,5)
    m1 = root%val(3,7)
    A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3 !/(xk- xk1)**3
    B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
    ffy(0,0,1) = 6*A*xk + 2*B
    ffy(0,1,1) = 6*A*xk1 + 2*B
    
    !fourth y edge
    
    
    f0 = root%val(1,6)
    f1 = root%val(1,8)
    m0 = root%val(3,6)
    m1 = root%val(3,8)
    A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3 !/(xk- xk1)**3
    B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
    ffy(1,0,1) = 6*A*xk + 2*B
    ffy(1,1,1) = 6*A*xk1 + 2*B
    
    
    !write(use_unit,*)"ffy", ffy(0:1,0:1,0:1)
    !write(use_unit,*)" "
    !z
    
    !first z edge
    
    xk = coord_v2(3,1)
    xk1 = coord_v2(3,2)
    f0 = root%val(1,1)
    f1 = root%val(1,5)
    m0 = root%val(4,1)
    m1 = root%val(4,5)
    
    
    A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3 !/(xk- xk1)**3
    B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
    ffz(0,0,0) = 6*A*xk + 2*B
    ffz(0,0,1) = 6*A*xk1 + 2*B
    
    !second z edge
    
    f0 = root%val(1,2)
    f1 = root%val(1,6)
    m0 = root%val(4,2)
    m1 = root%val(4,6)
    A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3 !/(xk- xk1)**3
    B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
    ffz(1,0,0) = 6*A*xk  + 2*B
    ffz(1,0,1) = 6*A*xk1 + 2*B
    
    !third z edge
    
    f0 = root%val(1,3)
    f1 = root%val(1,7)
    m0 = root%val(4,3)
    m1 = root%val(4,7)
    A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3 !/(xk- xk1)**3
    B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
    ffz(0,1,0) = 6*A*xk + 2*B
    ffz(0,1,1) = 6*A*xk1 + 2*B
    
    !fourth z edge
    
    f0 = root%val(1,4)
    f1 = root%val(1,8)
    m0 = root%val(4,4)
    m1 = root%val(4,8)
    A =   ( -2*f0 + 2*f1 + (m0 + m1)*(xk - xk1)   )*invdx3 !/(xk- xk1)**3
    B =   invdx3 * ( 3*f0*(xk + xk1) - 3*f1*(xk + xk1) - (xk -xk1)*(m0*xk + 2*m1*xk + 2*m0*xk1 + m1*xk1)    )
    ffz(1,1,0) = 6*A*xk + 2*B
    ffz(1,1,1) = 6*A*xk1 + 2*B
    
  
    !write(use_unit,*)"ffy",ffy(:,:,:)
    
    do k = 0,1
       do j = 0,1
          do i = 0,1
             density(0,i,j,k) = f(1,i,j,k)
             density(1,i,j,k) = ffx(i,j,k)
             density(2,i,j,k) = ffy(i,j,k)
             density(3,i,j,k) = ffz(i,j,k)
          enddo
       enddo
    enddo
    
    
    !integrated the prefactors;
    
    dx = coord_v2(1,2) - coord_v2(1,1)
    
    a1 = 0.5d0*dx
    a2 = 0.5d0*dx
    a3 = 0.5d0*dx
    b1 = 0.5d0*dx
    b2 = 0.5d0*dx
    b3 = 0.5d0*dx
    pref = 1d0/24d0
    c1 = -dx**3 * pref
    c2 = -dx**3 * pref
    c3 = -dx**3 * pref
    d1 = -dx**3 * pref
    d2 = -dx**3 * pref
    d3 = -dx**3 * pref
    
    
    
    klo1=0
    klo2=0
    klo3=0
    
    answer = 0d0
    
    answer=       a1*a2*a3*density(0,klo1,klo2,klo3)+b1*a2*a3*density(0,klo1+1,klo2,klo3)
    answer=answer+c1*a2*a3*density(1,klo1,klo2,klo3)+d1*a2*a3*density(1,klo1+1,klo2,klo3)
    answer=answer+a1*c2*a3*density(2,klo1,klo2,klo3)+b1*c2*a3*density(2,klo1+1,klo2,klo3)
    answer=answer+a1*a2*c3*density(3,klo1,klo2,klo3)+b1*a2*c3*density(3,klo1+1,klo2,klo3)
    
    answer=answer+a1*b2*a3*density(0,klo1,klo2+1,klo3)+b1*b2*a3*density(0,klo1+1,klo2+1,klo3)
    answer=answer+c1*b2*a3*density(1,klo1,klo2+1,klo3)+d1*b2*a3*density(1,klo1+1,klo2+1,klo3)
    answer=answer+a1*d2*a3*density(2,klo1,klo2+1,klo3)+b1*d2*a3*density(2,klo1+1,klo2+1,klo3)
    answer=answer+a1*b2*c3*density(3,klo1,klo2+1,klo3)+b1*b2*c3*density(3,klo1+1,klo2+1,klo3)
    
    answer=answer+a1*a2*b3*density(0,klo1,klo2,klo3+1)+b1*a2*b3*density(0,klo1+1,klo2,klo3+1)
    answer=answer+c1*a2*b3*density(1,klo1,klo2,klo3+1)+d1*a2*b3*density(1,klo1+1,klo2,klo3+1)
    answer=answer+a1*c2*b3*density(2,klo1,klo2,klo3+1)+b1*c2*b3*density(2,klo1+1,klo2,klo3+1)
    answer=answer+a1*a2*d3*density(3,klo1,klo2,klo3+1)+b1*a2*d3*density(3,klo1+1,klo2,klo3+1)
    
    answer=answer+a1*b2*b3*density(0,klo1,klo2+1,klo3+1)+b1*b2*b3*density(0,klo1+1,klo2+1,klo3+1)
    answer=answer+c1*b2*b3*density(1,klo1,klo2+1,klo3+1)+d1*b2*b3*density(1,klo1+1,klo2+1,klo3+1)
    answer=answer+a1*d2*b3*density(2,klo1,klo2+1,klo3+1)+b1*d2*b3*density(2,klo1+1,klo2+1,klo3+1)
    answer=answer+a1*b2*d3*density(3,klo1,klo2+1,klo3+1)+b1*b2*d3*density(3,klo1+1,klo2+1,klo3+1)
    
    q = answer
    
    
    
  end subroutine get_herm_charge
  











recursive subroutine add_child(eps,coord_v2,root,reclevel,deepest,count_f11)
  
  ! Octree for density interpolation. Charge inside node is calculated two or more different ways,
  ! and if these differ significantly, divide node. 
  ! Visualize tree by uncommenting write statements at end. 
  ! SAG
  
  use dimensions, only: n_centers_basis_integrals
  use pbc_lists, only: centers_basis_integrals

  implicit none 
  real*8 x1,x2,eps,max_dx,min_dx, q
  real*8  interval, error
  real*8, dimension(4) :: center_est, center_exact, topface, bottomface
  integer i, j,k,i_child,count,acc_not_achieved
  integer  acc_achieved,divide_overly_large, to_be_divided, exact_evals,count_f11
  real*8, dimension(18) ::  ydiff, ydiff_test
  real*8 diff, diff1, diff2, diff3
  real*8, dimension(4,18) :: yest, yest_test
  real*8, dimension(4,19) :: ym
  real*8, dimension(3,8) :: c_origin,coord, c_origin1
  real*8, dimension(3,2) :: coord_v2 
  real*8, dimension(3,3,3,3) :: coordinate
  real*8, dimension(3,3,3,4) :: yvalue
  real*8, dimension(3,8,8) :: c_coord,cc_coord
  real*8, dimension(3,2,8) :: ccc_coord   
  real*8, dimension(3) :: dx, dy, dz
  integer reclevel, interpol_evals, deepest, num
  type(treenode) :: root
  logical divide,writ, linear, hermite, store_est, shift
  parameter(max_dx=1.d0)
  parameter(min_dx=0.000d0)
  real*8 p0, p1, m0, m1, xk, x, xk1, t, charge1, chargesum1, charge2, charge2b
  real*8, dimension(4) :: chargesum2

  real*8, dimension(3) :: box_corner, old_coord, coord_diff
  logical, dimension(8) :: outside_uc
  integer :: i_corner, i_atom
  real*8, parameter :: tol = 1.0e-6
  integer, dimension(8) :: closest_atom
  real*8 :: dist_tab_sq(n_centers_basis_integrals)
  real*8 :: dir_tab(3,n_centers_basis_integrals)
  real*8 :: dist_tab(n_centers_basis_integrals)

  writ = .false.  !debugging output
  divide=.false.

  if(writ)write(use_unit,*)"coord_v2:"  
  if(writ)write(use_unit,*)coord_v2(:,:)
  if(writ)write(use_unit,*)"vals"
  if(writ)write(use_unit,*)root%val(:,:)
  
  interval = abs(coord_v2(1,2) - coord_v2(1,1))
 
  if(writ) write(use_unit,*)"dx ",interval
  if((reclevel).gt.deepest)then
     deepest = reclevel 
  endif
  
  
  !Get 27 cube point coordinates
  dx = abs(coord_v2(1,2)-coord_v2(1,1))*0.5d0*(/1d0,0d0,0d0/)
  dy = abs(coord_v2(2,2)-coord_v2(2,1))*0.5d0*(/0d0,1d0,0d0/)
  dz = abs(coord_v2(3,2)-coord_v2(3,1))*0.5d0*(/0d0,0d0,1d0/)
  
  
  do k = 1,3
     do j = 1,3
        do i = 1,3
           
           coordinate(:,i,j,k) = coord_v2(:,1) + (i-1)*dx + (j-1)*dy + (k-1)*dz           
           
        enddo
     enddo
  enddo
  
  center_exact(:) = fmp(coordinate(:,2,2,2) ,count_f11)   
  
 
  
  !Calculate charge in node. Various ways to do that: use No. 6. 
  
  if(.false.)then  !option 1    !not used.      
     !calculating charge in the cube two different ways
     
     !without center point
     chargesum1 = 0d0
     do i = 1,8
        chargesum1 = chargesum1 + root%val(1,i)
     enddo
     charge1 = (interval**3)*(0.125d0)*chargesum1
     
     !with center point
     chargesum2(:) = 0d0
     chargesum2(1) = (0.0625d0)*chargesum1 + 0.5d0*center_exact(1)
     charge2 = (interval**3)*chargesum2(1)
     
     diff = abs(charge1 - charge2)
     if(diff.gt.eps) divide = .true.   
     
     
  endif
   
  
  
  if(.false.)then !option 2     !not used. 
     
     !calculating charge in the cube two different ways
     
     !charge from vertices only
     chargesum1 = 0d0
     do i = 1,8
        chargesum1 = chargesum1 + root%val(1,i)
     enddo
     charge1 = (interval**3d0)*(0.125d0)*chargesum1
     
     !now with with interpolated center, face, edges.     
     chargesum2(:) = 0d0
     
     !add vertices
     do i = 1,8
        chargesum2(1) = chargesum2(1) + root%val(1,i)   
     enddo
     
     
     !add facepoints (interp)
     chargesum2(:) = chargesum2(:) + 2d0 * (  tree_interp_v3( coordinate(:,2,1,1) , root, coord_v2) + & 
          tree_interp_v3( coordinate(:,3,2,1) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,1,2,1) , root, coord_v2)  + &
          tree_interp_v3( coordinate(:,2,3,1) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,1,1,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,3,1,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,1,3,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,3,3,2) , root, coord_v2)+ & 
          tree_interp_v3( coordinate(:,2,1,3) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,3,2,3), root, coord_v2) + &
          tree_interp_v3( coordinate(:,1,2,3) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,2,3,3) , root, coord_v2)  )
     
     
     !add edges   (interp)
     chargesum2(:) = chargesum2(:) + 4d0 * ( tree_interp_v3( coordinate(:,2,2,1) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,2,1,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,3,2,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,1,2,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,2,3,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,2,2,3) , root, coord_v2) )
          
     !add center   (interp)
     chargesum2(:) = chargesum2(:) + 8d0* tree_interp_v3( coordinate(:,2,2,2) , root, coord_v2)
          
     !normalize, multiply by total volume:
     charge2 = chargesum2(1) * (interval**3d0) * (0.015625d0)
     
     diff = abs(charge1 - charge2)
     if(diff.gt.eps) divide = .true.
     
  endif
  
  
  
  
  
  if(.false.)then !option 3    !not used.
     
     !calculating charge in the cube two different ways
     
     !charge from vertices only
     chargesum1 = 0d0
     do i = 1,8
        chargesum1 = chargesum1 + root%val(1,i)
     enddo
     charge1 = (interval**3d0)*(0.125d0)*chargesum1
     
     
     !use center, face and edge points(with appropriate weights). 
     chargesum2(:) = 0d0
     
     !add vertices
     do i = 1,8
        chargesum2(1) = chargesum2(1) + root%val(1,i)   
     enddo
     
     !add facepoints (EXACT
     chargesum2(:) = chargesum2(:) + 2d0 * (  fmp( coordinate(:,2,1,1) ,count_f11) + &
          fmp( coordinate(:,3,2,1) ,count_f11)+ &
          fmp( coordinate(:,1,2,1) ,count_f11)+ &
          fmp( coordinate(:,2,3,1) ,count_f11)+ &
          fmp( coordinate(:,1,1,2) ,count_f11)+ &
          fmp( coordinate(:,3,1,2) ,count_f11)+ &
          fmp( coordinate(:,1,3,2) ,count_f11)+ &
          fmp( coordinate(:,3,3,2) ,count_f11)+ &
          fmp( coordinate(:,2,1,3) ,count_f11)+ &
          fmp( coordinate(:,3,2,3),count_f11)+ &
          fmp( coordinate(:,1,2,3) ,count_f11)+ &
          fmp( coordinate(:,2,3,3) ,count_f11) )
     
     !add edges   (EXACT
     chargesum2(:) = chargesum2(:) + 4d0 * ( fmp( coordinate(:,2,2,1) ,count_f11)+ &
          fmp( coordinate(:,2,1,2) ,count_f11)+ &
          fmp( coordinate(:,3,2,2) ,count_f11)+ &
          fmp( coordinate(:,1,2,2) ,count_f11)+ &
          fmp( coordinate(:,2,3,2) ,count_f11)+ &
          fmp( coordinate(:,2,2,3) ,count_f11)  )
     
     !add center   (interp)
     chargesum2(:) = chargesum2(:) + 8d0 * center_exact(1)
     
     !normalize, multiply by total volume:
     charge2 = chargesum2(1) * (interval**3d0) * (0.015625d0)
     
     diff = abs(charge1 - charge2)
     if(diff.gt.eps) divide = .true.
     
     
  endif
  
  
  
  
  
  if(.false.)then  !option 4.
     !combine options 1 and 2.  Not used. 
     
     !OPTION 1 
     
     !without center point
     chargesum1 = 0d0
     do i = 1,8
        chargesum1 = chargesum1 + root%val(1,i)
     enddo
     charge1 = (interval**3)*(0.125d0)*chargesum1
     
     !with center point
     chargesum2(:) = 0d0
     chargesum2(1) = (0.0625d0)*chargesum1 + 0.5d0*center_exact(1)
     charge2 = (interval**3)*chargesum2(1)
     
     
     !OPTION 2
     !interpolated edge, face, center
     chargesum2(:) = 0d0
     
     do i = 1,8
        chargesum2(1) = chargesum2(1) + root%val(1,i)   
     enddo
     
     !add facepoints (interp)
     chargesum2(:) = chargesum2(:) + 2d0 * (  tree_interp_v3( coordinate(:,2,1,1) , root, coord_v2) + & 
          tree_interp_v3( coordinate(:,3,2,1) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,1,2,1) , root, coord_v2)  + &
          tree_interp_v3( coordinate(:,2,3,1) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,1,1,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,3,1,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,1,3,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,3,3,2) , root, coord_v2)+ & 
          tree_interp_v3( coordinate(:,2,1,3) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,3,2,3), root, coord_v2) + &
          tree_interp_v3( coordinate(:,1,2,3) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,2,3,3) , root, coord_v2)  )
     
     !add edges   (interp)
     chargesum2(:) = chargesum2(:) + 4d0 * ( tree_interp_v3( coordinate(:,2,2,1) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,2,1,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,3,2,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,1,2,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,2,3,2) , root, coord_v2) + &
          tree_interp_v3( coordinate(:,2,2,3) , root, coord_v2) )
     
     !add center   (interp)
     chargesum2(:) = chargesum2(:) + 8d0* tree_interp_v3( coordinate(:,2,2,2) , root, coord_v2)
     
     !normalize, multiply by total volume:
     charge2b = chargesum2(1) * (interval**3d0) * (0.015625d0)
     
     
     diff1 = abs(charge1 - charge2)
     diff2 =  abs(charge1 - charge2b)
     
     if((diff1.gt.eps).or.(diff2.gt.eps)) divide = .true.
     
     
  endif
  
  
  
  
  
  if(.false.)then  !option 5, integration of Hermite Polynomial to get charge. Not used. 
     
     
     !first without center point
     chargesum1 = 0d0
     do i = 1,8
        chargesum1 = chargesum1 + root%val(1,i)
     enddo
     charge1 = (interval**3)*(0.125d0)*chargesum1   
     
     !now through integration over whole cube.
     charge2 = 0d0
     chargesum2(:) = 0d0
     call  get_herm_charge( root, coord_v2, q  )
     charge2 = q
     
     !compare
     diff = abs(charge1 - charge2)
     if(diff.gt.eps) divide = .true.
     
  endif
  
  
  
  if(.true.)then   !option 6. Used in production.
     !combines options 1,5. 
     
     charge1 = 0d0
     charge2b = 0d0
     charge2 = 0d0
     
   !without center point
     chargesum1 = 0d0
     do i = 1,8
        chargesum1 = chargesum1 + root%val(1,i)
     enddo
     charge1 = (interval**3)*(0.125d0)*chargesum1
     
     !with center point
     chargesum2(:) = 0d0
     chargesum2(1) = (0.0625d0)*chargesum1 + 0.5d0*center_exact(1)
     charge2 = (interval**3)*chargesum2(1)
     
     !with integration
     charge2b = 0d0
     call  get_herm_charge( root, coord_v2, q  )
     charge2b = q
     
     
     diff1 = abs(charge1 - charge2)
     diff2 =  abs(charge1 - charge2b)
     
     if((diff1.gt.eps).or.(diff2.gt.eps)) divide = .true.    
     
     
  endif

  if ((n_periodic > 0).or.(treeflag==4)) then
     outside_uc = .false.
     closest_atom = -1
     i_corner = 0
     do i = 1,2
        do j = 1,2
           do k = 1,2
              box_corner(1) = coord_v2(1,i)
              box_corner(2) = coord_v2(2,j)
              box_corner(3) = coord_v2(3,k)
              i_corner = i_corner + 1
              if (treeflag==4) then
                 call tab_atom_centered_coords_p0 &
                      ( box_corner, &
                      dist_tab_sq, &
                      dir_tab, &
                      n_centers_basis_integrals, centers_basis_integrals )
                 dist_tab = sqrt(dist_tab_sq)
                 do i_atom = 1, n_centers_basis_integrals
                    if (dist_tab(i_atom) < min_atom_atom_tab(i_atom)) closest_atom(i_corner) = i_atom
                 end do
              end if
!!$           if (myid==0.and.ANY(closest_atom > 0)) print *, closest_atom
!!$           if (myid==0) print *, 'dist_tab: ', dist_tab
!!$           if (myid==0) print *, 'min_atom_atom_tab: ', min_atom_atom_tab

              if (n_periodic > 0) then
              ! If, in the periodic case, the entire octree node is outside the unit cell we can leave it as it is
              ! since the interpolation will never get there.
              
              ! If the mapped coordinate is different from the original one we are ouside the unit cell.
              ! Be warned, if the coordinate is on the border of the unit cell it can actually end up
              ! begin mapped. For normal-sized unit cells the max_dx condition below will force the division if
              ! the test results in that the entire unit cell is inside the box.
                 old_coord(:) = box_corner(:)
                 call map_to_center_cell(box_corner(:))
                 ! if (myid == 0) print *, i_corner, box_corner(:), old_coord(:)
                 coord_diff(:) = box_corner(:) - old_coord(:)
                 if (dot_product(coord_diff(:),coord_diff(:)) > tol) outside_uc(i_corner) = .true.
              end if
           end do
        end do
     end do
  end if
  if (n_periodic > 0) then
     if (ALL(outside_uc)) then
        divide = .false.
        ! if (myid == 0) print *,'Outside UC.'
     end if
  end if
  if (treeflag==4) then
     do i_atom = 1, n_centers_basis_integrals
        if (ALL(closest_atom == i_atom)) then
           divide = .false.
!!$        if (myid==0) print *, closest_atom
        end if
     end do
  end if

  !counters
  if(.not.divide)then
     acc_achieved = acc_achieved +1
  endif
  
  if(divide) to_be_divided = to_be_divided + 1    
  
  if(divide.and.(interval.lt.min_dx))then
     acc_not_achieved = acc_not_achieved + 1
  endif
  
  if(.not.divide.and.(interval.gt.max_dx))then
     divide_overly_large = divide_overly_large +1
  endif
  
  
  
  
  
  if(divide.or.(interval.gt.max_dx))then    !Division condition.
     
     !Dividing the node:
     
     
     
     !calculate 12 edges, 6 faces exact vals.  
     
     ym(:,1) = fmp( coordinate(:,2,1,1) ,count_f11)  
     ym(:,2) = fmp( coordinate(:,3,2,1) ,count_f11)
     ym(:,3) = fmp( coordinate(:,1,2,1) ,count_f11)
     ym(:,4) = fmp( coordinate(:,2,3,1) ,count_f11)
     ym(:,5) = fmp( coordinate(:,1,1,2) ,count_f11)
     ym(:,6) = fmp( coordinate(:,3,1,2) ,count_f11)
     ym(:,7) = fmp( coordinate(:,1,3,2) ,count_f11)
     ym(:,8) = fmp( coordinate(:,3,3,2) ,count_f11)
     ym(:,9) = fmp( coordinate(:,2,1,3) ,count_f11)
     ym(:,10) = fmp( coordinate(:,3,2,3),count_f11)
     ym(:,11) = fmp( coordinate(:,1,2,3) ,count_f11)
     ym(:,12) = fmp( coordinate(:,2,3,3) ,count_f11)
     
     ym(:,13) = fmp( coordinate(:,2,2,1) ,count_f11)
     ym(:,14) = fmp( coordinate(:,2,1,2) ,count_f11)
     ym(:,15) = fmp( coordinate(:,3,2,2) ,count_f11)
     ym(:,16) = fmp( coordinate(:,1,2,2) ,count_f11)
     ym(:,17) = fmp( coordinate(:,2,3,2) ,count_f11)
     ym(:,18) = fmp( coordinate(:,2,2,3) ,count_f11)
     
     
     !27 cube points:
     
     yvalue(1,1,1,:) = root%val(:,1)
     yvalue(3,1,1,:) = root%val(:,2)
     yvalue(1,3,1,:) = root%val(:,3)
     yvalue(3,3,1,:) = root%val(:,4)
     yvalue(2,1,1,:) = ym(:,1)
     yvalue(1,2,1,:) = ym(:,3)
     yvalue(3,2,1,:) = ym(:,2)
     yvalue(2,3,1,:) = ym(:,4)
     yvalue(2,2,1,:) = ym(:,13)  
     yvalue(1,1,2,:) = ym(:,5)
     yvalue(3,1,2,:) = ym(:,6)
     yvalue(1,3,2,:) = ym(:,7)
     yvalue(3,3,2,:) = ym(:,8)
     yvalue(2,1,2,:) = ym(:,14)
     yvalue(3,2,2,:) = ym(:,15)
     yvalue(1,2,2,:) = ym(:,16)
     yvalue(2,3,2,:) = ym(:,17)
     yvalue(2,2,2,:) = center_exact(:)
     yvalue(1,1,3,:) = root%val(:,5)
     yvalue(3,1,3,:) = root%val(:,6)
     yvalue(1,3,3,:) = root%val(:,7)
     yvalue(3,3,3,:) = root%val(:,8)
     yvalue(2,1,3,:) = ym(:,9)
     yvalue(1,2,3,:) = ym(:,11)
     yvalue(3,2,3,:) = ym(:,10)
     yvalue(2,3,3,:) = ym(:,12)
     yvalue(2,2,3,:) = ym(:,18)
     

     
     !Treat Children
     
     allocate(root%branch(8))
     
     num =0
     do k = 1,2
        do j = 1,2
           do i = 1,2
              num = num+1
              c_origin(:,num) = coordinate(:,i,j,k)
           enddo
        enddo
     enddo
     
     
     !give values to Children
     
     count = 0
     do k = 1,2
        do j = 1,2
           do i = 1,2
              count = count+1
              root%branch(1)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,1) = c_origin(:,count) 
           enddo
        enddo
     enddo

     count = 0
     do k = 1,2
        do j = 1,2
           do i = 2,3
              count = count+1
              root%branch(2)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,2) = c_origin(:,count) + (coordinate(:,3,1,1) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo

      count = 0
     do k = 1,2
        do j = 2,3
           do i = 1,2
              count = count+1
              root%branch(3)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,3) = c_origin(:,count) + (coordinate(:,1,3,1) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo

      count = 0
     do k = 1,2
        do j = 2,3
           do i = 2,3
              count = count+1
              root%branch(4)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,4) = c_origin(:,count) + (coordinate(:,3,3,1) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo

      count = 0
     do k = 2,3
        do j = 1,2
           do i = 1,2
              count = count+1
              root%branch(5)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,5) = c_origin(:,count) + (coordinate(:,1,1,3) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo

     count = 0
     do k = 2,3
        do j = 1,2
           do i = 2,3
              count = count+1
              root%branch(6)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,6) = c_origin(:,count) + (coordinate(:,3,1,3) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo

      count = 0
      do k = 2,3
         do j = 2,3
            do i = 1,2
               count = count+1
               root%branch(7)%val(:,count) = yvalue(i,j,k,:)
               cc_coord(:,count,7) = c_origin(:,count) + (coordinate(:,1,3,3) - coordinate(:,1,1,1))*0.5 
            enddo
         enddo
      enddo

      count = 0
      do k = 2,3
        do j = 2,3
           do i = 2,3
              count = count+1
              root%branch(8)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,8) = c_origin(:,count) + (coordinate(:,3,3,3) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo
     


    do i_child = 1,8  
       
       nullify(root%branch(i_child)%branch)
       
       if(writ)    write(use_unit,*)"MOVING INTO CHILD NO.",i_child
       if(writ)    write(use_unit,*)"yvals(1) of this child"
       if(writ)    write(use_unit,*) root%branch(i_child)%val(1,:)
       
       ccc_coord(:,1,i_child)=cc_coord(:,1,i_child)
       ccc_coord(:,2,i_child)=cc_coord(:,8,i_child)

       call add_child(eps,ccc_coord(:,:,i_child),root%branch(i_child),&
            reclevel+1, deepest, count_f11)
    enddo
   
    
 else
    
  

!Uncomment to visualize tree:

  
!!$    if(   (abs(coord_v2(3,1)).lt.1d-5)  ) then   !XY
!!$              
!!$       write(144,*)coord_v2(1,1),coord_v2(2,1)!,coord(3,1)
!!$       write(144,*)coord_v2(1,2),coord_v2(2,1)!,coord(3,2)
!!$       write(144,*)coord_v2(1,2),coord_v2(2,2)!,coord(3,4)
!!$       write(144,*)coord_v2(1,1),coord_v2(2,2)!,coord(3,3)
!!$       write(144,*)coord_v2(1,1),coord_v2(2,1)!,coord(3,1)
!!$       write(144,*)" "
!!$       
!!$    endif
!!$
!!$
!!$    
!!$    if(   (abs(coord_v2(2,1)).lt.1d-5)  ) then    !XZ
!!$       
!!$       write(145,*)coord_v2(1,1),coord_v2(3,1)!,coord(3,1)
!!$       write(145,*)coord_v2(1,2),coord_v2(3,1)!,coord(3,2)
!!$       write(145,*)coord_v2(1,2),coord_v2(3,2)!,coord(3,4)
!!$       write(145,*)coord_v2(1,1),coord_v2(3,2)!,coord(3,3)
!!$       write(145,*)coord_v2(1,1),coord_v2(3,1)!,coord(3,1)
!!$       write(145,*)" "
!!$       
!!$    endif
!!$
!!$    if(   (abs(coord_v2(1,1)).lt.1d-1)  ) then   !YZ
!!$       
!!$       write(146,*)coord_v2(2,1),coord_v2(3,1)!,coord(3,1)
!!$       write(146,*)coord_v2(2,2),coord_v2(3,1)!,coord(3,2)
!!$       write(146,*)coord_v2(2,2),coord_v2(3,2)!,coord(3,4)
!!$       write(146,*)coord_v2(2,1),coord_v2(3,2)!,coord(3,3)
!!$       write(146,*)coord_v2(2,1),coord_v2(3,1)!,coord(3,1)
!!$       write(146,*)" "
!!$       
!!$    endif

 endif

end subroutine add_child






recursive subroutine add_child_interpolation(hermite,epsabs,eps,coord_v2,root,count_ended,count_started,& 
     acc_not_achieved, acc_achieved,divide_overly_large,&
     to_be_divided,exact_evals,reclevel,deepest,interpol_evals,count_f11)
  
  !  Octree for density interpolation.  Node division based on center point interpolation accuracy.
  !  Inferior to charge-based octrees. Unused but left here for development purposes.
  !  SAG    
  
  
  implicit none 
  real*8 x1,x2,eps,epsabs,max_dx,min_dx
  real*8  interval, error
  real*8, dimension(4) :: center_est, center_exact, topface, bottomface, center_est_linear
  integer i, j,k,i_child,count, count_ended,count_started,acc_not_achieved
  integer  acc_achieved,divide_overly_large, to_be_divided, exact_evals,count_f11
  real*8, dimension(18) ::  ydiff, ydiff_test
  real*8 diff
  real*8, dimension(4,18) :: yest, yest_test
  real*8, dimension(4,19) :: ym
  real*8, dimension(3,8) :: c_origin,coord, c_origin1
  real*8,dimension(3,2) ::coord_v2 
  real*8, dimension(3,3,3,3) :: coordinate
  real*8, dimension(3,3,3,4) :: yvalue
  real*8, dimension(3,8,8) :: c_coord,cc_coord
  real*8, dimension(3,2,8) :: ccc_coord 
  real*8,dimension(3):: dx, dy, dz
  integer reclevel, interpol_evals, deepest, num
  type(treenode) :: root
  logical divide,writ, linear, hermite, store_est, shift
  parameter(max_dx=1.d0)
  parameter(min_dx=0.000d0)
  real*8 p0, p1, m0, m1, xk, x, xk1, t

  
 

  writ = .false.
  count_started = count_started+1
  if(writ) write(use_unit,*)"total started:",count_started
  if(writ)write(use_unit,*)"coord_v2:"  
  if(writ)write(use_unit,*)coord_v2(:,:)
  if(writ)write(use_unit,*)"vals"
  if(writ)write(use_unit,*)root%val(:,:)
  
  interval = abs(coord_v2(1,2) - coord_v2(1,1))
  if(writ) write(use_unit,*)"dx ",interval
  if((reclevel).gt.deepest)then
     deepest = reclevel 
  endif
  
  

  !Open node coordinates. 
  dx = abs(coord_v2(1,2)-coord_v2(1,1))*0.5d0*(/1d0,0d0,0d0/)
  dy = abs(coord_v2(2,2)-coord_v2(2,1))*0.5d0*(/0d0,1d0,0d0/)
  dz = abs(coord_v2(3,2)-coord_v2(3,1))*0.5d0*(/0d0,0d0,1d0/)
  
  
  do k = 1,3
     do j = 1,3
        do i = 1,3
           coordinate(:,i,j,k) = coord_v2(:,1) + (i-1)*dx + (j-1)*dy + (k-1)*dz           
        enddo
     enddo
  enddo
  
  
  
  center_exact(:) = fmp(coordinate(:,2,2,2) ,count_f11)  !Center POINT 
  center_est(:) = tree_interp_v3( coordinate(:,2,2,2) , root, coord_v2) 
  
  divide=.false.
  
  if(writ)write(use_unit,*)"test3"
  if(writ)write(use_unit,*)" "
  
  
  !this is original interpolation tree condition. Doesn't work with constant rho at center for Xe. 
  diff = abs(center_est(1) - center_exact(1)) 
  if(   (diff.gt.(eps*center_exact(1))) .and. (diff.gt.epsabs) ) divide = .true.
  
  
  if(writ)write(use_unit,*)" "
  if(writ) write(use_unit,*)"test4"
  
  if(.not.divide)then
     acc_achieved = acc_achieved +1
  endif
  
  if(divide) to_be_divided = to_be_divided + 1    
  
  if(divide.and.(interval.lt.min_dx))then
     acc_not_achieved = acc_not_achieved + 1
  endif
  
  if(.not.divide.and.(interval.gt.max_dx))then
     divide_overly_large = divide_overly_large +1
  endif
   
  if(writ)write(use_unit,*)"decision factors:"
  if(writ) write(use_unit,*)center_est(1), center_exact(1)
  if(writ) write(use_unit,*) diff,"vs.", eps*center_exact(1),diff.gt.(eps*center_exact(1)) 
  if(writ) write(use_unit,*)"AND"
  if(writ) write(use_unit,*) diff,"vs.", epsabs, diff.gt.epsabs
  if(writ) write(use_unit,*)divide, interval.gt.max_dx
  
  
  if(divide.or.(interval.gt.max_dx))then     !Division condition.
     
     
     !Dividing the node:
     
     !calculate 12 edges, 6 faces exact vals. 
     
     ym(:,1) = fmp( coordinate(:,2,1,1) ,count_f11)  
     ym(:,2) = fmp( coordinate(:,3,2,1) ,count_f11)
     ym(:,3) = fmp( coordinate(:,1,2,1) ,count_f11)
     ym(:,4) = fmp( coordinate(:,2,3,1) ,count_f11)
     ym(:,5) = fmp( coordinate(:,1,1,2) ,count_f11)
     ym(:,6) = fmp( coordinate(:,3,1,2) ,count_f11)
     ym(:,7) = fmp( coordinate(:,1,3,2) ,count_f11)
     ym(:,8) = fmp( coordinate(:,3,3,2) ,count_f11)
     ym(:,9) = fmp( coordinate(:,2,1,3) ,count_f11)
     ym(:,10) = fmp( coordinate(:,3,2,3),count_f11)
     ym(:,11) = fmp( coordinate(:,1,2,3) ,count_f11)
     ym(:,12) = fmp( coordinate(:,2,3,3) ,count_f11)
     
     ym(:,13) = fmp( coordinate(:,2,2,1) ,count_f11)
     ym(:,14) = fmp( coordinate(:,2,1,2) ,count_f11)
     ym(:,15) = fmp( coordinate(:,3,2,2) ,count_f11)
     ym(:,16) = fmp( coordinate(:,1,2,2) ,count_f11)
     ym(:,17) = fmp( coordinate(:,2,3,2) ,count_f11)
     ym(:,18) = fmp( coordinate(:,2,2,3) ,count_f11)
     
     
     !27 cube points
     
     yvalue(1,1,1,:) = root%val(:,1)
     yvalue(3,1,1,:) = root%val(:,2)
     yvalue(1,3,1,:) = root%val(:,3)
     yvalue(3,3,1,:) = root%val(:,4)
     yvalue(2,1,1,:) = ym(:,1)
     yvalue(1,2,1,:) = ym(:,3)
     yvalue(3,2,1,:) = ym(:,2)
     yvalue(2,3,1,:) = ym(:,4)
     yvalue(2,2,1,:) = ym(:,13)     
     yvalue(1,1,2,:) = ym(:,5)
     yvalue(3,1,2,:) = ym(:,6)
     yvalue(1,3,2,:) = ym(:,7)
     yvalue(3,3,2,:) = ym(:,8)
     yvalue(2,1,2,:) = ym(:,14)
     yvalue(3,2,2,:) = ym(:,15)
     yvalue(1,2,2,:) = ym(:,16)
     yvalue(2,3,2,:) = ym(:,17)
     yvalue(2,2,2,:) = center_exact(:)
     yvalue(1,1,3,:) = root%val(:,5)
     yvalue(3,1,3,:) = root%val(:,6)
     yvalue(1,3,3,:) = root%val(:,7)
     yvalue(3,3,3,:) = root%val(:,8)
     yvalue(2,1,3,:) = ym(:,9)
     yvalue(1,2,3,:) = ym(:,11)
     yvalue(3,2,3,:) = ym(:,10)
     yvalue(2,3,3,:) = ym(:,12)
     yvalue(2,2,3,:) = ym(:,18)
     
     
     
     !Treat Children
     
     allocate(root%branch(8))
     
     num =0
     do k = 1,2
        do j = 1,2
           do i = 1,2
              num = num+1
              c_origin(:,num) = coordinate(:,i,j,k)
           enddo
        enddo
     enddo
     
     
     !yvalues to Children
     
     count = 0
     do k = 1,2
        do j = 1,2
           do i = 1,2
              count = count+1
              root%branch(1)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,1) = c_origin(:,count) 
           enddo
        enddo
     enddo

     count = 0
     do k = 1,2
        do j = 1,2
           do i = 2,3
              count = count+1
              root%branch(2)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,2) = c_origin(:,count) + (coordinate(:,3,1,1) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo

      count = 0
     do k = 1,2
        do j = 2,3
           do i = 1,2
              count = count+1
              root%branch(3)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,3) = c_origin(:,count) + (coordinate(:,1,3,1) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo

      count = 0
     do k = 1,2
        do j = 2,3
           do i = 2,3
              count = count+1
              root%branch(4)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,4) = c_origin(:,count) + (coordinate(:,3,3,1) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo

      count = 0
     do k = 2,3
        do j = 1,2
           do i = 1,2
              count = count+1
              root%branch(5)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,5) = c_origin(:,count) + (coordinate(:,1,1,3) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo
     
     count = 0
     do k = 2,3
        do j = 1,2
           do i = 2,3
              count = count+1
              root%branch(6)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,6) = c_origin(:,count) + (coordinate(:,3,1,3) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo
     
      count = 0
     do k = 2,3
        do j = 2,3
           do i = 1,2
              count = count+1
              root%branch(7)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,7) = c_origin(:,count) + (coordinate(:,1,3,3) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo
     
     count = 0
     do k = 2,3
        do j = 2,3
           do i = 2,3
              count = count+1
              root%branch(8)%val(:,count) = yvalue(i,j,k,:)
              cc_coord(:,count,8) = c_origin(:,count) + (coordinate(:,3,3,3) - coordinate(:,1,1,1))*0.5 
           enddo
        enddo
     enddo
     
     

    do i_child = 1,8  
             
       nullify(root%branch(i_child)%branch)
       
       if(writ)    write(use_unit,*)"MOVING INTO CHILD NO.",i_child
       if(writ)    write(use_unit,*)"yvals(1) of this child"
       if(writ)    write(use_unit,*) root%branch(i_child)%val(1,:)
       
       ccc_coord(:,1,i_child)=cc_coord(:,1,i_child)
       ccc_coord(:,2,i_child)=cc_coord(:,8,i_child)

       call add_child_interpolation(hermite,epsabs,eps,ccc_coord(:,:,i_child),root%branch(i_child),count_ended,count_started,&
            acc_not_achieved, acc_achieved,divide_overly_large,to_be_divided, exact_evals,reclevel+1,&
            deepest, interpol_evals,count_f11)
       if(writ)  write(use_unit,*)"test12"
       
    enddo
    
    !write(use_unit,*)"recursion loop ended"
    !write(11,*)"recursion loop ended"
    
    
 else
    
    



 

!Uncomment to visualize tree:

  
!!$    if(   (abs(coord_v2(3,1)).lt.1d-5)  ) then   !XY
!!$              
!!$       write(144,*)coord_v2(1,1),coord_v2(2,1)!,coord(3,1)
!!$       write(144,*)coord_v2(1,2),coord_v2(2,1)!,coord(3,2)
!!$       write(144,*)coord_v2(1,2),coord_v2(2,2)!,coord(3,4)
!!$       write(144,*)coord_v2(1,1),coord_v2(2,2)!,coord(3,3)
!!$       write(144,*)coord_v2(1,1),coord_v2(2,1)!,coord(3,1)
!!$       write(144,*)" "
!!$       
!!$    endif


!!$    
!!$    if(   (abs(coord_v2(2,1)).lt.1d-5)  ) then    !XZ
!!$       
!!$       write(145,*)coord_v2(1,1),coord_v2(3,1)!,coord(3,1)
!!$       write(145,*)coord_v2(1,2),coord_v2(3,1)!,coord(3,2)
!!$       write(145,*)coord_v2(1,2),coord_v2(3,2)!,coord(3,4)
!!$       write(145,*)coord_v2(1,1),coord_v2(3,2)!,coord(3,3)
!!$       write(145,*)coord_v2(1,1),coord_v2(3,1)!,coord(3,1)
!!$       write(145,*)" "
!!$       
!!$    endif
!!$
!!$    if(   (abs(coord_v2(1,1)).lt.1d-1)  ) then   !YZ
!!$       
!!$       write(146,*)coord_v2(2,1),coord_v2(3,1)!,coord(3,1)
!!$       write(146,*)coord_v2(2,2),coord_v2(3,1)!,coord(3,2)
!!$       write(146,*)coord_v2(2,2),coord_v2(3,2)!,coord(3,4)
!!$       write(146,*)coord_v2(2,1),coord_v2(3,2)!,coord(3,3)
!!$       write(146,*)coord_v2(2,1),coord_v2(3,1)!,coord(3,1)
!!$       write(146,*)" "
!!$       
!!$    endif

    
 endif
 
 
 
 count_ended = count_ended+1


end subroutine add_child_interpolation

 





end module octree_routines
