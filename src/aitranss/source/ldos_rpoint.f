c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           jul 2011
c     last revision:  jul 2012
c###########################################################

      module ldos_rpoint
c     ***************************************************
c     simplest possible simulation of the STM-images: 
c     computes r- and E-resolved local density of states 
c     on a 3d-grid and saves result in a binary .plt file
c     ***************************************************
       use globalvars
       use tools
       implicit none

c      temprorary arrays, seen within the module >>>
c      green's function matrix elements in non-orthogonal basis
       complex (8), allocatable :: gf_mn(:,:), 
     &                             gf_mn_a(:,:), gf_mn_b(:,:)

       contains

       subroutine make_rldos
c      **********************************************************
c      computes r- and E-resolved DOS: saves results in .plt file
c      **********************************************************
        implicit none

        integer, parameter :: plt1 = 3, plt2 = 2
        double precision, parameter :: offset = 1.0
        double precision :: dx = 0.25d0, dy = 0.25d0, dz = 0.25d0
c       double precision :: dx = 1.0d0, dy = 1.0d0, dz = 1.0d0
c       double precision :: dx = 2.5d0, dy = 2.5d0, dz = 2.5d0
        double precision xmin, xmax, ymin, ymax, rpoint(3), 
     &  	         rldos, rldos_a, rldos_b
        integer          ierr, ispin, iat, ofile(-1:2), m, n, 
     &	                 nx, ny, nz, xsegments, ysegments, zsegments,
     &                   nrecord, npoint, nfiles, monitor
	character(32) :: ofilename(-1:2)

        real(4)  stime, ftime, secs
        integer  mnts

        print '(/,a)', 
     &	  ' COMPUTING A(r,E): space & energy resolved spectral density >>>'

c       define a 3d box and a grid onto which the 
c       local spectral function is to be mapped >>>
        
        xmin = atom(1)%pos(1) ; xmax = atom(1)%pos(1)
        ymin = atom(1)%pos(2) ; ymax = atom(1)%pos(2)

        do iat = 2, num_atoms

	 if (atom(iat)%pos(1).lt.xmin) xmin=atom(iat)%pos(1)
	 if (atom(iat)%pos(1).gt.xmax) xmax=atom(iat)%pos(1)

	 if (atom(iat)%pos(2).lt.ymin) ymin=atom(iat)%pos(2)
	 if (atom(iat)%pos(2).gt.ymax) ymax=atom(iat)%pos(2)
	
	end do

        xmin = xmin - offset ; xmax = xmax + offset 
        xsegments = int((xmax - xmin)/dx)
	if (xsegments.eq.0) then
	  xsegments = 1 
	else 
	  dx = (xmax - xmin)/dble(xsegments)
        end if 

        ymin = ymin - offset ; ymax = ymax + offset  
        ysegments = int((ymax - ymin)/dy)
	if (ysegments.eq.0) then
	  ysegments = 1
	else
	  dy = (ymax - ymin)/dble(ysegments)
        end if
       
        zsegments = int(abs(zmax - zmin)/dz)
	if (zsegments.eq.0) then
	  zsegments = 1
	else
	  dz = (zmax - zmin)/dble(zsegments)
        end if

        print '(/,a,$)', '  -- spectral density will be computed in a box with '
	print '(i4,a,i4,a,i4,a)', 
     & 	       xsegments + 1, ' *', ysegments + 1, ' *', zsegments + 1, ' points' 
	print '(a,f9.6,a,f9.6,a,f9.6,a,/)', '  -- grid steps/bins are: dx =', 
     &	       dx, ' a.u.   dy =', dy, ' a.u.   dz =', dz, ' a.u.'

        monitor = int((xsegments+1)*(ysegments+1)/80)
        if (monitor.eq.0) monitor = 1

c       <<< done with definition of the grid
	
c       perform evaluation of A(r,E) 

c       allocate termporary arrays ...
        allocate(gf_mn(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
          print *
          stop '[SUBROUTINE make_rldos]: <gf_mn> allocation failure'
        end if

        if (nspin.eq.2) then 
         allocate(gf_mn_a(nsaos,nsaos),stat=ierr)
         if (ierr.ne.0) then
          print *
          stop '[SUBROUTINE make_rldos]: <gf_mn> allocation failure'
         end if
         allocate(gf_mn_b(nsaos,nsaos),stat=ierr)
         if (ierr.ne.0) then
          print *
          stop '[SUBROUTINE make_rldos]: <gf_mn> allocation failure'
         end if
        end if ! nspin

c       ! computing GF matrix elements
        do ispin = 1, nspin  

c        here we compute matrix elements of the green's function
         call cpu_time(stime)
         if (do_ldos3d_ewindow) then
c	   integrate green's function over the energy window
	   print '(2x,a,i1,a,$)', 'ispin = ', ispin, 
     &		 ' : computing energy integrated matrix elements of the Green''s function ... '
           if (.not.evunits) then 
	    call integrated_gf(ener,eend,.true.,ispin,heigvec(:,:,ispin),invheigvec(:,:,ispin),gf_mn)
	   else 
	    call integrated_gf(efermi+ener/hartree,efermi+eend/hartree,.true.,
     &                         ispin,heigvec(:,:,ispin),invheigvec(:,:,ispin),gf_mn)
	   end if
	 else if (do_ldos3d) then
c          use only the first energy point          
	   print '(2x,a,i1,a,$)', 'ispin = ', ispin, 
     & 	       ' : computing matrix elements of the Green''s function ... '
           if (.not.evunits) then
	    call integrated_gf(ener,0.0d0,.false.,ispin,heigvec(:,:,ispin),invheigvec(:,:,ispin),gf_mn)
           else
	    call integrated_gf(efermi+ener/hartree,0.0d0,.false.,
     &                         ispin,heigvec(:,:,ispin),invheigvec(:,:,ispin),gf_mn)
           end if
         else
	   print '(2x,a)', 'something went wrong here: internal confusion between '
	   print '(2x,a)', '"$ldos 3d" and "$ldos 3d_ewdindow" flags'
	   print '(2x,a)', 'transport module will be terminated now ...'
         end if
c        that is an old call, not used any more >>>
c        call gfmatrix(ener,ispin,tmp_b,tmp_binv,gf_mn)

         if (nspin.eq.2) then
	  select case (ispin)
	   case (1) 
	    forall(m=1:nsaos,n=1:nsaos)
	     gf_mn_a(m,n) = gf_mn(m,n)	
	    end forall    
	   case (2) 
	    forall(m=1:nsaos,n=1:nsaos)
	     gf_mn_b(m,n) = gf_mn(m,n)	
	    end forall    
           end select
	 end if

	 print *, 'done '
         call cpu_time(ftime)
         mnts = int((ftime-stime)/60.0d0,4)
         secs = ftime - stime - mnts*60.0d0
         print '(2x,a,i3,a,f5.2,a,/)',
     &         'time spent:', mnts, ' min ', secs, ' sec'
        
        end do 
c	! ispin : computation of GF matrix elements

        if (scaling_fctr < 0) then
c        ! estimate scaling factor 'scaling_fctr':
c        ! on output to .plt file spectral function A(r,E) will be 
c        ! multiplied by 'scaling_fctr', so that exponentially small 
c        ! density can be correctly visualized/plotted with gOpenMol 
c
c        ! take (xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2 as a reference point
c        ! get exponent, as int[-log10(A)] + 1

         rpoint = (/ (xmin+xmax)/2.0d0, (ymin+ymax)/2.0d0, (zmin+zmax)/2.0d0  /)
         rldos = a_re(rpoint,gf_mn)  
c        ! in case of open shell, gf_mn = gf_mn_b (see above)
         n = int(-log10(rldos)) + 1
         scaling_fctr = dexp( dble(n) * log(10.0d0) ) 
         print '(2x,a)', 'SCALING FACTOR is estimated automatically!'
	else 
c        ! otherwise, a value of 'scaling_fctr' will be taken from <tcontrol>
         print '(2x,a)', 'WARNING: scaling factor is taken from <tcontrol>!'
         print '(2x,a)', '         there is no garantee that gOpenMol can visualize your LDOS properly'
         print '(2x,a)', '         another option would be to comment out a keyword "scaling_fctr" in <tcontrol>'
         print '(2x,a)', '         and rely on the automatic estimation of the scaling factor'
        end if

        if (do_ldos3d) then
	 print '(/,2x,a,e12.6)', 
     &         'on output to .plt file, A(r,E) [in 1/Hartree units] will be multiplied by ', scaling_fctr 
         if (.not.evunits) then
          print '(/,2x,a,f14.10,a,/)', 'A(r,E) is computed at E = ',ener,' H'
         else
          print '(/,2x,a,f14.10,a,/)', 'A(r,E) is computed at E - EF = ',efermi+ener/hartree,' eV'
         end if
        else
	 print '(/,2x,a,e12.6)', 
     &         'on output to .plt file, \int A(r,E)dE will be multiplied by ', scaling_fctr 
         if (.not.evunits) then
          print '(/,2x,a,f14.10,a,f14.10,a/)', 
     &	        'A(r,E) is integrated from E1 = ',ener,' H till E2 = ', eend, ' H'
         else
          print '(/,2x,a,f14.10,a,f14.10,a/)', 
     &	        'A(r,E) is integrated from E1-EF = ',efermi+ener/hartree,
     &          ' eV till E2-EF = ', efermi+eend/hartree, ' eV'
         end if
        end if

c       open output binary files >>>
        if (nspin.eq.2) then ;  nfiles = 4
	                else ;  nfiles = 2
	end if  
        ofile(-1)= 49 
        ofile(0) = 50 ; ofile(1) = 51 ; ofile(2) = 52 
	ofilename(-1) = 'zmap.plt'     
	ofilename(0)  = 'r-resolved.ldos.plt'     
	ofilename(1)  = 'r-resolved.ldos.alpha.plt'     
	ofilename(2)  = 'r-resolved.ldos.beta.plt'     
      
        do n = -1, nfiles-2
          open(unit=ofile(n),
     &	       file=trim(ofilename(n)),access='direct',recl=4,iostat=ierr)	
	  if (ierr.ne.0) then
           stop '[make_rldos] : can not open output .plt-file!'
          end if

c        reference: calcgrd.f from turbomole moloch-library 
c        !!! caution: use ifort with "-assume byterecl" option 
c            to have record length in bytes	
c        header of the plt.file >>>
         write(ofile(n),rec=1) plt1
         write(ofile(n),rec=2) plt2

         write(ofile(n),rec=3) zsegments+1
         write(ofile(n),rec=4) ysegments+1
         write(ofile(n),rec=5) xsegments+1

         write(ofile(n),rec=6)  real(zmin*bohr_radius,4)
         write(ofile(n),rec=7)  real(zmax*bohr_radius,4)
         write(ofile(n),rec=8)  real(ymin*bohr_radius,4)
         write(ofile(n),rec=9)  real(ymax*bohr_radius,4)
         write(ofile(n),rec=10) real(xmin*bohr_radius,4)
         write(ofile(n),rec=11) real(xmax*bohr_radius,4)

c        << header of the plt.file 
c        set a current record number : 
         nrecord = 11

	end do ! files      
c       <<< done with file opening 
      
c       make cicles over grid points >>>
        npoint = 0
        do nz = 0, zsegments 
	  print '(2x,a,i5,a,$)', '### z-layer : ', nz+1, ' : | '          
	  rpoint(3) = zmin + nz * dz
	  do ny = 0, ysegments
	   rpoint(2) = ymin + ny * dy
	   do nx = 0, xsegments

	     rpoint(1) = xmin + nx * dx
             npoint = npoint + 1
             if ( mod(npoint,monitor)==0 ) print '(a,$)', '='  ! monitor

c            shift index of the record >>>
             nrecord = nrecord + 1 
c	     compute spectral function at a given grid point 
c            and save it to an external file 
             if (nspin.eq.1) then 

c             closed shell case >>>
              rldos = a_re(rpoint,gf_mn) 
	      if (rldos.lt.0.0d0) then 
	       write (80,'(a,f10.6,2x,f10.6,2x,f10.6,a,e20.14,a)')
     &          'r = ', rpoint(1), rpoint(2), rpoint(3), 
     &          ' :   A(r,E) = ', rldos, ' is NEGATIVE !!! ' 
              end if
	      write (90,'(a,f10.6,2x,f10.6,2x,f10.6,a,e20.14)')
     &         'r = ', rpoint(1), rpoint(2), rpoint(3), 
     &         ' :   A(r,E) = ', rldos * scaling_fctr 
	      write(ofile(0),rec=nrecord) real(2.0d0 * rldos * scaling_fctr,4)
	      write(ofile(-1),rec=nrecord) real(rpoint(3),4)

             else                              
c             open shell case >>>

c             --- alpha channel ---
              rldos_a = a_re(rpoint,gf_mn_a) 
	      if (rldos_a.lt.0.0d0) then 
	       write (81,'(a,f10.6,2x,f10.6,2x,f10.6,a,e20.14,a)')
     &           'r = ',rpoint(1),rpoint(2),rpoint(3),
     &           ' :   A(r,E) = ',rldos_a,' is NEGATIVE !!! ' 
              end if
	      write (91,'(a,f10.6,2x,f10.6,2x,f10.6,a,e20.14)')
     &           'r = ',rpoint(1),rpoint(2),rpoint(3),
     &            ' :   A(r,E) = ', rldos_a*scaling_fctr 
	      write(ofile(1),rec=nrecord) real(rldos_a * scaling_fctr,4)

c             --- beta channel ---
              rldos_b = a_re(rpoint,gf_mn_b) 
	      if (rldos_b.lt.0.0d0) then 
	       write (82,'(a,f10.6,2x,f10.6,2x,f10.6,a,e20.14,a)')
     &          'r = ',rpoint(1),rpoint(2),rpoint(3),
     &          ' :   A(r,E) = ',rldos_b,' is NEGATIVE !!! ' 
              end if
	      write (92,'(a,f10.6,2x,f10.6,2x,f10.6,a,e20.14)')
     &          'r = ',rpoint(1),rpoint(2),rpoint(3),
     &          ' :   A(r,E) = ',rldos_b*scaling_fctr 
	      write(ofile(2),rec=nrecord) real(rldos_b * scaling_fctr,4)

c             sum of two spin channels goes into a separate file >>>
	      write(ofile(0),rec=nrecord) real((rldos_a+rldos_b)*scaling_fctr,4)
	      write(ofile(-1),rec=nrecord) real(rpoint(3),4)

	     end if  ! closed vc open shell case

	   end do ! loop over x points
	  end do ! loop over y points
	  print '(a)', ' | done'
	end do ! loop over z points

        if (nspin.eq.2) then 
	 print '(/,2x,a,a,a)', 
     &	  'spectral function for a-channel is exported to file <',trim(ofilename(1)),'>'
	 print '(2x,a,a,a)', 
     &	  'spectral function for b-channel is exported to file <',trim(ofilename(2)),'>'
	 print '(2x,a,a,a)', 
     &	  'full spectral function (a+b channels) is exported to file <',trim(ofilename(0)),'>'
	else 
	 print '(/,2x,a,a,a)', 
     &	  'spectral function is exported to file <', trim(ofilename(0)),'>'
        end if
	print '(2x,a,a,a)', 
     &	  'z components of the 3d mesh are exported to file  <',trim(ofilename(-1)),'>'
        do n = -1, nfiles-2
	 close(ofile(n))
	end do

        print '(/,a)', ' <<< DONE WITH A(r,E) '
       
c       deallocate local arrays ...
        deallocate(gf_mn,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE make_rldos]: ',
     &           'impossible to deallocate temp. arrays'
         print *, 'nevertheless, proceed further ...'
        end if

        if (nspin.eq.2) then 
c        extra arrays were open in the open shell case >>>
	 deallocate(gf_mn_a,gf_mn_b,stat=ierr)
         if (ierr.ne.0) then
          print '(/,a,a)', ' [SUBROUTINE make_rldos]: ',
     &           'impossible to deallocate <gf_mn_a> or <gf_mn_b>'
          print *, 'nevertheless, proceed further ...'
	 end if 
        end if

       
       end subroutine make_rldos 

       double precision function basis_func_r(mu,fnorm,r_mu,r_point)
c      ***************************************
c      returns a value of the CGT orbital 'mu' 
c      (localized at r_mu) at given point r 
c      ***************************************       
       implicit none
       
       type(cgto)       mu          ! input: atomic orbital / basis function
       double precision fnorm       ! its norm
       double precision r_mu(3)     ! its localization   
       double precision r_point(3)  ! requested point in 3D    
       
       double precision r_part, omega_part, x, y, z, r2, tmp_coeff, tmp_xi,
     &                  sqrt12, sqrt24, sqrt40, sqrt60
       integer ncontr, n, lm
       
       x=r_point(1)-r_mu(1) ;  y=r_point(2)-r_mu(2) ;  z=r_point(3)-r_mu(3)
       r2 = x*x + y*y + z*z

c      compute radial part >>>
       ncontr = mu%ngto ;  r_part = 0.0d0
       do n = 1, ncontr
        tmp_coeff = mu%icoeff(n) ; tmp_xi = mu%xi(n)
        r_part = r_part + tmp_coeff * dexp(-tmp_xi*r2)
       end do

c      compute angular dependent part >>>

ccc    difine some constants ...
       sqrt12 = sqrt(12.0d0) ;  sqrt24 = sqrt(24.0d0)
       sqrt40 = sqrt(40.0d0) ;  sqrt60 = sqrt(60.0d0)

       lm = mu%lm
       select case (lm)
ccc    s-function
        case(1) ; omega_part= 1.0d0   ! s 
ccc    p-functions
	case(2) ; omega_part= x       ! px
        case(3) ; omega_part= y       ! py 
        case(4) ; omega_part= z       ! pz
ccc    d-functions
        case(5) ; omega_part=(-x*x-y*y+2.0d0*z*z)/sqrt12 !(-xx-yy+2zz)/sqrt(12) 
        case(6) ; omega_part= x*z                        ! xz
        case(7) ; omega_part= y*z                        ! yz 
        case(8) ; omega_part= x*y                        ! xy 
	case(9) ; omega_part=(x*x-y*y)/2.0d0             ! (xx-yy)/2 
ccc    f-functions
        case(10) ; omega_part=(2.0d0*z*z*z-3.0d0*x*x*z-3.0d0*y*y*z)/sqrt60
                                                    ! (2zzz-3xxz-3yyz)/sqrt(60)
        case(11) ; omega_part=(-x*x*x-x*y*y+4.0d0*x*z*z)/sqrt40
	                                            ! (-xxx-xyy+4xzz)/sqrt(40)
        case(12) ; omega_part=(-y*y*y-x*x*y+4.0d0*y*z*z)/sqrt40
	                                            ! (-yyy-xxy+4yzz)/sqrt(40)
        case(13) ; omega_part=x*y*z
	                                            ! xyz

       	case(14) ; omega_part=(x*x*z-y*y*z)/2.0d0
	   	                                    ! (xxz-yyz)/2 

	case(15) ; omega_part=(x*x*x-3.0d0*x*y*y)/sqrt24
	                                            ! (xxx-3xyy)/sqrt(24)

	case(16) ; omega_part=(y*y*y-3.0d0*x*x*y)/sqrt24
	                                            ! (yyy-3xxy)/sqrt(24)
        case default
	 print *
	 stop '<basis_func_r>: lm > 16 ???'
       end select ! case over lm	

       basis_func_r = omega_part * r_part / fnorm
       end function basis_func_r


       subroutine integrated_gf(e1,e2,flag,ispin,input_b,input_binv,gfout)
c       **********************************************************
c       computes green's function in the non-orthogonal CGTO basis
c       g(E) = s^{-1/2} * b * [E-e(p)]^{-1} * b^{-1} * s^{-1/2}
c       **********************************************************
        implicit none

c       energy window: [e1,e2]
        double precision, intent (in) :: e1, e2
c       if flag==.true. integral from e1 to e2 is computed
c       otherwise, green's function is evaluated at point e1            
	logical, intent (in)      :: flag 
        integer,      intent (in) :: ispin
        complex(8),   intent (in) :: input_b(:,:), input_binv(:,:)
c       output array: greens function (GF) for E=epoint
        complex(8), intent (out) :: gfout(:,:)

        complex(8), allocatable :: gftmp(:), bg(:,:), bgb(:,:), 
     &	                           sgf(:,:), s12inv(:,:)
	integer n, m, p, ierr 

        double precision, parameter :: tiny_epsilon = 1.0d-12
        double precision eps1, eps2, eta, tmp_atan1, tmp_atan2
	complex(8) zp

c       allocate temp. arrays: 
        allocate(gftmp(nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE gfmatrix]: <gftmp> allocation failure'
        end if

        allocate(bg(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE gfmatrix]: <bg> allocation failure'
        end if

        allocate(bgb(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE gfmatrix]: <bgb> allocation failure'
        end if

        allocate(sgf(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE gfmatrix]: <sgf> allocation failure'
        end if

        allocate(s12inv(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE gfmatrix]: <s12inv> allocation failure'
        end if

        if (flag) then
c        here we compute analitical formula for 
c        energy intergal: \int_e1^e2 de/(e-z_p)
         do p = 1, nsaos
	  zp   = zpoles(p,ispin)   
	  eps1 = e1 - dble(zp) ; eps2 = e2 - dble(zp)  
	  eta  = -dimag(zp)
	  tmp_atan1 = datan(eps1/eta)
	  tmp_atan2 = datan(eps2/eta)
          gftmp(p) = 0.5d0 * dlog( (eps2*eps2+eta*eta)/(eps1*eps1+eta*eta) ) 
     &             + ione * (tmp_atan1 - tmp_atan2)
	 end do	  
        else
c        evaluate green's function at point e1 only
         do p = 1, nsaos 
           gftmp(p) = cone/(cmplx(e1,0.0d0,8) - zpoles(p,ispin))
         end do
	end if 

        do n = 1, nsaos
	 do p = 1, nsaos
	  bg(n,p) = input_b(n,p) * gftmp(p) 
	 end do
	end do

c       compute gf in the orthogonal basis:  bgb = bg * input_binv >>>
        call dcmplx_ab('n','n',nsaos,bg,input_binv,bgb)

c       transform results to the non-orthogonal basis >>>

c       -- fill up temporary array with overlap elements
        forall(n=1:nsaos,m=1:nsaos)
         s12inv(n,m) = cmplx(smat12inv(n,m),0.0d0,8)
        end forall

c       -- step 1: sgf = smat12inv * bgb
        call dcmplx_ab('n','n',nsaos,s12inv,bgb,sgf)

c       -- step 2: gfout = sgf * smat12inv
        call dcmplx_ab('n','n',nsaos,sgf,s12inv,gfout)
 
        deallocate(gftmp,bg,bgb,sgf,s12inv,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE gfmatrix]: ',
     &           'impossible to deallocate temp. arrays'
         print *, 'nevertheless, proceed further ...'
        end if

       end subroutine integrated_gf

       subroutine gfmatrix(epoint,ispin,input_b,input_binv,gfout)
c       **********************************************************
c       computes green's function in the non-orthogonal CGTO basis
c       g(E) = s^{-1/2} * b * [E-e(p)]^{-1} * b^{-1} * s^{-1/2}
c       **********************************************************
        implicit none

        double precision, intent (in) :: epoint
        integer,      intent (in) :: ispin
        complex(8),   intent (in) :: input_b(:,:), input_binv(:,:)
c       output array: greens function (GF) for E=epoint
        complex(8), intent (out) :: gfout(:,:)

        complex(8), allocatable :: gftmp(:), bg(:,:), bgb(:,:), 
     &	                           sgf(:,:), s12inv(:,:)
	integer n, m, p, ierr 

c       allocate temp. arrays: 
        allocate(gftmp(nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE gfmatrix]: <gftmp> allocation failure'
        end if

        allocate(bg(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE gfmatrix]: <bg> allocation failure'
        end if

        allocate(bgb(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE gfmatrix]: <bgb> allocation failure'
        end if

        allocate(sgf(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE gfmatrix]: <sgf> allocation failure'
        end if

        allocate(s12inv(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE gfmatrix]: <s12inv> allocation failure'
        end if

        do p = 1, nsaos 
          gftmp(p) = cone/(cmplx(epoint,0.0d0,8) - zpoles(p,ispin))
        end do

        do n = 1, nsaos
	 do p = 1, nsaos
	  bg(n,p) = input_b(n,p) * gftmp(p) 
	 end do
	end do

c       compute gf in the orthogonal basis:  bgb = bg * input_binv >>>
        call dcmplx_ab('n','n',nsaos,bg,input_binv,bgb)

c       transform gf to the non-orthogonal basis >>>

c       -- fill up temporary array with overlap elements
        forall(n=1:nsaos,m=1:nsaos)
         s12inv(n,m) = cmplx(smat12inv(n,m),0.0d0,8)
        end forall

c       -- step 1: sgf = smat12inv * bgb
        call dcmplx_ab('n','n',nsaos,s12inv,bgb,sgf)

c       -- step 2: gfout = sgf * smat12inv
        call dcmplx_ab('n','n',nsaos,sgf,s12inv,gfout)
 
        deallocate(gftmp,bg,bgb,sgf,s12inv,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE gfmatrix]: ',
     &           'impossible to deallocate temp. arrays'
         print *, 'nevertheless, proceed further ...'
        end if

       end subroutine gfmatrix

       double precision function a_re(r_point,gf_elements)
c       *******************************************************
c       returns a(r,e) : r- & energy-resolved spectral function
c       *******************************************************      
c
c       -- request for matrix elements of the GF --> upper routine
c       call gfmatrix(e_point,ispin,gf_elements)
c
        implicit none
        double precision r_point(3)
	complex (8)      gf_elements(:,:)   ! gf matrix at e_point    
c       integer          ispin

        integer          iat, iorb, iatype, n, m, norb, ierr
        double precision imgf, r_nu(3), inorm
        type(cgto)       nu
	double precision, allocatable :: phi_r(:), tmpg(:,:), tmp_g_phi(:) 

c       BLAS level-one function to compute a dot product 
c       of the two double precision vectors
        double precision, external :: ddot

c       allocate temporary arrays
        allocate(tmp_g_phi(nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[FUNCTION a_re]: <tmpg> allocation failure'
        end if

        allocate(tmpg(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[FUNCTION a_re]: <tmpg> allocation failure'
        end if
        forall(n=1:nsaos,m=1:nsaos)
	 tmpg(n,m) = -imag(gf_elements(n,m))/pi
	end forall

c       values of CGT-orbital at the fixed r-point are  
c       to be computed here and to be put into temporary array 
        allocate(phi_r(nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[FUNCTION a_re]: <phi_r> allocation failure'
        end if
 
        phi_r = 0.0d0
	n = 0
        do iat = 1, num_atoms
         r_nu   = atom(iat)%pos         ! save atom coordinates
	 iatype = atom(iat)%atype       ! take atom-type index 
         norb   = n_basis_func(iatype) 	    
 	 do iorb = 1, norb
          n = n + 1
          nu = aos(iatype,iorb)
          inorm = orbnorm(n)
c         compute value of the 'nu' CGT-orbital at r-point 
          phi_r(n) = basis_func_r(nu,inorm,r_nu,r_point)  
         end do
	end do 

c       compute GF(x,x;E) * phi(x)  <-- tmp_g_phi
        call dgemv('n',nsaos,nsaos,1.0d0,tmpg,nsaos,phi_r,1,0.0d0,tmp_g_phi,1)
c       compute tmp_g_phi * phi_r --> LDOS(x)
        imgf = ddot(nsaos,tmp_g_phi,1,phi_r,1)

c       deallocate temporary arrays
        deallocate(phi_r,tmpg,tmp_g_phi,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [FUNCTION a_re]: ',
     &           'impossible to deallocate temp. array <phi_r>'
         print *, 'nevertheless, proceed further ...'
        end if

c       returns local spectral density ...
c       -- is it always positive even though gaussian 
c          like basis set is incomplete? 
c       -- if not, warning message will appear in a higher level routine
        a_re = imgf

       end function a_re

      end module ldos_rpoint
      