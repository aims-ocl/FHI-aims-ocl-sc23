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

      module cond_chann_wf
c     ***************************************************************************
c     computes r- and E-resolved "scattering wave functions"
c     defined as eigenvectors of the transmission matrix
c     tt^{+}(E) <-- (Gamma_L)^{1/2} G_LR(E) Gamma_R G^{+}_RL(E) (Gamma_L)^{1/2} ;
c     data are computed on the 3d grid is saved in a binary .plt file:
c     thus a visialization using gOpenMol is possible
c     ***************************************************************************
       use globalvars
       use ldos_rpoint
       use tools
       implicit none

       contains

       subroutine compute_scatt_wave_func(nchannels)
c      *********************************************
c      computes r-resolved scattering wave functions 
c      exports results in .plt files
c      *********************************************
        implicit none

        integer, intent(in) :: nchannels

        integer, parameter :: plt1 = 3, plt2 = 2
        double precision, parameter :: offset = 1.0
        double precision :: dx = 0.25d0, dy = 0.25d0, dz = 0.25d0
c       double precision :: dx = 2.5d0, dy = 2.5d0, dz = 2.5d0
        double precision xmin, xmax, ymin, ymax, rpoint(3) 
        complex(8)       rpsi
        integer          ichannel, ierr, iat, ofile(4), n, nx, ny, nz, 
     &                   xsegments, ysegments, zsegments, 
     &                   nrecord, npoint, nfiles, monitor
	character(32) :: ofilename(4)
	character        channel_index

c       real(4)  stime, ftime, secs
c       integer  mnts

        print '(/,a)', 
     &	  ' COMPUTING SCATTERING WAVE FUNCTIONS (WFs) >>>'

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
       
        zsegments = int(abs(wf_zmax - wf_zmin)/dz)
	if (zsegments.eq.0) then
	  zsegments = 1
	else
	  dz = (wf_zmax - wf_zmin)/dble(zsegments)
        end if

        print '(/,a,$)', '  -- wave functions will be computed in a box with '
	print '(i4,a,i4,a,i4,a)', 
     & 	       xsegments + 1, ' *', ysegments + 1, ' *', zsegments + 1, ' points' 
	print '(a,f9.6,a,f9.6,a,f9.6,a,/)', '  -- grid steps/bins are: dx =', 
     &	       dx, ' a.u.   dy =', dy, ' a.u.   dz =', dz, ' a.u.'

        monitor = int((xsegments+1)*(ysegments+1)/80)
        if (monitor.eq.0) monitor = 1

c       <<< done with definition of the grid
	
c       perform evaluation of A(r,E) 

        if (wf_scaling_fctr < 0) then
c        ! estimate scaling factor 'wf_scaling_fctr':
c        ! on output to .plt file scattering wave functions |Psi_n(x,E)| 
c        ! will be multiplied by 'wf_scaling_fctr', so that exponentially 
c        ! small densities can be correctly visualized/plotted with gOpenMol 
c
c        ! take x=0, y=0, z=wf_zmin as a reference point
c        ! get exponent, as int[-log10(|Psi|)] + 1

         rpoint = (/0.0d0,0.0d0,wf_zmin/)
         rpsi = compute_rpsi(rpoint,utt(:,nsaos,1))  
         n = int( -log10(abs(rpsi)) ) + 1
         wf_scaling_fctr = dexp( dble(n) * log(10.0d0) ) 
         print '(2x,a)', 'SCALING FACTOR is estimated automatically!'
	else 
c        ! otherwise, a value of 'wf_scaling_fctr' will be taken from <tcontrol>
         print '(2x,a)', 'WARNING: scaling factor is taken from <tcontrol>!'
         print '(2x,a)', '         there is no garantee that gOpenMol can visualize your WFs properly'
         print '(2x,a)', '         another option would be to comment out a keyword "wf_scaling_fctr" in <tcontrol>'
         print '(2x,a)', '         and rely on the automatic estimation of the scaling factor'
        end if

	print '(/,2x,a,e12.6)', 
     &         'on output to .plt file, WFs will be multiplied by ', wf_scaling_fctr 

c       cycling over few relevant conduction channels    
        do ichannel = 1, nchannels
         print '(/,2x,a,i2,/)', 'CONDUCTION EIGENCHANNEL # ', ichannel 
c        open output binary files >>>
         if (nspin.eq.2) then ;  nfiles = 4
	                 else ;  nfiles = 2
	 end if  
         ofile(1) = 51 ; ofile(2) = 52 
         ofile(3) = 53 ; ofile(4) = 54 

         write(channel_index,'(i1)') ichannel
         if (nspin.eq.1) then
  	  ofilename(1) = 'psi.'//trim(channel_index)//'.re.plt'
  	  ofilename(2) = 'psi.'//trim(channel_index)//'.im.plt'
         else
  	  ofilename(1) = 'psi.'//trim(channel_index)//'.re.alpha.plt'
  	  ofilename(2) = 'psi.'//trim(channel_index)//'.im.alpha.plt'
	  ofilename(3) = 'psi.'//trim(channel_index)//'.re.beta.plt'
  	  ofilename(4) = 'psi.'//trim(channel_index)//'.im.beta.plt'
         end if
      
         do n = 1, nfiles
          open(unit=ofile(n),
     &	       file=trim(ofilename(n)),access='direct',recl=4,iostat=ierr)	
	  if (ierr.ne.0) then
           stop '[compute_scatt_wave_func] : can not open output .plt-file!'
          end if

c         reference: calcgrd.f from turbomole moloch-library 
c         !!! caution: use ifort with "-assume byterecl" option 
c             to have record length in bytes	
c         header of the plt.file >>>
          write(ofile(n),rec=1) plt1
          write(ofile(n),rec=2) plt2

          write(ofile(n),rec=3) zsegments+1
          write(ofile(n),rec=4) ysegments+1
          write(ofile(n),rec=5) xsegments+1

          write(ofile(n),rec=6)  real(wf_zmin*bohr_radius,4)
          write(ofile(n),rec=7)  real(wf_zmax*bohr_radius,4)
          write(ofile(n),rec=8)  real(ymin*bohr_radius,4)
          write(ofile(n),rec=9)  real(ymax*bohr_radius,4)
          write(ofile(n),rec=10) real(xmin*bohr_radius,4)
          write(ofile(n),rec=11) real(xmax*bohr_radius,4)

c         << header of the plt.file 
	 end do ! files      
c        set a current record number : 
         nrecord = 11
c        <<< done with file opening 
      
c        make cicles over grid points >>>
         npoint = 0
         do nz = 0, zsegments 
	  print '(2x,a,i5,a,$)', '### z-layer : ', nz+1, ' : | '          
	  rpoint(3) = wf_zmin + nz * dz
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
              rpsi = compute_rpsi(rpoint,utt(:,nsaos+1-ichannel,1))
              if (testing) then
	       write (80+ichannel,'(a,f10.6,2x,f10.6,2x,f10.6,a,e16.10,4x,e16.10)')
     &          'r = ', rpoint(1), rpoint(2), rpoint(3), 
     &          ' :   Psi(r,E) = ', dble(rpsi * wf_scaling_fctr), dimag(rpsi * wf_scaling_fctr) 
              end if 
	      write(ofile(1),rec=nrecord) real( dble(rpsi * wf_scaling_fctr),4)
	      write(ofile(2),rec=nrecord) real(dimag(rpsi * wf_scaling_fctr),4)

             else                              
c             open shell case >>>
c             --- alpha channel ---
              rpsi = compute_rpsi(rpoint,utt(:,nsaos+1-ichannel,1))
              if (testing) then
	       write (80+ichannel,'(a,f10.6,2x,f10.6,2x,f10.6,a,e16.10,4x,e16.10)')
     &          'r = ', rpoint(1), rpoint(2), rpoint(3), 
     &          ' :   Psi.alpha(r,E) = ', dble(rpsi * wf_scaling_fctr), dimag(rpsi * wf_scaling_fctr) 
              end if 
	      write(ofile(1),rec=nrecord) real( dble(rpsi * wf_scaling_fctr),4)
	      write(ofile(2),rec=nrecord) real(dimag(rpsi * wf_scaling_fctr),4)

              rpsi = compute_rpsi(rpoint,utt(:,nsaos+1-ichannel,2))
              if (testing) then
	       write (90+ichannel,'(a,f10.6,2x,f10.6,2x,f10.6,a,e16.10,4x,e16.10)')
     &          'r = ', rpoint(1), rpoint(2), rpoint(3), 
     &          ' :   Psi.beta(r,E) = ', dble(rpsi * wf_scaling_fctr), dimag(rpsi * wf_scaling_fctr) 
              end if 
	      write(ofile(3),rec=nrecord) real( dble(rpsi * wf_scaling_fctr),4)
	      write(ofile(4),rec=nrecord) real(dimag(rpsi * wf_scaling_fctr),4)

             end if  ! closed vc open shell case

	   end do ! loop over x points
	  end do ! loop over y points
	  print '(a)', ' | done'
	 end do ! loop over z points

         do n = 1, nfiles
	  close(ofile(n))
	 end do

        end do ! ichannel

        if (nspin.eq.2) then 
	  print '(/,2x,a)', 
     &	   'scattering wave functions for alpha-spin electrons are exported to files <psi.nn.re/im.alpha.plt>'
	  print '(2x,a)', 
     &	   'scattering wave functions for beta-spin  electrons are exported to files <psi.nn.re/im.beta.plt>'
        else 
	  print '(/,2x,a)', 
     &	   'scattering wave functions are exported to files <psi.nn.re/im.plt>'
	end if 

        print '(/,a)', ' <<< DONE WITH SCATTERING WAVE FUNCTIONS '
       
       end subroutine compute_scatt_wave_func 

       complex(8) function compute_rpsi(r_point,eigenvector)
        double precision r_point(3)
        complex(8)   eigenvector(:) 

        double precision, allocatable :: re_vec(:), im_vec(:), s12vec(:)
        double precision  re_psi, im_psi, r_nu(3), inorm
        integer ierr, n, iat, iorb, iatype, norb

        type(cgto) nu
        double precision, allocatable :: phi_r(:)

c       BLAS level-one function to compute a dot product.
c       of the two double precision vectors
        double precision, external :: ddot

        allocate(s12vec(nsaos),stat=ierr)
        if (ierr.ne.0) then
          print *
          stop '[FUNCTION compute_rpsi]: <s12vec> allocation failure'
        end if

        allocate(re_vec(nsaos),stat=ierr)
        if (ierr.ne.0) then
          print *
          stop '[FUNCTION compute_rpsi]: <re_vec> allocation failure'
        end if
        allocate(im_vec(nsaos),stat=ierr)
        if (ierr.ne.0) then
          print *
          stop '[FUNCTION compute_rpsi]: <im_vec> allocation failure'
        end if

        forall(n=1:nsaos)
          re_vec(n) = dble(eigenvector(n))
          im_vec(n) = dimag(eigenvector(n))
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

c       compute  s12inv * re_vec  --> s12vec
        s12vec = 0.0d0
        call dgemv('n',nsaos,nsaos,1.0d0,smat12inv,nsaos,re_vec,1,0.0d0,s12vec,1)
c       compute   phi_r * s12vec  --> re_psi
        re_psi = ddot(nsaos,phi_r,1,s12vec,1)

c       compute  s12inv * im_vec  --> s12vec
        s12vec = 0.0d0
        call dgemv('n',nsaos,nsaos,1.0d0,smat12inv,nsaos,im_vec,1,0.0d0,s12vec,1)
c       compute   phi_r * s12vec  --> im_psi
        im_psi = ddot(nsaos,phi_r,1,s12vec,1)
        
        deallocate(s12vec,re_vec,im_vec,phi_r,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [FUNCTION compute_rpsi]: ',
     &           'impossible to deallocate temp. arrays'
         print *, 'nevertheless, proceed further ...'
        end if

c       re_psi =  1.0d-2
c       im_psi = -1.0d-2
        compute_rpsi = cmplx(re_psi,im_psi,8) 

       end function compute_rpsi

      end module cond_chann_wf
      
