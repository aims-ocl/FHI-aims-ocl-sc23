c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           march/april 2008
c     last revision:  jan 2012
c###########################################################

      module make_overlap
c      ***********************************
c      builds up overlap matrix: S=<mu|nu>
c      ***********************************
       use globalvars
       use math_functions
       implicit none

      contains

      subroutine init_overlap
       implicit none

       type(cgto) mu, nu   
       integer    iat,jat, iatype, jatype, 
     &            norb, morb, iorb, jorb, n,m, ierr
       double precision r_nu(3), r_mu(3)           
       logical, allocatable :: initdone(:,:)
       
       integer, parameter :: monitor = 25

c      external function
c       interface
c        double precision function overlap(mu,nu,r_mu,r_nu)
c         type(cgto)       mu,nu            ! input: atomic orbitals       
c         double precision r_mu(3), r_nu(3) !   & their localization   
c        end function overlap
c       end interface	 
       
       allocate(smat(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        stop '[SUBROUTINE init_overlap]: <smat> allocation failure'
       end if

       allocate(orbnorm(nsaos),stat=ierr)
       if (ierr.ne.0) then
        stop '[SUBROUTINE init_overlap]: <orbnorm> allocation failure'
       end if

       allocate(initdone(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        stop '[SUBROUTINE init_overlap]: <initdone> allocation failure'
       end if

       orbnorm = 1.0d0

       print '(/,a)', ' CALCULATING OVERLAP MATRIX ELEMENTS >>>'

       print '(a,$)', ' ' 	    
       initdone = .false.
       smat = 0.0d0
       n=0       
       do iat = 1, num_atoms
c       monitor
        if (mod(iat,2)==0) print '(a,$)', '*' 	    
    	r_nu = atom(iat)%pos  ! save position of atom iat
        iatype = atom(iat)%atype
	norb = n_basis_func(iatype)
	do iorb = 1, norb
         n = n+1 
         m = 0
         do jat = 1, num_atoms
          r_mu = atom(jat)%pos  ! save position of atom jat
          jatype = atom(jat)%atype
	  morb = n_basis_func(jatype)
          do jorb = 1, morb
	   m = m+1
ccc        here is a check of whether smat(m,n) is done already
	   if (.not.initdone(m,n)) then
            nu = aos(iatype,iorb)  ! take required
	    mu = aos(jatype,jorb)  ! orbitals
            smat(n,m) = overlap(nu,mu,r_nu,r_mu)
            initdone(n,m) = .true. 
c	    if (n==m) print *, ' n = ', n, 'smat(n,n) = ', smat(n,n) 
	   else
	    smat(n,m) = smat(m,n)    	  
 	    initdone(n,m) = .true. 
	   end if
c          saving normalization ...
           if (m==n) orbnorm(n) = sqrt(smat(n,n))

	  end do ! jorb
         end do ! jat  
        end do ! iorb
       end do ! iat

ccc    checking ...
       if (n.ne.nsaos) then
        print * 
        stop '[SUBROUTINE init_overlap]: n<>nsaos!'
       end if	

c      further we care about the normalization >>>
       initdone = .false.
       do n = 1, nsaos
        do m = 1, nsaos
	   if (.not.initdone(m,n)) then
	    smat(n,m) = smat(n,m)/(orbnorm(n)*orbnorm(m))
            initdone(n,m) = .true.
	   else 
            smat(n,m) = smat(m,n)    	  
            initdone(n,m) = .true.
           end if
	end do 
       end do 
c      <<< all done
       print '(/,a)', ' <<< DONE WITH THE OVERLAP MATRIX'

       deallocate(initdone,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE init_overlap]: ', 
     &           'impossible to deallocate <initdone>'
        print *, 'nevertheless, proceed further ...' 
       end if

c      if flag 'save_omat' is active, save the overlap matrix to a hard disk
       if (save_omat) call saveomat(outomat_file_name)

      end subroutine init_overlap

      double precision function overlap(mu,nu,r_mu,r_nu)
       implicit none
       
       type(cgto)       mu,nu            ! input: atomic orbitals       
       double precision r_mu(3), r_nu(3) !      & their localization   

       
       type(cart_cgto)  cart_mu(3), cart_nu(3) 

       integer          lm, i,j         
       double precision dr(3), r2, mxi, nxi,
     &                  sqrt12, sqrt24, sqrt40, sqrt60

       integer          ncart_mu, ncart_nu, m_go, n_go, m, n, p_a, p_b
       double precision x_a, x_b, icoeff_m, icoeff_n,
     &                  tmp_overlap, tmp_radial, 
     &                  tmp_sph(3), tmp_spherical 

ccc    external functions used : from module <math_functions>
c       double precision r00   ! normalization
c       double precision p     ! polynoms
       
ccc    transforming spherical harmonics to cartesian    

ccc    difine some constants ...
        sqrt12 = sqrt(12.0d0) ;  sqrt24 = sqrt(24.0d0)
        sqrt40 = sqrt(40.0d0) ;  sqrt60 = sqrt(60.0d0)

ccc    just to be on the safe side
       do i = 1, 3

        cart_mu(i)%icoeff = (/(0.0d0,j=1,contr)/)
	cart_mu(i)%pows   = (/0,0,0/)

        cart_nu(i)%icoeff = (/(0.0d0,j=1,contr)/)
	cart_nu(i)%pows   = (/0,0,0/)

       end do
     
ccc ==============  mu function =======================
ccc    first we deal with polynoms pows
       lm = mu%lm
       select case (lm)
ccc    s-function
        case(1) 
	 cart_mu(1)%pows = (/0,0,0/) ! all pows are zeros    
ccc    p-functions
	case(2)  
         cart_mu(1)%pows = (/1,0,0/) ! px
        case(3)  
	 cart_mu(1)%pows = (/0,1,0/) ! py 
        case(4)  
	 cart_mu(1)%pows = (/0,0,1/) ! pz 
ccc    d-functions
        case(5)  ! (-xx-yy+2zz)/sqrt(12) 
         cart_mu(1)%pows = (/2,0,0/) ! -xx          
	 cart_mu(2)%pows = (/0,2,0/) ! -yy 
         cart_mu(3)%pows = (/0,0,2/) ! +2zz 
        case(6)  ! xz
	 cart_mu(1)%pows = (/1,0,1/) ! xz 
        case(7)  ! yz 
	 cart_mu(1)%pows = (/0,1,1/) ! yz 
        case(8)  ! xy 
	 cart_mu(1)%pows = (/1,1,0/) ! xy
	case(9)  ! (xx-yy)/2 
	 cart_mu(1)%pows = (/2,0,0/) !  xx
	 cart_mu(2)%pows = (/0,2,0/) ! -yy
ccc    f-functions
        case(10) ! (2zzz-3xxz-3yyz)/sqrt(60)
	 cart_mu(1)%pows = (/0,0,3/) !  2zzz
	 cart_mu(2)%pows = (/2,0,1/) ! -3xxz
	 cart_mu(3)%pows = (/0,2,1/) ! -3yyz
        case(11) ! (-xxx-xyy+4xzz)/sqrt(40)
	 cart_mu(1)%pows = (/3,0,0/) !  -xxx
	 cart_mu(2)%pows = (/1,2,0/) !  -xyy
	 cart_mu(3)%pows = (/1,0,2/) !  4xzz
        case(12) ! (-yyy-xxy+4yzz)/sqrt(40)
	 cart_mu(1)%pows = (/0,3,0/) !  -yyy
	 cart_mu(2)%pows = (/2,1,0/) !  -xxy
	 cart_mu(3)%pows = (/0,1,2/) !  4yzz
        case(13) ! xyz
	 cart_mu(1)%pows = (/1,1,1/) !   xyz
       	case(14) ! (xxz-yyz)/2 
	 cart_mu(1)%pows = (/2,0,1/) !   xxz
	 cart_mu(2)%pows = (/0,2,1/) !  -yyz
	case(15) ! (xxx-3xyy)/sqrt(24)
	 cart_mu(1)%pows = (/3,0,0/) !   xxx 
	 cart_mu(2)%pows = (/1,2,0/) ! -3xyy
	case(16) ! (yyy-3xxy)/sqrt(24)
	 cart_mu(1)%pows = (/0,3,0/) !   yyy 
	 cart_mu(2)%pows = (/2,1,0/) ! -3xxy
        case default
	 print *
	 stop '<overlap>: lm>16 ???'
       end select ! case over lm	
	
ccc    then we deal with coefficients and exponents
       if (lm==5) then ! (-xx-yy+2zz)/sqrt(12)
        ncart_mu = 3
        cart_mu(1)%icoeff = -(1.0d0/sqrt12)*mu%icoeff
        cart_mu(2)%icoeff = -(1.0d0/sqrt12)*mu%icoeff
        cart_mu(3)%icoeff =  (2.0d0/sqrt12)*mu%icoeff
       else if (lm==9) then ! (xx-yy)/2
        ncart_mu = 2
        cart_mu(1)%icoeff =  0.5d0*mu%icoeff
        cart_mu(2)%icoeff = -0.5d0*mu%icoeff
       else if (lm==10) then ! (2zzz-3xxz-3yyz)/sqrt(60)
        ncart_mu = 3
        cart_mu(1)%icoeff =  (2.0d0/sqrt60)*mu%icoeff
        cart_mu(2)%icoeff = -(3.0d0/sqrt60)*mu%icoeff
        cart_mu(3)%icoeff = -(3.0d0/sqrt60)*mu%icoeff
       else if (lm==11.or.lm==12) then 
c       ! (-xxx-xyy+4xzz)/sqrt(40); (-yyy-xxy+4yzz)/sqrt(40)  
        ncart_mu = 3
        cart_mu(1)%icoeff = -(1.0d0/sqrt40)*mu%icoeff
        cart_mu(2)%icoeff = -(1.0d0/sqrt40)*mu%icoeff
        cart_mu(3)%icoeff =  (4.0d0/sqrt40)*mu%icoeff
       else if (lm==14) then ! (xxz-yyz)/2 
        ncart_mu = 2
        cart_mu(1)%icoeff =  0.5d0*mu%icoeff
        cart_mu(2)%icoeff = -0.5d0*mu%icoeff
       else if (lm==15.or.lm==16) then
c       ! (xxx-3xyy)/sqrt(24) ; (yyy-3xxy)/sqrt(24)
        ncart_mu = 2
        cart_mu(1)%icoeff =  (1.0d0/sqrt24)*mu%icoeff
        cart_mu(2)%icoeff = -(3.0d0/sqrt24)*mu%icoeff
       else
	ncart_mu = 1
        cart_mu(1)%icoeff = mu%icoeff 
       end if
ccc ============== done with mu function ==============

ccc ==============  nu function =======================
ccc    first we deal with polynoms pows
       lm = nu%lm
       select case (lm)
ccc    s-function
        case(1) 
	 cart_nu(1)%pows = (/0,0,0/) ! all pows are zeros    
ccc    p-functions
	case(2)  
         cart_nu(1)%pows = (/1,0,0/) ! px
        case(3)  
	 cart_nu(1)%pows = (/0,1,0/) ! py 
        case(4)  
	 cart_nu(1)%pows = (/0,0,1/) ! pz 
ccc    d-functions
        case(5)  ! (-xx-yy+2zz)/sqrt(12) 
         cart_nu(1)%pows = (/2,0,0/) ! -xx          
	 cart_nu(2)%pows = (/0,2,0/) ! -yy 
         cart_nu(3)%pows = (/0,0,2/) ! +2zz 
        case(6)  ! xz
	 cart_nu(1)%pows = (/1,0,1/) ! xz 
        case(7)  ! yz 
	 cart_nu(1)%pows = (/0,1,1/) ! yz 
        case(8)  ! xy 
	 cart_nu(1)%pows = (/1,1,0/) ! xy
	case(9)  ! (xx-yy)/2 
	 cart_nu(1)%pows = (/2,0,0/) !  xx
	 cart_nu(2)%pows = (/0,2,0/) ! -yy
ccc    f-functions
        case(10) ! (2zzz-3xxz-3yyz)/sqrt(60)
	 cart_nu(1)%pows = (/0,0,3/) !  2zzz
	 cart_nu(2)%pows = (/2,0,1/) ! -3xxz
	 cart_nu(3)%pows = (/0,2,1/) ! -3yyz
        case(11) ! (-xxx-xyy+4xzz)/sqrt(40)
	 cart_nu(1)%pows = (/3,0,0/) !  -xxx
	 cart_nu(2)%pows = (/1,2,0/) !  -xyy
	 cart_nu(3)%pows = (/1,0,2/) !  4xzz
        case(12) ! (-yyy-xxy+4yzz)/sqrt(40)
	 cart_nu(1)%pows = (/0,3,0/) !  -yyy
	 cart_nu(2)%pows = (/2,1,0/) !  -xxy
	 cart_nu(3)%pows = (/0,1,2/) !  4yzz
        case(13) ! xyz
	 cart_nu(1)%pows = (/1,1,1/) !   xyz
       	case(14) ! (xxz-yyz)/2 
	 cart_nu(1)%pows = (/2,0,1/) !   xxz
	 cart_nu(2)%pows = (/0,2,1/) !  -yyz
	case(15) ! (xxx-3xyy)/sqrt(24)
	 cart_nu(1)%pows = (/3,0,0/) !   xxx 
	 cart_nu(2)%pows = (/1,2,0/) ! -3xyy
	case(16) ! (yyy-3xxy)/sqrt(24)
	 cart_nu(1)%pows = (/0,3,0/) !   yyy 
	 cart_nu(2)%pows = (/2,1,0/) ! -3xxy
        case default
	 print *
	 stop '<overlap>: lm>16 ???'
       end select ! case over lm	
	
ccc    then we deal with coefficients and exponents
       if (lm==5) then ! (-xx-yy+2zz)/sqrt(12)
        ncart_nu = 3
        cart_nu(1)%icoeff = -(1.0d0/sqrt12)*nu%icoeff
        cart_nu(2)%icoeff = -(1.0d0/sqrt12)*nu%icoeff
        cart_nu(3)%icoeff =  (2.0d0/sqrt12)*nu%icoeff
       else if (lm==9) then ! (xx-yy)/2
        ncart_nu = 2
        cart_nu(1)%icoeff =  0.5d0*nu%icoeff
        cart_nu(2)%icoeff = -0.5d0*nu%icoeff
       else if (lm==10) then ! (2zzz-3xxz-3yyz)/sqrt(60)
        ncart_nu = 3
        cart_nu(1)%icoeff =  (2.0d0/sqrt60)*nu%icoeff
        cart_nu(2)%icoeff = -(3.0d0/sqrt60)*nu%icoeff
        cart_nu(3)%icoeff = -(3.0d0/sqrt60)*nu%icoeff
       else if (lm==11.or.lm==12) then 
c       ! (-xxx-xyy+4xzz)/sqrt(40); (-yyy-xxy+4yzz)/sqrt(40)  
        ncart_nu = 3
        cart_nu(1)%icoeff = -(1.0d0/sqrt40)*nu%icoeff
        cart_nu(2)%icoeff = -(1.0d0/sqrt40)*nu%icoeff
        cart_nu(3)%icoeff =  (4.0d0/sqrt40)*nu%icoeff
       else if (lm==14) then ! (xxz-yyz)/2 
        ncart_nu = 2
        cart_nu(1)%icoeff =  0.5d0*nu%icoeff
        cart_nu(2)%icoeff = -0.5d0*nu%icoeff
       else if (lm==15.or.lm==16) then
c       ! (xxx-3xyy)/sqrt(24) ; (yyy-3xxy)/sqrt(24)
        ncart_nu = 2
        cart_nu(1)%icoeff =  (1.0d0/sqrt24)*nu%icoeff
        cart_nu(2)%icoeff = -(3.0d0/sqrt24)*nu%icoeff
       else
	ncart_nu = 1
        cart_nu(1)%icoeff = nu%icoeff 
       end if
ccc ============== done with nu function ==============


ccc making all cycles: summing up overlap integral
       do i = 1, 3
        dr(i) = r_mu(i) - r_nu(i)  
       end do
       r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)        

       tmp_overlap = 0.0d0

       do m_go = 1, mu%ngto
        mxi = mu%xi(m_go)
	do n_go = 1, nu%ngto
         nxi = nu%xi(n_go)
	 tmp_radial = r00(mxi,nxi,r2) 

         tmp_spherical = 0.0d0
         do m=1, ncart_mu
          do n=1, ncart_nu
c           ! expansion coefficients which
c           ! take into account a spherical part 
            icoeff_m = cart_mu(m)%icoeff(m_go)
            icoeff_n = cart_nu(n)%icoeff(n_go) 
c            .... here polynoms weighted with 
c            .... gaussian functions are integrated 
c          cycle over x,y,z ! 
           do i = 1, 3
            x_a = r_mu(i)   ! 1st: mu-orbital   
            x_b = r_nu(i)   ! 2nd: nu-orbital   
            p_a = cart_mu(m)%pows(i)  ! 1st pow 	   
            p_b = cart_nu(n)%pows(i)  ! 2nd pow 	   
            tmp_sph(i) = p(mxi, nxi, x_a, x_b, p_a, p_b)
	   end do
           tmp_spherical = tmp_spherical
     &     + icoeff_m*icoeff_n*tmp_sph(1)*tmp_sph(2)*tmp_sph(3)
          end do ! n
         end do ! m
         tmp_overlap = tmp_overlap + tmp_radial*tmp_spherical

        end do ! n_go
       end do ! m_go

       overlap = tmp_overlap
      end function overlap


      subroutine saveomat(ofile)
c     *************************************      
c     saves overlap matrix to external file
c     *************************************
       character(*), intent (in) :: ofile
       integer, parameter        :: monitor = 100
       integer extfile, n,k,j, mincol, maxcol, ierr, ierr1
       character(128) syscall

       extfile = 29
       open(extfile,file=ofile,status='replace',action='write',iostat=ierr)
c      print *, ' ierr = ', ierr
       if (ierr.ne.0) then
        syscall = 'rm '//ofile ;  call system(trim(syscall))
        open(extfile,file=ofile,status='new',action='write',iostat=ierr1)
c       print *, ' ierr1 = ', ierr1
        if (ierr1.ne.0) then
         print '(/,a,a,a,/)',
     &            ' can not open file "', ofile, '" for writing'
         stop ': encounters a serious problem here!'
        end if
       end if

       write(extfile,fmt='(a,i5)') '$omat   format(4d20.14)   nsaos= ',nsaos
c       write(extfile1,fmt='(a)')
c    &  '#upper triangular of sqrt(overlap):      format(4d22.16))'

       print '(/,2x,a,$)', 'saving overlap matrix : | '
       k = 0
       do n = 1, nsaos
        write(extfile,fmt='(a,i5)') '#column ', n
        k = k + n
        if ( mod(k,monitor)==0 ) print '(a,$)', '='  ! monitor
c        only upper triangular part of matrices is stored
         maxcol = 0
         do while (maxcol<n)
          mincol = maxcol+1
          maxcol = min(maxcol+4,n)
          write(extfile,fmt='(4d20.14)') (smat(j,n),j=mincol,maxcol)
         end do
       end do
       print *, '| done'
       write(extfile,fmt='(a)') '$end'
       close(extfile)

      end subroutine saveomat


      subroutine readomat(omatfile)
c     ***************************************      
c     reads overlap matrix from external file
c     ***************************************      

        character(*), intent (in) :: omatfile
        integer, parameter        :: monitor = 100
        integer extfile, n,k,j, mincol, maxcol, ierr
        character(64) tmpstring

        character(80) headline, cut_headline
	character(16) fmt_omat
	integer ifmt, tmp_nsaos
				
        extfile = 28  ! external omat-file
        open(extfile,file=omatfile,status='old',action='read',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a,a)',
     &            ' can not open file "', trim(omatfile), '" for reading'
         print *,
     &    'please, check your directory content or access rights'
         print *
         stop ': transport module is terminated now'
        end if

        if (.not.allocated(smat)) then
         allocate(smat(nsaos,nsaos),stat=ierr)
         if (ierr.ne.0) then
          stop '[SUBROUTINE readomat]: <smat> allocation failure'
         end if
        end if

        print '(/,a,/)', ' === importing omat file ==='	
c       read header line: get i/o format of the omat-file
        read(extfile,'(a80)') headline
c       searching for substring like: format(4d20.14)
c       other choices of format could be : format(4d24.18), format(4d26.20), etc ...
        ifmt  = index(headline,'format')
        fmt_omat  = headline(ifmt+6:ifmt+14)  ! e.g. (4d20.14)
c       print out info message on found format :
        write(6,'(1x,a)') trim(omatfile)//' : ASCII file, format'//trim(fmt_omat)

c       after format is indentified, proceed with checking nsaos ...
        ifmt  = index(headline,'nsaos=')
        cut_headline  = headline(ifmt+6:len(headline))
        read(cut_headline,*) tmp_nsaos
        if (nsaos.ne.tmp_nsaos) then
         print '(/,1x,a,i6)',  'ERROR: internal conflict: you have specified nsaos= ', nsaos
         print '(1x,a,a,a)',   '       in your <tcontrol> file, while your <',trim(omatfile),'> file'
         print '(1x,a,i6,a)',  '       tells "nsaos" should be ', tmp_nsaos, '  : i am confused'
         print '(/,1x,a,/)',   '       please, check carefully what you are doing ...'
         stop ' : transport module will be terminated now'
        end if

c       after check is done, proceed with reading data ...
        print '(1x,a,$)', 'importing overlap matrix : | '
        k = 0
        do n = 1, nsaos
         read(extfile,fmt=*,iostat=ierr) tmpstring
         if (ierr.eq.-1) then
	    print '(/,/,1x,a)',   'ERROR: end of file is reached unexpectedly'
	    print '(1x,a,a,a,/)', '       your "',trim(omatfile),'" file seems to be corrupted'
	    stop ' : transport module is terminated now'
	 else if (ierr.gt.0) then
	    print '(/,/,1x,a)',   'reading of file resulted in error'
	    print '(1x,a,a,a,/)', 'your "',trim(omatfile),'" file seems to be corrupted'
	    stop ' : transport module is terminated now'
	 end if
         k = k + n
         if ( mod(k,monitor)==0 ) print '(a,$)', '='  ! monitor
c         only upper triangular part is saved to file
          maxcol = 0
          do while (maxcol<n)
           mincol = maxcol+1
           maxcol = min(maxcol+4,n)
           read(extfile,fmt=fmt_omat,iostat=ierr) (smat(j,n),j=mincol,maxcol)
           if (ierr.eq.-1) then
	     print '(/,/,1x,a)',   'ERROR: end of file is reached unexpectedly'
	     print '(1x,a,a,a,/)', '       your "',trim(omatfile),'" file seems to be corrupted'
	     stop ' : transport module is terminated now'
	   else if (ierr.gt.0) then
	     print '(/,/,1x,a)',   'reading of file resulted in error'
	     print '(1x,a,a,a,/)', 'your "',trim(omatfile),'" file seems to be corrupted'
	     stop ' : transport module is terminated now'
	   end if
c          read(extfile,fmt='(4d20.14)') (smat(j,n),j=mincol,maxcol)
c          take care about transposed elements :
           forall(j = mincol:maxcol) smat(n,j) = smat(j,n)
          end do
        end do
        print *, '| done'
        close(extfile)

      end subroutine readomat

      double precision function basis_func_r(mu,fnorm,r_mu,r_point)
c     *****************************************************
c     returns a value of the normalized CGTO basis 
c     function 'mu' (localized at r_mu) at given point r 
c     *****************************************************
       implicit none
       
       type(cgto)       mu          ! input: atomic orbital/basis function       
       double precision fnorm       ! its norm
       double precision r_mu(3)     ! its localization   
       double precision r_point(3)  ! required point in 3D    
       
       double precision r_part, omega_part, x, y, z, r2, tmp_coeff, tmp_xi,
     &                  sqrt12, sqrt24, sqrt40, sqrt60

       integer ncontr, n, lm
       
       x=r_point(1)-r_mu(1) ;  y=r_point(2)-r_mu(2) ;  z=r_point(3)-r_mu(3)
       r2 = x*x + y*y + z*z

c      compute radial part
       ncontr = mu%ngto ;  r_part = 0.0d0
       do n = 1, ncontr
        tmp_coeff = mu%icoeff(n) ; tmp_xi = mu%xi(n)
        r_part = r_part + tmp_coeff * dexp(-tmp_xi*r2)
       end do

c      compute angular dependent part

ccc    difine some constants ...
       sqrt12 = sqrt(12.0d0) ;  sqrt24 = sqrt(24.0d0)
       sqrt40 = sqrt(40.0d0) ;  sqrt60 = sqrt(60.0d0)

       lm = mu%lm
       select case (lm)
ccc    s-function
        case(1) ; omega_part = 1.0d0  ! s 
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

	case(15) ; omega_part=(x*x*x-3.0*x*y*y)/sqrt24
	                                            ! (xxx-3xyy)/sqrt(24)

	case(16) ; omega_part=(y*y*y-3.0d0*x*x*y)/sqrt24
	                                            ! (yyy-3xxy)/sqrt(24)
        case default
	 print *
	 stop '<basis_func_r>: lm > 16 ???'
       end select ! case over lm	

       basis_func_r = omega_part * r_part / fnorm
      end function basis_func_r

      end module make_overlap


