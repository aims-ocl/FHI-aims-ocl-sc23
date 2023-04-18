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

      module selfenergy
c     *****************************************************
c     computes the simplest possible self-energy
c     satisfying absorbing boundary conditions,
c     updates a hamiltonian: h = h0 + sigma
c     
c     refs: 
c     (a) Evers, F.; Arnold, A., arXiv:cond-mat/0611401v1
c
c     (b) Arnold, A., Weigend, F. & Evers, F. 
c         Quantum chemistry calculations for molecules 
c         coupled to reservoirs: Formalism, implementation, 
c         and application to benzenedithiol.
c         J. Chem. Phys. 126, 174101 (2007).
c     *****************************************************

       use globalvars
       use read_externals
       use tools
       implicit none    

       integer, private, parameter :: linelength = 120
       integer, private, parameter :: wordlength = 32
       integer, private, parameter :: wmark = 5
      
      contains

      subroutine rotate(fixpoints,inatom,outatom)
       implicit none

       integer, intent (in) :: fixpoints(3)
       type(type_atom), intent (in)  :: inatom(:)      
       type(type_atom), intent (out) :: outatom(:)
       
       double precision, parameter :: crdtiny = 1.0d-12
       
       double precision, allocatable :: crd0(:,:), 
     &                                  crd1(:,:)
       double precision catom(3), xatom(3), yatom(3),
     &                  v1(3), v2(3), vnormal(3), 
     &                  r, r2, rho, rho2, cost, sint  
       integer iat, i, ierr, jerr
       logical flag

       allocate(crd0(num_atoms,3),stat=ierr)
       allocate(crd1(num_atoms,3),stat=jerr)
       flag = (ierr.ne.0).or.(jerr.ne.0)       
       if (flag) then
        stop ' [SUBROUTINE rotate]: <crd> allocation failure'
       end if
       
c      first, we 'copy-paste' untouched information
       do iat = 1, num_atoms
        outatom(iat)%symbol = inatom(iat)%symbol 
	outatom(iat)%atype  = inatom(iat)%atype
	outatom(iat)%llead  = inatom(iat)%llead
	outatom(iat)%rlead  = inatom(iat)%rlead
       end do

c      define a normal to the surface plane
       iat = fixpoints(1) 
       catom = inatom(iat)%pos

       iat = fixpoints(2) 
       xatom = inatom(iat)%pos

       iat = fixpoints(3) 
       yatom = inatom(iat)%pos
       
       v1 = xatom - catom
       v2 = yatom - catom
       call crossprod(v1,v2,vnormal)

c      fist checking whether we've got correct normal
       r2 = vnormal(1)*vnormal(1) 
     &    + vnormal(2)*vnormal(2)
     &    + vnormal(3)*vnormal(3)
       if (r2 < crdtiny) then
        print '(/,/,2x,a)', 
     &      '[SUBROUTINE rotate]: normal vector is close to zero!'
        print '(2x,a,/)', 
     &	    '-- something seems to be wrong with electrode set-up'
        stop ': tranport module is terminated now'
       end if

c      move the structure to 'catom'
       do i = 1, 3  ! cycle over x,y,z  
        do iat = 1, num_atoms  
         crd0(iat,i) = inatom(iat)%pos(i) - catom(i)
	end do
       end do
   
       r2 = vnormal(1)*vnormal(1) + vnormal(2)*vnormal(2)
       r  = sqrt(r2)
       if (r < crdtiny) then
c       structure is kept as it is 
        do iat = 1, num_atoms
         outatom(iat)%pos(1) =  crd0(iat,1)          
	 outatom(iat)%pos(2) =  crd0(iat,2)
         outatom(iat)%pos(3) =  crd0(iat,3)
        end do

       else 
c       then we need to rotate the structure
c       performing the 1st rotation ...
        cost = vnormal(1)/r
        sint = vnormal(2)/r
        do iat = 1, num_atoms
         crd1(iat,1) =  cost * crd0(iat,1) + sint * crd0(iat,2)
         crd1(iat,2) = -sint * crd0(iat,1) + cost * crd0(iat,2)
         crd1(iat,3) =  crd0(iat,3)
        end do
c       performing the 2nd rotation ...
        rho2 = vnormal(1)*vnormal(1) + vnormal(2)*vnormal(2)
        rho = sqrt(rho2)
        r2  = rho2 + vnormal(3)*vnormal(3)  
        r   = sqrt(r2)
        cost = vnormal(3)/r
        sint = rho/r
        do iat = 1, num_atoms
         outatom(iat)%pos(3) =  cost*crd1(iat,3) + sint*crd1(iat,1)
         outatom(iat)%pos(1) = -sint*crd1(iat,3) + cost*crd1(iat,1)
         outatom(iat)%pos(2) =  crd1(iat,2)
        end do
       end if

       deallocate(crd0,crd1,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE rotate]: ',
     &        'impossible to deallocate <crd> arrays'
        print *, ' nevertheless, proceed further ...'  	
       end if
     
      end subroutine rotate      

      subroutine buildsigma(dz)
       implicit none
c      input: tolerance parameter to identify group
c             of atoms belonging to the same atomic plane
       double precision, intent(in) :: dz

       type aztype
        integer           type  ! index
	double precision  z     ! z-component of coords 
       end type 
       
       type(type_atom),  allocatable :: tmpatom(:)
       type(aztype),     allocatable :: ztype(:)
       double precision, allocatable :: z0comp(:), 
     &                                  zcomp(:), zsort(:)
       integer, allocatable          :: zindex(:)
    
       double precision  z0, ztmp, hex
       complex(8) :: snn, bsnn       
       integer i, iat, iatype, ilead, j, jmax, ierr, tmpfile, 
     &         nu, iorb, norb, itype, n1, n2, n3, layer_index_tmp
       logical atfound, zfound, flag, outinfo
       character(atsymbol_len) satname, tmpatname
    
c      for testing ...
       character(11), parameter :: lcrdfile = 'coord.l.tmp'
       character(11), parameter :: rcrdfile = 'coord.r.tmp'
       character(11)  crdfile
             
        
       allocate(tmpatom(1:num_atoms),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUBROUTINE buildsigma]: <tmpatom> allocation failure'
       end if

       allocate(z0comp(1:num_atoms),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUBROUTINE buildsigma]: <z0comp> allocation failure'
       end if

       allocate(ztype(num_atoms),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE buildsigma]: <ztype> allocation failure'
       end if

       if (.not.allocated(sigmann)) then
        if (.not.ecp) then
         allocate(sigmann(nsaos),stat=ierr)
        else
         allocate(sigmann(nvalence),stat=ierr)
        end if 
        if (ierr.ne.0) then
          print *
          stop '[SUBROUTINE buildsigma]: <sigmann> allocation failure'
        end if
       end if	

       if (.not.allocated(at_self_energy)) then
        allocate(at_self_energy(num_atoms),stat=ierr)
        if (ierr.ne.0) then
          print *
          stop '[SUBROUTINE buildsigma]: <sigmann> allocation failure'
        end if
       end if	

       if ( sp_leads.and.(.not.allocated(bsigmann)) ) then
        if (.not.ecp) then
          allocate(bsigmann(nsaos),stat=ierr)
        else
          allocate(bsigmann(nvalence),stat=ierr)
        end if   
        if (ierr.ne.0) then
          print *
          stop '[SUBROUTINE buildsigma]: <bsigmann> allocation failure'
        end if
       end if	

       if (.not.ecp) then
        allocate(gamma_left(nsaos),stat=ierr)
       else
        allocate(gamma_left(nvalence),stat=ierr)
       end if    
       if (ierr.ne.0) then
         print *
         stop 
     & 	 '[SUBROUTINE buildsigma]: <gamma_left> allocation failure'
       end if

       if (.not.ecp) then
        allocate(gamma_right(nsaos),stat=ierr)
       else
        allocate(gamma_right(nvalence),stat=ierr)
       end if    
       if (ierr.ne.0) then
         print *
         stop 
     & 	 '[SUBROUTINE buildsigma]: <gamma_right> allocation failure'
       end if

       at_self_energy = 0.0d0
       sigmann = czero
       if (sp_leads)  bsigmann = czero 
       gamma_left  = 0.0d0
       gamma_right = 0.0d0
       print '(/,a)', ' BUILD UP A SELF-ENERGY >>>'
c      ! cycle over left & right electrodes      
       do ilead = 1, 2  

	if (ilead == 1) then
          print '(/,a,$)', '  left electrode,  surface atoms: '
	  call rotate(left_satoms,atom,tmpatom)
	else 
          print '(/,a,$)', '  right electrode, surface atoms: '
          call rotate(right_satoms,atom,tmpatom)
        end if

ccc     just to be on the safe side
        z0comp = 0.0d0

        if (ilead == 1) then
	 iat = left_satoms(1)
	else
	 iat = right_satoms(1)
	end if
	
	satname = atom(iat)%symbol 
	print '(a2)', satname(1:2)
	z0comp(1) = abs(tmpatom(iat)%pos(3))

        jmax = 1
        do iat = 1, num_atoms
        
	 z0 = tmpatom(iat)%pos(3)
	 z0 = abs(z0)
         atfound = .false.
         j = 0
         do while ((.not.atfound).and.(j<jmax))
          j = j + 1
          ztmp = z0comp(j)
          if ( abs(ztmp-z0)<=dz .or. abs(ztmp+z0)<=dz ) then
	   atfound = .true.
	  end if
         end do

         if (.not.atfound) then
          jmax = jmax + 1
          z0comp(jmax) = z0
         end if
        
	end do ! iat 

        allocate(zcomp(jmax),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE buildsigma]: <zcomp> allocation failure'
        end if
        allocate(zsort(jmax),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE buildsigma]: <zcomp> allocation failure'
        end if
        allocate(zindex(jmax),stat=ierr)
        if (ierr.ne.0) then
          print *
          stop '[SUBROUTINE buildsigma]: <zindex> allocation failure'
        end if

c       save a copy ...
        forall (j=1:jmax) zcomp(j) = z0comp(j)
c       ... and sorting according to z-component
        call sortarray(zcomp, zsort, jmax, zindex)

c       fill up a structure <ztype>
        do iat = 1, num_atoms
	 z0 = tmpatom(iat)%pos(3)  
	 ztype(iat)%z   = z0
	 zfound = .false.
	 j = 0 
	 do while ((.not.zfound).and.(j<jmax))
	  j = j + 1 
	  ztmp = zsort(j)    
          if ( abs(ztmp-z0)<=dz .or. abs(ztmp+z0)<=dz ) then
	   zfound = .true.
	  end if
	 end do  
	 if (.not.zfound) then
	  print '(/,a)', 
     &	    ' [SUBROUTINE selfenergy]: something goes wrong here ...'
          stop ' : i am going to quit now'
	 end if
	 ztype(iat)%type = j
	end do

        deallocate(zcomp,zsort,zindex,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE buildsigma]: ',
     &           'impossible to deallocate temp. arrays'
         print *, 'nevertheless, proceed further ...'
        end if

c       testing ...
c       print out transformed coordinates and z-types
        if (testing) then

	 if (ilead == 1) then ;  crdfile = lcrdfile
         else ; crdfile = rcrdfile
	 end if
         tmpfile = 23
	 open(tmpfile,file=crdfile,status='unknown',iostat=ierr)
         if (ierr.ne.0) then
          stop ': can not open tmp-coord file!'
         end if

         write(tmpfile,'(a)') '$coord.tmp'
         if (ilead == 1) then
	  write(tmpfile,'(a)') '#checking left electrode'
	 else
	  write(tmpfile,'(a)') '#checking right electrode'
	 end if
	 write(tmpfile,'(a,i4)') '#number of atoms :', num_atoms
         write(tmpfile,'(a,i4)') '#atom z-types    :', jmax

         do iat = 1, num_atoms
c         for testing, keep all atomic positions in a.u.!
c         if (aims_input) then
c            write(tmpfile,'(3f22.14,6x,a2,4x,a,i4)')
c     &             (tmpatom(iat)%pos(i) * bohr_radius,i=1,3),
c     &             tmpatom(iat)%symbol(1:2), 'ztype =', ztype(iat)%type
c          else
           write(tmpfile,'(3f22.14,6x,a2,4x,a,i4)')
     &             (tmpatom(iat)%pos(i),i=1,3),
     &             tmpatom(iat)%symbol(1:2), 'ztype =', ztype(iat)%type
c	   end if ! aims_input
         end do
         write(tmpfile,'(a)') '$end'
         close(tmpfile)

        end if ! testing

c       build up self-energy
        n1 = 0 ; n2 = 0 ; n3 = 0;  
        nu = 0
	do iat = 1, num_atoms

c        checking whether given atom belongs to coupling region
c        if yes, ascribe self-energy contribution to it
         snn  = czero 
         bsnn = czero
c	 if (sp_leads) bsnn = czero
         ztmp = ztype(iat)%z
         layer_index_tmp = ztype(iat)%type
         tmpatname = atom(iat)%symbol 
	   
c        >>> modified on jan 2012: 
c            allow to use either of keywords, $dist or $nlayers     
         flag = .false.
	 if (dist_found) then
          flag = (abs(ztmp)<=dist .and. satname(1:2)==tmpatname(1:2) )	 
	 else if (nlayers_found) then
	  flag = (layer_index_tmp<=nlayers .and. satname(1:2)==tmpatname(1:2) )
	 else
	  stop ' [SUBROUTINE buildsigma]: neither "dist" nor "nlayers" found ??? '
	 end if

         if (flag) then
c         we cought atom from coupling region 
c         if (ilead == 1) then ; atom(iat)%llead = .true.
c	                  else ; atom(iat)%rlead = .true.
c	  end if		  

c         we cought atom from coupling region 
          if (ilead == 1) then
	    if (atom(iat)%rlead) then 
              print * 
              print *, ' your LEFT and RIGHT electrodes share the same atoms: '
	      print *, ' are you shure you are doing it right?  '    
	      print *, ' please, check your <tcontrol> file ...'
              print * 
	      stop ': transport module will terminate now!'
 	    else 
              atom(iat)%llead = .true. 	    
	    end if
	  else
	    if (atom(iat)%llead) then 
              print * 
              print *, ' your LEFT and RIGHT electrodes share the same atoms: '
	      print *, ' are you shure you are doing it right?  '    
	      print *, ' please, check your <tcontrol> file ...'
              print * 
	      stop ': transport module will terminate now!'
 	    else 
              atom(iat)%rlead = .true. 	    
	    end if
	  end if

	  itype = ztype(iat)%type     

          if (.not.sp_leads) then
c           non-spin-polarized case          
	    select case (itype)
              case(1)      ;  snn = sigma(1) ; n1 = n1 + 1
	      case(2)      ;  snn = sigma(2) ; n2 = n2 + 1 
	      case default ;  snn = sigma(3) ; n3 = n3 + 1 
	    end select
          else         
c           spin-polarized case  
            if (ilead == 1) then 
	      hex = 0.5d0*hex_left
	    else 
	      hex = 0.5d0*hex_right            
	    end if
	    select case (itype)
              case(1)      
	        snn=sigma(1)-hex ; bsnn=sigma(1)+hex ; n1=n1+1
	      case(2)      
	        snn=sigma(2)-hex ; bsnn=sigma(2)+hex ; n2=n2+1 
	      case default  
	        snn=sigma(3)-hex ; bsnn=sigma(3)+hex ; n3=n3+1 
	    end select
          end if  ! sp_leads
	 
	 end if  ! atom catch
c	 otherwise, contributions snn/bsnn remain 0  

c        add self-energy contributions	 
         at_self_energy(iat) = at_self_energy(iat) - dimag(snn) 
	 iatype = atom(iat)%atype
         if (.not.ecp) then
	  norb = n_basis_func(iatype)
	 else
	  norb = n_val_basis_func(iatype)
	 end if 
	 do iorb = 1, norb
          nu = nu + 1
c         important: here we remember about two electrodes	  
c         that's why we add correction snn
          sigmann(nu) = sigmann(nu) + snn
	  if (sp_leads) bsigmann(nu) = bsigmann(nu) + bsnn
          if (ilead == 1) then ! left electrode
	   gamma_left(nu)  = -2.0d0 * dimag(snn)
	  else                 ! right electrode 
	   gamma_right(nu) = -2.0d0 * dimag(snn)
          end if	  
	 end do

        end do  ! iat
c       print *, ' check: nu = ', nu, '   nsaos = ', nsaos

	print '(a,i3)', 
     & 	  '  number of atoms coupled to a reservoir:', n1+n2+n3
	print '(a,i3)', '  -- layer 1:', n1
	print '(a,i3)', '  -- layer 2:', n2
	print '(a,i3)', '  -- layer 3:', n3

       end do ! ilead = left, right

c      in case of $landauer flag, output 
c      self energy to external file and update a <tcontrol> file
c      adding a line $self_energy file=<...filename...>
       if (landauer) then
        call output_self_energy(trim(self_energy_file), outinfo) 
        if (outinfo) then 
         print '(/,a,a,a)',
     &	 ' information on the self-energy has been written to a file "',trim(self_energy_file),'"'       
         call updatetcontrol(trim(tcntrl_file_name),'$self_energy',0,0.0d0)
        end if
       end if	
    
       print '(/,a)', ' <<< DONE WITH THE SELF-ENERGY'

       deallocate(tmpatom,z0comp,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE buildsigma]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
       end if

      end subroutine buildsigma

      subroutine output_self_energy(se_file_name,outinfo)
c      ************************************
c      outputs imaginary piece of the model 
c      self-energy matrix to external file
c      ************************************
       implicit none
       character(*), intent (in) :: se_file_name
       logical, intent (out) :: outinfo

       integer     :: sefile = 33, ierr, iat, i
       character(5)   atommark

       open(sefile,file=se_file_name,status='unknown',iostat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a,a,/,a)', 
     &    ' WARNING: was not able to open file "',trim(se_file_name),'" to save a model ',
     &    '          self-energy; will proceed working further anyway ! '
         outinfo = .false.
         close(sefile)
       else
c        do the rest         
 	 outinfo = .true.
         write(sefile,'(a)') '$self.energy: imaginary piece per atom'
         do iat = 1, num_atoms 
	   if (atom(iat)%llead) then 
	    atommark = 'left '
	   else if (atom(iat)%rlead) then
	    atommark = 'right'  
           else
	    atommark = 'empty'         	   
	   end if    
           if (aims_input) then
            write(sefile,'(i5,1x,3f15.7,4x,a2,4x,a5,4x,d20.14)')
     &              iat, (atom(iat)%pos(i) * bohr_radius,i=1,3),
     &              atom(iat)%symbol(1:2), atommark, at_self_energy(iat)
           else
            write(sefile,'(i5,1x,3f15.7,4x,a2,4x,a5,4x,d20.14)')
     &              iat, (atom(iat)%pos(i),i=1,3),
     &              atom(iat)%symbol(1:2), atommark, at_self_energy(iat)
	   end if
 	 end do 
         write(sefile,'(a)') '$end'
	 close(sefile)
       end if

      end subroutine output_self_energy
  
c     ********** reading self-energy-file ************* 
      subroutine get_self_energy_data(st,iat,r,atname,atommark,imse)
c       ***************************************************************
c       given a line 'st' from coord-file,
c       extract a position 'r' of atom, and its full name 'atname', 
c       and value '|Im_Sigma| [Hartree]' of the local self-energy 'imse'
c       ***************************************************************  
        character(*), intent(in)              :: st
        integer, intent (out)                 :: iat
        double precision, intent (out)        :: r(3)
        character(atsymbol_len), intent (out) :: atname
        character(wmark), intent (out)        :: atommark
        double precision, intent (out)        :: imse

        integer i, nw
        character(wordlength)   words(10)
        character(2*wordlength) tmp_name, tmp_atommark

        call split_line(st,words,nw)
        if (nw.ne.7) then
c        a non-consistent format of hsource-file
         print *
         stop '[SUB. get_self_energy_data]: wrong format of the self-energy file'
        else 
         read(st,*) iat, (r(i),i=1,3)
         tmp_name = words(5)
         atname = tmp_name(1:atsymbol_len)
         read(words(6),*) tmp_atommark     
	 atommark = tmp_atommark(1:wmark)
         read(words(7),*) imse     
         write(*,fmt='(i5,2x,3f15.7,4x,a2,4x,a5,4x,d20.14)') 
     &                 iat, (r(i),i=1,3),atname,atommark,imse
	end if
	
       end subroutine get_self_energy_data

       subroutine read_self_energy(sefile,outinfo)
        implicit none
        character(*), intent (in) :: sefile 
        logical, intent (out)     :: outinfo

        integer  ifile, ierr, iat, tmpiat
        logical  start_reading, flag        
	double precision  rxyz(3), extse
	character(atsymbol_len) atom_name 
        character             ch
	character(linelength) st
        character(wmark)      atommark

        logical, allocatable :: atominfo(:)

        allocate(atominfo(num_atoms),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE read_self_energy]: <atominfo> allocation failure'
        end if
        atominfo = .false.

        print '(/,a,a,a,/)', ' === reading file <', sefile, '> ==='

        ifile = 25  ! external self-energy file
        open(ifile,file=sefile,status='old',action='read',iostat=ierr)
        if (ierr.ne.0) then
         print '(a,a,a/,a,/,a)', 
     & 	        ' WARNING: i was not able to open file "',trim(sefile),'" for reading : ',
     &	        '          a self-energy matrix could not be imported; instead, it will be',
     &	        '          constructed based on parameters found in your <tcontrol> file'
         outinfo = .false.
c	 print '(/,a,a,a)', 
c     &	           ' can not open file "', trim(sefile), '" for reading'
c         print '(a,/)',
c     &    ' please, check your directory content or access rights'
c        stop ': transport module is terminated now'
	else
         outinfo = .true.
         start_reading = .false.
         do while (.not.start_reading)
c         read(ifile,'(a1)') ch   ! read first symbol
          read(ifile,*) ch   ! read first symbol
          flag = (ch=='$').or.(ch=='#')
          if (.not.flag)  start_reading = .true.
         end do
         backspace(ifile)

         do tmpiat = 1, num_atoms
          read(ifile,fmt='(a)',iostat=ierr) st
          if (ierr.eq.-1) then
            print '(/,1x,a)',     'ERROR: end of file is reached unexpectedly'
            print '(1x,a,a,a,/)', '       your "',trim(sefile),'" file seems to be corrupted'
            stop ' : transport module is terminated now'
          else if (ierr.gt.0) then
            print '(/,1x,a)',     'reading of file resulted in error'
            print '(1x,a,a,a,/)', 'your "',trim(sefile),'" file seems to be corrupted'
            stop ' : transport module is terminated now'
          end if
          call get_self_energy_data(st,iat,rxyz,atom_name,atommark,extse)
c         print *, atommark, extse
          atom(iat)%imse = dabs(extse)
	  select case (trim(atommark))
	   case('left')  
	    atom(iat)%llead = .true.
	    atom(iat)%rlead = .false.
	   case('right')  
	    atom(iat)%rlead = .true.
	    atom(iat)%llead = .false.
	   case('empty')  
	    atom(iat)%rlead = .false.
	    atom(iat)%llead = .false.
           case default 
	    atom(iat)%rlead = .false.
	    atom(iat)%llead = .false.
          end select
c         print *, atom(iat)%llead, atom(iat)%rlead, atom(iat)%imse 
          atominfo(iat) = .true.
	 end do ! tmpiat

c        check whether all diagonal matrix elements have been read ...
         do iat = 1, num_atoms
	  if (.not.atominfo(iat)) then
	   print '(/,1x,a)',   'WARNING: not all matrix elements have been read from your'
	   print '(1x,a,a,a)', '         <',trim(sefile),'> file: it seems to be corrupted ;'
	   print '(1x,a)',     '         a self-energy will be constructed based on parameters'
	   print '(1x,a)',     '         found in your <tcontrol> file'
           outinfo = .false.
	  end if
  	 end do ! iat

	end if  ! succesfull file reading

        close(ifile)

        deallocate(atominfo,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE read_self_energy]: ',
     &         'impossible to deallocate <atominfo> array'
         print *, ' nevertheless, proceed further ...'  	
        end if

       end subroutine read_self_energy

      subroutine fill_sigma_matrix
c      *************************************************
c      importing a self-energy matrix from external file      
c      *************************************************
       implicit none

       integer    iat, iatype, iorb, norb, nu, ierr
       complex(8) snn

       if (.not.allocated(sigmann)) then
        if (.not.ecp) then
         allocate(sigmann(nsaos),stat=ierr)
        else
         allocate(sigmann(nvalence),stat=ierr)
        end if 
        if (ierr.ne.0) then
          print *
          stop '[SUBROUTINE import_sigma]: <sigmann> allocation failure'
        end if
       end if    

       if (.not.allocated(gamma_left)) then
        if (.not.ecp) then
         allocate(gamma_left(nsaos),stat=ierr)
        else
         allocate(gamma_left(nvalence),stat=ierr)
        end if 
        if (ierr.ne.0) then
         print *
         stop 
     & 	 '[SUBROUTINE buildsigma]: <gamma_left> allocation failure'
        end if
       end if    

       if (.not.allocated(gamma_right)) then
        if (.not.ecp) then
         allocate(gamma_right(nsaos),stat=ierr)
        else
         allocate(gamma_right(nvalence),stat=ierr)
        end if 
        if (ierr.ne.0) then
         print *
         stop 
     & 	 '[SUBROUTINE buildsigma]: <gamma_right> allocation failure'
        end if
       end if    
 
       sigmann     = czero  
       gamma_left  = 0.0d0
       gamma_right = 0.0d0

c      initialize self-energy & gamma matrices >>>
       nu = 0
       do iat = 1, num_atoms
         snn = -ione * dabs(atom(iat)%imse)
	 iatype = atom(iat)%atype
         if (.not.ecp) then
	  norb = n_basis_func(iatype)
	 else
	  norb = n_val_basis_func(iatype)
	 end if 
	 do iorb = 1, norb
          nu = nu + 1
c         important: here we remember about two electrodes
c         that's why we add a correction snn
          sigmann(nu) = sigmann(nu) + snn
          if ( atom(iat)%llead ) then       !  left electrode
	   gamma_left(nu)  = 2.0d0 * dabs(atom(iat)%imse)
	  else if ( atom(iat)%rlead ) then  !  right electrode 
	   gamma_right(nu) = 2.0d0 * dabs(atom(iat)%imse)
          end if	  
	 end do  ! iorb
       end do ! iat

      end subroutine fill_sigma_matrix   

      end module selfenergy    
