c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           march-may 2008
c     last revision:  march 2012
c###########################################################

      module read_externals
c     *************************************
c     read data from external files:
c       -- tcontrol, coord/geometry.in 
c       -- basis/basis-indices.out 
c       -- mos/alpha/beta
c       -- (optionally): smat12, smat12inv
c     *************************************
       use globalvars
       use periodic_table
       implicit none

       integer, private, parameter :: linelength = 128
       integer, private, parameter :: wordlength = 32
       
       double precision, private, parameter :: def_isigma = 0.1d0
       
       logical, private :: nelfound = .false.
       
       contains

c      ********** reading tcontrol-file ********** 
       subroutine readtcontrol(tcntrl)
        implicit none
        character(*), intent (in) :: tcntrl 
        
        double precision tmpz
        integer          tfile, ierr, tmplen, i, j, 
     &                   tmp_occ, ra, tmp_val_electrons
c  	                 tmplo, nw 
	character(64)    keystring, st1, st2
c       character(linelength) st
c       character(wordlength) words(64)

	logical :: readnext, flag, tmpf
	
	logical :: mos_found    = .false. ,  
     &             alpha_found  = .false. ,
     &             beta_found   = .false. ,
     &             basis_found  = .false. ,
     &             coord_found  = .false. ,
     &             natoms_found = .false. 
c    &             dist_found   = .false. 
c    >>> remark added, jan 2012: dist_found is moved to 
c        globalvars.f : becomes a global variable
c        together with a new variable nlayers_found 

        logical :: iterfound       = .false. ,
     &             biasfound       = .false. ,
     &             dmixfound       = .false. ,
     &             hexfound        = .false. ,
     &             locexfound      = .false. , 
     &             z_range_found   = .false. ,
     &             wf_z_range_found= .false. ,
c   >>> as default, 'neqflag' is modified to .true. from now: jan 2012 
c       if keyword '$non-equilibrium off' is found then neqflag=.false.
     &	           neqflag         = .true.  

c    "hubbard-u" related variables
        logical :: ldau_found        = .false. ,
     &             num_ratoms_found  = .false. ,
     &             ratom_types_found = .false. ,
     &             u_j_found         = .false. 

c    "mol_subsystem" related variables
        logical :: mol_subsystem_found  = .false. ,
     &             num_res_atoms_found  = .false. ,
     &             res_atom_types_found = .false. 

c    energy mesh 
        logical :: ener_found  = .false. ,
     &             estep_found = .false. ,
     &             eend_found  = .false.
	
	logical    found_rsigma(deep), found_isigma(deep)

c    effective core potential flag
        logical :: ecp_found = .false.
	
        tfile = 10  ! tcontrol file
        open(tfile,file=tcntrl,status='old',action='read',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,1x,a,a,a)', 
     &	       'ERROR: can not open file "', trim(tcntrl), '" for reading'
         print '(1x,a,/)',
     &         '       please, check your directory content or access rights'
         stop ': transport module is terminated now'
        end if

c       to be on the safe side
        left_satoms  = 0
	right_satoms = 0
	num_atoms    = 0

c       to be on the safe side, initialize self-energy parameters: 
        do i = 1, deep
         rsigma(i) = 0.0d0
         isigma(i) = def_isigma/dble(2**(i-1))
         found_rsigma(i) = .false.
         found_isigma(i) = .false.
        end do

	print '(/,a,/)', ' === reading <tcontrol> ==='
	readnext = .true.
	do while (readnext)

         read(tfile,*) keystring

	 select case (trim(keystring))

          case ('$aims_input')
           backspace(tfile)
           read(tfile,*) st1, st2
           if (trim(st2).eq.'on') then
            aims_input = .true.
            print *, 'INPUT DATA are provided by FHI-aims!'
c           else
c           ! CAUTION : this line is to protect a default fhi-aims user; 
c           !           flag aims_input should be active anyway!
c            aims_input = .true.
c            print *, 'INPUT DATA are provided by FHI-aims anyway!'
           end if

          case ('$landauer')
           backspace(tfile)
           read(tfile,*) st1, st2
           if (trim(st2).eq.'on') then
            landauer    = .true.
            do_landauer = .true.
            print *, 'Landauer ballistic transmission will be computed'
c           activate "effective core potential" flag: 
c           default user should profit from it anyway !
            ecp = .true.
           else
c           ! CAUTION : this line is to protect a default fhi-aims user ;
c           !           in the first publicly available release, it 
c           !           prevents using any advanced options except of a basic
c           !           option "$landauer on" to compute transmission function
c           !           and (optionally, "$landauer off" and "$ldos on") 
c           !           the local density of states
c            print *
c            print *, '"$landauer" flag in NOT active ... but if later on "$ldos" flag'
c            print *, 'is NOT active as well, Landauer ballistic transmission function'
c            print *, 'will be computed anyway as a default choice'
c            print *
            landauer    = .true.
            do_landauer = .false.
            ecp = .true.
           end if

c         >>> next keyword is analized for compatibility reasons with
c             previous tcontrol-files & versions of the transport code  
c             default value from now (jan 2012) is 'neqflag = .true.'
          case ('$non-equilibrium')
	   backspace(tfile)
	   read(tfile,*) st1, st2
	   if (trim(st2).eq.'off') neqflag = .false.
c         <<< in all other cases, neqflag is kept .true. 

          case ('$no_leads')
	   backspace(tfile)
	   read(tfile,*) st1, st2
	   if (trim(st2).eq.'on') then
	    no_leads = .true.
	    print *, 'NO LEADS flag is active !'
	   else
	    no_leads = .false.
	    print *, 'NO LEADS flag is switched OFF'
           end if	       

          case('$control') ;

          case ('$coord')
           backspace(tfile)
           read(tfile,*) st1, st2
	   tmplen = len_trim(st2) 
           crd_file_name = st2(6:tmplen) 
	   print *, 'coord-file name = ', trim(crd_file_name)
           coord_found = .true.

          case ('$natoms')
           backspace(tfile)
           read(tfile,*) st1, num_atoms
           print '(a,i4)', ' number of atoms = ', num_atoms
           natoms_found = .true.

          case ('$basis')
           backspace(tfile)
           read(tfile,*) st1, st2
	   tmplen = len_trim(st2) 
           basis_file_name = st2(6:tmplen) 
	   print *, 'basis-file name = ', trim(basis_file_name)
           basis_found = .true.

          case ('$ecp')
           backspace(tfile)
           read(tfile,*) st1, st2
           if (trim(st2).eq.'on') then
             ecp = .true.
             print *, 'ecp flag is switched on: core states will be integrated out'
           else
             print *, 'ecp flag is switched off'
             ecp = .false.
           end if
           ecp_found = .true.

          case ('$scfmo')
           backspace(tfile)
           read(tfile,*) st1, st2
	   tmplen = len_trim(st2) 
           mos_file_name = st2(6:tmplen) 
	   nspin = 1
	   print *, 'mos-file name   = ', trim(mos_file_name)
	   mos_found = .true.

	  case ('$uhfmo_alpha')
           backspace(tfile)
           read(tfile,*) st1, st2
	   tmplen = len_trim(st2) 
           alpha_file_name = st2(6:tmplen) 
	   nspin = 2
	   print *, 'alpha-file name = ', trim(alpha_file_name)
           alpha_found = .true.

	  case ('$uhfmo_beta')
           backspace(tfile)
           read(tfile,*) st1, st2
	   tmplen = len_trim(st2) 
           beta_file_name = st2(6:tmplen) 
	   nspin = 2
	   print *, 'beta-file name  = ', trim(beta_file_name)
           beta_found = .true.

          case ('$self_energy')
           backspace(tfile)
           read(tfile,*) st1, st2
	   tmplen = len_trim(st2) 
           self_energy_file = st2(6:tmplen) 
	   print *, 'self-energy file name = ', trim(self_energy_file)
	   import_self_energy = .true.

          case ('$nsaos')
           backspace(tfile)
           read(tfile,*) st1, nsaos
           print 100, ' nsaos = ', nsaos

          case ('$lsurc')
           backspace(tfile)
           read(tfile,*) st1, left_satoms(1)
           print 100, ' left  surface center = ', left_satoms(1)

          case ('$lsurx')
           backspace(tfile)
           read(tfile,*) st1, left_satoms(2)
           print 100, ' left  surface x      = ', left_satoms(2)

          case ('$lsury')
           backspace(tfile)
           read(tfile,*) st1, left_satoms(3)
           print 100, ' left  surface y      = ', left_satoms(3)

          case ('$rsurc')
           backspace(tfile)
           read(tfile,*) st1, right_satoms(1)
           print 100, ' right surface center = ', right_satoms(1)

          case ('$rsurx')
           backspace(tfile)
           read(tfile,*) st1, right_satoms(2)
           print 100, ' right surface x      = ', right_satoms(2)

          case ('$rsury')
           backspace(tfile)
           read(tfile,*) st1, right_satoms(3)
           print 100, ' right surface y      = ', right_satoms(3)

          case ('$dist')
           backspace(tfile)
           read(tfile,*) st1, dist
           print 200, ' coulping region range= ', dist, 'a.u.'
           dist_found = .true.

          case ('$nlayers')
           backspace(tfile)
           read(tfile,*) st1, nlayers
           print '(a,i3,a)', ' coulping region range=', nlayers, ' atomic layers'
           nlayers_found = .true.

          case ('$s1r')
           backspace(tfile)
           read(tfile,*) st1, rsigma(1)
           print 211, ' r(sigma) 1st layer = ', rsigma(1), 'H'
           found_rsigma(1)= .true.

          case ('$s1i')
           backspace(tfile)
           read(tfile,*) st1, isigma(1)
           print 211, ' i(sigma) 1st layer = ', isigma(1), 'H'
           found_isigma(1)= .true.

          case ('$s2r')
           backspace(tfile)
           read(tfile,*) st1, rsigma(2)
           print 211, ' r(sigma) 2nd layer = ', rsigma(2), 'H'
           found_rsigma(2)= .true.

          case ('$s2i')
           backspace(tfile)
           read(tfile,*) st1, isigma(2)
           print 211, ' i(sigma) 2nd layer = ', isigma(2), 'H'
           found_isigma(2)= .true.

          case ('$s3r')
           backspace(tfile)
           read(tfile,*) st1, rsigma(3)
           print 211, ' r(sigma) 3rd layer = ', rsigma(3), 'H'
           found_rsigma(3)= .true.

          case ('$s3i')
           backspace(tfile)
           read(tfile,*) st1, isigma(3)
           print 211, ' i(sigma) 3rd layer = ', isigma(3), 'H'
           found_isigma(3)= .true.

          case ('$dz_tolerance')
           backspace(tfile)
           read(tfile,*) st1, dz_tolerance
           print 200, ' dz_tolerance = ', dz_tolerance, 'a.u.'

          case('$hexchange')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            sp_leads = .true. ;  hexfound = .true.
            print '(/,a)', ' electrodes are SPIN POLARIZED'
           end if

          case('$hsource')
	   locexfound = .true.
	   localexf   = .true.
           backspace(tfile)
           read(tfile,*) st1, st2
	   tmplen = len_trim(st2) 
           hsrc_file_name = st2(6:tmplen) 
           write(*,'(a,a)') 
     &	        ' local ex-fields source file = ', trim(hsrc_file_name)

          case ('llead') 
           backspace(tfile)
           read(tfile,*) st1, hex_left
           print 200, ' left lead : exchange field = ', hex_left, 'H'

          case ('rlead') 
           backspace(tfile)
           read(tfile,*) st1, hex_right
           print 200, ' right lead: exchange field = ', hex_right, 'H'

          case('$symmetrize')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            symmetrize  = .true.
            print *, 'symmetrization is switched ON'
           else
            symmetrize  = .false.
            print *, 'symmetrization is switched OFF'
           end if

          case('$z-symmetrize')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='ap') then
            z_symmetrize  = .true.
            ap_symmetrize = .true.
            p_symmetrize  = .false.
            print *, 'ap symmetrization is switched ON'
           else if (st2=='p') then
            z_symmetrize    = .true.
            p_symmetrize  = .true.
            ap_symmetrize = .false.
            print *, 'p symmetrization is switched ON'
           else
            z_symmetrize  = .false.
            p_symmetrize  = .false.
            ap_symmetrize = .false.
            print *, 'symmetrization is switched OFF'
           end if

          case('$x-symmetrize')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            x_symmetrize  = .true.
            print *, 'x-symmetrization is switched ON'
           else
            x_symmetrize  = .false.
            print *, 'x-symmetrization is switched OFF'
           end if

          case('$y-symmetrize')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            y_symmetrize  = .true.
            print *, 'y-symmetrization is switched ON'
           else
            y_symmetrize  = .false.
            print *, 'y-symmetrization is switched OFF'
           end if

          case('$restart_mo')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            mo_restart  = .true.
            print *, 'MOs restart information will be appended to temporary mos file(s)'
           else
            mo_restart  = .false.
            print *, 'MOs restart information will NOT be appended to temporary mos file(s)'
           end if

          case('$expert')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            expert = .true.
            print *, 'CAUTION! "expert" flag is switched ON !'
           else
            expert = .false.
            print *, '"expert" flag is switched off'
	   end if

          case('$adjust_rsigma')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            tune_rsigma = .true.
            print *, 'tuning real piece of the self-energy: ON'
           else
            tune_rsigma = .false.
            print *, 'tuning real piece of the self-energy: OFF'
	   end if

          case ('rfactor')
           backspace(tfile)
           read(tfile,*) st1, rfactor
           print 220, ' # rfactor= ', rfactor

          case ('rguess')
           backspace(tfile)
           read(tfile,*) st1, rguess
           print 220, ' # rguess = ', rguess

          case ('nelcnv') 
           backspace(tfile)
           read(tfile,*) st1, nelcnv
           print 220, ' # nelcnv = ', nelcnv

          case ('iterlimit') 
           backspace(tfile)
           read(tfile,*) st1, iter_limit
           print 100, ' # iterlimit =', iter_limit

          case('$densmat_cycle')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            do_dmat_cycle = .true.
            print *, 'self-consistency cycle: ON'
           else
c           ! default value defined in module <globalvars.f>
c           do_dmat_cycle = .false.
            print *, 'self-consistency cycle: OFF'
           end if

          case('$densmat') ; 
          case('charge.dens')
           print *, 'charge-density file= dmat'
          case('spin.dens')
           print *, 'spin-density file  = smat'

          case ('$nelectr')
           backspace(tfile)
           read(tfile,*) st1, nelectr
           print 100, ' number of electrons= ', nelectr
           if (nelectr.gt.0) then ; nelfound = .true.
           else ; stop ' number of electrons is 0 or negative ??? '
	   end if

          case ('$valence_electrons')
           backspace(tfile)
           read(tfile,*) st1, tmp_val_electrons
           print 100, ' number of valence electrons= ', abs(tmp_val_electrons)

          case ('$bias')
           backspace(tfile)
           read(tfile,*) st1, bias
           print 200, ' bias-voltage = ', bias, 'eV'
           biasfound = .true.

          case('$dmix') 
           backspace(tfile)
           read(tfile,*) st1, dmix
           print 220, ' admixing factor = ', dmix
           dmixfound = .true.
c          if keyword $dmix is not found, default value (0.1) 
c          is used, as defined in <globalvars.f>

          case('$densconv') 
           backspace(tfile)
           read(tfile,*) st1, densconv
           print 220, ' dens.convergence criterium   = ', densconv

c         == hubbard-u subgroup ==
          case('$hubbard-u')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            ldau_found = .true.
            print '(/,a)', ' <Hubbard-U> flag is ACTIVE'
           else if (st2=='old') then
            ldau_found = .true.
            ldau_old   = .true.
            print '(/,a)', ' <Hubbard-U> is ACTIVE: CAUTION! old formulation will be used'
           else
            print '(/,a)', ' <Hubbard-U> flag is switched OFF'
           end if

          case ('num_ratom_types')
           backspace(tfile)
           read(tfile,*) st1, num_ratoms
           print 100, ' <Hubbard-U> number of reservoir atom types: ', num_ratoms        
	   if (num_ratoms.ge.1) num_ratoms_found = .true.

          case ('ratom_types')
           backspace(tfile)
	   if (num_ratoms_found) then
	    if (.not.allocated(ratom_types)) 
     &	      allocate(ratom_types(num_ratoms),stat=ierr)
            if (ierr.ne.0) 
     &        stop '[SUBR. externals]: <ratoms_types> allocation failure ' 
            read(tfile,*) st1, (ratom_types(ra),ra=1,num_ratoms)
            print '(a,25(i3))', ' <Hubbard-U> reservoir atom types: ', 
     &	      (ratom_types(ra),ra=1,num_ratoms)        
            ratom_types_found = .true.
	   else
	    print '(/,a,/,a,/,a,/)', 
     &	     ' ERROR : Hubbard-U: number <n> of reservoir atom types is not found!',
     &	     '         please, check your <tcontrol> file:',
     &	     '         "num_ratom_types <n>" line should open "$hubbard-u on/off" group'
	    stop ': transport module is terminated now'
	   end if

          case('u-j')
           backspace(tfile)
           read(tfile,*) st1, u_j
           u_j_found = .true.
           print 200, ' <Hubbard-U> screened coulomb interaction =', u_j, 'eV'

          case('spectral_function')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            spectr_func = .true.
            print '(a)', ' <Hubbard-U> spectral function will be output'
           else
            spectr_func = .false.
            print '(a)', ' <Hubbard-U> spectral function will NOT be output'
           end if

          case('read_dmat')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            read_dmat = .true.
            print '(a)', ' <Hubbard-U> density matrix will be imported'
           else
            read_dmat = .false.
            print '(a)', ' <Hubbard-U> density matrix will NOT be imported'
           end if
           print *

          case ('$hion') ;           
          case ('action')            
	   backspace(tfile)
           read(tfile,*) st1, st2
	   if (trim(st2).eq.'save') then 
	     save_hion = .true.
             print *, 
     &	       '<Hubbard-U>: hamiltonian of the 3d-ion and its dens.matrix will be saved'
           else if (trim(st2).eq.'read') then 
	     read_hion = .true.
             print *, 
     &	       '<Hubbard-U>: hamiltonian of the 3d-ion and its dens.matrix will be imported'
	   end if
 
          case ('source')
           backspace(tfile)
           read(tfile,*) st1, st2
	   tmplen = len_trim(st2) 
           hion_file_name = st2(6:tmplen) 
	   print '(1x,a,a)', '<Hubbard-U>: 3d-ion hamiltonian file name = ', trim(hion_file_name)
           print *

c         == end of hubbard-u subgroup ==

c         == "mol_subsystem" subgroup ==
          case('$mol_subsystem')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            mol_subsystem_found = .true.
            print '(/,a)', ' <mol_sybsystem> flag is ACTIVE'
           else
            print '(/,a)', ' <mol_sybsystem> flag is switched OFF'
           end if

          case ('num_res_atom_types')
           backspace(tfile)
           read(tfile,*) st1, num_res_atoms
           print 100, ' <mol_subsystem>: number of reservoir atom types: ', num_res_atoms        
	   if (num_res_atoms.ge.1) num_res_atoms_found = .true.

          case ('res_atom_types')
           backspace(tfile)
	   if (num_res_atoms_found) then
	    if (.not.allocated(res_atom_types)) 
     &	      allocate(res_atom_types(num_res_atoms),stat=ierr)
            if (ierr.ne.0) 
     &        stop '[SUBR. externals]: <res_atoms_types> allocation failure ' 
            read(tfile,*) st1, (res_atom_types(ra),ra=1,num_res_atoms)
            print '(a,25(i3))', ' <mol_subsystem>: reservoir atom types: ', 
     &	      (res_atom_types(ra),ra=1,num_res_atoms)        
            res_atom_types_found = .true.
	   else
	    print '(/,a,/,a,/,a,/)', 
     &	     ' ERROR: [mol_subsystem]: number <n> of reservoir atom types is not found!',
     &	     '        please, check your <tcontrol> file:',
     &	     '        "num_res_atom_types <n>" line should open "$mol_subsystem on/off" group'
	    stop ': transport module is terminated now'
	   end if

          case ('read_hmat')
           backspace(tfile)
           read(tfile,*) st1, st2
           tmplen = len_trim(st2)
           inhmat_file_name = st2(6:tmplen)
           print '(a,/,a,a)', ' <mol_subsystem>: KS-hamiltonian of the molecule ', 
     &	         '                  to be read from file: ', trim(inhmat_file_name)
           read_hmat = .true.

          case('mol_spectral_function')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            mol_spectr_func = .true.
            print '(a)', ' <mol_subsystem>: spectral function will be output'
           else
            mol_spectr_func = .false.
            print '(a)', ' <mol_subsystem>: spectral function will NOT be output'
           end if

          case ('occ_mol_orbitals')
           backspace(tfile)
           read(tfile,*) st1, occ_mol_orbitals
	   occ_mol_orbitals = max(1,occ_mol_orbitals)
           print '(a,i3,a)', ' <mol_subsystem>: local DOS will be projected to ', 
     &	                      occ_mol_orbitals, ' occupied molecular orbital(s)'        

          case ('empty_mol_orbitals')
           backspace(tfile)
           read(tfile,*) st1, empty_mol_orbitals
	   empty_mol_orbitals = max(1,empty_mol_orbitals)
           print '(a,i3,a,/)', ' <mol_subsystem>: local DOS will be projected to ', 
     &	                      empty_mol_orbitals, ' empty molecular orbital(s)'        

c         == end of "mol_subsystem" subgroup ==

          case ('$save_hmat')
           backspace(tfile)
           read(tfile,*) st1, st2
           tmplen = len_trim(st2)
           outhmat_file_name = st2(6:tmplen)
           print *, 'KS-hamiltonian to be save to file: ', trim(outhmat_file_name)
           save_hmat = .true.

          case('$dnorm') 
           backspace(tfile)
           read(tfile,*) st1, dnorm
           dnormfound = .true.
	  
          case ('$efsearch') ;
          case ('nelconv') 
           backspace(tfile)
           read(tfile,*) st1, nelconv
           print 220, ' conv.criterium for el.number = ',nelconv

          case ('nclstmo') 
           backspace(tfile)
           read(tfile,*) st1, nclstmo
c	   print 100, ' nclstmo = ', nclstmo

           case ('$saveoccn') 
            print *, '"$saveoccn" is an obsolete keyword: it will not have any effect'
c	   backspace(tfile)
c           read(tfile,*) st1, st2
c	   if (trim(st2).eq.'on') then 
c	     store_occn = .true.
c             print *, 'pseudo-occupation numbers to be saved: yes'
c           else
c	     store_occn = .false.
c             print *, 'pseudo-occupation numbers to be saved: no'
c	   end if     

          case ('$output_charges') 
           backspace(tfile)
           read(tfile,*) st1, st2
	   if (trim(st2).eq.'on') then 
	     qoutput = .true.
             print *, 'loewdin charges/magn.moments will be printed out'
           else
	     qoutput = .false.
             print *, 'loewdin charges/magn.moments will NOT be printed out'
	   end if     

          case('$population')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            do_pop = .true.
            print *, 'population analysis is switched on'
           else if (st2=='lm-projected') then
            do_pop   = .false.
            do_lmpop = .true.
            print *, 'lm-decomposed population analysis is switched on'
           end if

          case('$ldos')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            do_ldos = .true.
            print *, 'LDOS calculation is switched on'
           else if (st2=='lm-projected') then
            do_ldos  = .true.
	    do_lmdos = .true.
            print *, 'lm-projected LDOS calculation is switched on'
           else if (st2=='3d') then
            do_ldos3d  = .true.
	    print *, '-- computation of the space resolved spectral function is switched on'
           else if (st2=='3d_ewindow') then
            do_ldos3d_ewindow  = .true.
	    print *, '-- computation of the space resolved spectral function is switched on'
	    print *, '-- spectral function will be integrated over given energy window'
           else 
            do_ldos = .false.
            print *, 'calculation of LDOS is switched off'
           end if

          case ('scaling_factor')
           backspace(tfile)
           read(tfile,*) st1, scaling_fctr 	   
	   print '(a,e12.6)', 
     &	    ' -- spectral function will be multiplied by a scaling factor = ', scaling_fctr
c          if not found, scaling factor will be estimated automatically later on, 
c          based on the value of A(r,EF) at r=(0,0,zmin)

          case ('z_range')
           backspace(tfile)
           read(tfile,*) st1, zmin, zmax
	   print '(a,f8.4,a,f8.4,a)', ' -- zmin = ', zmin, 
     &	                              ' a.u. ;  zmax = ', zmax, ' a.u.'
           z_range_found = .true.    
           if (zmin > zmax) then
            tmpz = zmin ; zmin = zmax ; zmax = tmpz
           end if

          case('$transmission')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            do_trans_calc = .true.
            print *, 'transmission calculation: ON'
c          >>> modified by RK and AB :
           else if (st2=='channels') then
            do_trans_calc    = .true.
            do_cond_channels = .true.
            print *, 'transmission calculation: ON'
            print *, '-- conduction eigenchannels are requested'
c         <<< done with modification by RK and AB
           else if (st2=='wave_functions') then
            do_trans_calc      = .true.
            do_cond_channels   = .true.
            do_scatt_wave_func = .true.
            print *, 'transmission calculation: ON'
            print *, '-- scattering wave functions are requested'
          else
c         ! default valued defined in module <globalvars.f>
c           do_trans_calc = .false.
            print *, 'transmission calculation: OFF'
          end if

          case ('wf_scaling_factor')
           backspace(tfile)
           read(tfile,*) st1, wf_scaling_fctr 	   
	   print '(a,e12.6)', 
     &	    ' -- scattering wave functions will be multiplied by a scaling factor = ', wf_scaling_fctr
c          if not found, scaling factor will be estimated automatically later on, 
c          based on the value of |Psi_n(r,E)| at r=(0,0,zmin)

          case ('wf_z_range')
           backspace(tfile)
           read(tfile,*) st1, wf_zmin, wf_zmax
	   print '(a,f8.4,a,f8.4,a)', ' -- wf_zmin = ', wf_zmin, 
     &	                              ' a.u. ;  wf_zmax = ', wf_zmax, ' a.u.'
           wf_z_range_found = .true.    
           if (wf_zmin > wf_zmax) then
            tmpz = wf_zmin ; wf_zmin = wf_zmax ; wf_zmax = tmpz
           end if

          case('$evunits')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            evunits = .true.
            print '(/,1x,a)', 'CAUTION: you have requested an energy mesh to be defined in eV vs Fermi'
            print '(1x,a,/)', '         level, not in (default) Hartree units as may be written above!'
           end if

          case ('$ener')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
	   if (st2=='undefined') then
	    ener = 0.0d0 ! just to be on the safe side
	    print *, 'first E-point = ', trim(st2)
	    conductance = .true.
	   else
            ener_found = .true.
	    backspace(tfile)
      	    read(tfile,*) st1, ener
            if (.not.evunits) then
             print 210, ' first E-point = ', ener, 'H'
	    else
             print 210, ' first E-point = ', ener, 'eV vs EF'
            end if
	   end if

          case ('$estep')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
	   if (st2=='undefined') then
	    estep = 0.0d0 ! just to be on the safe side
	    print *, 'dE step       = ', trim(st2)
	    conductance = .true.
	   else
	    estep_found = .true.
	    backspace(tfile)
      	    read(tfile,*) st1, estep
            if (.not.evunits) then
             print 210, ' dE step       = ', estep, 'H'
	    else
	     print 210, ' dE step       = ', estep, 'eV'
	    end if
	   end if

          case ('$eend')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
	   if (st2=='undefined') then
	    eend = 0.0d0 ! just to be on the safe side
	    print *, 'last E-point  = ', trim(st2)
            conductance = .true.
	   else
	    eend_found = .true.
	    backspace(tfile)
      	    read(tfile,*) st1, eend
      	    if (.not.evunits) then
             print 210, ' last E-point  = ', eend, 'H'
	    else
	     print 210, ' last E-point  = ', eend, 'eV vs EF'
	    end if
	   end if

          case ('$eta')
           backspace(tfile)
           read(tfile,*) st1, eta
           eta = dabs(eta)
	   print 300, ' infinit. eta = ',  eta, 'H'

          case ('$output')
           backspace(tfile)
           read(tfile,*) st1, st2
	   tmplen = len_trim(st2) 
           output_file_name = st2(6:tmplen) 
	   print *, 'output-file name = ', trim(output_file_name)

          case ('$save_omat')
           backspace(tfile)
           read(tfile,*) st1, st2
           tmplen = len_trim(st2)
           outomat_file_name = st2(6:tmplen)
           print *, 'overlap matrix to be save to file: ', trim(outomat_file_name)
           save_omat = .true.

          case ('$read_omat')
           backspace(tfile)
           read(tfile,*) st1, st2
           tmplen = len_trim(st2)
           inomat_file_name = st2(6:tmplen)
           print *, 'overlap matrix to be read from file: ', trim(inomat_file_name)
           read_omat = .true.
        
          case ('$dens_matrix')
	   backspace(tfile)
	   read(tfile,*) st1, st2
	   if (trim(st2).eq.'on') then 
	    calc_densmat0 = .true.
	    print *, 'equlibrium density matrix calculation = on'
	   else 
	    calc_densmat0 = .false.
	    print *, 'equilibrium density matrix calculation = off'
           end if
          
	  case ('num.of.occupied')
           backspace(tfile)
           read(tfile,*) st1, st2, tmp_occ
           select case (trim(st2))
	    case ('orbitals') ; occ = tmp_occ 
	      print 100, ' number of occ. orbitals  = ', occ
            case ('alpha-orbitals') ; a_occ = tmp_occ
             print 100, 
     &	      ' alpha spin: number of occ. orbitals = ', a_occ
	    case ('beta-orbitals') ; b_occ = tmp_occ 
             print 100, 
     &	      ' beta spin:  number of occ. orbitals = ', b_occ
    	    case default ; 
	   end select    

          case ('$testing')
           backspace(tfile)
           read(tfile,*) st1, st2
           if (trim(st2).eq.'on') then
             testing = .true.
             print *, 'TESTING is switched on '
           else
             testing = .false.
           end if

          case ('$zspectrum')
           backspace(tfile)
           read(tfile,*) st1, st2
           if (trim(st2).eq.'on') then
             zspectrum = .true.
             print *, 'complex-valued spectrum to be output: yes'
           else
             zspectrum = .false.
             print *, 'complex-valued spectrum to be output: no'
           end if

          case('$storesmat')
           backspace(tfile)
           read(tfile,*) st1, st2
           st2 = trim(st2)
           if (st2=='on') then
            storesmat = .true.
            print *, 'allow to use stored overlap matrices: yes'
           else
            storesmat = .false.
            print *, 'allow to use stored overlap matrices: no'
           end if

          case ('$overlapint')
           importsmat = .true.
           print *, 'overlap matrices to be read from hard drive'

          case ('smat12')
           backspace(tfile)
           read(tfile,*) st1, st2
	   tmplen = len_trim(st2) 
           smat12_file_name = st2(6:tmplen) 
	   print *, 'smat12-file name    = ', trim(smat12_file_name)

          case ('smat12inv')
           backspace(tfile)
           read(tfile,*) st1, st2
	   tmplen = len_trim(st2) 
           smat12inv_file_name = st2(6:tmplen) 
	   print *, 'smat12inv-file name = ', trim(smat12inv_file_name)

          case ('$efermi')
	   backspace(tfile)
           read(tfile,*) st1, efermi
	   print 211, ' Fermi energy guess = ', efermi, 'H'
	   efermi_guess = .true.

          case ('$fixed_efermi')
           backspace(tfile)
           read(tfile,*) st1, st2
           if (trim(st2).eq.'on') then
             fixed_efermi = .true.
             print *, 'CAUTION! fixed_efermi flag is switched ON !'
           else
             fixed_efermi = .false.
           end if

          case ('$iter')
	   iterfound = .true.
	   backspace(tfile)
           read(tfile,*) st1, iter  	     
	   
	  case('$u_energy') ; 
	   
	  case ('$end') 
	   readnext = .false.

          case default ;
           if (keystring(1:1)=='$') then
            print '(/,1x,a,a,a)', 'ERROR: "',trim(keystring),'" is unknown keyword'
            print '(1x,a,/)',     '       may be it is a typo? ... please, check'
            stop ': transport module is terminated'
           end if

         end select
	end do  ! readnext
	
	close (tfile)

        if (.not.aims_input) then
c       checking for 'coord' file in case of turbomole input
         if (.not.coord_found) then
	  print '(/,a,/,a,/,/,a,/)', 
     &      ' ERROR : you have not specified file with atomic coordinates',
     &      '         required keyword is $coord file=<my_coord_file>',
     &      '         please, check your <tcontrol> file'    
	  stop ': transport module is terminated now'
	 end if
        else
c       checking for 'geometry.in' file in case of fhi-aims input
         if (.not.coord_found) then
	  print '(/,a,/,a,/,a,/,/,a,/)', 
     &      ' ERROR : you have not specified a file with atomic coordinates ',
     &      '         this file is called "geometry.in" in case of FHI-aims',
     &      '         required keyword is $coord file=<my_geo_file>',
     &      '         please, check your <tcontrol> file'    
	  stop ': transport module is terminated now'
         end if
	end if ! .not.aims_input

c       # check whether number of atoms is read properly
        if (.not.natoms_found) then 
         print '(/,a,/,a,/,a,/)', 
     &	  ' ERROR : number of atoms has not been specified:',
     &    '         this number should follow a keyword',
     &    '         $natoms in your <tcontrol> file'
         stop ': transport module is terminated now'
        end if
      
c       # take care about left & right electrodes 
        flag = left_satoms(1)<=0.or.left_satoms(1)>num_atoms.or.
     &         left_satoms(2)<=0.or.left_satoms(2)>num_atoms.or.
     &         left_satoms(3)<=0.or.left_satoms(3)>num_atoms
        if (flag) then
         print '(/,a,/,a,/)', 
     &	  ' ERROR : atoms fixing LEFT electrode are out of range',
     &    '         please, check your <tcontrol> file'
         stop ': transport module is terminated now'
        end if

        flag = right_satoms(1)<=0.or.right_satoms(1)>num_atoms.or.
     &         right_satoms(2)<=0.or.right_satoms(2)>num_atoms.or.
     &         right_satoms(3)<=0.or.right_satoms(3)>num_atoms
        if (flag) then
         print '(/,a,/,a,/,a,/)', 
     &	  ' ERROR : atoms fixing RIGHT electrode are out of range',
     &    '         please, check your <tcontrol> file'
         stop ': transport module is terminated now'
        end if

        flag = left_satoms(1)==left_satoms(2) .or.
     &         left_satoms(2)==left_satoms(3) .or.
     &         left_satoms(3)==left_satoms(1)
        if (flag) then
         print '(/,a,/,a,/,a,/)', 
     &	  ' ERROR : your LEFT electrode is not specified correctly:',
     &    '         atoms fixing the surface plane should be different',
     &    '         please, check your <tcontrol> file'
         stop ': transport module is terminated now'
        end if

        flag = right_satoms(1)==right_satoms(2) .or.
     &         right_satoms(2)==right_satoms(3) .or.
     &         right_satoms(3)==right_satoms(1) 
        if (flag) then
         print '(/,a,/,a,/,a,/)', 
     &	  ' ERROR: your RIGHT electrode is not specified correctly:',
     &    '        atoms fixing the surface plane should be different',
     &    '        please, check your <tcontrol> file'
         stop ': transport module is terminated now'
        end if         

        if ((.not.dist_found).and.(.not.nlayers_found)) then
c        self-energy region is not specfied!	
         print '(/,a,/,a,/,a,/)', 
     &	  ' ERROR: amount of atomic layers defining a "coupling region" should follow',
     &	  '        a keyword $nlayers ; this keyword is missing in <tcontrol> file: ', 
     &    '        a self-energy matrix can not be initialized'
         stop ': transport module is terminated now'
        end if

c       if both keywords $dist and $nlayers are found, alarm 
c       a conflict between the new ($nlayers) and old-style ($dist) keywords: 
        if (dist_found.and.nlayers_found) then
         print '(/,a,/,a,/,a,/)', 
     &	  ' ERROR: you have used both keywords, $dist and $nlayers to define a',
     &	  '        "coupling region" for the self-enerfy matrix: i am confused ;',      
     &    '        please, use either of these keywords, but only one !'
         stop ': transport module is terminated now'
        end if
        
        flag = .false.
        do i = 1, 3
         do j = 1, 3 
	  tmpf = flag
          flag = tmpf.or.(left_satoms(i) == right_satoms(j))
         end do
        end do
        if (flag) then
         print '(/,a,/,a,/)', 
     &    ' ERROR : it seems, your LEFT and RIGHT electrodes share same atoms',
     &    '         please, check your <tcontrol> file'
         stop ': transport module is terminated now'
        end if
c       # done with self-energy and elecrodes

c       checking for mos/alpha/beta files >>>
        if (.not.( mos_found.or.(alpha_found.and.beta_found) ) ) then
	  print '(/,a,/,a,/)', 
     &      ' ERROR : you have not specified file(s) with molecular orbitals',
     &      '         please, check your <tcontrol> file' 
	  stop ': transport module is terminated now'
	end if 

c       checking for overlap matrix in case of aims-input  ...
        if (aims_input.and.(.not.read_omat)) then
	  print '(/,a,/,a,/,a,/)', 
     &      ' ERROR : you have not specified file with overlap matrix elements',
     &      '         required keyword is $read_omat file=<your_overlap_file>',
     &      '         please, check your <tcontrol> file'    
	  stop ': transport module is terminated now'
	end if 

        if (.not.aims_input) then
c       checking for 'basis' file in case of turbomole input
         if (.not.basis_found) then
	  print '(/,a,/,a,/,/,a,/)', 
     &      ' ERROR : you have not specified a "basis" file with information on',
     &      '         CGT orbitals: required keyword is $basis file=<my_basis_file>',
     &      '         please, check your <tcontrol> file'    
	  stop ': transport module is terminated now'
	 end if
        else
c       checking for 'basis-indices.out' file in case of fhi-aims input
         if (.not.basis_found) then
	  print '(/,a,/,a,/,a,/,/,a,/)', 
     &      ' ERROR : you have not specified a file indexing your basis functions; this file',
     &      '         is written out by FHI-aims and usually has a name "basis-indices.out" ;',
     &      '         required keyword is $basis file=<my_basis_file>',
     &      '         please, check your <tcontrol> file'    
	  stop ': transport module is terminated now'
         end if
	end if ! .not.aims_input

        IF (landauer) THEN
c        a flag '$landauer' is introduced ...
c        (i)   for testing purposes;..
c        (ii)  it fixes the only choice available so far for
c              the fhi-aims user in the first common release.
c        (iii) it means:
c              take a self-energy with imaginarty piece only,
c              build up a Green's function and compute ballistic.
c              Landauer transmission function......
          
c        #caution! '$landauer' flag overrides all other possible flags
c                   and options ; switch it off or remove it from <tcontrol> 
c                   if you are interested in advanced tasks ... 

c         flag = tune_rsigma.or.do_dmat_cycle.or.
c     &          do_ldos.or.do_ldos3d.or.do_ldos3d_ewindow.or.
c     &          do_cond_channels.or.do_scatt_wave_func.or.
c     &          mol_subsystem_found.or.ldau_found.or.no_leads.or.do_pop
c
         flag = tune_rsigma.or.do_dmat_cycle.or.
     &          do_ldos3d.or.do_ldos3d_ewindow.or.
     &          do_cond_channels.or.do_scatt_wave_func.or.
     &          mol_subsystem_found.or.ldau_found.or.no_leads.or.do_pop
    
         if (flag) then
          print '(/,1x,a)', 
c    &       'WARNING! "$landauer on" option will override all other flags'
     &       'WARNING! "$landauer on/off" option will override all other flags (except of "$ldos")'
          print '(1x,a)',   
     &       '         introduced by you, including "$tune_rsigma", "$densmat_cycle",'
          print '(1x,a)',   
c    &       '         "$transmission", "$ldos", "mol_subsystem", "$ldau",  '
     &       '         "$transmission", "mol_subsystem", "$ldau",  '
          print '(1x,a)',   
     &       '         "no_leads", "$population", etc. ...  '
          print '(/,1x,a)',   
     &       '         if you wish to proceed with any of advanced tasks, please'
          print '(1x,a)',   
     &       '         comment out or remove $landauer flag in your <tcontrol> file '

         end if

c        >>> checking, whether transsmion calculation or ldos calculation is required ...
         if ( (.not.do_landauer).and.(.not.do_ldos) ) then
          print '(/,a)',
     &      ' ERROR : both job flags "$landauer" and "$ldos" are switched off '
          print '(a,/)',
     &      '         please, specify correctly your choice'
          stop ': transport module is terminated now'
         end if
 
c        >>> checking, whether self-energy parameters have been read ...
         if ( .not.( (found_isigma(1)).and.(found_isigma(2)).and.(found_isigma(3)) ) ) then
          print '(/,1x,a)',
     &       'WARNING: a model self-energy has not been specified properly,'
          print '(1x,a)',
     &       '         a default set of parameters (0.1H, 0.05H and 0.025H) will'
          print '(1x,a)',
     &       '         be used to parametrize an imaginary piece of the self-energy'

          do i = 1, deep
           isigma(i) = def_isigma/dble(2**(i-1))
          end do

         end if ! check for isigma

c        in case of $landauer flag, set a real piece of self-energy
c        to zero in any case, independent of the input parameters
         do i = 1, deep
          rsigma(i) = 0.0d0
         end do

c        >>> checking for electron number in case of turbomole input ;
c        >>> in case of all electron treatment (fhi-aims) electron number 
c        >>> is computed based on periodic table
c        >>> comment added: jan 2012 
         if ((.not.nelfound).and.(.not.aims_input)) then
	  print '(/,a,/)',' number of electrons is not set!'
	  stop ': transport module is terminated now'
         end if

c        "effective core potentials" can be active in case of fhi-aims only
         if (ecp.and.ecp_found.and.(.not.aims_input)) then
          print '(/,1x,a)',
     &       'WARNING: $ecp flag has no effect in case of TURBOMOLE input'
         end if
         if (.not.aims_input) ecp = .false.

c        output warning message if ecp flag is explicitly switched off
         if (aims_input.and.(.not.ecp)) then
          print '(/,1x,a,/,1x,a,/,1x,a,/,1x,a)',
     &       'WARNING: you have switched off flag $ecp! however, in case of "fhi-aims" based input',
     &       '         usage of "effective core potentials" ($ecp on) flag in strongly recommended to',
     &       '         speed up calculations: in that case core states will be effectively "integrated out"',
     &       '         and dimension of the hilbert space will be given by the amount of valence states".'
         end if

c        if flag $ecp has not been found in <tcontrol>, 
c        put it there in case of fhi-aims input
         if (aims_input.and.(.not.ecp_found)) then
          call updatetcontrol(trim(tcntrl_file_name),'$ecp',1,0.0d0)
         end if

        ELSE
c        check all other advanced options >>>>        

c        for compatiblity with previous versions, checking for 
c        "$non-equilibrium" keyword >>>
         if (.not.neqflag) then
	  print '(/,a,/,a,/,a,/,a,/)', 
     &	        ' $non-equilibrium keyword is explicitly switched off: ' , 
     &          ' to proceed working with current version of the code, ' ,
     &          ' comment out a line "$non-equilibrium off" in your '    ,
     &          ' <tcontrol> file or modify it to "$non-equilibrium on" '
	  stop  ': transport module is terminated now'
	 end if

c        do not use self-energy input file if no $landauer flag is found
         if (import_self_energy.and.(.not.landauer)) then
	   print '(/,a)',' WARNING: self-energy parameters will be imported from your <tcontrol> file'
	   print '(a)',  '          but not from the external file, since no keyword $landauer is found'
           import_self_energy = .false.
         end if

c        >>> motivated by a development of the transport code based on fhi-aims
c            from jan 2012, '$tune_rsigma' flag does not conflict any more
c            with a 'transmission' or 'ldos' type of calculation 
c        !!! the only exception is a conflict with the self-consistency 
c            dens.matrix cycle

         if (tune_rsigma.and.do_dmat_cycle.and.(.not.expert)) then
	  print '(/,a,/,a,/)',
     &     ' ERROR : job flag "$adjust_rsigma" conflicts with "$densmat_cycle" flag',
     &     '         only one of above flags should be switched "on" ' 
	  print '(a,/)',  
     &	   ' please, specify correctly which type of calculation you are willing to do ...'
	  stop ': transport module is terminated now'
         end if

c        checking for consistency of job flags
         if (do_dmat_cycle.and.do_trans_calc) then
	  print '(/,a,/,a)',
     &	        ' ERROR : job flags "$densmat_cycle" & "$transmission"',
     &	        '         are not specified correctly! '
	  print '(a,/)',  
     &	        '         only one of above flags should be "on"'
	  stop ': transport module is terminated now'
         end if

c        further checking for consistency of job flags
         if (do_dmat_cycle.and.(do_ldos.or.do_ldos3d.or.do_ldos3d_ewindow)) then
	  print '(/,a,/,a)',' ERROR : job flags "$densmat_cycle" & "$ldos"',
     &	                    '         are not specified correctly! '
	  print '(a,/)',    '         only one of above flags should be "on"'
	  stop ': transport module is terminated now'
         end if

c        further checking for consistency of job flags >>> removed on jan 2012 !
c        flag = (do_dmat_cycle.or.do_ldos.or.do_ldos3d.or.do_ldos3d_ewindow.or.do_trans_calc)
c        if (tune_rsigma.and.flag) then
c	 print '(/,a,/,a)',
c     &    ' job flag "$adjust_rsigma" conflicts with other flags:',
c     &	  ' "$densmat_cycle" or "$ldos" or "$transmission"'
c	 print '(a,/)',  
c     &	  ' please, switch them off'
c	 stop ': transport module is terminated now'
c        end if

c        check for $bias flag: if not found, assume bias-voltage=0.0eV
         if (.not.biasfound) then
c         for now (jan 2012) black-box-like calculations with transport code & fhi-aims 
c         are performed only for bias-voltage=0eV: so, alarm only in case of TURBOMOLE calcs:
	  if (.not.aims_input) then
	   print '(/,a)',' WARNING: bias-voltage ($bias) is not specified: assume bias = 0.0 eV'
	  end if
	  bias = 0.0d0 
         end if

c        checking for energy points if ldos-calculation is required
         if (do_ldos.or.do_lmdos) then
          flag = (ener_found.and.estep_found.and.eend_found)
          if (.not.flag) then
           print '(/,a,/,a,/)',
     &     ' ERROR : you have not specified energy points for LDOS calculation',
     &     '         please, check your <tcontrol> file'
           stop ': transport module is terminated now'
          end if
         end if
         
c        # if energy mesh is not defined, transmission will be 
c        # evaluated at the Fermi energy 
c        # next check is commented out on jan 2012 >>>
c        checking for energy points if transmission calculation is required
c         if (do_trans_calc.and.(.not.ener_found)) then
c           print '(/,a,/,a,/,a,/)',
c     &     ' you have not specified the first energy point or',
c     &     ' energy mesh required for transmission calculation ...',
c     &     ' please, check your <tcontrol> file'
c           stop ': transport module is terminated now'
c         end if

c        if aims_input = .true. do not allow 3d/2d visualization of LDOS
         if (aims_input) then
          if (do_ldos3d.or.do_ldos3d_ewindow) then	 
           print '(/,a,/,a,/)',
     &     ' ERROR : visualization of the space resolved spectral function',
     &     '         is not yet implemented in case of FHI-aims input'
           stop ': transport module is terminated now'
	  end if
	 end if

c        checking for at least one energy point if the 
c        space-resolved spectral function has to be evaluated 
         if (do_ldos3d.and.(.not.ener_found)) then
           print '(/,a,/,a,/,a,/)',
     &     ' ERROR : you have not specified at least one energy point required',
     &     '         for calculation of the space resolved spectral function ...',
     &     '         please, check your <tcontrol> file'
           stop ': transport module is terminated now'
         end if

c        checking for energy window if the energy integrated 
c        space-resolved spectral function has to be evaluated 
         if (do_ldos3d_ewindow) then
          flag = (ener_found.and.eend_found)
          if (.not.flag) then
           print '(/,a,/,a,/)',
     &     ' ERROR : energy window is not specified for computing the r-resolved LDOS',
     &     '         please, set up values $ener and $eend [in Hartrees] in <tcontrol> file'
           stop ': transport module is terminated now'
          end if
         end if

c        checking for z-window if the space-resolved 
c        spectral function has to be evaluated 
         if ((do_ldos3d.or.do_ldos3d_ewindow).and.(.not.z_range_found)) then
           print '(/,a,/,a,/,a,/)',
     &     ' ERROR : you have not specified an interval along z-axis required',
     &     '         for computation of the space resolved spectral function ...',
     &     '         introduce a line "z_range <z_min> <z_max>" in <tcontrol> file'
           stop ': transport module is terminated now'
         end if

c        if aims_input = .true. do not allow visualization 
c        of the scattering wave functions
         if (aims_input.and.do_scatt_wave_func) then
           print '(/,a,/,a,/)',
     &     ' ERROR : visualization of the scattering wave functions',
     &     '         is not yet implemented in case of FHI-aims input'
           stop ': transport module is terminated now'
	 end if

c        checking for at least one energy point if 
c        scattering wave function are to be evaluated 
         if (do_scatt_wave_func.and.(.not.ener_found)) then
           print '(/,a,/,a,/,a,/)',
     &     ' ERROR : you have not specified an energy point "$ener" required',
     &     '         for the calculation of the scattering wave functions ...',
     &     '         please, check your <tcontrol> file'
           stop ': transport module is terminated now'
         end if

c        checking for z-window if scattering wave function are to be evaluated
         if (do_scatt_wave_func.and.(.not.wf_z_range_found)) then
           print '(/,a,/,a,/,a,/)',
     &     ' ERROR : you have not specified an interval along z-axis required',
     &     '         for computation of the scattering wave functions ...',
     &     '         introduce a line "wf_z_range <z_min> <z_max>" in <tcontrol> file'
           stop ': transport module is terminated now'
         end if

c        check for $nelectr, in case of fhi-aims input, 
c        (all electron treatment) electron number can be computed
c        >>> comment added: jan 2012 
         if ((.not.nelfound).and.(.not.aims_input)) then
	  print '(/,a,/)',' ERROR : number of electrons is not set!'
	  stop ': transport module is terminated now'
         end if

c        checking dmix:
         if (dmixfound) then
          if (dmix<0.0d0 .or. dmix>1.0d0) then
           print '(/,a)', ' found admixing factor is out of range'
	   print '(a,f4.2,a)', 
     &	     ' default value = ', defdmix, ' will be used'
	   dmix = defdmix
          end if
	 else
c          use a default value
	   dmix = defdmix
c          print a warning message, if dmix is really needed 
	   if (do_dmat_cycle) then
            print '(/,a)', 
     &	     ' admixing factor "$dmix" for the density matrix update is not specified'
	    print '(a,f4.2,a)', 
     &	     ' default value = ', defdmix, ' will be used'
	   end if 
	 end if 

c        spin-polarized electrodes: check for consistency of flags  
	 if (hexfound.and.nspin==1) then
	  sp_leads = .false.
	  hex_left = 0.0d0 ; hex_right = 0.0d0
	  print '(/,a,/,a,/)',
     &	   ' WARNING : exchange field parameters will not have effect',
     &	   '           in case of non-spin-polarized calculation!'
	 end if
	 if (locexfound.and.nspin==1) then
	  localexf = .false.
	  print '(/,a,/,a,/)',
     &	   ' WARNING : local exchange fields will not be active',
     &	   '           in case of non-spin-polarized calculation!'
	 end if
	 if (symmetrize.and.nspin==1) then
	  symmetrize = .false.
	  print '(/,a,/,a,/)',
     &	   ' WARNING : hamiltonian symmetrization flag does not have ',
     &	   '           effect in case of non-spin-polarized calculation!'
	 end if	
c	 if (sp_leads.and.tune_rsigma) then	
c	   print '(/,a,/,a,/)',
c     &	       ' active flag "$adjust_rsigma on" is not ',
c     &        ' supported in case of spin-polarized leads '
c	   stop ': transport module is terminated now'
c	 end if

c       ##########################################################################
c       check consistency of input in case of a "molecular subsystem" analysis ...
        if (mol_subsystem_found) then

         if ((.not.num_res_atoms_found).or.(.not.res_atom_types_found)) then
          print '(/,a,/,a,/)', 
     &     ' ERROR: <mol_subsystem> : reservoir is not specified !',
     &     '        please, check your <tcontrol> file'
          stop ': transport module is terminated now'
	 end if

c        checking for energy points if spectral function is required
         if (mol_spectr_func) then
          flag = (ener_found.and.estep_found.and.eend_found)
          if (.not.flag) then
           print '(/,a,/,a,/,a,/)',
     &     ' ERROR : <mol_subsystem>: you have not specified energy',
     &     '         points for the spectral function output',
     &     '         please, check your <tcontrol> file'
           stop ': transport module is terminated now'
          end if
         end if

c        further checking for consistency of job flags
         flag = (do_dmat_cycle.or.do_ldos.or.do_trans_calc.or.tune_rsigma)
         if (flag) then
	  print '(/,a,/,a)',
     &     ' ERROR : job flag "mol_subsystem" conflicts with other flags:',
     &	   '         "$densmat_cycle" or "$ldos" or "$transmission" or "$adjust_rsigma"'
	  print '(a,/)',  
     &	   '         please, switch them off'
	  stop ': transport module is terminated now'
         end if

c        when all checks are got through, set up the mol_subsystem flag >>>
         mol_subsystem = .true.
        end if ! mol_subsystem_found 
c       ##########################################################################

c       further checking for consistency of job flags
        if ((.not.do_dmat_cycle).and.(.not.do_trans_calc).and.(.not.do_ldos).and.(.not.do_ldos3d)
     & 	   .and.(.not.do_ldos3d_ewindow).and.(.not.tune_rsigma).and.(.not.mol_subsystem)) then
	 print '(/,a,/,a)',
     &      ' ERROR : all job flags "$densmat_cycle", "$ldos", $transmission", ',
     &	    '         "$adjust_rsigma", "mol_subsystem" are switched off'
	 print '(a,/)',  
     &	    '         please, specify correctly your choice'
	 stop ': transport module is terminated now'
        end if

c       ####################################################################
c       check consistency of input in case of ldau+u type of calculation >>>
        if (ldau_found) then

         if ((.not.num_ratoms_found).or.(.not.ratom_types_found)) then
          print '(/,a,/,a,/)', 
     &     ' <HUBBARD-U> : reservoir is not specified !',
     &     ' please, check your <tcontrol> file'
          stop ': transport module is terminated now'
         end if

         if (.not.u_j_found) then
          print '(/,a,/,a,/)', 
     &     ' <HUBBARD-U> : sreened Coulomb interaction parameter U-J is not specified !',
     &     ' please, check your <tcontrol> file'
          stop ': transport module is terminated now'
	 end if

c        checking for energy points if spectral function is required
         if (spectr_func) then
          flag = (ener_found.and.estep_found.and.eend_found)
          if (.not.flag) then
           print '(/,a,/,a,/)',
     &     ' <HUBBARD-U>: you have not specified energy points for spectral function output',
     &     ' please, check your <tcontrol> file'
           stop ': transport module is terminated now'
          end if
         end if

         if (no_leads) then

	  if (.not.calc_densmat0) then
           print '(/,a,/,a,/,a,/)', 
     &      ' <HUBBARD-U> : flag "$no_leads" is active',
     &      ' but orbitals occupation numbers are not found !',
     &      ' please, specify properly a group "$dens_matrix" in your <tcontrol> file'
           stop ': transport module is terminated now'

	  else if ((nspin.eq.1).and.(occ.ne.nelectr)) then
           print '(/,a,/,a,/)', 
     &      ' <HUBBARD-U> : number of occupied orbitals does not match number of electrons !',
     &      ' please, check your <tcontrol> file'
           stop ': transport module is terminated now'

          else if ((nspin.eq.1).and.( (nelectr/2)*2.ne.nelectr)) then 
           print '(/,a,/,a,/)', 
     &      ' <HUBBARD-U> : you have CLOSED SHELL calculation and ODD number of electrons !',
     &      ' are you sure about that ?  please, check your <tcontrol> file'
           stop ': transport module is terminated now'

          else if ((nspin.eq.2).and.((a_occ+b_occ).ne.nelectr)) then
           print '(/,a,a,/,a,/)', 
     &      ' <HUBBARD-U> : sum of occupied alpha & beta orbitals',
     &      ' do not match number of electrons !',
     &      ' please, check your <tcontrol> file'
           stop ': transport module is terminated now'
          end if
	 
 	 end if ! no_leads

c        when all checks are got through, set up the ldau-flag >>>
         ldau = .true.
        end if ! ldau_found
c       ###################################################################        

c       "no_leads" flag is active only together with hubbard-u
        if (no_leads.and.(.not.ldau)) then
         no_leads = .false.
         print '(/,1x,a)', 
     &	       'WARNING: flag "no_leads" is meaningless without <HUBBARD-U>'
        end if	
	
c       check for applicability of "save_ion" & "read_hion" flags
        if (save_hion.and.(.not.ldau)) then
	 save_hion = .false.
         print '(/,a,/,a,/)', 
     &     ' WARNINIG: "save_hion" flag will be ignored: ',
     &     '           it can be only active in case of LDA+U type of calculation' 
        end if

        if (read_hion.and.(.not.ldau)) then
	 read_hion = .false.
         print '(/,a,/,a,/)', 
     &     ' WARNINIG: "read_hion" flag will be ignored: ',
     &     '           it can be only active in case of LDA+U type of calculation' 
        end if

c       "effective core potentials" can be active in case of fhi-aims only
        if (ecp.and.ecp_found.and.(.not.aims_input)) then
          print '(/,1x,a)',
     &       'WARNING: $ecp flag has no effect in case of TURBOMOLE input'
        end if
        if (.not.aims_input) ecp = .false.

c       ecp flag is not supported in case on non-eq-density update, ldau & no_leads
        if (ecp.and.do_dmat_cycle) then
          print '(/,1x,a)',
     &       'WARNING: $ecp flag is yet supported if $densmat_cycle is active' 
          ecp = .false.
        end if
        if (ecp.and.no_leads) then
          print '(/,1x,a)',
     &       'WARNING: $ecp flag is not supported if $no_leads is active' 
          ecp = .false.
        end if
        if (ecp.and.ldau) then
          print '(/,1x,a)',
     &       'WARNING: $ecp flag is yet supported if $ldau is active' 
          ecp = .false.
        end if

c       updating iterations count and tcontrol
        if (do_dmat_cycle) then 

         if (iterfound) then
          iter = iter + 1
         else
          iter = 1
         end if
         call updatetcontrol(tcntrl,'$iter',iter,0.0d0)
        end if	

c       >>> checking, whether self-energy parameters have been read ...
        if (tune_rsigma) then
c        in this case, check for imaginary piece only
         if ( .not.( (found_isigma(1)).and.(found_isigma(2)).and.(found_isigma(3)) ) ) then
          print '(/,1x,a)',
     &       'WARNING: a model self-energy has not been specified properly,'
          print '(1x,a)',
     &       '         a default set of parameters (0.1H, 0.05H and 0.025H) will'
          print '(1x,a)',
     &       '         be used to parametrize an imaginary piece of the self-energy'

          do i = 1, deep
           rsigma(i) = 0.0d0
           isigma(i) = def_isigma/dble(2**(i-1))
          end do

         end if
        
        else ! flag $tune_rsigma is not active
c        in all other cases, check for all parameters
         if ( .not.( (found_isigma(1)).and.(found_isigma(2)).and.(found_isigma(3)).and. 
     &               (found_rsigma(1)).and.(found_rsigma(2)).and.(found_rsigma(3)) ) ) then
          print '(/,1x,a)',
     &       'ERROR : a model self-energy has not been specified properly'
          print '(1x,a,/)',
     &       '        please, check your <tcontrol> file'
          stop ': transport module will be terminated'
         end if

        end if ! tune_rsigma

c       check for efermi_flag ....
        if (fixed_efermi) then
c         if (.not.do_dmat_cycle) then
c          fixed_efermi = .false.
c          print '(/,1x,a)',
c     &       'WARNING: fixed_efermi flag is active only in case of'
c          print '(1x,a)',
c     &       '         iterative cycle for the density matrix'
c         else
          if (.not.efermi_guess) then
           print '(/,1x,a)',
     &       'ERROR : fixed_efermi flag is active but fermi energy guess is not given'
           print '(1x,a,/)',
     &       '        please, check your <tcontrol> file'
           stop ': transport module will be terminated'
          end if
c         end if
        end if
       
       END IF ! landauer

c      check for energy mesh if ballistic transmission is required
       if (landauer.or.do_trans_calc.and.(.not.do_scatt_wave_func)) then
        flag = ener_found.and.estep_found.and.eend_found
        if (.not.flag) then
         print '(/,1x,a)', 'WARNING : you have requested to output transmission function T(E)' 
         print '(1x,a)',   '          but have not specified an energy mesh properly (keywords ' 
         print '(1x,a)',   '          $ener, $estep, and $eend in the <tcontrol> file) ; ' 
         print '(1x,a)',   '          instead, only transmission at the Fermi energy will be computed' 
c        to be on the safe side, confirm conductance = .true.
         conductance = .true.
        end if
       end if

c      remove 'conductance' flag if scattering wave functions are required
       if (conductance.and.do_scatt_wave_func) then
        conductance = .false.
       end if

c      # take care about complex self-energy:
c      # it has NEGATIVE imaginary piece
       forall (i=1:deep) isigma(i) = -abs(isigma(i))
       forall (i=1:deep)  sigma(i) =  cmplx(rsigma(i),isigma(i),8)

c      take care about efermi (just for safety reasons)       
       if (.not.efermi_guess) efermi = 0.0d0
 
100    format(a,i4)

200    format(a,f10.6,1x,a)
c 201    format(a,f18.14,1x,a)
210    format(a,f12.8,1x,a)
211    format(a,f16.12,1x,a)
220    format(a,e12.6)
c 230  format(a,f7.3,a,f7.3,a)
c 240  format(a,f6.4)
300    format(a,e12.6,1x,a)
       end subroutine readtcontrol

c      ********** reading coord-file ********** 
       subroutine split_line(st,words,nw)
c       ***********************************************************
c       split a line 'st' into different words separated by spaces,
c       'nw' counts number of words found
c       ***********************************************************
        character(*) st
        character(*) words(*)

        integer j, wb, we, nw, stlength
        stlength = len_trim(st)
        nw = 0 ; wb = 0 ; we = 0
        do j = 1, stlength
          if (st(j:j)==' ') then
             if (wb > 0) then
c             we have found a space, store the word
              nw = nw + 1
              words(nw) = st(wb:we)
              wb = 0
             end if
          else if (wb==0) then
c           beginning of the new word is found
            wb = j ; we = j
          else
c           just go further along current word
            we = we + 1
          end if
        end do
c       adding last word to the array >>>
        if (wb > 0) then
          nw = nw + 1
          words(nw) = st(wb:we)
        end if
       end subroutine split_line

       subroutine get_atom_data(st,r,atname)
c       ************************************************************
c       given line 'st' from coord-file,
c       extract a position 'r' of atom, and its full name 'atname',
c       including possible sorts (like 'c chain', 'c ring', etc ...)
c       ************************************************************
        character(*) st
        double precision  r(3)
        character(atsymbol_len) atname

        integer i, nw
        character(wordlength) words(10)
        character(2*wordlength) tmp_name

        read(st,*) (r(i),i=1,3)
        call split_line(st,words,nw)
        if (nw.eq.4) then
c        no extra information is specified after atomic name
         tmp_name = words(4)
         atname = tmp_name(1:atsymbol_len)
        else
c        take care about additional atom "sort" following atomic symbol
         tmp_name = trim(words(4))//' '//trim(words(5))
         atname = tmp_name(1:atsymbol_len)
        end if
       end subroutine get_atom_data

       subroutine get_aims_atom_data(st,r,atname)
c       ************************************************************
c       given line 'st' from coord-file,
c       extract a position 'r' of atom, and its full name 'atname',
c       including possible sorts (like 'c chain', 'c ring', etc ...)
c       ************************************************************
        character(*) st
        double precision  r(3)
        character(atsymbol_len) atname

        integer i, nw
        character(wordlength) words(10)
        character(2*wordlength) tmp_name

        read(st,*) tmp_name, (r(i),i=1,3)
	forall (i=1:3)
c        be careful here: in case of fhi-aims, 
c        positions of atoms in geometry.in a given in Angstroms
c        while we are using atomic units through the code 	
	 r(i) = r(i)/bohr_radius
	end forall
        call split_line(st,words,nw)
c       !!! cation: here slight modifications are done to
c                   accept fhi-aims format of the coord-file
        if (nw.eq.5) then
c        no extra information is specified after atomic name
         tmp_name = words(5)
         atname = tmp_name(1:atsymbol_len)
c        covert atom symbol to a lower case >>>
         call lower_case(atname)
	else
c        take care about additional atom "sort" following atomic symbol
c        covert atom symbol to a lower case >>>
         call lower_case(words(5))
         tmp_name = trim(words(5))//' '//trim(words(6))
         atname = tmp_name(1:atsymbol_len)
        end if
       end subroutine get_aims_atom_data

       subroutine grepatoms(infilename,outfilename)
c       *****************************************************
c       relevant to fhi-aims only !
c       searches for lines with 'atom' keyword in geometry.in
c       and exports these lines to a temporary file
c       *****************************************************
        character (len = *) :: infilename, outfilename
        integer       infile, outfile, ierr
c       character*256 chline
c       character*32  st1
        character*256 syscall

        infile = 41  ! geometry.in file
        open(infile,file=trim(infilename),status='old',action='read',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a,a)',
     &             ' can not open file "', trim(infilename), '" for reading'
         print '(a,/)',
     &    ' please, check your directory content or access rights'
         stop ': transport module is terminated now'
        end if

        outfile = 42  ! geometry.tmp file
        open(outfile,file=trim(outfilename),status='unknown',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)',
     &         ' can not open temporary file "', trim(outfilename)
         print *, 'please, check access rights of your directory'
         print *
         stop  ': transport module is terminated now'
        end if
        write(outfile,'(a)') '#aims.geometry.tmp'
        close(infile) ; close(outfile)

c       do while (.true.)
c         read(infile,fmt='(a)',end=999) chline
c         read(chline,*) st1
c         if (trim(st1)=='atom') then
c          write(outfile,'(a)') trim(chline)
c         end if
c       enddo
c999    continue
        syscall = 'grep atom '//trim(infilename)//' >> '//trim(outfilename)
        call system(trim(syscall))
        syscall = 'echo '//'"#end"'//' >> '//trim(outfilename)
        call system(trim(syscall))

       end subroutine grepatoms

       subroutine lower_case(word)
c      ****************************
c      convert a word to lower case
c      ****************************
        character (len = *) :: word
        integer :: i, ic, nlen, icdiff, ica, icz

        ica = ichar('A')
        icz = ichar('Z')
        icdiff = ichar('a') - ichar('A')
        nlen = len(word)
        do i = 1, nlen
         ic = ichar(word(i:i))
         if (ic >= ica .and. ic <= icz) word(i:i) = char(ic+icdiff)
	end do
       
       end subroutine lower_case 

       subroutine readcoord(crdfile)
c       *****************************
c       reads coord/geometry.in files
c       *****************************
        implicit none
        character(*), intent (in) :: crdfile 

        integer  ifile, tmpfile, ierr, iat, i, j, jmax, nw, isort
        logical  start_reading, flag, atfound
        double precision  rxyz(3)
	character(atsymbol_len), allocatable :: atlist(:)
	character(atsymbol_len)  atom_name, myatom, tmp_atom 
	character ch
        character(linelength) st
        character(wordlength) words(10)
	
        ifile = 20  ! coord file
        open(ifile,file=crdfile,status='old',action='read',iostat=ierr)
        if (ierr.ne.0) then
	  print '(/,a,a,a)', 
     &	           ' can not open file "', trim(crdfile), '" for reading'
          print '(a,/)',
     &     ' please, check your directory content or access rights'
          stop ': transport module is terminated now'
	end if 

        allocate(atom(1:num_atoms),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE readcoord]: <atom> allocation failure'
        end if

        allocate(atlist(1:num_atoms),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE readcoord]: <atlist> allocation failure'
        end if

        start_reading = .false.
        do while (.not.start_reading)
c        read(ifile,'(a1)') ch   ! read first symbol
         read(ifile,*) ch   ! read first symbol
         flag = (ch=='$').or.(ch=='#')
         if (.not.flag)  start_reading = .true.
        end do
        backspace(ifile)

        if (aims_input) then
	 print '(/,a)', ' === reading geometry file ==='
	else
         print '(/,a,a,a)', ' === reading file <', trim(crdfile), '> ==='
	end if 
	 
        do iat = 1, num_atoms
c        read(ifile,*) (rxyz(i),i=1,3), atom_name
         read(ifile,fmt='(a)',iostat=ierr) st
	 if (ierr.eq.-1) then
	   print '(/,1x,a)',     'ERROR: end of file is reached unexpectedly' 
	   print '(1x,a,a,a,/)', '       your "',trim(crdfile),'" file seems to be corrupted' 
	   stop ' : transport module is terminated now'
	 else if (ierr.gt.0) then
	   print '(/,1x,a)',     'reading of file resulted in error' 
	   print '(1x,a,a,a,/)', 'your "',trim(crdfile),'" file seems to be corrupted' 
	   stop ' : transport module is terminated now'
	 end if
         if (aims_input) then
           call get_aims_atom_data(st,rxyz,atom_name)
	 else
           call get_atom_data(st,rxyz,atom_name)
	 end if  
         atom(iat)%pos    = rxyz
         atom(iat)%symbol = atom_name
         atom(iat)%llead  = .false.
         atom(iat)%rlead  = .false.
	 atom(iat)%charge = get_charge(atom_name(1:2))
        end do
        close(ifile)

ccc     next we build up reference arrays <atlist> and <at_slist>
ccc     with different atomic types and sorts being found

	atom(1)%atype = 1 
ccc     just to be on the safe side
        forall (iat=1:num_atoms) atlist(iat) = 'xxxx'

	atlist(1) = atom(1)%symbol
        jmax = 1
	do iat = 2, num_atoms
         myatom = atom(iat)%symbol
         atfound = .false.
	 j = 0
	 do while ((.not.atfound).and.(j<jmax))
          j = j + 1
          tmp_atom = atlist(j)
	  if (myatom(1:2)==tmp_atom(1:2)) atfound=.true.	
	 end do
	 if (.not.atfound) then
	  jmax = jmax + 1
	  atlist(jmax) = myatom(1:2)
	  atom(iat)%atype = jmax 
	 else
	  atom(iat)%atype = j
	 end if
        end do

ccc     put a local array <atlist> to the global one <atom_type>
        num_atom_types = jmax
        allocate(atom_type(num_atom_types),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE readcoord]: <atom_type> allocation failure'
        end if
	atom_type(1:num_atom_types) = atlist(1:num_atom_types)

c       then make search over atom sorts, where we distinguish
c       between atomic attributes, e.g. "c chain" and "c ring"

        atom(1)%asort = 1
ccc     just to be on the safe side
        forall (iat=1:num_atoms) atlist(iat) = 'xxxx'

        atlist(1) = atom(1)%symbol
        jmax = 1
        do iat = 2, num_atoms
         myatom = atom(iat)%symbol
         atfound = .false.
         j = 0
         do while ((.not.atfound).and.(j<jmax))
          j = j + 1
          tmp_atom = atlist(j)
          if (myatom(1:atsymbol_len)==tmp_atom(1:atsymbol_len))
     &      atfound=.true.
         end do
         if (.not.atfound) then
          jmax = jmax + 1
          atlist(jmax) = myatom
          atom(iat)%asort = jmax
         else
          atom(iat)%asort = j
         end if
        end do

ccc     put a local array <atlist> to the global one <atom_type>
        num_atom_sorts = jmax
        allocate(atom_sort(num_atom_sorts),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE readcoord]: <atom_sort> allocation failure'
        end if
        atom_sort(1:num_atom_sorts) = atlist(1:num_atom_sorts)

        deallocate(atlist,stat=ierr)

c       build up an array with "reference" atomic names
        allocate(ref_sort(num_atom_sorts),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE readcoord]: <ref_sort> allocation failure'
        end if

        do isort = 1, num_atom_sorts
         tmp_atom = atom_sort(isort)
         call split_line(tmp_atom,words,nw)
         if (nw.eq.1) then
          ref_sort(isort) = tmp_atom
         else ! nw == 2
          ref_sort(isort) = trim(words(1))//'_'//trim(words(2))
         end if
c        print *, ' - isort ', isort, ': ', ref_sort(isort)
        end do

ccc     print out info >>>
        print '(/,a,i3)', 
     &	      ' number of atom types found: ', num_atom_types 

        print '(5x,a)', '-------------'     
        print '(5x,a)', ' atom   type '     
        print '(5x,a)', '-------------'     
        do j = 1, num_atom_types
         print '(7x,a2,2x,i4)', atom_type(j)(1:2), j
        end do     	
        print '(5x,a)', '-------------'     

ccc     print out info >>>
        if (num_atom_types.ne.num_atom_sorts) then
          print '(/,a,i3)',
     &       ' number of atom SORTS found: ', num_atom_sorts

          print '(5x,a)', '---------------------------'
          print '(5x,a)', '   ATOM             SORT '
          print '(5x,a)', '---------------------------'
          do j = 1, num_atom_sorts
           print '(7x,a,2x,i4)', atom_sort(j), j
          end do
          print '(5x,a)', '---------------------------'
        end if

        if (testing) then
         tmpfile = 21
         open(tmpfile,file='coord.tmp',status='unknown',iostat=ierr)
         if (ierr.ne.0) then
          stop ': can not open coord.tmp!'
         end if

         write(tmpfile,'(a)') '$coord.tmp'
         write(tmpfile,'(a,i4)') '#number of atoms :', num_atoms
         write(tmpfile,'(a,i4)') '#atom types      :', num_atom_types
         write(tmpfile,'(a,i4)') '#atom sorts      :', num_atom_sorts

	 do iat = 1, num_atoms
c          for testing, keep all atomic positions in a.u. !!!
c          if (aims_input) then
c	   write(tmpfile,'(3f22.14,6x,a2,i5,i5)') 
c     &	    (atom(iat)%pos(i)*bohr_radius,i=1,3), 
c     &      atom(iat)%symbol(1:2), atom(iat)%atype, atom(iat)%asort
c          else
 	   write(tmpfile,'(3f22.14,6x,a2,i5,i5)') 
     &	        (atom(iat)%pos(i),i=1,3), 
     &          atom(iat)%symbol(1:2), atom(iat)%atype, atom(iat)%asort
c          end if ! aims_input
	 end do
         write(tmpfile,'(a)') '$end'
         close(tmpfile)	
	end if
	
       end subroutine readcoord

       subroutine check_electron_number
c       ############################################################## 
c       -- this routine is introduced to avoid possible user's mistake
c          in setting overall electron number in <tcontrol> file
c       -- in case of full-electron treatment, it counts amount
c          of electrons based on periodic table of elements
c       -- relevant only for the fhi-aims based input 
c          not for TURBOMOLE, where pseudo-potentials are used !
c       ############################################################## 
        integer :: iat, out_charge = 0
	logical :: term = .true.

	do iat = 1, num_atoms
	  out_charge = out_charge + atom(iat)%charge
	  if (atom(iat)%charge .eq. -1) then
	   term = .false.
	   print '(/,1x,a,a,a)', 
     &	     'WARNING: can not define a charge for element "',trim(atom(iat)%symbol(1:2)),'"'
	  end if 
	end do
       
        if (term) then
c         >>> compare a found number with the input value 
          if (.not.nelfound) then
            nelectr = out_charge
            print '(/,1x,a,i5)', 'Overall amount of electrons in your system: ', nelectr
          else
c          >>> compare a found number with the input value 
	   if (out_charge == nelectr) then
            print '(/,1x,a,i5)', 'Overall amount of electrons in your system: ', nelectr
	   else
	    print '(/,1x,a,i5,a)', 'found ', out_charge, ' electrons in your system,'
	    print '(1x,a,i5,a)',   'while you have introduced ', nelectr, ' electrons'
	    print '(1x,a)',        'using a keyword "$nelectr" in <tcontrol>'
	    print '(/,1x,a)',      'are you sure you are doing it right ???'
	    print '(/,1x,a)',      'please, correct electron number or remove'
	    print '(1x,a,/)',      'the keyword "$nelectr" from your <tcontrol> file'
	    stop ': transport module will be terminated now ...'
	   end if
          end if  ! nelfound
	else
c 	   >>> here term = .false.
           if (nelfound) then
c	     >>> will use the input value, if exists, but a warning message will appear  
	     print '(/,1x,a)',    'WARNING: i was not able to count amount of electrons !'
             print '(1x,a,i5,a)', '         assume, your system contains ', nelectr, ' electrons'
	     print '(1x,a)',      '         as has been set using a keyword "$nelectr" in <tcontrol>'	 
	     print '(1x,a)',      '         anyway, check again carefully what you are doing ...'
	   else
c            >>> neither electron number is introduced, nor automatic procedure failed	   
             if (.not.landauer) then
c             here it is assumed that a user is doing advanced tasks, 
c             so that keyword $nelectr is recommended
              print '(/,1x,a)', 'ERROR :  i was not able to count amount of electrons !'
	      print '(1x,a)',   '         please, check your <geometry.in> file or, at least'
	      print '(1x,a,/)', '         use a keyword "$nelectr" in your <tcontrol> file'
	      stop ': transport module will be terminated now ...'
	     else
c             in case of basic $landauer flag, print out just an error message
              print '(/,1x,a)', 'ERROR :  i was not able to count amount of electrons !'
	      print '(1x,a,/)', '         please, check your <geometry.in> file'
	      stop ': transport module will be terminated now ...'
	     end if ! landauer
	   end if ! nelfound     
	end if ! term
	
       end subroutine check_electron_number

c      ********** reading hsource-file ************* 

       subroutine get_hsource_data(st,r,atname,hsrc)
c       **************************************************************
c       given a line 'st' from coord-file,
c       extract a position 'r' of atom, and its full name 'atname', 
c       and value 'hsrc [Hartree]' of the local exchange field on atom
c       **************************************************************  
        character(*), intent(in)              :: st
        double precision, intent (out)        :: r(3)
        character(atsymbol_len), intent (out) :: atname
        double precision, intent (out)        :: hsrc

        integer i, nw
        character(wordlength) words(10)
        character(2*wordlength) tmp_name

        call split_line(st,words,nw)
        if (nw.ne.5) then
c        a non-consistent format of hsource-file
         print *
         stop '[SUB. get_hsource_data]: wrong format of hsource-file'
        else 
         read(st,*) (r(i),i=1,3)
         tmp_name = words(4)
         atname = tmp_name(1:atsymbol_len)
         read(words(5),*) hsrc     
         write(*,fmt='(3f22.14,6x,a2,f12.6)') (r(i),i=1,3),atname,hsrc
	end if
       end subroutine get_hsource_data

       subroutine readhsource(hsrcfile)
        implicit none
        character(*), intent (in) :: hsrcfile 

        integer  ifile, ierr, iat
        logical  start_reading, flag        
	double precision  rxyz(3), exf
	character(atsymbol_len) atom_name 
        character ch
	character(linelength)   st

        ifile = 21  ! hsource file
        open(ifile,file=hsrcfile,status='old',action='read',iostat=ierr)
        if (ierr.ne.0) then
	 print '(/,a,a,a)', 
     &	           ' can not open file "', trim(hsrcfile), '" for reading'
         print '(a,/)',
     &    ' please, check your directory content or access rights'
         stop ': transport module is terminated now'
        end if

        start_reading = .false.
        do while (.not.start_reading)
c        read(ifile,'(a1)') ch   ! read first symbol
         read(ifile,*) ch   ! read first symbol
         flag = (ch=='$').or.(ch=='#')
         if (.not.flag)  start_reading = .true.
        end do
        backspace(ifile)

        print '(/,a,a,a,/)', ' === reading file <', hsrcfile, '> ==='
        do iat = 1, num_atoms
         read(ifile,'(a)') st
         call get_hsource_data(st,rxyz,atom_name,exf)
c        atom(iat)%pos    = rxyz
c        atom(iat)%symbol = atom_name
c        atom(iat)%llead = .false.
c        atom(iat)%rlead = .false.
         atom(iat)%exfield = exf
        end do
        close(ifile)
       end subroutine readhsource

c      ********** reading mos-file ********** 
       subroutine readmos(ispin, mosfile)
        implicit none
        integer, intent (in) ::     ispin  ! =1 or 2
        character(*), intent (in) :: mosfile

c       length of string in the mos-file until the eigenvalue is met
	integer, parameter :: stlen=26  
	character(stlen) evstr ! all characters before the eigenvalue 

	integer ifile, tmpfile, ierr, n, j, mincol, maxcol, tmp_nsaos
        integer, parameter :: monitor = 50          
        double precision en
        logical   start_reading, flag
	character ch
	
        character(80) headline, cut_headline
	character(16) fmt_mos, fmt_ener
	integer       ifmt 

        ifile = 30  ! external mos-file
	open(ifile,file=mosfile,status='old',action='read',iostat=ierr)
        if (ierr.ne.0) then
	 print '(/,a,a,a)',
     &	          ' can not open file "', trim(mosfile), '" for reading'
	 print *, 
     &	  'please, check your directory content or access rights'
         print *
	 stop ': transport module is terminated now'
        end if

c       read header line: get i/o format of the mos-file
	start_reading = .false.
        do while (.not.start_reading)
c        read(ifile,'(a1)') ch   ! read first symbol
	 read(ifile,*) ch   ! read first symbol
	 flag = (ch=='#')
	 if (.not.flag)  start_reading = .true.
	end do
	backspace(ifile)
        read(ifile,'(a80)') headline

c       searching for substring like: format(4d20.14) 
c       other choices of format could be : format(4d24.18), format(4d26.20), etc ...
        ifmt  = index(headline,'format')
	fmt_mos  = headline(ifmt+6:ifmt+14)  ! e.g. (4d20.14)
	fmt_ener = headline(ifmt+8:ifmt+13)  ! e.g.    20.14

c       print out info message on found format :
        write(6,'(/,1x,a)') trim(mosfile)//' : ASCII file, format'//trim(fmt_mos) 

c       after format is indentified, proceed with checking nsaos ... 

c       searching for substring like: nsaos= ... 
c       other choices of format could be : format(4d24.18), format(4d26.20), etc ...
        rewind(ifile)
        start_reading = .false.
        do while (.not.start_reading)
c        read(ifile,'(a1)') ch   ! read first symbol
	 read(ifile,*) ch   ! read first symbol
	 flag = (ch=='$').or.(ch=='#')
	 if (.not.flag)  start_reading = .true.
	end do
	backspace(ifile)
        read(ifile,'(a80)') headline

        ifmt  = index(headline,'nsaos=')
	cut_headline  = headline(ifmt+6:len(headline)) 
        read(cut_headline,*) tmp_nsaos
        if (nsaos.ne.tmp_nsaos) then
         print '(/,1x,a,i6)', 'ERROR: internal conflict: you have specified nsaos= ', nsaos
         print '(1x,a,a,a)',  '       in your <tcontrol> file, while your <',trim(mosfile),'> file'
         print '(1x,a,i6,a)', '       tells "nsaos" should be ', tmp_nsaos, '  : i am confused'
         print '(/,1x,a,/)',  '       please, check carefully what you are doing ...'
         stop ' : transport module will be terminated now'
	end if

c       after check is done, proceed with reading data ... 
        rewind(ifile)
        start_reading = .false.
        do while (.not.start_reading)
c        read(ifile,'(a1)') ch   ! read first symbol
	 read(ifile,*) ch   ! read first symbol
	 flag = (ch=='$').or.(ch=='#')
	 if (.not.flag)  start_reading = .true.
	end do
	backspace(ifile)

        if (testing) then
         tmpfile = 31
         if (nspin == 1) then
	  open(tmpfile,file='mos.tmp',status='unknown',iostat=ierr)
          if (ierr.ne.0) then
	   stop ': can not open mos.tmp!'
          end if
	  write(tmpfile,'(a)') '$mos.tmp'
	 else  ! spin-polarized case
          if (ispin == 1) then	 
	   open(tmpfile,file='alpha.tmp',status='unknown',iostat=ierr)
           if (ierr.ne.0) then
	    stop ': can not open alpha.tmp!'
           end if
	   write(tmpfile,'(a)') '$alpha.tmp'
	  else
	   open(tmpfile,file='beta.tmp',status='unknown',iostat=ierr)
           if (ierr.ne.0) then
	    stop ': can not open beta.tmp!'
           end if
	   write(tmpfile,'(a)') '$beta.tmp'
	  end if 
	 end if
	 write(tmpfile,'(a,i5)') '#nsaos=', nsaos
        end if ! testing

        print '(a,a,a,$)', ' importing file <', mosfile, '>  : | '

c       reading nsaos molecular orbitals
	do n = 1, nsaos

         if ( mod(n,monitor)==0 ) print '(a,$)', '='  ! monitor
c        read(ifile,'(a26,d20.14)') evstr, en
c        print '(i4,a,d20.14)', n, ' : ', en 
         read(ifile,fmt='(a26,'//trim(fmt_ener)//')',iostat=ierr) evstr, en
	 if (ierr.eq.-1) then
	   print '(/,/,1x,a)',   'ERROR: end of file is reached unexpectedly' 
	   print '(1x,a,a,a,/)', '       your "',trim(mosfile),'" file seems to be corrupted' 
	   stop ' : transport module is terminated now'
	 else if (ierr.gt.0) then
	   print '(/,/,1x,a)',   'reading of file resulted in error' 
	   print '(1x,a,a,a,/)', 'your "',trim(mosfile),'" file seems to be corrupted' 
	   stop ' : transport module is terminated now'
	 end if
	 mo_en(n,ispin)=en

c        if (testing) then
c         write(tmpfile,'(a26,d20.14)') evstr, mo_en(n,ispin)
         if (testing) then  
	   write(tmpfile,'(i6,4x,e'//fmt_ener(2:6)//')') n, mo_en(n,ispin)
c          write(tmpfile,'(i6,4x,e20.14)') n, mo_en(n,ispin)
         end if
	 maxcol = 0
         do while (maxcol<nsaos)
          mincol = maxcol+1
          maxcol = min(maxcol+4,nsaos)
c          read(ifile,'(4d20.14)') (mo_coeff(j,n,ispin),j=mincol,maxcol)
cccc       write(tmpfile,'(4(d20.14,1x))') 
c          write(tmpfile,'(4(d20.14))') 
c     &	                          (mo_coeff(j,n,ispin),j=mincol,maxcol)
           read(ifile,fmt=trim(fmt_mos),iostat=ierr) (mo_coeff(j,n,ispin),j=mincol,maxcol)
	   if (ierr.eq.-1) then
	     print '(/,/,1x,a)',   'ERROR: end of file is reached unexpectedly' 
	     print '(1x,a,a,a,/)', '       your "',trim(mosfile),'" file seems to be corrupted' 
	     stop ' : transport module is terminated now'
	   else if (ierr.gt.0) then
	     print '(/,/,1x,a)',   'reading of file resulted in error' 
	     print '(1x,a,a,a,/)', 'your "',trim(mosfile),'" file seems to be corrupted' 
	     stop ' : transport module is terminated now'
	   end if
         end do ! do_while
        end do ! nsaos
	
	if (testing) then
 	 write(tmpfile,'(a)') '$end'
	 close(tmpfile)
        end if
	
        print *, '| done'
        close(ifile)

       end subroutine readmos

c      importing mos/alpha/beta files ...
       subroutine import_mos
        use globalvars
     	integer       ierr

        print '(/,a,/)', ' === importing mos/alpha/beta files ==='
        print '(a,i1)',  ' nspin = ', nspin
        print '(a,i4)',  ' nsaos = ', nsaos
	
        allocate(mo_en(nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
	 stop ': impossible to allocalte <mo_en>'
	end if
c 	print '(/,a)', ' <mo_en>:    successful allocation'

        allocate(mo_coeff(nsaos,nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
	 stop ': impossible to allocalte <mo_coeff>'
	end if
c 	print *, '<mo_coeff>: successfull allocation'

ccc     call readmos
        if (nspin==1) then 
         call readmos(1,trim(mos_file_name))
        else     
         call readmos(1,trim(alpha_file_name)) ! ispin=1 : alpha
         call readmos(2,trim(beta_file_name))  ! ispin=2 : beta
        end if
   
       end subroutine import_mos

c      ********** reading basis-file **********
       subroutine err_message(tmp_dim,dim0)
        implicit none
	integer, intent (in) :: tmp_dim, dim0
	
	print *
	print 5050, 
     &	 ' uups! ... used parameter "orbdim" =', dim0, 
     &   ' is too small'
        print 5050, 
     &    ' please, recompile the code, using "orbdim" >=', tmp_dim, 
     &	  ' in module <globalvars.f>'
        print *
	stop ': transport module is terminated now'  
5050   format(a,i2,a)	
       end subroutine err_message

       subroutine readbasis(bssfile)
        implicit none
        character(*), intent (in) :: bssfile

        type tmp_cgto
	 integer           ngto
	 double precision  icoeff(contr)
	 double precision  xi(contr)
	end type

        integer bfile, tmpfile, ierr, itype, ncontr, 
     &         	n, i, iorb, j, lm, iat, iatype
c       ! lmax0 (=3) is defined in <globalvars.f>
        integer  nfunc(lmax0+1), n_ve_func(lmax0+1), tmp_count
c	! for current atom type: amount of s, p, d, f basis functions
c       ! in case of fhi-aims code, formally g & h functions are needed
	integer, allocatable :: s_count(:), p_count(:), 
     &	                        d_count(:), f_count(:),
     &	                        g_count(:), h_count(:)
	logical, allocatable :: attype_found(:)
	logical    exitnow, flag, findstar, readall
	character  myatom*2, symbols*4, ch, ltype

        type(tmp_cgto), allocatable :: s_basfunc(:,:),
     &  	                       p_basfunc(:,:), 
     &                                 d_basfunc(:,:), 
     &                                 f_basfunc(:,:),
     &                                 g_basfunc(:,:), 
     &                                 h_basfunc(:,:) 
        double precision icoeff(contr), xi(contr),
     &                   tmp_coeff, tmp_exp, tmp_norm, tmp_norm3	

        character(128) syscall
     
        allocate(s_basfunc(num_atom_types,orbdim), stat=ierr) 
        if (ierr.ne.0) then
          print *
          stop '[SUBR. readbasis]: <s_basfunc> allocation failure'
        end if

        allocate(p_basfunc(num_atom_types,orbdim), stat=ierr) 
        if (ierr.ne.0) then
          print *
          stop '[SUBR. readbasis]: <p_basfunc> allocation failure'
        end if

        allocate(d_basfunc(num_atom_types,orbdim), stat=ierr) 
        if (ierr.ne.0) then
          print *
          stop '[SUBR. readbasis]: <d_basfunc> allocation failure'
        end if

        allocate(f_basfunc(num_atom_types,orbdim), stat=ierr) 
        if (ierr.ne.0) then
          print *
          stop '[SUBR. readbasis]: <f_basfunc> allocation failure'
        end if

c       >>> fhi-aims related modifications  
        allocate(g_basfunc(num_atom_types,orbdim), stat=ierr) 
        if (ierr.ne.0) then
          print *
          stop '[SUBR. readbasis]: <g_basfunc> allocation failure'
        end if

        allocate(h_basfunc(num_atom_types,orbdim), stat=ierr) 
        if (ierr.ne.0) then
          print *
          stop '[SUBR. readbasis]: <h_basfunc> allocation failure'
        end if
c       <<< fhi-aims related modifications 
		
        bfile = 40  ! basis file
        open(bfile,file=bssfile,status='old',action='read',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a,a)',
     &	          ' can not open file "', trim(bssfile), '" for reading'
         print *,
     &    'please, check your directory content or access rights'
         print *
         stop ': transport module is terminated now'
        end if

        if (.not.aims_input) then
	 print '(/,a,a,a,/)', ' === reading file <', bssfile, '> ==='
	else
	 print '(a,/)', '  collecting information on basis functions ...'
	end if 

c       search cycle over all different atom types
        allocate(attype_found(num_atom_types),stat=ierr)
        if (ierr.ne.0) then
          print *
          stop '[SUBR. readbasis]: <attype_found> allocation failure'
        end if

        attype_found = .false.
c        forall (itype=1:num_atom_types) 
c	 attype_found(itype) = .false.
c	end forall
	 
        do itype = 1, num_atom_types
	 myatom = atom_type(itype)(1:2)
         rewind(bfile)
 	 flag = .false.
         do while (.not.flag)
          read(bfile,*) symbols
          flag = (symbols(1:2)==myatom).or.
     &	   (symbols=='$end').or.(symbols=='$ecp')
	 end do
          
	 if (symbols(1:2)==myatom) then 
           print *, ' -- found atom type: ', myatom
           attype_found(itype) = .true.
c         else 
c           attype_found(itype) = .false.
	 end if
	end do	

c       check whether we have found all atom types in a basis file
        exitnow = .false.
        do itype = 1, num_atom_types
         if (.not.attype_found(itype)) then
	  myatom = atom_type(itype)(1:2) 
          print *, 'info for atom "', trim(myatom), '" is not found!'
          exitnow = .true.
	 end if
        end do
        
	if (exitnow) then
         close(bfile)
	 print '(a,/)', ' please, check your basis-file'
         stop ': transport module is terminated now'
	end if
        deallocate(attype_found,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,/,a)',
     &   ' [SUBR. readbasis]: <attype_found> deallocation failure ',
     &   ' ... will proceed further anyway'
        end if
ccc     all checks are OK
ccc     we are ready to read exponents and coeffs 

        allocate(n_basis_func(1:num_atom_types),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBR. readbasis]: <n_basis_func> allocation failure'
        end if

c       to manage "effective core potentials" and fhi-aims input, 
c       we need to allocate additional array, counting only 
c       basis functions associated with valence electrons
        allocate(n_val_basis_func(1:num_atom_types),stat=ierr)
        if (ierr.ne.0) then
          print *
          stop '[SUBR. readbasis]: <n_val_basis_func> allocation failure'
        end if
        allocate(ncore(1:num_atom_types),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBR. readbasis]: <ncore> allocation failure'
        end if

        allocate(s_count(1:num_atom_types),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBR. readbasis]: <s_count> allocation failure'
        end if
        allocate(p_count(1:num_atom_types),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBR. readbasis]: <p_count> allocation failure'
        end if
        allocate(d_count(1:num_atom_types),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBR. readbasis]: <d_count> allocation failure'
        end if
        allocate(f_count(1:num_atom_types),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBR. readbasis]: <f_count> allocation failure'
        end if

c       >>> fhi-aims related modifications 
        allocate(g_count(1:num_atom_types),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBR. readbasis]: <g_count> allocation failure'
        end if
        allocate(h_count(1:num_atom_types),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBR. readbasis]: <h_count> allocation failure'
        end if
c       <<< fhi-aims related modifications 

	s_count = 0 
	p_count = 0
        d_count = 0 
	f_count = 0
c       >>> fhi-aims related modifications 
        g_count = 0 
	h_count = 0
c       <<< fhi-aims related modifications 

        if (testing) then
         tmpfile = 32
	 open(tmpfile,file='basis.tmp',status='unknown',iostat=ierr)
         if (ierr.ne.0) then
          stop ': can not open basis.tmp!'
         else
          write(tmpfile,fmt='(a)') '$basis.tmp'
         end if
	end if

!       in case of "ecp" fill up the array 'ncore'
        do itype = 1, num_atom_types
         myatom = atom_type(itype)(1:2)
         ncore(itype) = get_ncore(get_charge(myatom)) 
        end do
	
        do itype = 1, num_atom_types

	 myatom = atom_type(itype)(1:2)
         rewind(bfile)
 	 flag = .false.
         do while (.not.flag)
          read(bfile,*) symbols
          flag = (symbols(1:2)==myatom)
	 end do
	 findstar = .false.
	 do while (.not.findstar) 
	  read(bfile,*) ch
	  if (ch=='*')  findstar = .true.
	 end do
         
	 forall (i=1:lmax0+1) nfunc(i)=0
c        ! case of "ecp" and fhi-aims input
	 forall (i=1:lmax0+1) n_ve_func(i)=0
	 
         readall = .false.	
c        ! now we are reading basis of given atom (myatom)

	 print '(/,a,a,a)', ' basis of atom "', trim(myatom), '" '
	 if (testing) then
	  write(tmpfile,fmt='(a,a,a)') 
     &	             '$basis of atom "', trim(myatom), '" '
         end if 

         do while (.not.readall)    

	  read(bfile,*) ch
	  if (ch=='*') then 
	   readall = .true. ; cycle
	  end if 

	  backspace(bfile)
	  read(bfile,*) ncontr, ltype
	  if (ncontr>contr) then
	   print '(/,a,i2,a,i2)', ' uups!... it happend that ncontr =',
     &	                        ncontr, ' > contr =', contr
	   print '(a)',' please, recompile the code, by setting up' 
	   print '(a,i2,a,/)', ' "contr" >=', ncontr, 
     &	                       ' in module <globalvars.f>'
	   stop ': transport module is terminated now'
	  end if
	   
	  if (testing) then
c	   print *, ncontr, ltype
 	   write(tmpfile,*) ncontr, ltype
 	  end if
	  
c	  flag = (ltype=='s').or.(ltype=='p').or.
c     &	         (ltype=='d').or.(ltype=='f')
c
c       >>> fhi-aims related modifications 
	  flag = (ltype=='s').or.(ltype=='p').or.(ltype=='d').or.
     &	         (ltype=='f').or.(ltype=='g').or.(ltype=='h')
c       <<< fhi-aims related modifications 
	  if (.not.flag) then
	    close(bfile)
	    print '(/,a,a)', ' [readbasis]: ', 
     &            'only s, p, d, f, g and h-type functions are accepted'
            print '(a,a,a,/)', ' atom "', trim(myatom), 
     &	                       '" has caused a problem'
	    stop ': transport module is terminated now'
	  end if

	  forall(n=1:contr)
	   icoeff(n)= 0.0d0; xi(n) = 0.0d0
	  end forall 
	  do n = 1, ncontr
           read(bfile,*) tmp_exp, tmp_coeff
	   if (testing) then
c           print *, tmp_exp, tmp_coeff
 	    write(tmpfile,*) tmp_exp, tmp_coeff
 	   end if
           if (.not.aims_input) then
ccc         next we rescale coefficients ... (ref: ferdinand evers)
ccc         approved to be correct (alexej bagrets) 
c <old>	    tmp_norm = sqrt(0.5d0*pi/tmp_exp)
	    tmp_norm = sqrt(1.0d0/tmp_exp)
	    tmp_norm3 = tmp_norm*tmp_norm*tmp_norm 
            if (ltype == 'p') then
c <old>	     tmp_norm3 = tmp_norm3/(2.0d0*tmp_exp)
	     tmp_norm3 = tmp_norm3/(tmp_exp)
	    else if (ltype=='d') then
c <old>	     tmp_norm3 = tmp_norm3*3.0d0/(4.0d0*tmp_exp*tmp_exp)  	 
	     tmp_norm3 = tmp_norm3/(tmp_exp*tmp_exp)  	 
	    else if (ltype=='f') then
	     tmp_norm3 = tmp_norm3/(tmp_exp*tmp_exp*tmp_exp)  	 
	    else if (ltype=='g') then
	     tmp_norm3 = tmp_norm3/(tmp_exp*tmp_exp*tmp_exp*tmp_exp)  	 
	    else if (ltype=='h') then
	     tmp_norm3 = tmp_norm3/(tmp_exp*tmp_exp*tmp_exp*tmp_exp*tmp_exp)  	 
	    end if
	   else
	     tmp_norm3 = 1.0d0
	   end if    
	   xi(n) = tmp_exp
	   icoeff(n) = tmp_coeff/sqrt(tmp_norm3)
ccc        end testing	 
	  end do  ! ncontr

          select case (ltype)
	    case ('s')
	      nfunc(1) = nfunc(1) + 1
	      if ( ecp .and.(.not.( (int(xi(1))<=ncore(itype)).and.(int(icoeff(1))==2) )) ) then
	       n_ve_func(1) = n_ve_func(1) + 1
	      end if 
	      s_count(itype) = s_count(itype) + 1
	      tmp_count = s_count(itype)
	      if (tmp_count>orbdim) call err_message(tmp_count,orbdim)
	      s_basfunc(itype,tmp_count)%ngto = ncontr
	      s_basfunc(itype,tmp_count)%icoeff = icoeff
	      s_basfunc(itype,tmp_count)%xi = xi
	    case ('p')
	      nfunc(2) = nfunc(2) + 3
	      if ( ecp .and.(.not.( (int(xi(1))<=ncore(itype)).and.(int(icoeff(1))==2) )) ) then
	       n_ve_func(2) = n_ve_func(2) + 3
	      end if 
	      p_count(itype) = p_count(itype) + 1
	      tmp_count = p_count(itype)
	      if (tmp_count>orbdim) call err_message(tmp_count,orbdim)
	      p_basfunc(itype,tmp_count)%ngto = ncontr
	      p_basfunc(itype,tmp_count)%icoeff = icoeff
	      p_basfunc(itype,tmp_count)%xi = xi
	    case ('d')
	      nfunc(3) = nfunc(3) + 5
	      if ( ecp .and.(.not.( (int(xi(1))<=ncore(itype)).and.(int(icoeff(1))==2) )) ) then
	       n_ve_func(3) = n_ve_func(3) + 5
	      end if 
	      d_count(itype) = d_count(itype) + 1
	      tmp_count = d_count(itype)
	      if (tmp_count>orbdim) call err_message(tmp_count,orbdim)
	      d_basfunc(itype,tmp_count)%ngto = ncontr
	      d_basfunc(itype,tmp_count)%icoeff = icoeff
	      d_basfunc(itype,tmp_count)%xi = xi
	    case ('f')
	      nfunc(4) = nfunc(4) + 7
	      if ( ecp .and.(.not.( (int(xi(1))<=ncore(itype)).and.(int(icoeff(1))==2) )) ) then
	       n_ve_func(4) = n_ve_func(4) + 7
	      end if 
	      f_count(itype) = f_count(itype) + 1
	      tmp_count = f_count(itype)
	      if (tmp_count>orbdim) call err_message(tmp_count,orbdim)
	      f_basfunc(itype,tmp_count)%ngto = ncontr
	      f_basfunc(itype,tmp_count)%icoeff = icoeff
	      f_basfunc(itype,tmp_count)%xi = xi
c           >>> fhi-aims related modifications
	    case ('g')
	      nfunc(5) = nfunc(5) + 9
	      if ( ecp .and.(.not.( (int(xi(1))<=ncore(itype)).and.(int(icoeff(1))==2) )) ) then
	       n_ve_func(5) = n_ve_func(5) + 9
	      end if 
	      g_count(itype) = g_count(itype) + 1
	      tmp_count = g_count(itype)
	      if (tmp_count>orbdim) call err_message(tmp_count,orbdim)
	      g_basfunc(itype,tmp_count)%ngto = ncontr
	      g_basfunc(itype,tmp_count)%icoeff = icoeff
	      g_basfunc(itype,tmp_count)%xi = xi
	    case ('h')
	      nfunc(6) = nfunc(6) + 11
	      if ( ecp .and.(.not.( (int(xi(1))<=ncore(itype)).and.(int(icoeff(1))==2) )) ) then
	       n_ve_func(6) = n_ve_func(6) + 11
	      end if 
	      h_count(itype) = h_count(itype) + 1
	      tmp_count = h_count(itype)
	      if (tmp_count>orbdim) call err_message(tmp_count,orbdim)
	      h_basfunc(itype,tmp_count)%ngto = ncontr
	      h_basfunc(itype,tmp_count)%icoeff = icoeff
	      h_basfunc(itype,tmp_count)%xi = xi
c           <<< fhi-aims related modifications
	  end select

	 end do  ! readall: finished reading basis of given atom

         n_basis_func(itype) = 0
         do i = 1, lmax0+1 
          n_basis_func(itype) = n_basis_func(itype) + nfunc(i)
         end do
         print '(a,i2)', ' number of basis functions =  ',
     &            n_basis_func(itype)
         if (testing) then
          write(tmpfile,fmt='(a,i2)') 
     &            '#number of basis functions =  ', n_basis_func(itype)
	 end if

         if (ecp) then
          n_val_basis_func(itype) = 0
          do i = 1, lmax0+1 
           n_val_basis_func(itype) = n_val_basis_func(itype) + n_ve_func(i)
          end do
          print '(a,i2)', ' number of valence orbitals = ',
     &            n_val_basis_func(itype)
          if (testing) then
           write(tmpfile,fmt='(a,i2)') 
     &            '#number of valence orbitals = ', n_val_basis_func(itype)
	  end if
	 end if

	end do  ! itype	
        print *
	close(bfile)
        
        num_at_gto = 0
        do itype = 1, num_atom_types
         if (n_basis_func(itype) > num_at_gto) then
          num_at_gto = n_basis_func(itype)
         end if
        end do

        if (ecp) then
         num_val_gto = 0
         do itype = 1, num_atom_types
          if (n_val_basis_func(itype) > num_val_gto) then
           num_val_gto = n_val_basis_func(itype)
          end if
         end do
        end if
         
	if (testing) then
         write(tmpfile,fmt='(a,i2)') '$num_at_gto = ', num_at_gto
         if (ecp) write(tmpfile,fmt='(a,i2)') '$num_val_gto = ', num_val_gto
	 write(tmpfile,fmt='(a)') '$end'
	 close(tmpfile)
        end if

c       building up atomic orbitals array :
c       aos(num_atom_type,num_at_gto)
        allocate(aos(num_atom_types,num_at_gto), stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE readbasis]: <aos> allocation failure'
        end if
       
ccc      contracted GTOs info.......
ccc       type cgto.
ccc        integer           ngto           ! degree of contraction == number of GTO..
ccc        double precision  icoeff(contr)  ! expansion coeffs.
ccc        double precision  xi(contr)      ! exponents.
ccc        integer           lm             ! lm = 1 2  3  4   5 6 7 8 9
ccc                                         !      s px py pz  d  ...  d
ccc        logical           valence        ! .false. --> core, .true. --> valence state.
ccc       end type
ccc     >>> to start with, set all numbers/values to 0
        do itype = 1, num_atom_types
         do iorb = 1, num_at_gto
          aos(itype,iorb)%ngto = 0
          aos(itype,iorb)%icoeff = 0.0d0
          aos(itype,iorb)%xi = 0.0d0
          aos(itype,iorb)%lm = 0
          aos(itype,iorb)%valence = .false.
         end do  
        end do
   
        do itype = 1, num_atom_types
         iorb = 0 
c        s-functions 
         tmp_count = s_count(itype)
         do j = 1, tmp_count
	  iorb = iorb + 1
	  aos(itype,iorb)%ngto   = s_basfunc(itype,j)%ngto
	  aos(itype,iorb)%icoeff = s_basfunc(itype,j)%icoeff
	  aos(itype,iorb)%xi     = s_basfunc(itype,j)%xi
	  aos(itype,iorb)%lm     = 1    ! s type
         end do
c        p-functions
         tmp_count = p_count(itype)
         do j = 1, tmp_count
          do lm = 2, 4
	   iorb = iorb + 1
	   aos(itype,iorb)%ngto   = p_basfunc(itype,j)%ngto
	   aos(itype,iorb)%icoeff = p_basfunc(itype,j)%icoeff
	   aos(itype,iorb)%xi     = p_basfunc(itype,j)%xi
	   aos(itype,iorb)%lm     = lm  ! p type
	  end do
	 end do
c        d-functions
         tmp_count = d_count(itype)
         do j = 1, tmp_count
          do lm = 5, 9
	   iorb = iorb + 1
	   aos(itype,iorb)%ngto   = d_basfunc(itype,j)%ngto
	   aos(itype,iorb)%icoeff = d_basfunc(itype,j)%icoeff
	   aos(itype,iorb)%xi     = d_basfunc(itype,j)%xi
	   aos(itype,iorb)%lm     = lm  ! d type
	  end do
         end do        
c        f-functions
         tmp_count = f_count(itype)
         do j = 1, tmp_count
          do lm = 10, 16
	   iorb = iorb + 1
	   aos(itype,iorb)%ngto   = f_basfunc(itype,j)%ngto
	   aos(itype,iorb)%icoeff = f_basfunc(itype,j)%icoeff
	   aos(itype,iorb)%xi     = f_basfunc(itype,j)%xi
	   aos(itype,iorb)%lm     = lm  ! f type
	  end do
	 end do

c        >>> fhi-aims related modifications
c        g-functions
         tmp_count = g_count(itype)
         do j = 1, tmp_count
          do lm = 17, 25
	   iorb = iorb + 1
	   aos(itype,iorb)%ngto   = g_basfunc(itype,j)%ngto
	   aos(itype,iorb)%icoeff = g_basfunc(itype,j)%icoeff
	   aos(itype,iorb)%xi     = g_basfunc(itype,j)%xi
	   aos(itype,iorb)%lm     = lm  ! g type
	  end do
	 end do
c        h-functions
         tmp_count = h_count(itype)
         do j = 1, tmp_count
          do lm = 26, 36
	   iorb = iorb + 1
	   aos(itype,iorb)%ngto   = h_basfunc(itype,j)%ngto
	   aos(itype,iorb)%icoeff = h_basfunc(itype,j)%icoeff
	   aos(itype,iorb)%xi     = h_basfunc(itype,j)%xi
	   aos(itype,iorb)%lm     = lm  ! h type
	  end do
	 end do
c        <<< fhi-aims related modifications

c        internal check >>>
         if (iorb==n_basis_func(itype)) then
	  print '(a,i3,a)', ' atom type =', itype, ': consistency check is OK'
         else
	  print '(a,i3,a,/)', ' atom type =', itype, ': consistency check FAILS!'
          stop 'something goes wrong here; stop now '
	 end if	
	 
	end do  ! itype; array aos

c       in case of fhi-aims input and "ecp" 
c       split orbitals to core and valence ones >>>
        do itype = 1, num_atom_types
         do iorb = 1, n_basis_func(itype)
          if (.not.ecp) then
            aos(itype,iorb)%valence = .true.
          else if ( (aos(itype,iorb)%icoeff(1)==2).and.(aos(itype,iorb)%xi(1)<=ncore(itype)) ) then
c           we picked up core orbital
            aos(itype,iorb)%valence = .false.
c           print *, 'itype = ', itype, ' iorb =', iorb, ': core orbital'
          else
c           we picked up valence orbital
            aos(itype,iorb)%valence = .true.    
          end if
         end do  
        end do

c       counting number of valence orbitals in the system, 
c       i.e. reduced dimension of hilbert space 
        nvalence = 0
        do iat = 1, num_atoms
         iatype = atom(iat)%atype
         do iorb = 1, n_basis_func(iatype)
          if (aos(iatype,iorb)%valence) then
           nvalence = nvalence + 1
          end if 
         end do
        end do
        
        ncore_states = nsaos - nvalence 
        if (ecp) then
         print '(/,a)', ' (default) ecp flag is active: core states will be integrated out'
         print '(1x,a,i5)', 'reduced dimension of the hilbert space =', nvalence 
        else
         print '(/,1x,a,i5)', 'dimension of the hilbert space =', nvalence 
        end if 

c       deallocate(s_count,p_count,d_count,f_count,
c     &   	   s_basfunc,p_basfunc,d_basfunc,f_basfunc,stat=ierr)

c       >>> fhi-aims related modifications
        deallocate(s_count,p_count,d_count,f_count,g_count,h_count,
     &   	   s_basfunc,p_basfunc,d_basfunc,f_basfunc,g_basfunc,h_basfunc,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,/,a)',
     &     ' [SUBROUTINE readbasis]: deallocation of tmp-arrays failed',
     &     ' ... will proceed further anyway'
	end if

c       in case of aims-input, remove temporary basis file >>>
        if (aims_input.and.(.not.testing)) then
         syscall = 'rm '//tmpbasis_file_name ;  call system(trim(syscall))
        end if
       
       end subroutine readbasis

       subroutine read_aims_basis(bssfile,tmpbssfile)
        implicit none
        character(*), intent (in) :: bssfile, tmpbssfile

        integer bfile, tmpfile, ierr, iat
        logical start_reading, flag, continue_read_atom
        logical, allocatable :: basis_info(:)

        character(16)  tmp_word
        character(32)  tmp_contr
	character      ch
	character(2)   tmp_ch_number
        integer        tmp_index, tmp_iat, tmp_n, tmp_l, tmp_m,
     &                 current_iat, ii, indexcount, tmp_atype 
        integer        l_count(0:lmax0) ! lmax0=5: up to h-functions

        logical        atomic_like(0:lmax0,1:orbdim)    ! .true. if "atomic" like
        logical        atomic_icount(0:lmax0)           ! just counter
        integer        n_number(0:lmax0,1:orbdim)

        allocate(basis_info(1:num_atom_types),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE read_aims_basis]: <basis_info> allocation failure'
        end if
        basis_info = .false.

        bfile = 40  ! basis file
        open(bfile,file=bssfile,status='old',action='read',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a,a)',
     &         ' can not open file "', trim(bssfile), '" for reading'
         print *, 'please, check your directory content or access rights'
         print *
         stop  ': transport module is terminated now'
        end if

        tmpfile = 41
        open(tmpfile,file=tmpbssfile,status='unknown',iostat=ierr)
        if (ierr.ne.0) then
          print '(/,a,a,a)',
     &          ' can not create a temporary file "', trim(tmpbssfile), '" !'
          print *, 'please, check your directory content or access rights'
          print *
          stop  ': transport module is terminated now'
        end if

        print '(/,a,a,a,/)', ' === reading file <', bssfile, '> ==='

        write(tmpfile,fmt='(a)') '$aims.pseudo.basis'
        write(tmpfile,fmt='(a)') '#caution! this file is created for testing purposes only'
	write(tmpfile,fmt='(a)') '#         gaussian exponents are dummy numbers '
	write(tmpfile,fmt='(a)') '*'

c       read comment lines >>>
        start_reading = .false.
        do while (.not.start_reading)
c        read(bfile,'(a1)') ch   ! read first symbol
         read(bfile,*) ch   ! read first symbol
         flag = (ch=='$').or.(ch=='#').or.(ch=='f')
         if (.not.flag)  start_reading = .true.
        end do
        backspace(bfile)

c       make a loop over atoms >>>
        iat = 0
        indexcount = 0
        do while (iat<num_atoms)

         iat = iat + 1
c 	 print *, ' iat = ', iat 
         read(bfile,fmt=*,iostat=ierr) tmp_index, tmp_word, current_iat, tmp_n, tmp_l, tmp_m

	 if (ierr.eq.-1) then
	   print '(1x,a)',       'ERROR: end of file is reached unexpectedly' 
	   print '(1x,a,a,a,/)', '       your "',trim(bssfile),'" file seems to be corrupted' 
	   stop ' : transport module is terminated now'
	 else if (ierr.gt.0) then
	   print '(1x,a)',       'reading of file resulted in error' 
	   print '(1x,a,a,a,/)', 'your "',trim(bssfile),'" file seems to be corrupted' 
	   stop ' : transport module is terminated now'
	 end if

         continue_read_atom = .true.
         backspace(bfile)

         l_count = 0
         atomic_like = .false.
         n_number = 0
         atomic_icount = 0
         do while (continue_read_atom.and.(tmp_index<nsaos))

          read(bfile,fmt=*,iostat=ierr) tmp_index, tmp_word, tmp_iat, tmp_n, tmp_l, tmp_m
c         print *, tmp_index, tmp_word, tmp_iat, tmp_n, tmp_l, tmp_m

	  if (ierr.eq.-1) then
	   print '(1x,a)',       'ERROR: end of file is reached unexpectedly' 
	   print '(1x,a,a,a,/)', '       your "',trim(bssfile),'" file seems to be corrupted' 
	   stop ' : transport module is terminated now'
	  else if (ierr.gt.0) then
	   print '(1x,a)',       'reading of file resulted in error' 
	   print '(1x,a,a,a,/)', 'your "',trim(bssfile),'" file seems to be corrupted' 
	   stop ' : transport module is terminated now'
	  end if

          if (tmp_iat.ne.current_iat) then
           continue_read_atom = .false.
           backspace(bfile)
          else

           atomic_icount(tmp_l) = atomic_icount(tmp_l) + 1
           n_number(tmp_l,atomic_icount(tmp_l)) = tmp_n
           if (trim(tmp_word)=='atomic') atomic_like(tmp_l,atomic_icount(tmp_l)) = .true.
           
           l_count(tmp_l) = l_count(tmp_l) + 1
           indexcount = indexcount + (2*tmp_l+1)
           do ii = 1, 2*tmp_l
            read(bfile,fmt=*,iostat=ierr) tmp_index
c           print *, tmp_index
	    if (ierr.eq.-1) then
	      print '(/,1x,a)',     'ERROR: end of file is reached unexpectedly' 
	      print '(1x,a,a,a,/)', '       your "',trim(bssfile),'" file seems to be corrupted' 
	      stop ' : transport module is terminated now'
	    else if (ierr.gt.0) then
	      print '(/,1x,a)',     'reading of file resulted in error' 
	      print '(1x,a,a,a,/)', 'your "',trim(bssfile),'" file seems to be corrupted' 
	      stop ' : transport module is terminated now'
	    end if
	   end do
	   
          end if
	  
         end do ! ... continue_read_atom 

         tmp_atype = atom(iat)%atype
         if (.not.basis_info(tmp_atype)) then
          basis_info(tmp_atype) = .true.
c         put dummy info into temporary basis file >>> 
	  write(tmpfile,fmt='(a,2x,a)') trim(atom(iat)%symbol(1:2)), 'fhi-aims-basis'
          
	  tmp_contr = '['     	  
	  
	  if (l_count(0).ne.0) then
	   if (l_count(0).lt.10) then ; write(tmp_ch_number,'(i1)') l_count(0)
	   else  ; write(tmp_ch_number,'(i2)') l_count(0) ; end if
	   tmp_contr = trim(tmp_contr)//trim(tmp_ch_number)//'s'
	  end if

	  if (l_count(1).ne.0) then
	   if (l_count(1).lt.10) then ; write(tmp_ch_number,'(i1)') l_count(1)
	   else  ; write(tmp_ch_number,'(i2)') l_count(1) ; end if
	   tmp_contr = trim(tmp_contr)//'/'//trim(tmp_ch_number)//'p'
	  end if

	  if (l_count(2).ne.0) then
	   if (l_count(2).lt.10) then ; write(tmp_ch_number,'(i1)') l_count(2)
	   else ; write(tmp_ch_number,'(i2)') l_count(2) ; end if     
	   tmp_contr = trim(tmp_contr)//'/'//trim(tmp_ch_number)//'d'
	  end if

	  if (l_count(3).ne.0) then
	   if (l_count(3).lt.10) then ; write(tmp_ch_number,'(i1)') l_count(3)
	   else ; write(tmp_ch_number,'(i2)') l_count(3) ; end if     
	   tmp_contr = trim(tmp_contr)//'/'//trim(tmp_ch_number)//'f'
	  end if

	  if (l_count(4).ne.0) then
	   if (l_count(4).lt.10) then ; write(tmp_ch_number,'(i1)') l_count(4)
	   else ; write(tmp_ch_number,'(i2)') l_count(4) ; end if     
	   tmp_contr = trim(tmp_contr)//'/'//trim(tmp_ch_number)//'g'
	  end if

	  if (l_count(5).ne.0) then
	   if (l_count(5).lt.10) then ; write(tmp_ch_number,'(i1)') l_count(5)
	   else ; write(tmp_ch_number,'(i2)') l_count(5); end if     
	   tmp_contr = trim(tmp_contr)//'/'//trim(tmp_ch_number)//'h'
	  end if
          tmp_contr = trim(tmp_contr)//']'

	  write(tmpfile,fmt='(a,1x,a,2x,a)') 
     &	        '#', trim(atom(iat)%symbol(1:2)), trim(tmp_contr)
          write(tmpfile,fmt='(a)') '*'
	  do tmp_l = 0, lmax0
           if (l_count(tmp_l).gt.0) then
            do ii = 1, l_count(tmp_l)
	     select case(tmp_l)
              case(0) ; write(tmpfile,'(i4,2x,a)') 1, 's'
              case(1) ; write(tmpfile,'(i4,2x,a)') 1, 'p'
              case(2) ; write(tmpfile,'(i4,2x,a)') 1, 'd'
              case(3) ; write(tmpfile,'(i4,2x,a)') 1, 'f'
              case(4) ; write(tmpfile,'(i4,2x,a)') 1, 'g'
              case(5) ; write(tmpfile,'(i4,2x,a)') 1, 'h'
             end select
              if ( atomic_like(tmp_l,ii) ) then
c              ! "atomic" like function gets formal coeffiecient 2.0
c              write(tmpfile,'(f14.8,f14.8)') 1.0d0+tmp_l+0.01d0*ii, 2.0d0
               write(tmpfile,'(f14.8,f14.8)') dble(n_number(tmp_l,ii)), 2.0d0
              else
c              ! othes functions get formal coeffiecients 1.0
c              write(tmpfile,'(f14.8,f14.8)') 1.0d0+tmp_l+0.01d0*ii, 1.0d0
               write(tmpfile,'(f14.8,f14.8)') dble(n_number(tmp_l,ii)), 1.0d0
              end if
	    end do ! ii	   
	   end if
  	  end do ! tmp_l
	  
          write(tmpfile,fmt='(a)') '*'
         end if

        end do ! loop over atoms

        write(tmpfile,fmt='(a)') '$end'
        write(tmpfile,'(a,i6)') '# degrees of freedom =', indexcount
        close(tmpfile)
        close(bfile)
c       stop 'stop with reading aims-basis file!'

       end subroutine read_aims_basis

       subroutine readsmat(ifile1,ifile2,outstatus)
c      **************************************************
c      reads square-root (and its inverse) of overlap
c      matrix from external files:
c      sqrt(smat) <-- ifile1 ; inv(sqrt(smat)) <-- ifile2
c      **************************************************
       implicit none

        character(*), intent (in)  :: ifile1, ifile2
        logical,      intent (out) :: outstatus
        integer, parameter         :: monitor = 100
        integer extfile1, extfile2, n,k,j, 
     &	        mincol, maxcol, ierr, tmpnsaos
c	character(10)   str1, str2

c       print '(/,a)', ' === reading <smat12/smat12inv> files ==='

        extfile1 = 29
c       open(extfile1,file=ifile1,
c    &	     status='old',action='read',iostat=ierr)
        open(extfile1,file=ifile1,status='old',
     &	              action='read',form='unformatted',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,1x,a,a,a)',
     &            'can not open file "', trim(ifile1), '" for reading'
         print '(1x,a)', 'overlap matrices will be recalculated'
         close(extfile1)
	 outstatus = .false.
	 return
        end if

        extfile2 = 28
c       open(extfile2,file=ifile2,
c    &	     status='old',action='read',iostat=ierr)
        open(extfile2,file=ifile2,status='old',
     &	              action='read',form='unformatted',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,1x,a,a,a)',
     &            'can not open file "', trim(ifile2), '" for reading'
         print '(1x,a)', 'overlap matrices will be recalculated'
         close(extfile1) ;  close(extfile2)
	 outstatus = .false.
	 return
        end if

c        allocate(smat12(nsaos,nsaos),stat=ierr)
c        if (ierr.ne.0) then
c         print *
c         stop '[SUBROUTINE readsmat]: <smat12> allocation failure'
c        end if
c
c        allocate(smat12inv(nsaos,nsaos),stat=ierr)
c        if (ierr.ne.0) then
c         print *
c         stop '[SUBROUTINE readsmat]: <smat12inv> allocation failure'
c        end if

c       reading comment lines
c       1st file
c       read(extfile1,*) str1, str2, tmpnsaos
        read(extfile1) tmpnsaos
	if (tmpnsaos.ne.nsaos) then
	 print '(/,a,i4,a,i4)', 
     &	  ' tmpnsaos= ', tmpnsaos ,
     &	  ' read from <smat12> does not match nsaos=', nsaos
         print '(1x,a)', 'overlap matrices will be recalculated'
         close(extfile1) ;  close(extfile2)
	 deallocate(smat12,smat12inv,stat=ierr)
         outstatus = .false.
	 return
	end if
c	read(extfile1,*) str1

c       reading comment lines
c       2nd file
c       read(extfile2,*) str1, str2, tmpnsaos
        read(extfile2) tmpnsaos
	if (tmpnsaos.ne.nsaos) then
	 print '(/,a,i4,a,i4)', 
     &	  ' tmpnsaos= ', tmpnsaos ,
     &	  ' read from <smat12inv> does not match nsaos=', nsaos
         print '(1x,a)', 'overlap matrices will be recalculated'
         close(extfile1) ;  close(extfile2)
	 deallocate(smat12,smat12inv,stat=ierr)
         outstatus = .false.
	 return
	end if
c	read(extfile2,*) str1
        
        print '(/,1x,a,$)', 'importing overlap matrices : | '
        k = 0
        do n = 1, nsaos
c        read comment line	
c        read(extfile1,*) str1
c        read(extfile2,*) str1
         k = k + n
         if ( mod(k,monitor)==0 ) print '(a,$)', '='  ! monitor
c         only upper triangular part of matrices is stored
          maxcol = 0
          do while (maxcol<n)
           mincol = maxcol+1
           maxcol = min(maxcol+4,n)
c           read(extfile1,fmt='(4d22.16)')
c     &                   (smat12(j,n),j=mincol,maxcol)
           read(extfile1) (smat12(j,n),j=mincol,maxcol)
           do j=mincol,maxcol 
            if (j.ne.n) smat12(n,j) = smat12(j,n)
	   end do
c           read(extfile2,fmt='(4d22.16)')
c     &                  (smat12inv(j,n),j=mincol,maxcol)
           read(extfile2) (smat12inv(j,n),j=mincol,maxcol)
           do j=mincol,maxcol 
            if (j.ne.n) smat12inv(n,j) = smat12inv(j,n)
	   end do     
	  end do
        end do
        print *, '| done'
        outstatus = .true.

        close(extfile1) ; close(extfile2)
       
       end subroutine readsmat

       subroutine readdnsmtrx(outputdens,nspin,info)
c       ************************************
c       reads charge & spin density matrices        
c       ************************************
c       output array                   
        double precision,intent(out) :: outputdens(:,:,:)
        integer, intent (in)         :: nspin
        integer, intent (out)        :: info
ccc                         info = 0 : successful termination
ccc                              = 1 : charge_dens_file failed      
ccc                              = 2 : spin_dens_file failed      
        integer dmatfile, smatfile, ierr, ierr1
        double precision tmp_chdens, tmp_spdens
	integer n, m  
       
        dmatfile = 31
        open(dmatfile,file=charge_dens_file,status='old', 
     &	                       action='read',iostat=ierr)
        smatfile = 32
        open(smatfile,file=spin_dens_file,status='old', 
     &	                       action='read',iostat=ierr1)
        if (ierr.ne.0) then
          print '(/,a,a,a)',
     &     ' [readdnsmtrx]: can not open file "', trim(charge_dens_file), '" for reading'
          close(dmatfile)
          info = 1
	else if (ierr1.ne.0) then
          print '(/,a,a,a)',
     &     ' [readdnsmtrx]: can not open file "', trim(spin_dens_file), '" for reading'
          close(smatfile)
          info = 2
        else
c         read data from external files 	
          do n = 1, nsaos
           do m = 1, n
            if (nspin.eq.1) then
             read(dmatfile,*) tmp_chdens
             outputdens(n,m,1) = tmp_chdens/2.0d0
             if (n.ne.m) outputdens(m,n,1) = tmp_chdens/2.0d0
            else
             read(dmatfile,*) tmp_chdens
             read(smatfile,*) tmp_spdens
             outputdens(n,m,1) = (tmp_chdens + tmp_spdens)/2.0d0
             outputdens(n,m,2) = (tmp_chdens - tmp_spdens)/2.0d0
             if (n.ne.m) then
              outputdens(m,n,1) = outputdens(n,m,1)
              outputdens(m,n,2) = outputdens(n,m,2)
             end if
            end if
           end do
          end do
          close(dmatfile) ; close(smatfile)
          info = 0
        end if

       end subroutine readdnsmtrx

       subroutine updatetcontrol(tcntrl,keyword,ipar,dblpar)
c       ********************* 
c       updates tcontrol-file
c       ********************* 
        implicit none
        character(*), intent(in)     :: tcntrl
        character(*), intent(in)     :: keyword
        integer, intent(in)          :: ipar
	double precision, intent(in) :: dblpar	
	
	integer, parameter :: linelength = 256
	
	integer tfile, ierr, linecount, count, i
	character(linelength) chline, keystring*25
	logical  readnext
 
        character(linelength), allocatable :: tcline(:)
	
c	print '(/,1x,a,$)', 'updating <tcontrol> ...'
        tfile = 10  ! tcontrol file
        open(tfile,file=tcntrl,status='old',action='read',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a,a)', 
     &	       ' can not open file "', trim(tcntrl), '" for reading'
         print *,
     &    'please, check your directory content or access rights'
         print *
         stop ': transport module is terminated now'
        end if
        
c       first we count number of lines >>>	
        readnext = .true.
	linecount = 0
	do while (readnext)
	 read(tfile,*) keystring
	 linecount = linecount + 1
	 if (trim(keystring)=='$end') readnext = .false.
	end do	

c       next, we update tcontrol file >>>
        allocate(tcline(linecount-1),stat=ierr)
	if (ierr.ne.0) then
         print *
         stop '[SUBR. updatetcontrol]: <tcline> allocation failure'
        end if
   
        rewind(tfile)
	count = 0
	i = 1
        do while (i < linecount)
	 read(tfile,'(a)') chline ; read(chline,*) keystring
c        here we save only those strings which are not modified 
 	 if ( trim(keystring) == trim(keyword) ) then
	   if (trim(keyword)=='$overlapint') then
	    read(tfile,'(a)') chline 
	    read(tfile,'(a)') chline 
	    i = i + 3
           else if ( (trim(keyword)=='$dens_matrix').and.(nspin.eq.2)) then
	    read(tfile,'(a)') chline 
	    read(tfile,'(a)') chline 
	    i = i + 3
           else if ( (trim(keyword)=='$dens_matrix').and.(nspin.eq.1)) then
	    read(tfile,'(a)') chline 
	    i = i + 2
	   else 
            i = i + 1
	   end if
	 else
            count = count + 1
	    tcline(count) = chline
            i = i + 1
	 end if
	end do
	close(tfile)

c       open again <tcontrol> and replace it with a modified file
        tfile = 10  ! tcontrol file
        open(tfile,file=tcntrl,status='replace',
     &	                       action='write',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a,a)', 
     &	       ' can not open file "', trim(tcntrl), '" for writing'
         print *,
     &    'please, check your directory content or access rights'
         print *
         stop ': transport module is terminated now'
        end if
       
        do i = 1, count
	 write(tfile,'(a)') trim(tcline(i))
	end do

        select case (trim(keyword))
          case ('$overlapint')	   
	   write(tfile,'(a)') '$overlapint'
	   write(tfile,'(2x,a,2x,a)') 
     &	         'smat12   ', 'file='//trim(smat12_file_name)
	   write(tfile,'(2x,a,2x,a)') 
     &	         'smat12inv', 'file='//trim(smat12inv_file_name)
 	  case ('$iter')
	   if (ipar.lt.10) then
     	    write(tfile,'(a,6x,i1)') '$iter', ipar
           else if (ipar.lt.100) then
     	    write(tfile,'(a,6x,i2)') '$iter', ipar
           else
     	    write(tfile,'(a,4x,i5)') '$iter', ipar
	   end if	    
	  case ('$efermi')
     	   write(tfile,'(a,f18.12)') '$efermi', dblpar  
	  case ('$valence_electrons')
     	   write(tfile,'(a,i5)') '$valence_electrons ', ipar  
c	  case ('#$efermi0')
c     	   write(tfile,'(a,f18.12)') '#$efermi0', dblpar  
	  case ('$u_energy')
     	   write(tfile,'(a,f16.10)') '$u_energy', dblpar  
	  case ('$dnorm')
     	   write(tfile,'(a,5x,e10.4)') '$dnorm', dblpar  
	  case ('$s1r')
     	   write(tfile,'(a,4x,d26.20)') '$s1r', dblpar  
	  case ('$s1i')
     	   write(tfile,'(a,4x,d26.20)') '$s1i', dblpar  
	  case ('$s2r')
     	   write(tfile,'(a,4x,d26.20)') '$s2r', dblpar  
	  case ('$s2i')
     	   write(tfile,'(a,4x,d26.20)') '$s2i', dblpar  
	  case ('$s3r')
     	   write(tfile,'(a,4x,d26.20)') '$s3r', dblpar  
	  case ('$s3i')
     	   write(tfile,'(a,4x,d26.20)') '$s3i', dblpar  
          case ('$dens_matrix')
	   write(tfile,'(a)') '$dens_matrix on'
	   if (nspin.eq.1) then
            write(tfile,'(2x,a,i5)') 'num.of.occupied orbitals ', ipar	   
	   else
c           nspin == 2
            write(tfile,'(2x,a,i5)') 'num.of.occupied alpha-orbitals ', ipar	   
            write(tfile,'(2x,a,i5)') 'num.of.occupied beta-orbitals  ', int(dblpar)	   
	   end if
	  case ('$adjust_rsigma')
	   if (ipar > 0) then
	     write(tfile,'(a)') '$adjust_rsigma on'
	   else 
	     write(tfile,'(a)') '$adjust_rsigma off'	   
	   end if 
	  case ('rfactor') 
	    write(tfile,'(2x,a,3x,f5.2)') 'rfactor', dblpar 
	  case ('rguess') 
	    write(tfile,'(2x,a,4x,f5.2)') 'rguess', dblpar 
	  case ('nelcnv') 
	    write(tfile,'(2x,a,5x,d8.2)') 'nelcnv', dblpar 
	  case ('iterlimit') 
	    write(tfile,'(2x,a,1x,i3)')   'iterlimit', ipar 
          case ('$self_energy') 
	    write(tfile,'(a,a)') '$self_energy file=', trim(self_energy_file) 
	  case ('$ecp')
	   if (ipar > 0) then
	     write(tfile,'(a)') '$ecp on'
	   else 
	     write(tfile,'(a)') '$ecp off'
	   end if  
	  case default ; 
	end select
	
	write(tfile,'(a)') '$end'
        close(tfile)
        deallocate(tcline)
c       print *, 'done'

       end subroutine updatetcontrol

      end module read_externals


