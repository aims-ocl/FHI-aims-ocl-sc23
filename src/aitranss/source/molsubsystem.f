c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           april 2011
c     last revision:  jan 2012
c###########################################################

      module molsubsystem
c     ***************************************      
c     contains set of routines to perform 
c     an analysis of the molecular sybsystem  
c     ***************************************      

       use globalvars
       use read_externals 
c      use density_matrix0
       use hamiltonian
       use tools
       use rmatrix
       
       implicit none

       real(4), private ::  stime, ftime, secs
       integer, private ::  mnts  

c      local array to store the green's function for given energy
       complex(8), private, allocatable :: gf(:)

c      molecule-reservoir index-correlation matrix
       double precision, allocatable :: molecule_res_index(:)

c      hamiltonian of the molecule
       double precision, allocatable :: h_molecule(:,:,:)    
ccc      (external) hamiltonian of the bare molecule       
ccc       double precision, allocatable :: h_mol_ext(:,:,:)
c      spectrum of <h_molecule> or <h_mol_ext>
       double precision, allocatable :: ep_molecule(:,:)    
c      and its eigenvectors 
       double precision, allocatable :: umat_molecule(:,:,:)    

c      density matrix of the molecule 
       double precision, allocatable :: rhomat_molecule(:,:,:)    
c      density matrix of the whole system 
       double precision, allocatable :: rhomat(:,:,:)    

c      overlap matrix within the molecular subspace
       double precision, allocatable :: omatrix(:,:)
c      its square-root 
       double precision, allocatable :: omatrix12(:,:)
c      and inverse of <omatrix12>
       double precision, allocatable :: omatrix12inv(:,:)

c      projector operators
       double precision, allocatable :: osmatrix(:,:)
c      complex(8), allocatable       :: cmplx_osmatrix(:,:)
 
c      fermi level for the taget (reference) system
       double precision ref_efermi 

c      complex eigenvectors of the extended hamiltonian
c      projected onmolecular orbtitals
       complex(8), private, allocatable :: ub(:,:,:), bu(:,:,:)
 
      contains

      subroutine execute_molecular_subsystem
c      ****************************************************************
c      higher level routine for the <mol_subsystem> type of calculation
c      ****************************************************************
       implicit none
       integer bottom_orb, top_orb, ierr

       print '(/,a)', ' ANALYZING MOLECULAR SUBSYSTEM >>>'

c      -- get dimensions
       call get_dimensions
c      -- initialize molecule-reservoir index-correlation matrix
       call initialize_indexes
c      -- overlap integrals for the ion
       call init_omatrix

c      -- compute hamiltonian of the molecule
       call init_h_molecule       
c      -- compute spectrum of the molecular hamiltonian
       call h_molecule_spectrum(1)

c       if (read_hmat) then
c        -- read external molecular hamiltonian if required       
c         call readhmat(inhmat_file_name,h_molecule,molecule_size)
c	 .. and evaluates its spectrum 
c         call h_molecule_spectrum(2)
c       end if  

c      -- evaluate density matrix of the molecule
       call molecule_densmatrix(bottom_orb,top_orb)
       
c      -- build up projector operators on molecular states      
       call build_mol_orb_projectors

c      -- evaluates local DOS projected to molecular orbitals
       call projected_ldos(bottom_orb,top_orb)
        
c      -- deallocate all module related arrays
       deallocate(gf,molecule_res_index,h_molecule,ep_molecule,
     &            umat_molecule,rhomat_molecule,rhomat,omatrix,
     &		  omatrix12,omatrix12inv,ub,bu,osmatrix,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE execute_molecular_subsystem]: ',
     &           'impossible to deallocate arrays'
        print *, 'nevertheless, proceed further ...'
       end if
					    
       print '(/,a)', ' <<< DONE WITH ANALYSIS OF THE MOLECULAR SUBSYSTEM '
  
      end subroutine execute_molecular_subsystem

      subroutine get_dimensions      
c      **********************************************
c      computes dimensions of the "molecule" and the
c      "reservoir" blocks of the hamiltonian matrix         
c      **********************************************
       implicit none    
       integer iat, is, iatype, iasort ! , iorb, norb, lm, dcount
       logical flag

       molecule_size = 0 
       reservoir_size = 0

       do iat = 1, num_atoms

        iatype = atom(iat)%atype
        iasort = atom(iat)%asort
        flag = .false.
	do is = 1, num_res_atoms 
	 flag = ( flag.or.(iasort == res_atom_types(is)) )
	end do

        if (flag) then
c        given atom <iat> belongs to the "reservoir"
         reservoir_size = reservoir_size + n_basis_func(iatype) 	
	else
c        given atom <iat> belongs to the "molecule": 
         molecule_size = molecule_size + n_basis_func(iatype) 	
        end if

       end do ! iat

       print '(/,a,i4,a,i4)',
     &   ' BLOCK DIMENSIONS:  molecule --> ', molecule_size,
     &   '  reservoir(s) --> ', reservoir_size

       if (molecule_size.le.0) then
        print '(/,a)',
     &             ' dimension of the "molecule" related subspace could not be found !'
        print '(a,/)', ' please, check your <tcontrol> file ... i will quit now'
        stop ': transport module is terminated'
       end if

c      consistency check:
       if ((reservoir_size + molecule_size).ne.nsaos) then
        print '(/,a)',
     &    ' !!! WARNING: sum of the molecule"s and reservoir"s dimensions does not much <nsaos>!'
c       print '(a)', ' something went wrong here ... "reservoir_size" will be corrected'
        print '(a)', ' something went wrong here ... '
        stop ': transport module is terminated'
c        reservoir_size = nsaos - molecule_size
c        print '(/,a,i4,a,i4)',
c     &   ' !!! WARNING: CORRECTED BLOCK DECOMPOSITION:  molecule --> ', molecule_size,
c     &   '  reservoir(s) --> ', reservoir_size
       end if

      end subroutine get_dimensions

      subroutine initialize_indexes
c      ********************************************************
c      for the molecular subsystem finds its counterpart within
c      the larger set of orbitals of the whole system
c      ********************************************************
       implicit none
       integer ma, ra, na, iorb, iat, norb,
     &         iatype, iasort, is, ierr
       logical aflag

c      allocate index-correlation array
       allocate(molecule_res_index(molecule_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE initialize_indexes]: <molecule_res_index> allocation failure'
       end if

       na = 0 ; ma = 0 ; ra = 0
       do iat = 1, num_atoms
        iatype = atom(iat)%atype
        iasort = atom(iat)%asort
        norb = n_basis_func(iatype)
        aflag = .false.
        do is = 1, num_res_atoms
         aflag = ( aflag.or.(iasort == res_atom_types(is)) )
        end do

        do iorb = 1, norb
          na = na + 1
          if (aflag) then
            ra = ra + 1
          else
	    ma = ma + 1
            molecule_res_index(ma) = na
          end if  ! aflag

        end do ! iorb
       end do ! iat

c      check dimensions consistency
       print '(/,a,i4,a,i4,a,i4)', ' ma =', ma, ';  ra = ', ra, ';  na = ', na
       if ( (ma.ne.molecule_size) .or. (ra.ne.reservoir_size) .or. (na.ne.nsaos) ) then
         print '(/,a)',
     &     ' [initialize_indexes]: dimensions are out of range!'
         print '(a,/)', ' something went wrong here ... will quit now'
         stop ': transport module is terminated'
       end if
       print '(1x,a)', 'molecule-reservoir index-correlation matrix is initialized'

      end subroutine initialize_indexes

      subroutine init_omatrix
c     ****************************************************************
c     initializes the overlap matrix for a molecule: <omatrix>,
c     computes its square-root <omatrix12>, and inverse of <omatrix12>
c     ****************************************************************

       implicit none
       integer ierr, ia, ib, na, nb

       allocate(omatrix(molecule_size,molecule_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_omatrix]: <omatrix> allocation failure'
       end if

       allocate(omatrix12(molecule_size,molecule_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_omatrix]: <omatrix12> allocation failure'
       end if

       allocate(omatrix12inv(molecule_size,molecule_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_omatrix]: <omatrix12inv> allocation failure'
       end if

c      fill up <omatrix> matrix
       print '(/,1x,a,$)', 'molecular orbitals, initializing overlap matrix ...'
       do ib = 1, molecule_size
        nb = molecule_res_index(ib)
        do ia = 1, molecule_size
         na = molecule_res_index(ia)
         omatrix(ia,ib) = smat(na,nb)
        end do
       end do
       print *, 'done'

c      computes <omatrix12> and <omatrix12inv>
       call sqrt_smat(omatrix,omatrix12,omatrix12inv,molecule_size,3)

      end subroutine init_omatrix

      subroutine init_h_molecule
c      ****************************************
c      computes hamiltonian of the molecule by
c      projecting hamiltonian of the system 
c      onto molecular degrees of freedom    
c      ****************************************
       implicit none
       integer ierr, ia, ib, na, ma, ispin
       double precision, allocatable :: tmph0(:,:), xsmat12(:,:),
     &                                  ph0(:,:), tmph0ion(:,:) 
 
       logical :: loewdin_sub = .true. 

       allocate(h_molecule(molecule_size,molecule_size,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE correlate_indexes]: <h_molecule> allocation failure'
       end if

       allocate(tmph0ion(molecule_size,molecule_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE correlate_indexes]: <tmph0ion> allocation failure'
       end if

       allocate(osmatrix(molecule_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_h_molecule]: <osmatrix> allocation failure'
       end if
       osmatrix = 0.0d0

c       if (mol_spectr_func) then
c        allocate(cmplx_osmatrix(molecule_size,nsaos),stat=ierr)
c        if (ierr.ne.0) then
c        print *
c        stop
c     &   '[SUBROUTINE init_h_molecule]: <cmplx_osmatrix> allocation failure'
c        end if
c        cmplx_osmatrix = czero
c       end if

       allocate(xsmat12(molecule_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_h_molecule]: <xsmat12> allocation failure'
       end if
       xsmat12 = 0.0d0

       allocate(tmph0(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_h_molecule]: <tmph0> allocation failure'
       end if

       allocate(ph0(molecule_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_h_molecule]: <ph0> allocation failure'
       end if

       print '(/,a)', ' COMPUTING HAMILTONIAN OF THE MOLECULE >>>'
       print '(/,2x,a,$)', 'builds up projector operator ...'
c      take required block of the smat12
       do ia = 1, molecule_size
        ma = molecule_res_index(ia)
        do na = 1, nsaos
	  xsmat12(ia,na) = smat12(ma,na)
	end do
       end do

c      biuld up projector operator: osmatrix = omatrix12inv * xsmat12
       call xdble_ab('n','n',molecule_size,molecule_size,nsaos,omatrix12inv,xsmat12,osmatrix)
       print *, 'done'
c       do ia = 1, molecule_size
c         write(66,*) ia, osmatrix(ia,ia)
c       end do 

c       if (mol_spectr_func) then
c        forall(n=1:molecule_size,m=1:nsaos)
c	 cmplx_osmatrix(n,m) = cmplx(osmatrix(n,m),0.0d0,8)
c	end forall
c       end if

c      project hamiltonian onto ionic degrees of freedom
c      h0ion = osmatrix * h0mat * osmatrix^T

       if (loewdin_sub) then 
        
        print '(2x,a,$)', 'takes a molecular hamiltonian block ...'
        do ispin = 1, nspin
         do ia = 1, molecule_size
          ma = molecule_res_index(ia)
          do ib = 1, molecule_size 
	   na = molecule_res_index(ib)
	   h_molecule(ia,ib,ispin) = h0mat(na,ma,ispin)
	  end do
	 end do
	end do  
	print *, 'done'
       
       else 

        print '(2x,a,$)', 'projecting hamiltonian onto molecular orbitals ...'
        do ispin = 1, nspin
c        -- 1st step:  ph0 = osmatrix * h0mat
         call xdble_ab('n','n',molecule_size,nsaos,nsaos,osmatrix,h0mat(:,:,ispin),ph0)
c        -- 2nd step:  h_molecule = ph0 * osmatrix^T
         call xdble_ab('n','t',molecule_size,nsaos,molecule_size,ph0,osmatrix,h_molecule(:,:,ispin))
        end do ! ispin
        print *, 'done'
       
       end if ! loewdin subspace
 
       print '(/,a)', ' <<< DONE WITH MOLECULAR HAMILTONIAN'

c      deallocate(xsmat12,osmatrix,tmph0,ph0,tmph0ion,stat=ierr)
       deallocate(xsmat12,tmph0,ph0,tmph0ion,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE init_h_molecule]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
       end if
       
      end subroutine init_h_molecule

      subroutine h_molecule_spectrum(icall)
c      ******************************************
c      computes spectrum of the hamiltonian
c      projected on molecule's subspace
c      ******************************************
       implicit none 
       
       integer, intent(in) :: icall 
       integer ierr, ispin, n, m, ofile

c      temporary arrays 
       double precision, allocatable :: tmp_h_molecule(:,:)

       character(32) :: afile = 'mo.levels.alpha.dat',
     &                  bfile = 'mo.levels.beta.dat',
     &                  mfile = 'mo.levels.dat',
     &                  ofilename

       if (.not.allocated(ep_molecule)) then
        allocate(ep_molecule(molecule_size,nspin),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop
     &    '[SUBROUTINE h_molecule_spectrum]: <ep_molecule> allocation failure'
        end if
       end if

       if (.not.allocated(umat_molecule)) then
        allocate(umat_molecule(molecule_size,molecule_size,nspin),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop
     &    '[SUBROUTINE h_molecule_spectrum]: <umat_molecule> allocation failure'
        end if
       end if

       allocate(tmp_h_molecule(molecule_size,molecule_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE h_molecule_spectrum]: <tmp_h_molecule> allocation failure'
       end if

       print '(/,a,/)', ' SOLVING EIGENVALUE PROBLEM FOR THE MOLECULE >>>'
       do ispin = 1, nspin
 
        print '(2x,a,i1,a,$)', 'ispin = ', ispin, ': ... '

c       save local copy of the hamiltonian
        forall (n=1:molecule_size,m=1:molecule_size)
         tmp_h_molecule(n,m) = h_molecule(n,m,ispin)
        end forall

c       solving eigenvalue problem 
        call realspectrum(tmp_h_molecule,ep_molecule(:,ispin),umat_molecule(:,:,ispin),molecule_size)
        print *, 'done'

	if (icall.eq.1) then 
	   ofile = 90 + ispin
	   if (nspin.eq.1) then 
             ofilename = trim(mfile)
	   else if (ispin.eq.1) then 
	     ofilename = trim(afile)
	   else
 	     ofilename = trim(bfile)
	   end if
	   open(ofile,file=trim(ofilename),status='unknown',iostat=ierr)
	   if (ierr.ne.0) then
	     stop '[h_molecule_spectrum] : can not open output file!'
	   end if
	else  
	   ofile = 80 + ispin
	end if 

        if (nspin.eq.1) then 
	 write(ofile,'(a)') '$molecular.levels'
	else if (ispin.eq.1) then 
	 write(ofile,'(a)') '$molecular.levels: alpha channel'	
	else
	 write(ofile,'(a)') '$molecular.levels: beta channel'		
	end if
	write(ofile,'(a,f14.10,a,f14.10,a)') 
     &              '# E_Fermi = ', efermi, ' H  =', efermi * hartree, ' eV'
        write(ofile,'(a,9x,a,15x,a,12x,a)') '# orb','E [Hartree]','E [eV]','E - EF [eV]'
 	  
        do n = molecule_size, 1, -1
c        do n=1, molecule_size 
	  write(ofile,'(i5,4x,e20.14,4x,f16.8,4x,f16.8)') 
     &	        n, ep_molecule(n,ispin), ep_molecule(n,ispin) * hartree,  
     &          (ep_molecule(n,ispin) - efermi) * hartree        
        end do
        write(ofile,'(a)') '$end'

c        do n = 1, molecule_size
c         write(80+ispin,*) n, umat_molecule(n,n,ispin)
c        end do

       end do ! ispin
       print '(/,a)', ' <<< DONE WITH THE EIGENVALUE PROBLEM'
       
       deallocate(tmp_h_molecule,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE h_molecule_spectrum]: ',
     &           'impossible to deallocate a temporary array'
        print *, 'nevertheless, proceed further ...'
       end if


      end subroutine h_molecule_spectrum

      subroutine molecule_densmatrix(bottom_orb,top_orb)
c      **************************************
c      evaluates density matrix projected to 
c      molecule's orbitals
c      **************************************
       implicit none 
    
       integer, intent(out) :: bottom_orb, top_orb
        
       integer           ia, ja, na, ma, n, m, ispin, 
     &                   rhofile, ierr ! info  
       double precision, allocatable :: xsmat(:,:), 
     &                   tmp_osmat(:,:), pmat(:,:),
     &                   tmp_rho_molecule(:,:), p_rho(:,:)
       double precision  tmp_ch(2)
       character(64) ::  rhofilename = 'molecule.occ.numbers'

       allocate(rhomat(nsaos,nsaos,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop
     &    '[SUBROUTINE molecule_densmatrix]: <rho> allocation failure'
       end if  

c      read density matrices from external files
c       print '(/,1x,a,$)', '<MOLECULE''s SUBSPACE>: reading density matrices ... '
c       call readdnsmtrx(rhomat,nspin,info)
c       if (info.eq.0) then
c         print *, 'done'
c       else
c         print '(/,1x,a,i2,/)',
c     &         '[dion_densmatrix]: can not read external files, info =', info
c        stop ': transport module is terminated now'
c       end if

       forall(n=1:nsaos,m=1:nsaos,ispin=1:nspin)
        rhomat(n,m,ispin) = neqdmat(n,m,ispin)
       end forall

       allocate(xsmat(molecule_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE molecule_densmatrix]: <xsmat> allocation failure'
       end if

       allocate(tmp_osmat(molecule_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE molecule_densmatrix]: <tmp_osmat> allocation failure'
       end if

       allocate(pmat(molecule_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE molecule_densmatrix]: <pmat> allocation failure'
       end if

       allocate(rhomat_molecule(molecule_size,molecule_size,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE molecule_densmatrix]: <rhomat_molecule> allocation failure'
       end if  

       allocate(tmp_rho_molecule(molecule_size,molecule_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE molecule_densmatrix]: <tmp_rho_molecule> allocation failure'
       end if  

       allocate(p_rho(molecule_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE molecule_densmatrix]: <p_rho> allocation failure'
       end if  

       print '(/,1x,a/)', '<MOLECULE''s SUBSPACE>: EVALUATING DENSITY MATRIX FOR THE MOLECULE >>>'

       print '(2x,a,$)', 'builds up projector operator ...'

c      take required block of the smat
c       do ia = 1, molecule_size
c        ma = molecule_res_index(ia)
c        do na = 1, nsaos
c          xsmat(ia,na) = smat(ma,na)
c        end do
c       end do
c      biuld up first operator: tmp_osmat = omatrix12inv * xsmat
c       call xdble_ab('n','n',molecule_size,molecule_size,nsaos,omatrix12inv,xsmat,tmp_osmat)
c       print *, 'done'

       do ia = 1, molecule_size
        ma = molecule_res_index(ia) 
	do na = 1, nsaos 
	  tmp_osmat(ia,na) = smat12(ma,na)
	end do
       end do
       print *, 'done'

       do ispin = 1, nspin
        print '(2x,a,i1,a,$)', 'spin = ', ispin, ' : ... '

c       make a unitary transform of <tmp_osmat>:  pmat <-- umat_molecule^{t} * tmp_osmat
        call xdble_ab('t','n',molecule_size,molecule_size,nsaos,umat_molecule(:,:,ispin),tmp_osmat,pmat)
c       print '(a,$)', 'pmat ... '

c       project the density matrix onto the subspace of molecular orbitals: 
c       rhomat_molecule = pmat * rhomat * pmat^{t} 

c       -- 1st step: p_rho <-- pmat * rhomat 
         call xdble_ab('n','n',molecule_size,nsaos,nsaos,pmat,rhomat(:,:,ispin),p_rho)
c        print '(a,$)', 'pmat * rho ... '

c       -- 2nd step: tmp_rho_molecule <-- p_rho * pmat^{t}
         call xdble_ab('n','t',molecule_size,nsaos,molecule_size,p_rho,pmat,tmp_rho_molecule)
c        print '(a,$)', 'pmat * rho * pmat^{t} ... '

c        fill up array rhomat_molecule(ispin) <-- tmp_rho_molecule 
         forall(ia=1:molecule_size,ja=1:molecule_size)
           rhomat_molecule(ia,ja,ispin) = 
     &	           0.5d0*(tmp_rho_molecule(ia,ja)+tmp_rho_molecule(ja,ia))
         end forall
         print *, 'done'

       end do ! ispin

       print '(/,1x,a)', 
     &       '<MOLECULE''s SUBSPACE>: <<< DONE WITH DENSITY MATRIX FOR THE MOLECULE'

c      print out rho(m,m) into external file
       rhofile = 50
       open(rhofile,file=trim(rhofilename),status='unknown',iostat=ierr)
       if (ierr.ne.0) then
           stop '[molecule_densmatrix] : can not open output file!'
       end if
       
       if (nspin.eq.2) then 
        write(rhofile,'(a)') '$molecule.dens.matrix'
        write(rhofile,'(a)') '#occupation numbers, alpha & beta channels'
	write(rhofile,'(a,f14.10,a,f14.10,a)') 
     &              '# E_Fermi = ', efermi, ' H  =', efermi * hartree, ' eV'
        write(rhofile,'(a,10x,a)') 
     &                '# orb      E.a - EF [eV]      alpha_occ',
     &                'E.b - EF [eV]      beta_occ'
       else
        write(rhofile,'(a)') '$molecule.dens.matrix'
        write(rhofile,'(a)') '#occupation numbers'
	write(rhofile,'(a,f14.10,a,f14.10,a)') 
     &              '# E_Fermi = ', efermi, ' H  =', efermi * hartree, ' eV'
        write(rhofile,'(a)') '# orb       E - EF [eV]        orb_occ'
       end if	

       tmp_ch = 0.0d0
       do ia = molecule_size, 1, -1
        if (nspin.eq.2) then
          write(rhofile,'(i5,4x,f14.8,4x,f14.10,6x,f14.8,4x,f14.10)')	  
     &	               ia, (ep_molecule(ia,1)-efermi)*hartree, rhomat_molecule(ia,ia,1), 
     &                     (ep_molecule(ia,2)-efermi)*hartree, rhomat_molecule(ia,ia,2)  	 
         tmp_ch(1) = tmp_ch(1) + rhomat_molecule(ia,ia,1)
         tmp_ch(2) = tmp_ch(2) + rhomat_molecule(ia,ia,2)
	else
         write(rhofile,'(i5,4x,f14.8,4x,f14.10)') 
     &	               ia, (ep_molecule(ia,1)-efermi)*hartree, rhomat_molecule(ia,ia,1) 	 
         tmp_ch(1) = tmp_ch(1) + rhomat_molecule(ia,ia,1)
	end if
       end do  
       if (nspin.eq.2) then
         write(rhofile,'(a,f14.10,a,f14.10,a,f14.10,a,f14.10)') 
     &	  '# na =', tmp_ch(1), '  nb =', tmp_ch(2), 
     &    '  na+nb =', tmp_ch(1)+tmp_ch(2), '  na-nb =', tmp_ch(1)-tmp_ch(2)
         bottom_orb = int( (tmp_ch(1)+tmp_ch(2))/2.0d0 ) - occ_mol_orbitals
         top_orb    = int( (tmp_ch(1)+tmp_ch(2))/2.0d0 ) + empty_mol_orbitals               
       else 
         write(rhofile,'(a,f10.6)') '# total occupation per spin       : ', tmp_ch(1)
         write(rhofile,'(a,f10.6)') '# occupation in both spin channels: ', 2.0d0 * tmp_ch(1)
         bottom_orb = int(tmp_ch(1)) - occ_mol_orbitals
	 top_orb    = int(tmp_ch(1)) + empty_mol_orbitals  
       end if
       bottom_orb = max(bottom_orb,0)        
       top_orb    = min(top_orb,molecule_size)

ccc    optional >>>
ccc    print out density matrix
c       do ispin = 1, nspin
c        write(rhofile,'(a,i1)') '#density matrix: ispin = ', ispin
c        do ia = 1, molecule_size
c         do ja = 1, molecule_size
c	  write(rhofile,'(2x,e12.6,$)'), rhomat_molecule(ia,ja,ispin) 
c	 end do
c	 write(rhofile,*) 
c        end do
c      end do ! ispin

       write(rhofile,'(a)') '$end'
       close(rhofile)

       deallocate(xsmat,tmp_osmat,pmat,tmp_rho_molecule,p_rho,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE molecule_densmatrix]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
       end if

      end subroutine molecule_densmatrix

      subroutine build_mol_orb_projectors
c     *************************************************
c     builds projector operators on molecular orbitals 
c     *************************************************
       integer     ispin, ia, na, ma, ierr 
       complex(8), allocatable :: tmp_b(:,:),  tmp_binv(:,:), tmp_umat(:,:)

c      allocate global arrays
       allocate(ub(molecule_size,nsaos,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop
     &   '[SUBROUTINE build_mol_orb_projectors]: <ub> allocation failure'
       end if
c      print *, 'target 1: OK'
       
       allocate(bu(nsaos,molecule_size,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop
     &   '[SUBROUTINE build_mol_orb_projectors]: <bu> allocation failure'
       end if
c      print *, 'target 2: OK'

c      allocate temporary arrays
       allocate(tmp_b(molecule_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop
     &   '[SUBROUTINE build_mol_orb_projectors]: <tmp_b> allocation failure'
       end if
c      print *, 'target 5: OK'

       allocate(tmp_binv(nsaos,molecule_size),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop
     &   '[SUBROUTINE build_mol_orb_projectors]: <tmp_binv> allocation failure'
       end if
c      print *, 'target 6: OK'

       allocate(tmp_umat(molecule_size,molecule_size),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop
     &   '[SUBROUTINE build_mol_orb_projectors]: <tmp_umat> allocation failure'
       end if
c      print *, 'target 7: OK'

       print '(/,a,$)', ' projecting eigenvectors of H_full on molecular orbitals ... ' 	
c      evaluate projector operators >>>
       do ispin = 1, nspin

c       take selected block of heigvec & invheigvec
        do ia = 1, molecule_size
         ma = molecule_res_index(ia)
         do na = 1, nsaos
          tmp_b(ia,na) = heigvec(ma,na,ispin)
          tmp_binv(na,ia) = invheigvec(na,ma,ispin)
         end do
        end do
c       print *, 'target 8: OK'

c       takes a local copy of umat_molecule for a given spin
        forall(na=1:molecule_size,ma=1:molecule_size)
         tmp_umat(na,ma) = cmplx(umat_molecule(na,ma,ispin),0.0d0,8)
        end forall
c       print *, 'target 9: OK'

c       compute ub <-- tmp_umat^T * tmp_b
        call xcmplx_ab('t','n',molecule_size,molecule_size,nsaos,tmp_umat,tmp_b,ub(:,:,ispin))
c       print *, 'target 10: OK'

c       compute bu <-- tmp_binv * tmp_umat
        call xcmplx_ab('n','n',nsaos,molecule_size,molecule_size,tmp_binv,tmp_umat,bu(:,:,ispin))
c       print *, 'target 11: OK'

       end do ! ispin
       print *, 'done'    

c      deallocate temporary arrays
       deallocate(tmp_b,tmp_binv,tmp_umat,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE build_mol_orb_projectors]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'proceed further anyway ...'
       end if
c      print *, 'target 13: OK'

      end subroutine build_mol_orb_projectors
 
      subroutine projected_ldos(minorb, maxorb)
c     **********************************************
c     evaluates LDOS projected to molecular orbitals
c     **********************************************      
      
       integer, intent(in)        :: minorb, maxorb

       character(32), allocatable :: ofilename(:)
       character(6)     tmpstring
       character(64)    outfmt

       integer          iorb, ierr, ifile, ispin, n, wsize
       double precision energy, tmp_dos, tmp_mdos, btmp_dos, btmp_mdos,
     &                  ener_save, estep_save, eend_save

       complex(8), allocatable :: orbdos(:), b_orbdos(:) 
       complex(8)                 mol_dos, b_mol_dos        

c      allocation of array with GF:
       allocate(gf(nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE projected_ldos]: <gf> allocation failure'
       end if

       wsize = maxorb-minorb+1
c      allocate ldos arrays
       allocate(orbdos(wsize),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE projected_ldos]: <orbdos> allocation failure'
       end if

       if (nspin.eq.2) then
        allocate(b_orbdos(wsize),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop
     &    '[SUBROUTINE projected_ldos]: <b_orbdos> allocation failure'
        end if
       end if

c      take care about output files
       allocate(ofilename(wsize),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE projected_ldos]: <ofilename> allocation failure'
       end if
       
       do iorb = minorb, maxorb
        write(tmpstring,'(i5)') iorb 
        ifile = iorb-minorb+1
	ofilename(ifile) = 'ldos.orb_'//trim(adjustl(tmpstring))//'.dat'
	open(ifile,file=trim(ofilename(iorb-minorb+1)),status='unknown',iostat=ierr)
	if (ierr.ne.0) then
	  print '(a,i5)', ' [projected_ldos] : can not open output file, iorb = ', iorb
          stop  ': transport module is terminated '
	else 
         write(ifile,fmt='(a,i5)')  '#local dos projected to molecular orbital no.', iorb 
         if (nspin ==  1) then
          write(ifile,fmt='(a)') '#non-spin-polarized calculation'
         else
          write(ifile,fmt='(a)') '#spin-polarized calculation'
         end if ! nspin
         write(ifile,'(a,f10.6,a,f17.12,a)')
     &                  '#bias =',bias,' V    efermi =', efermi, ' H'
         write(ifile,fmt='(a,3x,a,7x,a)') '#', 'E [Hartree]', 'E-EF [eV]'
	end if
       end do
       
c      choose a format of the output line
       if (nspin.eq.1) then
         outfmt = '(1x,f14.10,f16.8,2x,2(E16.6))'
       else
         outfmt = '(1x,f14.10,f16.8,2x,4(E16.6))'
       end if

       print '(/,a,/)', ' <MOLECULE''s SUBSYSTEM>: OUTPUT ORBITAL PROJECTED LDOS >>>'

c      take care about energy mesh, if input is in eV
       if (evunits) then
        ener_save = ener ; estep_save = estep ; eend_save = eend
        ener  = efermi + ener/hartree
        estep = estep/hartree
        eend =  efermi + eend/hartree
       end if

c      cycle over energy points
       n = 0 ; energy = ener
       do while (energy <= eend)
        
        if (evunits) then
         print '(2x,a,f14.10,a,$)', 'E - EF = ',(energy-efermi)*hartree,' eV    '
        else
         print '(2x,a,f14.10,a,$)', 'E = ',energy,' H    '
        end if

c       cycle over spin channels
        do ispin = 1, nspin
c        evaluating GF at given energy
         call greensf(energy,ispin,gf)
c        compute local dos projected on molecular orbitals
         if (ispin.eq.1) then
           call get_projected_ldos(gf,1,orbdos,mol_dos,minorb,maxorb)
         else
           call get_projected_ldos(gf,2,b_orbdos,b_mol_dos,minorb,maxorb)
         end if
        end do ! spin channels

c       output computed data to external files >>>
        do iorb = minorb, maxorb
          ifile = iorb - minorb + 1
          if (nspin.eq.1) then
            tmp_dos  = dble(orbdos(ifile))/hartree
            tmp_mdos = dble(mol_dos)/hartree
	    write(ifile,trim(outfmt))
     &            energy,(energy-efermi)*hartree, tmp_dos, tmp_mdos
          else
            tmp_dos   = dble(orbdos(ifile))/hartree
            btmp_dos  = dble(b_orbdos(ifile))/hartree
	    tmp_mdos  = dble(mol_dos)/hartree
	    btmp_mdos = dble(b_mol_dos)/hartree
	    write(ifile,trim(outfmt))
     &            energy,(energy-efermi)*hartree,tmp_dos,tmp_mdos,btmp_dos,btmp_mdos
          end if
        end do ! iorb
        print *, ' done'

        n = n + 1
        energy = ener + n * estep
       end do !  cycle over E-points

c      close files:
       do iorb = minorb, maxorb
        ifile = iorb - minorb + 1
        write(ifile,'(a)') '#end'
	close(ifile)
       end do

       print '(/,a,/)', ' <MOLECULE''s SUBSYSTEM>: DONE WITH DOS OUTPUT <<<'

c      restore parameters for energy mesh
       if (evunits) then
        ener = ener_save ; estep = estep_save ; eend = eend_save
       end if

c      deallocate temporary arrays
       deallocate(orbdos,ofilename,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE projected_ldos]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'proceed further anyway ...'
       end if

       if (nspin.eq.2) then 
         deallocate(b_orbdos,stat=ierr)
         if (ierr.ne.0) then
           print '(/,a,a)', ' [SUBROUTINE projected_ldos]: ',
     &                      'impossible to deallocate temporary arrays'
           print *, 'proceed further anyway ...'
         end if
       end if

      end subroutine projected_ldos


      subroutine get_projected_ldos(gf,ispin,orbdos,mol_dos,minorb,maxorb)
c     ****************************************************************        
c     evaluates LDOS projected on molecular orbitals for given energy:
c     array with values of the GF is used as an input
c     ****************************************************************        
       complex(8), intent(in)  :: gf(:)  
       integer,    intent(in)  :: ispin, minorb, maxorb
       complex(8), intent(out) :: orbdos(:), mol_dos
    
       integer iorb, wsize, p, m
       complex(8) tmpgr, tmpdos

       wsize = maxorb - minorb + 1
       do iorb = 1, wsize
        orbdos(iorb) = czero
       end do
       mol_dos = czero 
       
c      computes full local dos at molecule
       do m = 1, molecule_size
        do p = 1, nsaos
         tmpgr =  ub(m,p,ispin) * gf(p) * bu(p,m,ispin)
         mol_dos  = mol_dos + ione*(tmpgr - dconjg(tmpgr))/(2.0d0*pi)
	end do
       end do

c      computes contributions projected on selected orbitals
       do m = minorb, maxorb
        tmpdos = czero
        do p = 1, nsaos
         tmpgr  = ub(m,p,ispin) * gf(p) * bu(p,m,ispin)
         tmpdos = tmpdos + ione*(tmpgr - dconjg(tmpgr))/(2.0d0*pi)
	end do
	orbdos(m-minorb+1) = tmpdos 
       end do

      end subroutine get_projected_ldos

      end module molsubsystem       
