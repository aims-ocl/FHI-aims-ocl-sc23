c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           may 2010
c     last revision:  jan 2012
c###########################################################

      module hubbardu
c     **************************************      
c     contains a set of routines to perform 
c     an lda+u type of calculation 
c     **************************************      

       use globalvars
       use read_externals 
       use density_matrix0
       use hamiltonian
       use tools
       use rmatrix
       
       implicit none

       real(4), private ::  stime, ftime, secs
       integer, private ::  mnts  

       integer, private, parameter :: ldim = 5  ! d-functions only
  
c      ion-reservoir index-correlation matrix
       double precision, allocatable :: ion_res_index(:)

c      hamiltonian of the 3d-ion(s)
       double precision, allocatable :: h_ion(:,:,:)    
c      hamiltonian of the 3d-ion(s) for the target system
       double precision, allocatable :: h_trgt(:,:,:)    
c      spectrum of <h_ion> or <h_trgt>
       double precision, allocatable :: ep_ion(:,:)    
c      and its eigenvectors 
       double precision, allocatable :: umat_ion(:,:,:)    

c      lda+u correction to the hamiltonian 
       double precision, allocatable :: dh_xc(:,:,:)

c      density matrix of the ion 
       double precision, allocatable :: rho_ion(:,:,:)    
c      density matrix of the whole system 
       double precision, allocatable :: rho(:,:,:)    
c      density matrix of the ion for the target system
c      double precision, allocatable :: rho_trgt(:,:,:)    

c      overlap matrix for ion(s)
       double precision, allocatable :: omat(:,:)
c      its square-root 
       double precision, allocatable :: omat12(:,:)
c      and inverse of <omat12>
       double precision, allocatable :: omat12inv(:,:)
 
c      projector operators  
c       double precision, allocatable :: osmat(:,:)
c       complex(8), allocatable       :: cmplx_osmat(:,:)

c      lm-weights  
       double precision, allocatable :: d_weight(:,:),
     &                                  dlm_weight(:,:,:)
 
c      fermi level for the taget (reference) system
       double precision ref_efermi 

c      infinitesimal small value
       double precision, parameter :: ueps = 1.0d-8
 
      contains

      subroutine make_lda_u
c      *******************************************************
c      higher level routine for the ldau+u type of calculation
c      *******************************************************
       implicit none
       integer ierr
       
c      -- get dimensions
       call dimdef
c      -- initialize ion-reservoir index-correlation matrix
       call correlate_indexes
c      -- overlap integrals for the ion
       call init_omat

c      -- compute hamiltonian of the ion(s)
       call init_h_ion
       
c      -- if read_hion is 'on', read hamiltonian and the
c         density matrix of the target system
       
       if (read_hion) then

        allocate(h_trgt(ion_size,ion_size,nspin),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop
     &   '[SUBROUTINE make_lda_u]: <h_trgt> allocation failure'
        end if

        allocate(rho_ion(ion_size,ion_size,nspin),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop
     &   '[SUBROUTINE make_lda_u]: <rho_trgt> allocation failure'
        end if
       
        call read_hion_data(h_trgt,rho_ion,ion_size,hion_file_name)  
c       stop 'end of testing'
	
       end if ! read_hion	

c      -- compute spectrum of the ionic hamiltonian
       call h_ion_spectrum

c      -- evaluate density matrix of the 3d-ion(s)
       if (.not.read_hion) then   
         call dion_densmatrix(1)
       else
         call dion_densmatrix(2)
       end if
       
c      save 3d-ion(s) related data to external file
       if (save_hion) 
     &   call save_hion_data(h_ion,rho_ion,ion_size,hion_file_name)  

c      -- compute correction to the xc-potential 
c         & update hamiltonian of the system
       if (u_j > ueps) then
        call dv_xc 
       else 
        print '(/,1x,a)', '<HUBBARD-U>: screened Coulomb repulsion is inifinesimally small'
        print '(1x,a)',   '             hamiltonian of the system will not be updated'
       end if

c      compute spectrum of the modified hamiltonian
       if (no_leads) call update_hspectrum

c      -- deallocate all temporary arrays
       call dealloc_ldau       
  
      end subroutine make_lda_u

      subroutine dimdef      
c      ********************************************
c      computes dimensions of the "ion(s)" and the
c      "reservoir" blocks of the hamiltonian matrix         
c      ********************************************
       implicit none    
       integer iat, is, iatype, iasort, iorb, norb, lm, dcount
       logical flag

       ion_size = 0 
       res_size = 0
       do iat = 1, num_atoms
        iatype = atom(iat)%atype
        iasort = atom(iat)%asort

        flag = .false.
	do is = 1, num_ratoms 
	 flag = ( flag.or.(iasort == ratom_types(is)) )
	end do

        if (flag) then
c        given atom <iat> belongs to a "reservoir"
         res_size = res_size + n_basis_func(iatype) 	
	else
c        given atom <iat> belongs to an "ion": 
c        here we pick up only the first d-type CGTOs
c        with large exponents but not the diffuse d-orbitals
         norb = n_basis_func(iatype) 
	 dcount = 0
         do iorb = 1, norb
          lm = aos(iatype,iorb)%lm
	  if ((lm.ge.5).and.(lm.le.9).and.(dcount.lt.5)) then
c	  if ((lm.ge.5).and.(lm.le.9)) then
	   dcount = dcount + 1
	  end if
	 end do ! orb
         ion_size = ion_size + dcount 	
         res_size = res_size + n_basis_func(iatype) - dcount
	end if
       end do ! iat

       print '(/,a,i4,a,i4)',
     &   ' U-term, DIMENSIONS:  3d-ion(s) --> ', ion_size,
     &   '  reservoir --> ', res_size

       if (ion_size.le.0) then
        print '(/,a)',
     &    ' dimension of the correlated subspace could not be found !'
        print '(/,a)',
     &    ' your "correlated" sites most likely do not have d-orbitals'
        print '(a,/)', ' please, check your <tcontrol> file ... i will quit now'
        stop ': transport module is terminated'
       end if

c      consistency check:
       if ((res_size + ion_size).ne.nsaos) then
        print '(/,a)',
     &    ' !!! WARNING: sum of the molecule/ion and reservoir dimensions does not much nsaos!'
        print '(a)', ' something went wrong here ... "res_size" will be corrected'
c       stop ': transport module is terminated'
        res_size = nsaos - ion_size
        print '(/,a,i4,a,i4)',
     &   ' !!! WARNING: U-term, CORRECTED DIMENSIONS:  3d-ion(s) --> ', ion_size,
     &   '  reservoir --> ', res_size
       end if

      end subroutine dimdef
      
      subroutine correlate_indexes
c      **********************************************
c      for ion d-orbital finds its counterpart within  
c      a larger set of orbitals of the whole system       
c      **********************************************
       implicit none
       integer ia, ra, na, iorb, iat, norb, 
     &         iatype, iasort, is, dcount, lm, ierr
       logical aflag
        
c      allocate an index-correlation array
       allocate(ion_res_index(ion_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE correlate_indexes]: <ion_res_index> allocation failure'
       end if

       na = 0 ; ia = 0 ; ra = 0       
       do iat = 1, num_atoms
        iatype = atom(iat)%atype
        iasort = atom(iat)%asort
	norb = n_basis_func(iatype)
        aflag = .false.
	do is = 1, num_ratoms 
	 aflag = ( aflag.or.(iasort == ratom_types(is)) )
        end do

        dcount = 0
        do iorb = 1, norb

          na = na + 1
          if (aflag) then  
	    ra = ra + 1
          else  
c           take care about a given site and its d-orbitals
	    lm = aos(iatype,iorb)%lm
	    if ((lm.ge.5).and.(lm.le.9).and.(dcount.lt.5)) then
c	    if ((lm.ge.5).and.(lm.le.9)) then
	     dcount = dcount + 1
	     ia = ia + 1
	     ion_res_index(ia) = na   
c            if (testing) write(99,*) ia, na
	    else
	     ra = ra + 1 
	    end if
  	  end if  ! aflag
	  
        end do ! iorb
       end do ! iat      

c      check dimensions consistency
       print '(/,a,i4,a,i4,a,i4)', ' ia =', ia, ';  ra = ', ra, ';  na = ', na
       if ( (ia.ne.ion_size) .or. (ra.ne.res_size) .or. (na.ne.nsaos) ) then
         print '(/,a)',
     &     ' [correlate_indexes]: dimensions are out of range!'
         print '(a,/)', ' something went wrong here ... will quit now'
         stop ': transport module is terminated'
       end if
       print '(1x,a)', 'ion-reservoir index-correlation matrix is initialized'

      end subroutine correlate_indexes

      subroutine init_omat
c     ***********************************************
c     initializes overlap matrix for 3d-ions: <omat>,
c     computes its square-root <omat12>, 
c     and inverse of <omat12>       
c     ***********************************************
       implicit none
       integer ierr, ia, ib, na, nb

       allocate(omat(ion_size,ion_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_omat]: <omat> allocation failure'
       end if

       allocate(omat12(ion_size,ion_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_omat]: <omat12> allocation failure'
       end if

       allocate(omat12inv(ion_size,ion_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_omat]: <omat12inv> allocation failure'
       end if

c      fill up <omat> matrix
       print '(/,1x,a,$)', 'ionic orbitals, initializing overlap matrix ...'
       do ib = 1, ion_size
        nb = ion_res_index(ib)
        do ia = 1, ion_size
         na = ion_res_index(ia)
         omat(ia,ib) = smat(na,nb) 
        end do
       end do
       print *, 'done'       

c      computes <omat12> and <omat12inv>
       call sqrt_smat(omat,omat12,omat12inv,ion_size,2)
     
      end subroutine init_omat


      subroutine init_h_ion
c      ****************************************
c      computes hamiltonian of the 3d-ion(s) by
c      projecting hamiltonian of the system 
c      onto ionic degrees of freedom    
c      ****************************************
       implicit none
       integer ierr, ia, na, ma, n, m, ispin
       double precision, allocatable :: xsmat12(:,:), ph0(:,:)

       allocate(h_ion(ion_size,ion_size,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE correlate_indexes]: <h_ion> allocation failure'
       end if

       allocate(osmat(ion_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_h_ion]: <osmat> allocation failure'
       end if
       osmat = 0.0d0

       if (spectr_func) then
        allocate(cmplx_osmat(ion_size,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_h_ion]: <cmplx_osmat> allocation failure'
        end if
        cmplx_osmat = czero
       end if

       allocate(xsmat12(ion_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_h_ion]: <xsmat12> allocation failure'
       end if
       xsmat12 = 0.0d0

       allocate(ph0(ion_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE init_h_ion]: <ph0> allocation failure'
       end if

       print '(/,a)', ' <HUBBARD-U>: COMPUTING HAMILTONIAN OF THE 3d-ION(S) >>>'
       print '(/,2x,a,$)', 'builds up projector operator ...'
c      take required block of the smat12
       do ia = 1, ion_size
        ma = ion_res_index(ia)
        do na = 1, nsaos
	  xsmat12(ia,na) = smat12(ma,na)
	end do
       end do

c      biuld up the projector operator: osmat = omat12inv * xsmat12
       call xdble_ab('n','n',ion_size,ion_size,nsaos,omat12inv,xsmat12,osmat)
       print *, 'done'
c       do ia = 1, ion_size
c         write(66,*) ia, osmat(ia,ia)
c       end do 

       if (spectr_func) then
        forall(n=1:ion_size,m=1:nsaos)
	 cmplx_osmat(n,m) = cmplx(osmat(n,m),0.0d0,8)
	end forall
       end if

c      project the hamiltonian onto ionic degrees of freedom
c      h0ion = osmat * h0mat * osmat^T
       print '(2x,a,$)', 'projecting hamiltonian onto 3d-ion(s) orbitals ...'
       do ispin = 1, nspin

c       -- 1st step:  ph0 = osmat * h0mat
        call xdble_ab('n','n',ion_size,nsaos,nsaos,osmat,h0mat(:,:,ispin),ph0)

c       -- 2nd step:  h_ion = ph0 * osmat^T
        call xdble_ab('n','t',ion_size,nsaos,ion_size,ph0,osmat,h_ion(:,:,ispin))

       end do ! ispin
       print *, 'done'

       print '(/,a)', ' <HUBBARD-U>: <<< DONE WITH THE IONIC HAMILTONIAN'

c      deallocate(xsmat12,osmat,tmph0,ph0,tmph0ion,stat=ierr)
       deallocate(xsmat12,ph0,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE init_h_mat]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
       end if
       
      end subroutine init_h_ion

      subroutine h_ion_spectrum
c      ******************************************
c      computes spectrum of the ionic hamiltonian
c      ******************************************
       implicit none 
       integer ierr, ispin, n, m, p, lm, ofile

c      temporary arrays 
       double precision, allocatable :: tmp_h_ion(:,:)
     
        character(32) :: afile = 'z_ion.alpha.tmp',
     &                   bfile = 'z_ion.beta.tmp',
     &                   mfile = 'z_ion.mos.tmp',
     &                   ofilename

       allocate(ep_ion(ion_size,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE h_ion_spectrum]: <ep_ion> allocation failure'
       end if

       allocate(umat_ion(ion_size,ion_size,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE h_ion_spectrum]: <umat_ion> allocation failure'
       end if

       allocate(tmp_h_ion(ion_size,ion_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE h_ion_spectrum]: <tmp_h_ion> allocation failure'
       end if

       print '(/,a,/)', ' <HUBBARD-U>: SOLVING EIGENVALUE PROBLEM FOR THE 3d-ION(S) >>>'
       do ispin = 1, nspin
 
        print '(2x,a,i1,a,$)', 'ispin = ', ispin, ': ... '

c       save local copy of the hamiltonian
        if (read_hion) then
c         take hamiltonian of the reference system
          forall (n=1:ion_size,m=1:ion_size)
	   tmp_h_ion(n,m) = h_trgt(n,m,ispin)
 	  end forall
        else 
          forall (n=1:ion_size,m=1:ion_size)
	   tmp_h_ion(n,m) = h_ion(n,m,ispin)
 	  end forall
	end if 

c       solving eigenvalue problem 
        call realspectrum(tmp_h_ion,ep_ion(:,ispin),umat_ion(:,:,ispin),ion_size)
        print *, 'done'

c        do n = 1, ion_size
c         write(80+ispin,*) n, umat_ion(n,n,ispin)
c        end do

       end do ! ispin
       print '(/,a)', ' <HUBBARD-U>: <<< DONE WITH THE EIGENVALUE PROBLEM'

       deallocate(tmp_h_ion,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE h_ion_spectrum]: ',
     &           'impossible to deallocate temporary arrays'
        print *, 'nevertheless, proceed further ...'
       end if

c      computes lm-weights for each ionic orbital >>>
       call weights
       print '(/,1x,a)', 'localization weights are assigned to ionic orbitals ...'
c      stop

c      print out spectrum to output file 
       do ispin = 1, nspin 

        if (testing.or.zspectrum.or.ldau) then
         ofile = 64-ispin
         if (nspin.eq.1) then
          ofilename = mfile
         else
          if (ispin.eq.1) then ; ofilename = afile
          else ; ofilename = bfile
          end if
         end if
         
	 open(ofile,file=trim(ofilename),status='unknown',iostat=ierr)
         if (ierr.ne.0) then
           stop '[h_ion_spectrum] : can not open output file!'
         end if

         write(ofile,'(a)') '$'//trim(ofilename)
	 write(ofile,'(a,5x,a,12x,a,5x,a,11x,a,10x,a,15x,a,15x,a,13x,a)') 
     &	             '# orbital','E [H]','E-EF [eV]', 'd & f-weight',
     &               '2zz-xx-yy','xz','yz','xy','xx-yy'
         do p=1, ion_size 
         write(ofile,'(i5,4x,e14.8,4x,e14.8,3x,e14.8,4x,5(e17.8))')
     &                p, ep_ion(p,ispin), (ep_ion(p,ispin)-efermi)*hartree, 
     &                d_weight(p,ispin), (dlm_weight(p,lm,ispin), lm = 1,5)
c                                         ! d-orbital contributions
         end do
         write(ofile,'(a)') '$end'
	 close(ofile)

        end if  ! testing or zspectrum or ldau
	 
       end do ! ispin
       print '(1x,a)', '... and stored to external files'

      end subroutine h_ion_spectrum

      subroutine weights
      implicit none 
c     *********************************************
c     for each eigenvector of the ionic hamiltonian
c     computes its partial weights after projection 
c     on the d-orbitals
c     *********************************************
       integer          ierr, na, ma, ra, p, ispin, norb,
     &                  is, iat, iorb, iatype, iasort, lm, dcount
       logical          aflag  
       double precision tmpcoeff  

       allocate(d_weight(ion_size,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE weights]: <d_weight> allocation failure'
       end if
       d_weight = 0.0d0

       allocate(dlm_weight(ion_size,ldim,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE weights]: <dlm_weight> allocation failure'
       end if
       dlm_weight = 0.0d0

       na = 0; ma = 0; ra = 0;
       do iat = 1, num_atoms

        iatype = atom(iat)%atype
        iasort = atom(iat)%asort
        norb = n_basis_func(iatype)

        aflag = .false.
        do is = 1, num_ratoms
          aflag = ( aflag.or.(iasort==ratom_types(is)) )
        end do

        dcount = 0 
        do iorb = 1, norb
         na = na + 1
         if (aflag) then
          ra = ra + 1
c         nothing to do anymore ...
         else
c         ion(s) degree of freedom is found here:
c         then we should count contributions from d-orbitals
c         and add them to the "localization" weights
	  lm = aos(iatype,iorb)%lm
	  if ((lm.ge.5).and.(lm.le.9).and.(dcount.lt.5)) then
c	  if ((lm.ge.5).and.(lm.le.9)) then
	   dcount = dcount + 1
	   ma = ma + 1
           do ispin = 1, nspin
            do p = 1, ion_size         
	     tmpcoeff = umat_ion(ma,p,ispin)
             d_weight(p,ispin)=d_weight(p,ispin)+tmpcoeff*tmpcoeff
c            collect partial m-weights:
             if (lm.le.9) then
              dlm_weight(p,lm-ldim+1,ispin) = 
     &	               dlm_weight(p,lm-ldim+1,ispin)+tmpcoeff*tmpcoeff
	     end if
            end do  ! p
           end do ! ispin 
	  else
	     ra = ra + 1 
	  end if
         end if ! aflag

        end do ! iorb

       end do ! iat

       if ( (ma.ne.ion_size) .or. (ra.ne.res_size) .or. (na.ne.nsaos) ) then
         print '(/,a)',
     &     ' [weigths]: dimensions are out of range!'
         print '(a,/)', ' something went wrong here ... will quit now'
         stop ': transport module is terminated'
       end if

      end subroutine weights

      subroutine dion_densmatrix(icall)
c      **************************************
c      evaluates density matrix of the 3d-ion
c      for the d-orbital subspace
*      **************************************
       implicit none 
       integer, intent (in) :: icall
       
       integer           ia, ja, na, ma, ispin, rhofile, info, ierr
       double precision, allocatable :: xsmat(:,:), osmat(:,:), 
     &                   pmat(:,:), tmp_rho_ion(:,:), p_rho(:,:)
       double precision  nel, nel_a, nel_b, uenergy, tmp_ch(2)
       character(64) ::  rhofilename = 'z_rho.tmp'

       select case (icall)

       case(1)
c      main call starts here >>>

       allocate(rho(nsaos,nsaos,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop
     &    '[SUBROUTINE dion_densmatrix]: <rho> allocation failure'
       end if  
 
       if (no_leads.and.(.not.read_dmat)) then
c       compute the equilibrium density matrix
        print '(/,1x,a)', '<HUBBARD-U>: CALCULATING EQUILIBRIUM DENSITY MATRIX >>> '
        if (nspin.eq.1) then
         call dens_matrix0(mo_coeff,rho,1,occ)        
         nel = mtrace(rho,smat,nsaos,1)
         print '(/,2x,a,f12.6)', 'number of electrons = ', 2.0d0*nel
        else
c        alpha channel 
         call dens_matrix0(mo_coeff,rho,1,a_occ)    
         nel_a = mtrace(rho,smat,nsaos,1)
         print '(/,2x,a,f12.6)', 'number of alpha-electrons = ', nel_a
c        beta channel
         call dens_matrix0(mo_coeff,rho,2,b_occ)        
         nel_b = mtrace(rho,smat,nsaos,2)
         print '(2x,a,f12.6)', 'number of beta-electrons  = ', nel_b
	end if
       print '(/,1x,a)', '<HUBBARD-U>: <<< DONE WITH EQUILIBRIUM DENSITY MATRIX'
       else 
c       read density matrices from external files
        print '(/,1x,a,$)', '<HUBBARD-U>: reading density matrices ... '
        call readdnsmtrx(rho,nspin,info) 
        if (info.eq.0) then
         print *, 'done'
        else
         print '(/,1x,a,i2,/)', 
     &	       '[dion_densmatrix]: can not read external files, info =', info
 	 stop ': transport module is terminated now'
        end if
c      stop ': stop after reading'
       end if ! no_leads 

       allocate(xsmat(ion_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dion_densmatrix]: <xsmat> allocation failure'
       end if

       allocate(osmat(ion_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dion_densmatrix]: <osmat> allocation failure'
       end if

       allocate(pmat(ion_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dion_densmatrix]: <pmat> allocation failure'
       end if

       allocate(rho_ion(ion_size,ion_size,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dion_densmatrix]: <rho_ion> allocation failure'
       end if  

       allocate(tmp_rho_ion(ion_size,ion_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dion_densmatrix]: <tmp_rho_ion> allocation failure'
       end if  

       allocate(p_rho(ion_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dion_densmatrix]: <p_rho> allocation failure'
       end if  

       print '(/,1x,a/)', '<HUBBARD-U>: EVALUATING A DENSITY MATRIX FOR THE 3d-ION(S) >>>'

       print '(2x,a,$)', 'builds up a projector operator ...'
c      take required block of the smat
       do ia = 1, ion_size
        ma = ion_res_index(ia)
        do na = 1, nsaos
          xsmat(ia,na) = smat(ma,na)
        end do
       end do
c      biuld up first operator: osmat = omat12inv * xsmat
       call xdble_ab('n','n',ion_size,ion_size,nsaos,omat12inv,xsmat,osmat)
       print *, 'done'

       do ispin = 1, nspin
        print '(2x,a,i1,a,$)', 'spin = ', ispin, ' : ... '

c       make a unitary transform of <osmat>:  pmat <-- umat_ion^{t} * osmat
        call xdble_ab('t','n',ion_size,ion_size,nsaos,umat_ion(:,:,ispin),osmat,pmat)
c       print '(a,$)', 'pmat ... '

c       project density matrix onto 3d-ion orbital subspace: 
c       rho_ion = pmat * rho * pmat^{t} 

c       -- 1st step: p_rho <-- pmat * rho 
         call xdble_ab('n','n',ion_size,nsaos,nsaos,pmat,rho(:,:,ispin),p_rho)
c        print '(a,$)', 'pmat * rho ... '

c       -- 2nd step: tmp_rho_ion <-- p_rho * pmat^{t}
         call xdble_ab('n','t',ion_size,nsaos,ion_size,p_rho,pmat,tmp_rho_ion)
c        print '(a,$)', 'pmat * rho * pmat^{t} ... '

c        fill up array rho_ion(ispin) <-- tmp_rho_ion 
         forall(ia=1:ion_size,ja=1:ion_size)
           rho_ion(ia,ja,ispin) = 
     &	           0.5d0*(tmp_rho_ion(ia,ja)+tmp_rho_ion(ja,ia))
         end forall
         print *, 'done'

       end do ! ispin

       print '(/,1x,a)', 
     &       '<HUBBARD-U>: <<< DONE WITH THE DENSITY MATRIX FOR THE 3d-ION(S)'

       case default 
c       in case of icall == 2, rho_ion matrix 
c       has been read from external file      
        ;
	 
       end select

c      print out rho(m,m) into external file
       rhofile = 50
       open(rhofile,file=trim(rhofilename),status='unknown',iostat=ierr)
       if (ierr.ne.0) then
           stop '[dion_densmatrix] : can not open output file!'
       end if
       
       if (nspin.eq.2) then 
        write(rhofile,'(a)') '$z_rho.tmp'
        write(rhofile,'(a)') '#occupation numbers, alpha & beta channels'
c        write(rhofile+2,'(a)') '$z_rho.tmp'
c        write(rhofile+2,'(a)') '#occupation numbers, beta channel'
       else
        write(rhofile,'(a)') '$z_rho.tmp'
        write(rhofile,'(a)') '#occupation numbers'
       end if	

       tmp_ch = 0.0d0
       do ia = 1, ion_size
        if (nspin.eq.2) then
c         write(rhofile+1,'(i4,4x,e16.10,4x,f14.10)') 
c     &	               ia, ep_ion(ia,1), rho_ion(ia,ia,1) 	 
c         write(rhofile+2,'(i4,4x,e16.10,4x,f14.10)') 
c     &	               ia, ep_ion(ia,2), rho_ion(ia,ia,2) 	 
          write(rhofile,'(i4,4x,e16.10,4x,f14.10,6x,e16.10,4x,f14.10)') 
     &	               ia, ep_ion(ia,1), rho_ion(ia,ia,1), 
     &                     ep_ion(ia,2), rho_ion(ia,ia,2)  	 
         tmp_ch(1) = tmp_ch(1) + rho_ion(ia,ia,1)
         tmp_ch(2) = tmp_ch(2) + rho_ion(ia,ia,2)
	else
         write(rhofile,'(i4,4x,e16.10,4x,f14.10)') 
     &	               ia, ep_ion(ia,1), rho_ion(ia,ia,1) 	 
         tmp_ch(1) = tmp_ch(1) + rho_ion(ia,ia,1)
	end if
       end do  
       if (nspin.eq.2) then
         write(rhofile,'(a,f14.10,a,f14.10,a,f14.10,a,f14.10)') 
     &	  '# na =', tmp_ch(1), '  nb =', tmp_ch(2), 
     &    '  na+nb =', tmp_ch(1)+tmp_ch(2), '  na-nb =', tmp_ch(1)-tmp_ch(2)
       else 
         write(rhofile,'(a,f10.6)') '#  ', tmp_ch(1)
       end if

ccc    optional >>>
ccc    print out density matrix
       do ispin = 1, nspin
        write(rhofile,'(a,i1)') '#density matrix: ispin = ', ispin
        do ia = 1, ion_size
         do ja = 1, ion_size
	  write(rhofile,'(2x,e12.6,$)'), rho_ion(ia,ja,ispin)
 	 end do
	 write(rhofile,*) 
        end do
       end do ! ispin

c      compute correction to the lda/gga energy
       uenergy = u_energy(u_j,rho_ion,ion_size)
       write(rhofile,'(a,f16.10)') '#u_energy ', uenergy
       call updatetcontrol(trim(tcntrl_file_name),'$u_energy',0,uenergy)
       
       write(rhofile,'(a)') '$end'
       close(rhofile)

       if (icall.eq.1) then
        deallocate(xsmat,osmat,pmat,tmp_rho_ion,p_rho,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE dion_densmatrix]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
        end if
       end if ! icall

      end subroutine dion_densmatrix

      subroutine dv_xc
c     **************************************************
c     -- computes an orbital dependent contribution to 
c        the xc-potential according to the lda+u recipe:
c        S.L.Dudarev et al. Phys.Rev.B 57, 1505 (1998)    
c     -- updates the hamiltonian <h0mat> of a system
c     **************************************************
       implicit none 
       double precision, allocatable :: dh0mm(:,:), umat_em(:,:),
     &                   xsmat12(:,:), o12s12(:,:),  
     &                   dh0ab(:,:), dhp(:,:), delta_h(:,:)
       integer ispin, n, m, m1, ma, na, nb, ia, ierr

       allocate(dh0mm(ion_size,ion_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dv_xc]: <dh0mm> allocation failure'
       end if

       allocate(dh0ab(ion_size,ion_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dv_xc]: <dh0ab> allocation failure'
       end if

       allocate(o12s12(ion_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dv_xc]: <o12s12> allocation failure'
       end if

       allocate(xsmat12(ion_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dv_xc]: <xsmat12inv> allocation failure'
       end if

       allocate(umat_em(ion_size,ion_size),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dv_xc]: <umat_em> allocation failure'
       end if

       allocate(dhp(ion_size,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dv_xc]: <dhp> allocation failure'
       end if

       allocate(delta_h(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dv_xc]: <delta_h> allocation failure'
       end if

       allocate(dh_xc(nsaos,nsaos,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE dv_xc]: <dh_xc> allocation failure'
       end if

c      take required block of the <smat12inv> matrix
       do ia = 1, ion_size
        ma = ion_res_index(ia)
        do na = 1, nsaos
          if (ldau_old) then
	   xsmat12(ia,na) = smat12inv(ma,na)
	  else
	   xsmat12(ia,na) = smat12(ma,na)
	  end if 
	end do
       end do
 
       if (ldau_old) then
c       build up a projector operator : o12s12 = omat12 * xsmat12inv
        call xdble_ab('n','n',ion_size,ion_size,nsaos,omat12,xsmat12,o12s12)
       else
c       build up a projector operator : o12s12 = omat12inv * xsmat12
        call xdble_ab('n','n',ion_size,ion_size,nsaos,omat12inv,xsmat12,o12s12)
       end if

c      in case of "read_hion" flag, take care about modifications for dh0ab
       if (read_hion) then
        print '(/,1x,a,/,1x,a)', 
     &    	'WARNING: hamiltonian of the 3d-ion(s) will be substituted',
     &          '         by the one of the reference system'	
       end if
 
       do ispin = 1, nspin
        
c       compute correction to the xc-potential (see dudarev et al.)
c       !! pay attention: u_j is in eV !!!
        forall(m=1:ion_size,m1=1:ion_size)
	  dh0mm(m,m1) = -u_j * rho_ion(m,m1,ispin) / hartree
	end forall               
	forall(m=1:ion_size)
	  dh0mm(m,m) = 0.5d0 * u_j/hartree + dh0mm(m,m)
	end forall
       
c       build up matrix dh0ab = umat_hion * dh0mm * umat_hion^t 
	call dble_ab('n','n',ion_size,umat_ion(:,:,ispin),dh0mm,umat_em)
	call dble_ab('n','t',ion_size,umat_em,umat_ion(:,:,ispin),dh0ab)

c       in case of "read_hion" flag, take care about modifications for dh0ab
        if (read_hion) then
         do n = 1, ion_size
	  do m = 1, ion_size
	   dh0ab(n,m) = dh0ab(n,m) - h_ion(m,n,ispin) + h_trgt(m,n,ispin)
	  end do
	 end do 
	 if (efermi_guess) then
          forall (n=1:ion_size)
	   dh0ab(n,n) = dh0ab(n,n) - ref_efermi + efermi
	  end forall
	 end if 
	end if ! read_hion 

c       build up matrix dh0 = o12s12^t * dh0ab * o12s12
	call xdble_ab('n','n',ion_size,ion_size,nsaos,dh0ab,o12s12,dhp)
	call xdble_ab('t','n',nsaos,ion_size,nsaos,o12s12,dhp,delta_h)
	
c       copy a local array to the global one
        forall (na=1:nsaos,nb=1:nsaos)
	  dh_xc(na,nb,ispin) = delta_h(na,nb)
	end forall	
       
       end do ! ispin

c      add an lda+u xc-contribution to the hamiltonian
       forall (na=1:nsaos,nb=1:nsaos,ispin=1:nspin)
	 h0mat(na,nb,ispin) = h0mat(na,nb,ispin) + dh_xc(na,nb,ispin)
       end forall	

       print '(/,1x,a,$)', '<HUBBARD-U>: correction to the xc-potential'
       print '(a)',        ' is added to the hamiltonian'

       deallocate(dh0mm,dh0ab,umat_em,xsmat12,o12s12,dhp,delta_h,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE dv_xc]: ',
     &           'impossible to deallocate temporary arrays'
        print *, 'nevertheless, proceed further ...'
       end if
      
      end subroutine dv_xc

      subroutine update_hspectrum
c      *********************************************
c      computes spectrum of the modified hamiltonian
c      *********************************************
       
       double precision, allocatable :: tmp_h0mat(:,:), tmp_utrans(:,:) 
       integer ispin, n, m, p, tmpfile, ierr
   
       allocate(tmp_h0mat(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE update_hspectrum]: <tmp_h0mat> allocation failure'
       end if

       allocate(tmp_utrans(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE update_hspectrum]: <tmp_utrans> allocation failure'
       end if
   
       print '(/,1x,a,/)', 
     &       '<HUBBARD-U>: COMPUTING SPECTRUM OF THE MODIFIED HAMILTONIAN >>>'        
            
       do ispin = 1, nspin

c       solving eigenvalue problem
        forall(n=1:nsaos,m=1:nsaos)
         tmp_h0mat(n,m) = h0mat(n,m,ispin)
        end forall

        print '(2x,a,i1,a,$)', 'ispin = ', ispin, ' ...'
        call cpu_time(stime) 
        call realspectrum(tmp_h0mat,mo_en(:,ispin),tmp_utrans,nsaos) 

        call cpu_time(ftime)
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts * 60.0d0
	print *, 'done'
        print '(2x,a,i3,a,f5.2,a,/)',
     &       'time spent:', mnts, ' min ', secs, ' sec'   

c       save a copy of <mo_ort> to the global array <mo_orth>
        if (do_pop.or.do_lmpop.or.no_leads) then
         if (.not.allocated(mo_orth)) then
           allocate(mo_orth(nsaos,nsaos,nspin),stat=ierr)
           if (ierr.ne.0) then
           print *
           stop
     &     '[SUBROUTINE update_hspectrum]: <mo_orth> allocation failure'
           end if
	 end if 
         forall (n=1:nsaos,p=1:nsaos)
     &    mo_orth(n,p,ispin) = tmp_utrans(n,p)

c        back-transform to the non-orthogonal basis >>>
c        mo_coeff = smat12inv * mo_ort.
         call dble_ab('n','n',nsaos,smat12inv,mo_orth(:,:,ispin),mo_coeff(:,:,ispin))

c        testing new routine >>>>.
         call savemos(mo_en(:,ispin),mo_coeff(:,:,ispin),ispin)
c
        end if ! do_pop .or. do_lmpop .or. no_leads
 
       end do ! ispin 
       print '(1x,a)', '<HUBBARD-U>: <<< DONE WITH MODIFIED SPECTRUM'
     
c      print out spectrum to external file, if required
c      if (testing.or.do_pop.or.do_lmpop) then
       if (testing.or.zspectrum) then
        do ispin = 1, nspin

         tmpfile = 31
         if (nspin == 1) then
          open(tmpfile,file='mos.corr.tmp',status='unknown',iostat=ierr)
          if (ierr.ne.0) then
           stop ': can not open mos.corr.tmp!'
          end if
          write(tmpfile,'(a)') '$mos.corr.tmp'
         else  ! spin-polarized case
          if (ispin == 1) then
           open(tmpfile,file='alpha.corr.tmp',status='unknown',iostat=ierr)
           if (ierr.ne.0) then
            stop ': can not open alpha.corr.tmp!'
           end if
           write(tmpfile,'(a)') '$alpha.corr.tmp'
          else
           open(tmpfile,file='beta.corr.tmp',status='unknown',iostat=ierr)
           if (ierr.ne.0) then
            stop ': can not open beta.corr.tmp!'
           end if
           write(tmpfile,'(a)') '$beta.corr.tmp'
          end if
         end if
         write(tmpfile,'(a,i5)') '#nsaos=', nsaos

         do n = 1, nsaos
          write(tmpfile,'(i6,4x,e20.14)') n, mo_en(n,ispin) 	
 	 end do

         write(tmpfile,'(a)') '$end'        
	 close(tmpfile)

        end do ! ispin
	
       end if ! testing
       
       deallocate(tmp_utrans,tmp_h0mat,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE dv_xc]: ',
     &           'impossible to deallocate temporary arrays'
        print *, 'nevertheless, proceed further ...'
       end if
     
      end subroutine update_hspectrum

      double precision function u_energy(u,rho_mm,isize)
c      *************************************************
c      computes correction to the lda(gga) energy ;
c      see S.L.Dudarev et al. Phys.Rev.B 57, 1505 (1998)    
c      *************************************************
       double precision, intent (in) :: u  
       double precision, intent (in) :: rho_mm(:,:,:)
       integer, intent (in) :: isize
      
       integer          ispin, m, m1
       double precision tmpe1, tmpe2
      
       tmpe1 = 0.0d0
       tmpe2 = 0.0d0
       do ispin = 1, nspin
        do m = 1, isize
         tmpe1 = tmpe1 + rho_mm(m,m,ispin) 
	 do m1 = 1, isize
 	  tmpe2 = tmpe2 + rho_mm(m,m1,ispin) * rho_mm(m1,m,ispin)
	 end do
        end do      
       end do ! ispin
       u_energy = 0.5d0*u*(tmpe1 - tmpe2)/hartree

      end function u_energy

      subroutine save_hion_data(hmatrix,dmatrix,msize,ofilename)
c     ***************************************
c     saves ion's hamiltonian and its density
c     matrix to external file      
c     ***************************************
       double precision, intent (in) :: hmatrix(:,:,:), dmatrix(:,:,:)
       integer, intent (in)          :: msize
       character (50), intent (in)  ::  ofilename
      
       integer ofile, ierr, ispin, m, n
      
       ofile = 80  
       open(ofile,file=trim(ofilename),status='unknown',iostat=ierr)
       if (ierr.ne.0) then
         print '(/,1x,a,a,a)',
     &         'can not open file "', trim(ofilename), '" for writing'
         print '(a,1x)',
     &    'please, check your directory content or access rights'
         print '(a,1x,/)',
     &    '3d-ion(s) hamiltonian and its density matrix will not be saved'
       else
c       do all job here >>>
        print '(/,1x,a,$)', ' 3d-ion(s) data ... '
        write(ofile,'(a)') '$3d-ion(s).data'
        write(ofile,'(a,i2)') '#msize: ',msize 
        write(ofile,'(a,i2)') '#nspin: ',nspin 
        write(ofile,'(a,f16.10,a)') '#efermi: ', efermi, ' H' 
        write(ofile,'(a)') '#hamiltonian'
c       saving hamiltonian
        do ispin = 1, nspin
         write(ofile,'(a,i1)'), '#ispin: ', ispin
	 do n = 1, msize
	  do m = 1, n         
	   write(ofile,'(1x,D20.14,$)') 
     &	               0.5d0*(hmatrix(n,m,ispin)+hmatrix(m,n,ispin))	 
	  end do
	  write(ofile,*) 
	 end do 
        end do 
c       saving density matrix
        write(ofile,'(a)') '#density-matrix'
        do ispin = 1, nspin
         write(ofile,'(a,i1)') '#ispin: ', ispin
	 do n = 1, msize
	  do m = 1, n
           write(ofile,'(1x,D20.14,$)') 
     &	               0.5d0*(dmatrix(n,m,ispin)+dmatrix(m,n,ispin))	 
	  end do
	  write(ofile,*) 
	 end do 
        end do 
        write(ofile,'(a)') '$end'
        close(ofile)
        print '(a)', 'are saved to external file'
       end if

      end subroutine save_hion_data

      subroutine read_hion_data(hmatrix,dmatrix,msize,ofilename)
c     ***********************************************
c     saves ions' hamiltonian and its density
c     matrix for the target system from external file      
c     ***********************************************
       double precision, intent(out) :: hmatrix(:,:,:), dmatrix(:,:,:)
       integer, intent (in)          :: msize
       character(50), intent (in)    :: ofilename
       
       character(64) tmpword
       character(1024) tmpfmt
       integer ofile, ierr, ispin, jspin, m, n, tmpsize, tmpnspin
      
       print '(/,1x,a,/)', 
     &       '<HUBBARD-U>: 3d-ION(S) DATA ARE READ FROM EXTERNAL FILE >>> '
       ofile = 80  
       open(ofile,file=trim(ofilename),status='old',action='read',iostat=ierr)
       if (ierr.ne.0) then
         print '(1x,a,a,a)',
     &         'can not open file "', trim(ofilename), '" for reading'
         print '(1x,a)',
     &         'please, check your directory content or access rights'
         stop  ' : i will stop now'
       else
c       do all job here >>>
        read(ofile,*)      tmpword
        print '(1x,a)',    trim(tmpword)
        read(ofile,*)      tmpword, tmpsize
        print '(1x,a,i4)', trim(tmpword), tmpsize
	if (tmpsize.ne.msize) then
         print '(/,1x,a,/,a)',
     &         'dimension of the correlated subspace for the target system',
     &         'is different from yours: it seems, you"ve made a mistake somewhere'       
         stop ': i will quit now.'
	end if
        read(ofile,*)      tmpword, tmpnspin
        print '(1x,a,i4)', trim(tmpword), tmpnspin
        if (tmpnspin.ne.nspin) then
         print '(/,1x,a,/,a)',
     &         'number of different spin channels for the target system',
     &         'is different from yours: it seems, you"ve made a mistake somewhere'       
         stop ': i will quit now.'
	end if

        read(ofile,*) tmpword, ref_efermi
        print '(1x,a,f16.10,a)', trim(tmpword), ref_efermi, ' H' 
	
c       reading the hamiltonian of the target system
        read(ofile,*)   tmpword
        print '(1x,a)', tmpword
        do jspin = 1, nspin
         read(ofile,*), tmpword, ispin
         print '(1x,a,i2)', trim(tmpword), ispin
	 do n = 1, msize
           tmpfmt = '(1x,D20.14'
           do m = 2, n 
	    tmpfmt = trim(tmpfmt)//',1x,D20.14'
	   end do
	   tmpfmt = trim(tmpfmt)//')'
           read(ofile,trim(tmpfmt)) (hmatrix(n,m,ispin),m=1,n)	 
c          read(ofile,*) (hmatrix(n,m,ispin),m=1,n)	 
	   do m = 1, n
	    hmatrix(m,n,ispin) = hmatrix(n,m,ispin)
	    write(6,'(1x,D20.14,$)') 
     &	               0.5d0*(hmatrix(n,m,ispin)+hmatrix(m,n,ispin))	 
	   end do
	   print *
	 end do ! n 
        end do  ! jspin

c       reading density matrix of the target system
        read(ofile,*)   tmpword
        print '(1x,a)', trim(tmpword)
        do jspin = 1, nspin
         read(ofile,*),     tmpword, ispin
         print '(1x,a,i2)', trim(tmpword), ispin
	 do n = 1, msize
           tmpfmt = '(1x,D20.14'
           do m = 2, n 
 	    tmpfmt = trim(tmpfmt)//',1x,D20.14'
 	   end do
	   tmpfmt = trim(tmpfmt)//')'
           read(ofile,trim(tmpfmt)) (dmatrix(n,m,ispin),m=1,n)	 
c          read(ofile,*) (dmatrix(n,m,ispin),m=1,n)	 
	   do m = 1, n
	    dmatrix(m,n,ispin) = dmatrix(n,m,ispin)
	    write(6,'(1x,D20.14,$)') 
     &	               0.5d0*(dmatrix(n,m,ispin)+dmatrix(m,n,ispin))	 
	   end do
	   print *
	 end do  ! n 
        end do  ! jspin
	
        print '(/,1x,a)', '<HUBBARD-U>: <<< DONE WITH THE ION(S) DATA'
       end if

      end subroutine read_hion_data

      subroutine dealloc_ldau
c      ***************************      
c      deallocate temporary arrays
c      ***************************      
       implicit none 
       integer ierr

c       print '(/,a)', ' deallocating temporary arrays ...'
c       print '(a,$)', ' >>> ion_res_index ... '
       deallocate(ion_res_index,stat=ierr)
c       print '(a,i4)', 'stat =', ierr
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE dealloc_ldau]: ',
     &           'impossible to deallocate <ion_res_index>'
         print *, 'nevertheless, proceed further ...'
       end if

       if (allocated(dh_xc)) deallocate(dh_xc,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE dealloc_ldau]: ',
     &      'impossible to deallocate <h_ion>'
         print *, 'nevertheless, proceed further ...'
       end if

       deallocate(h_ion,ep_ion,umat_ion,rho_ion,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE dealloc_ldau]: ',
     &      'impossible to deallocate <h_ion>, or <ep_ion>, or <umat_ion>, or <rho_ion>'
         print *, 'nevertheless, proceed further ...'
       end if

       if (.not.read_hion) then
c       print '(a,$)', ' >>> rho ... '
        deallocate(rho,stat=ierr)
c       print '(a,i4)', 'stat =', ierr
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE dealloc_ldau]: ',
     &         'impossible to deallocate <rho> '
         print *, 'nevertheless, proceed further ...'
        end if
       end if	

c      print '(a,$)', ' >>> omat, omat12, omat12inv, ... '
       deallocate(omat,omat12,omat12inv,stat=ierr)
c      print '(a,i4)', 'stat =', ierr
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE dealloc_ldau]: ',
     &         'impossible to deallocate overlap matricies'
         print *, 'nevertheless, proceed further ...'
       end if

c      print '(a,$)', ' >>> d_weight, dlm_weight ... '
       deallocate(d_weight,dlm_weight,stat=ierr)
c      print '(a,i4)', 'stat =', ierr
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE dealloc_ldau]: ',
     &         'impossible to deallocate <d_weight> or <dlm_weight>'
         print *, 'nevertheless, proceed further ...'
       end if
c      print *

       if (read_hion) then
c       print '(a,$)', ' >>> h_trgt ... '
        deallocate(h_trgt,stat=ierr)
c       print '(a,i4)', 'stat =', ierr
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE dealloc_ldau]: ',
     &         'impossible to deallocate <h_trgt> or <rho_trgt> '
         print *, 'nevertheless, proceed further ...'
        end if
c       print * 
       end if
 
      end subroutine dealloc_ldau

      end module hubbardu       

