c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           jun 2008
c     last revision:  jan 2012
c###########################################################

      module neq_density_matrix
c      ********************************************      
c      contains a set of subroutines for evaluation  
c      of the non-equilibrium density matrix
c      ********************************************
       use globalvars
       use read_externals
       use tools
       use domainwall
       
       implicit none
       
       double precision, private, parameter :: ef_small_diff = 0.1d0
       
c      a full set of molecular orbitals (alpha & beta channels); dim = 2*nsaos       
       type(spectrum), allocatable :: fullsetmo(:)
   
       complex(8), private, allocatable :: 
     &                         left_mpp(:,:),  right_mpp(:,:)

c      pseudo-occupation numbers:  occ(iorb,ispin) 
       double precision, allocatable :: occn(:,:)        
       
       real(4), private :: stime, ftime, secs
       integer, private :: mnts

c      double precision, private, parameter :: myprec = 1.0d-14

       contains

       subroutine make_dmat
       implicit none
c      **********************************************************
c      a higher level routine, which controls calculation 
c      of the non-equilibrium density matrix :
c      -- calls efermi search --> finds dmat in the orthogonal basis;
c      -- transforms calculated density-matrix to the non-orthogonal 
c         basis set representation;
c      -- saves result to external files "dmat" & "smat"
c      **********************************************************

c       temporary arrays           
        double precision, allocatable :: tmp_dmat(:,:), tmp_sd(:,:)
      	integer ispin, ierr,  n, m

      	double precision tmp_nelectr
	
        allocate(neqdmat(nsaos,nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &   '[SUBROUTINE make_dmat]: <neqdmat> allocation failure'
        end if
       
c       searching for the fermi-energy; evaluate dens.matrix

        if (.not.fixed_efermi) then
         call searchef(efermi,neqdmat) 
        else

c        allocate an array with occupation numbers.
         if (.not.allocated(occn)) then
          allocate(occn(nsaos,nspin),stat=ierr)
          if (ierr.ne.0) then
           stop
     &     '[SUBROUTINE searchef]: <occn> allocation failure'
          end if
         end if 

         print '(/,2x,a)',
     &         '-------------------------------------------------------'
         print '(2x,a)', 'assuming efermi is fixed: no iterations !'
         print '(2x,a,f16.12,a)', 'EFermi = ', efermi, ' H'
         call compute_densmat(efermi,bias/hartree,neqdmat,tmp_nelectr,.true.)
         print '(/,2x,a,f16.12,a,$)', 'EFermi = ', efermi, ' H'
         print '(6x,a,f17.12)', 'Ne = ', tmp_nelectr
         print '(2x,a)',
     &         '-------------------------------------------------------'
        end if 

c       testing >>>
c       in case of domain-wall solution, symmetrize the density matrix
ccc        if (symmetrize) then
ccc         call updatedmat(neqdmat,1)
ccc        end if ! symmetrize
c       <<< testing

c       prints out atomic loewdin charges/magn.moments
        if (qoutput) call qprint(neqdmat)

c       neqdmat --> smat12inv * neqdmat * smat12inv

c       a local copy of dens.matrix for a given spin
        allocate(tmp_dmat(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &   '[SUBROUTINE make_dmat]: <tmp_dmat> allocation failure'
        end if

        allocate(tmp_sd(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &   '[SUBROUTINE make_dmat]: <tmp_sd> allocation failure'
        end if

        print '(/,1x,a,$)', 
     &	'transform density-matrix to non-orthogonal basis ...'
        call cpu_time(stime)
        do ispin = 1, nspin
c        put a local copy of dens.matrix to temporary array  >>>
         forall (n=1:nsaos,m=1:nsaos) 
	   tmp_dmat(n,m) = neqdmat(n,m,ispin)
	 end forall	
	
c        transform dens.matrix to the non-orthogonal basis
c        tmp_dmat --> smat12inv * tmp_dmat * smat12inv
c        1st step: tmp_sd = sma12inv * tmp_dmat         
         call dble_ab('n','n',nsaos,smat12inv,tmp_dmat,tmp_sd)
c        2nd step: tmp_dmat = tmp_sd * smat12         
         call dble_ab('n','n',nsaos,tmp_sd,smat12inv,tmp_dmat)

c        save a local copy of dens.matrix >>>
         forall (n=1:nsaos,m=1:nsaos)
	   neqdmat(n,m,ispin) = tmp_dmat(n,m)
	 end forall	
	
	end do ! ispin   
        print *, 'done'
        call cpu_time(ftime)      
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts*60.0d0
        print '(1x,a,i3,a,f5.2,a)',
     &        'time spent:', mnts, ' min ', secs, ' sec'

c       in case of domain-wall solution, symmetrize the density matrix
c       if (symmetrize) then
c	 call updatedmat(neqdmat,1)
c       end if ! symmetrize

c       admix new density with the previous one, 
c       save result to external file
        if (.not.(mol_subsystem.or.do_ldos3d.or.do_ldos3d_ewindow)) then
         call mixdensmat(neqdmat)
        end if

c       all done: deallocate arrays >>>
c       deallocate(tmp_dmat,tmp_sd,neqdmat,stat=ierr)
        deallocate(tmp_dmat,tmp_sd,stat=ierr)
        if (ierr.ne.0) then
          print '(/,a,a)', ' [SUBROUTINE make_dmat]: ',
     &             ' impossible to deallocate temp. arrays'
          print *, ' anyway, proceed further ...'
        end if
      
       end subroutine make_dmat

       subroutine mixdensmat(cdmat)
c       **********************************************
c       admixes the new density with the previous one, 
c       saves result to external file
c       **********************************************

c       input:  density matrix <cdmat> taken from current (c) iteration
        double precision, allocatable :: cdmat(:,:,:)	

c       density matrix from previous iteration     
        double precision, allocatable ::  prvdmat(:,:,:)	

        double precision tmp_chdens, tmp_spdens, ! tmp_trace, 
     &	                 diffnorm, tmpdiff    
	integer ispin, ierr, dmatfile, smatfile, n, m, k 
        integer, parameter :: monitor = 100
	logical  :: readdens, convflag = .false. 

c       next, the new density matrix is admixed with the one from 
c       the previous iteration; case of the 1st iteration is excluded
        readdens = .true.
        if (iter > 1) then
c        reading old density matrix
         dmatfile = 31
         open(dmatfile,file=charge_dens_file,status='old',
     &                            action='read',iostat=ierr)
         if (ierr.ne.0) then
          print '(/,a,a,a)',
     &     ' can not open file "', trim(charge_dens_file), '" for reading'
          print *,
     &     'density matrices will not be admixed'
          readdens = .false.
	  close(dmatfile)
	 end if

         if (readdens) then
	  smatfile = 32
          open(smatfile,file=spin_dens_file,status='old',
     &                          action='read',iostat=ierr)
          if (ierr.ne.0) then
           print '(/,a,a,a)',
     &      ' can not open file "',trim(spin_dens_file), '" for reading'
           print *,
     &      'density matrices will not be admixed'
           print *
           readdens = .false.
	   close(smatfile) ; close(dmatfile)
          end if
	 end if
     
         if (readdens) then
          allocate(prvdmat(nsaos,nsaos,nspin),stat=ierr)
          if (ierr.ne.0) then
          print * 
          stop
     &     '[SUBROUTINE mixdensmat]: <prvdmat> allocation failure'
          end if
          prvdmat = 0.0d0

          print '(/,1x,a,$)', 'admixing density-matrices ...'
          do n = 1, nsaos
	   do m = 1, n
	    if (nspin.eq.1) then
             read(dmatfile,*) tmp_chdens 
	     prvdmat(n,m,1) = tmp_chdens/2.0d0	   
             if (n.ne.m) prvdmat(m,n,1) = tmp_chdens/2.0d0
	    else
             read(dmatfile,*) tmp_chdens	   
             read(smatfile,*) tmp_spdens
	     prvdmat(n,m,1) = (tmp_chdens + tmp_spdens)/2.0d0	   
	     prvdmat(n,m,2) = (tmp_chdens - tmp_spdens)/2.0d0	   
             if (n.ne.m) then
	      prvdmat(m,n,1) = prvdmat(n,m,1)
	      prvdmat(m,n,2) = prvdmat(n,m,2)
	     end if
	    end if
	   end do 
	  end do

          close(dmatfile) ; close(smatfile)

c         admixing density matrices ;
c         computing norm of the deviation matrix ...
          diffnorm = 0.0d0
          do ispin = 1, nspin
           do m = 1, nsaos
	    do n = 1, nsaos
	     tmpdiff = cdmat(n,m,ispin) - prvdmat(n,m,ispin)
	     diffnorm = diffnorm + tmpdiff*tmpdiff
	     cdmat(n,m,ispin) =  
     &	                  dmix * cdmat(n,m,ispin) + 
     &                    (1.0d0-dmix) * prvdmat(n,m,ispin)    
	    end do
	   end do
	  end do
          diffnorm = dsqrt(diffnorm)/nsaos

          print *, 'done' 
          print '(/,a,e10.4)', 
     &    ' |newdens-dens| : (deviation norm)/nsaos = ',diffnorm

c         put a new value of density norm to <tcontrol>
          call updatetcontrol(trim(tcntrl_file_name),
     &	                      '$dnorm',0,diffnorm)

c         checking whether density matrix is converged or not >>>
c         -- here 'dnorm' is a deviation norm of the dens.matrix found at 
c            previous iteration and having been read from file <tcontrol>
c         -- 'densconv' is conv.criterium specified by keyword $densconv
          if (dnormfound) then
	   convflag = ( diffnorm < densconv ) 
c     optional >>>
c	   convflag = ( diffnorm < densconv ) .and. 
c     &	              ( dabs(dnorm-diffnorm) < fctr*densconv )
c     where >>>    -- last condition means that a difference between two 
c                     dens.matrix deviation norms should not be quite large 
c                     so that abrupt jumps are avoided (which are likely not 
c                     consistent with convergence)      
c                  -- 'fctr' is defined in <globalvars.f>
           if (convflag) then
	    print '(/,a)',' DENSITY CONVERGENCE CRITERIUM IS SATISFIED !'
	   end if	  
     	  end if	 

          deallocate(prvdmat,stat=ierr)
          if (ierr.ne.0) then
           print '(/,a,a)', ' [SUBROUTINE mixdensmat]: ',
     &             ' impossible to deallocate array <prvdmat>'
           print *, ' anyway, proceed further ...'
          end if
	 
	 end if ! readdens	
	end if ! iter > 1
             
c       store density matrix on the hard drive >>>
        dmatfile = 31
        open(dmatfile,file=charge_dens_file,status='replace',
     &                         	  action='write',iostat=ierr)
        if (ierr.ne.0) then
          print '(/,a,a,a)',
     &     ' can not open file "',trim(charge_dens_file), '" for writing'
          print *,
     &     'please, check your directory content or access rights'
          print *
          stop ': transport module is terminated now'
	end if

        smatfile = 32
        open(smatfile,file=spin_dens_file,status='replace',
     &                          action='write',iostat=ierr)
        if (ierr.ne.0) then
          print '(/,a,a,a)',
     &     ' can not open file "',trim(spin_dens_file), '" for writing'
          print *,
     &     'please, check your directory content or access rights'
          print *
          stop ': transport module is terminated now'
        end if

        print '(/,1x,a,$)', 'saving density-matrices : | '
        if (nspin.eq.2) then
c	 spin-polarized case
         k = 0 
         do n = 1, nsaos
	  k = k + 1 
          if ( mod(k,monitor)==0 ) print '(a,$)', '='  ! monitor
	  do m = 1, n
c          charge density matrix >>>
           write(dmatfile,1000) cdmat(n,m,1)+cdmat(n,m,2)
c          spin density matrix >>>
           write(smatfile,1000) cdmat(n,m,1)-cdmat(n,m,2)
	  end do
	 end do
	else
c        non-spin-polarized case
         k = 0 
         do n = 1, nsaos
	  k = k + 1 
          if ( mod(k,monitor)==0 ) print '(a,$)', '='  ! monitor
	  do m = 1, n
c          charge density matrix >>>
           write(dmatfile,1000) 2.0d0*cdmat(n,m,1)
c          spin density matrix >>>
           write(smatfile,'(a)') ' 0.0D+00'
	  end do
	 end do
	end if ! nspin
        print *, '| done' 
	
        close(dmatfile) ; close(smatfile)
       
c1000   format(1x,d22.16)       
1000   format(1x,d26.20)       
       end subroutine mixdensmat

       subroutine searchef(efio,densmtrx)
c       ************************************
c       searches for efermi(efio), so that 
c       trace[D(efio)] = number of electrons      
c       ************************************
        implicit none
        double precision              :: efio    ! input/ouput
        double precision, intent(out) :: densmtrx(:,:,:)

        double precision efguess0, nelectr0, efguess1, nelectr1,
     &	                 diff_mu, dmu, first_ef_guess 
        integer :: n, nn, ierr, iterno = 0
c    &             tmpfile
c       variables for the self-consistent cycle   
        double precision tmpef0, tmpef1, tmpef2, 
     &	                 n0, n1, n2, dn0, dn1, deriv 
	logical :: convflag, tflag = .false., dmatdone = .false.
	double precision, parameter :: preconv = 1.0d+2 
        
c       temporary set of molecular orbitals (just alpha + beta channels)
c       which are not properly ordered 
        type(spectrum), allocatable :: tmpset(:)
        integer, allocatable :: indarr(:)
c       character chann
   
        allocate(fullsetmo(2*nsaos),stat=ierr)
        if (ierr.ne.0) then
        stop
     &   '[SUBROUTINE searchef]: <fullsetmo> allocation failure'
        end if
c       prepair full set of molecular orbitals
        print '(/,1x,a)', 'SEARCHING FOR CHEMICAL POTENTIAL(S) >>> '
        print '(/,2x,a,$)', 'merge mo-sets of two spin channels ...'
        if (nspin.eq.1) then
	 nn = 0
	 do n = 1, nsaos
	  nn = nn + 1
          fullsetmo(nn)%en      = mo_en(n,1)
          fullsetmo(nn)%spin    = 1  
          fullsetmo(nn)%moindex = n  
	  nn = nn + 1
          fullsetmo(nn)%en      = mo_en(n,1)
          fullsetmo(nn)%spin    = 2  
          fullsetmo(nn)%moindex = n  
	 end do
	else
c        nspin=2: 
c        first, we merge two spin channels
         allocate(tmpset(2*nsaos),stat=ierr)
         if (ierr.ne.0) then
         stop
     &    '[SUBROUTINE searchef]: <tmpset> allocation failure'
         end if
         allocate(indarr(2*nsaos),stat=ierr)
         if (ierr.ne.0) then
         stop
     &    '[SUBROUTINE searchef]: <indarr> allocation failure'
         end if
         nn = 0
         do n = 1, nsaos
	  nn = nn + 1
          tmpset(nn)%en      = mo_en(n,1)
          tmpset(nn)%spin    = 1  
          tmpset(nn)%moindex = n  
	  nn = nn + 1
          tmpset(nn)%en      = mo_en(n,2)
          tmpset(nn)%spin    = 2  
          tmpset(nn)%moindex = n  
         end do
c        second, we are sorting tmp-array vs mo-energies: 
c        output --> fullsetmo
	 call sort_moset(tmpset,fullsetmo,2*nsaos,indarr)

c        testing
c         tmpfile=55
c         write(tmpfile,'(a)') '$merged spectrum'
c         do n = 1, 2*nsaos
c	  if (fullsetmo(n)%spin.eq.1) then ; chann = 'a'
c	                              else ; chann = 'b'
c	  end if  
c	  write(55,'(d24.14,4x,a1,4x,i5)') 
c     &     fullsetmo(n)%en, chann, fullsetmo(n)%moindex	  
c    	 end do
c         write(tmpfile,'(a)') '$end'

         deallocate(tmpset,indarr,stat=ierr)
         if (ierr.ne.0) then
          print '(/,a,a)', '  [SUBROUTINE searchef]: ',
     &             ' impossible to deallocate temporary arrays'
          print *, ' anyway, proceed further ...'
         end if
       
	end if ! nspin
        print *, 'done'
        if (.not.ecp) then
	 print '(2x,a,i5)', 'Overall number of electrons = ', nelectr
	else
	 print '(2x,a,i5)', 'Number of valence electrons = ', nelectr
	end if 
	print '(2x,a,e12.6)', 'Convergence criterium = ', nelconv

c       allocate an array with occupation numbers 
        if (.not.allocated(occn)) then
         allocate(occn(nsaos,nspin),stat=ierr)
         if (ierr.ne.0) then
          stop
     &     '[SUBROUTINE searchef]: <occn> allocation failure'
         end if
        end if 
		
c       *** first guess for fermi energy ***	

c       >>> caution! an fhi-aims user should be protected:
c                    she/he is NOT supposed to use $efermi keyword and
c                    have a freedom to modify a fermi level guess 'ef0'.
c                    default choice is always ef0 = (ehomo+elumo)/2

c       >>> currently, testing of the negf implementation based on 
c           fhi-aims is under way: to reach EF with less iterations, 
c           a constrain below (.not.aims_input) can be eliminated  

        if (efermi_guess .and. iter>1 .and. (.not.aims_input)) then
c        ... use for efguess0 a value from tcontrol 
	 efguess0 = efio
	 first_ef_guess = efguess0
	else
c        ... take for efguess0 = (ehomo+elumo)/2
         efguess0 = 0.5d0*(fullsetmo(nelectr)%en + 
     &	                   fullsetmo(nelectr+1)%en)
         first_ef_guess = efguess0
	end if	
	print '(/,2x,a)', 
     &	  '-------------------------------------------------------'
        iterno = iterno + 1
        print '(2x,a,i3)', 'iteration :',iterno
        print '(2x,a,f16.12,a)', 'EFermi = ', efguess0, ' H'
c       convert bias-voltage from ev to hartree      
        diff_mu = bias/hartree 
        call compute_densmat(efguess0,diff_mu,densmtrx,nelectr0,tflag)
        print '(/,2x,a,f16.12,a,$)', 'EFermi = ', efguess0, ' H'
	print '(6x,a,f17.12)', 'Ne = ', nelectr0         
        print '(2x,a)', 
     &	  '-------------------------------------------------------'

c       searching for better (second) fermi-energy guess >>> 
        call second_guess(efguess0,nelectr0,nelectr,nclstmo,efguess1)

c       take second guess for ef, and find the electron number >>>   
        print '(/,2x,a)', 
     &	  '-------------------------------------------------------'
        iterno = iterno + 1
        print '(2x,a,i3)', 'iteration :',iterno
        print '(2x,a,f16.12,a)', 'EFermi = ', efguess1, ' H'
        call compute_densmat(efguess1,diff_mu,densmtrx,nelectr1,tflag)
        print '(/,2x,a,f16.12,a,$)', 'EFermi = ', efguess1, ' H'
	print '(6x,a,f17.12)', 'Ne = ', nelectr1         
        print '(2x,a)', 
     &	  '-------------------------------------------------------'
c       if we are close to nelectr, dens.matrix will be 
c       calculated at the next iteration 
        if ( dabs(nelectr1-nelectr)< preconv*nelconv ) tflag = .true.

ccc     ### self-consistent cycle for fermi energy ###
ccc     ... not too much work has to be done here ...
        convflag = .false.
	do while (.not.convflag)
c        put left/right boundaries for search
         if (nelectr0 < nelectr1) then
          tmpef0 = efguess0 ;  tmpef1 = efguess1 
	  n0     = nelectr0 ;  n1     = nelectr1   
 	 else
          tmpef0 = efguess1 ;  tmpef1 = efguess0 
	  n0     = nelectr1 ;  n1     = nelectr0   
	 end if

         deriv = (n1-n0)/(tmpef1-tmpef0)
         dmu = (nelectr-n0)/deriv 
c        now, a better guess for efermi is >>>
         tmpef2 = tmpef0 + dmu
c        searching for the electron number with better efermi >>>	 
         print '(/,2x,a)', 
     &	   '-------------------------------------------------------'
         iterno = iterno + 1
         print '(2x,a,i3)', 'iteration :',iterno
         print '(2x,a,f16.12,a)', 'EFermi = ', tmpef2, ' H' 
         call compute_densmat(tmpef2,diff_mu,densmtrx,n2,tflag)
         print '(/,2x,a,f16.12,a,$)', 'EFermi = ', tmpef2, ' H'
	 print '(6x,a,f17.12)', 'Ne = ', n2         
         print '(2x,a)', 
     &	   '-------------------------------------------------------'
         if (tflag) dmatdone = .true.

c        testing >>>
c         print *, ' >>> input <<<'
c         print *, ' ef0 = ', tmpef0, ' n0 = ', n0
c         print *, ' ef1 = ', tmpef1, ' n1 = ', n1
c         print *, ' ef2 = ', tmpef2, ' n2 = ', n2

c        new guess is taken as one boundary for the next iteration
         efguess0 = tmpef2 ;  nelectr0 = n2
c        while the second one is the most closest 
c        among two previous boundaries
         dn0 = dabs(nelectr-n0)
         dn1 = dabs(nelectr-n1)
         if (dn0 < dn1) then
           efguess1 = tmpef0 ;  nelectr1 = n0
	 else
           efguess1 = tmpef1 ;  nelectr1 = n1
	 end if

c        testing
c         print *, ' >>> output <<<'
c         print *, ' efguess0 = ', efguess0, 
c     &	          ' nelectr0 = ', nelectr0
c         print *, ' efguess1 = ', efguess1, 
c     &	          ' nelectr1 = ', nelectr1

c        if we are close to convergence, change tflag to .true.
c        then full dens.matrix (not just trace) will be calculated
c        within the next iteration  
         if ( dabs(n2-nelectr)< preconv*nelconv ) tflag = .true.
c        checking for the next iteration 
         if ( dabs(n2-nelectr)< nelconv ) convflag = .true.
        end do
ccc     ### end of the self-consistent cycle ###
	
        print '(/,2x,a)','cycle converged!'

        print '(/,2x,a)', 
     &	   '======================================================='
        print '(2x,a,f16.12,a,$)', 'E_Fermi= ', efguess0, ' H'
        print '(6x,a,f17.12)', 'Ne = ', nelectr0         
        print '(2x,a)', 
     &	   '======================================================='

c       if density matrix is still not calculated, then build it >> 
        if (.not.dmatdone) then
         print '(/,2x,a)', 
     &	  '-------------------------------------------------------'
         print '(2x,a)', 'final iteration'
         print '(2x,a,f16.12,a)', 'EFermi = ', efguess0, ' H'
	 tflag = .true.
         call compute_densmat(efguess0,diff_mu,densmtrx,nelectr0,tflag)
         print '(/,2x,a,f16.12,a,$)', 'EFermi = ', efguess0, ' H'
	 print '(6x,a,f17.12)', 'Ne = ', nelectr0         
         print '(2x,a)', 
     &	  '-------------------------------------------------------'
        end if
 	
c       in case of fhi-aims: check, how far a way is found EF value
c       from its first guess; this check is only there for the non-sc
c       calculation 
 	
 	if (aims_input) then
 	 if (dabs(efguess0 - first_ef_guess) .ge. (ef_small_diff/hartree) ) then
c         print a warning message >>>
          print '(/,2x,a)',       
     &          'WARNING : a difference between the first guess for the'
          print '(2x,a,f10.6,a)',  
     &          '          fermi energy and its found value is ', dabs(efguess0 - first_ef_guess)*hartree, ' eV'
          print '(2x,a)',        
     &          '          which is more than 0.1 eV ! '
          print '(/,2x,a)', '          it seems, something goes wrong here. '
          print '(/,2x,a)', '          most likely, you are using incorrect parameters '
          print '(2x,a)',   '          for the real part of the self-energy ; to resolve '
          print '(2x,a)',   '          this problem, switch a flag $tune_rsigma" to "on" '
          print '(2x,a)',   '          and find a new set of parameters for your system '
  	 end if
 	end if
 	
c       put a new value of the fermi-energy to <tcontrol>
        call updatetcontrol(trim(tcntrl_file_name),'$efermi',0,efguess0)
c       ... and update output for efio >>> 
        efio = efguess0
     
        print '(/,1x,a)', '<<< DONE WITH CHEMICAL POTENTIALS'

c       >>> "pseudo-occupation numbers" were needed for testing only ...
c       >>> commented on jan 2012
c       if (store_occn) then
c       saving pseudo-occupation numbers to external file "occ" 
c         print '(/,1x,a,$)', 'saving occupation numbers ...'
c         call save_occn(trim(occ_nmbrs_file))
c         print *, 'done' 
c	end if 

        deallocate(fullsetmo,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', '  [SUBROUTINE searchef]: ',
     &            ' impossible to deallocate temp. array <fullsetmo>'
         print *, ' anyway, proceed further ...'
        end if

        deallocate(occn,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', '  [SUBROUTINE searchef]: ',
     &            ' impossible to deallocate array <occn>'
         print *, ' anyway, proceed further ...'
        end if

       end subroutine searchef


       subroutine second_guess(ef0,n0,nel,norb,ef1)
c       ***********************************************
c       estimate density of states (dos) ;
c       knowing dos, extrapolate \int_E (dos) and finds 
c       the second guess for the fermi-energy 
c       ***********************************************
        implicit none
c       !                 input guess for the 'efermi' and 
c       !                 corresponding electrons number
        double precision, intent(in)  :: ef0, n0  
c       !                 number of electrons in the system
        integer, intent(in)           :: nel     
c       ! number of orbitals taken for the dos-estimation
        integer, intent(in)           :: norb     
        double precision, intent(out) :: ef1

        integer          n_up, n_down, dnmo
        double precision dn, e_up, e_down, dos0, dmu
	
	n_up   = nel + 2*norb 
	n_down = nel + 1 - 2*norb 

        dnmo = n_up - n_down + 1  
c              that means, dnmo = 4*norb, which is ok 
        e_up  =(fullsetmo(n_up)%en + fullsetmo(n_up+1)%en)/2.0d0
        e_down=(fullsetmo(n_down)%en + fullsetmo(n_down-1)%en)/2.0d0
    
c       now, an estimate for dos is dn/de >>>   
        dos0 = dnmo/(e_up-e_down)        
        dn  = nel-n0
	dmu = dn/dos0
c       finally, output >>>         
        ef1 = ef0 + dmu

c       testing	
c       just print out info >>>
c        print *, ' === testing ==='
c        print *, ' E_up = ', e_up, '   E_down = ', e_down  
c        print *, ' dnmo = ', dnmo
c	 print *, ' dos0 = ', dos0
c	 print *, ' dn   = ', dn
c        print *, ' output >>>'
c	 print *, ' dmu = ', dmu
c	 print *, ' ef1 = ', ef1
c        print *, ' === testing ==='

c       all done
        return
       end subroutine second_guess

       subroutine compute_densmat(av_mu,diff_mu,densmtrx,nel,fintask)
        implicit none
c       **************************************************************
c       computes a non-equilibrium density matrix 
c       (in the orthogonal basis) assuming as
c       -- input efermi=av_mu & bias-voltage = diff_mu
c          be careful!  units are "hartree" 
c          while global variable "bias" is in ev ! 
c       -- output : array <desmat> = density matrix
c                   nel = overall number of electrons
c
c       optimization flag added >>>
c       -- input: 
c          fintask = .false.: only trace (electron number) is computed;
c                  = .true. : full density matrix is computed
c       *************************************************************
        double precision, intent(in)  :: av_mu, diff_mu
        double precision, intent(out) :: densmtrx(:,:,:)
	double precision, intent(out) :: nel 
	logical, intent (in)          :: fintask	    
	
        integer ierr, ispin, n, m, p
        double precision mu_left, mu_right
        complex(8), allocatable :: left_dmat(:,:), right_dmat(:,:),
     &	                           left_bm(:,:),   right_bm(:,:),
     &                             tmp_b(:,:)
        complex(8)  ctrace
	double precision tmptrace,tmpocc,tmpsum,occsum(2),tmp_nel(2) 
	
c       arrays to be used in [fppmat]
        allocate(left_mpp(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &   '[SUBROUTINE compute_densmat]: <left_mpp> allocation failure'
        end if
        allocate(right_mpp(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &   '[SUBROUTINE compute_densmat]: <right_mpp> allocation failure'
        end if

c       allocation of the temporary arrays
        allocate(left_dmat(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &  '[SUBROUTINE compute_densmat]: <left_dmat> allocation failure'
        end if
        allocate(right_dmat(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &  '[SUBROUTINE compute_densmat]: <right_dmat> allocation failure'
        end if

        allocate(left_bm(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &    '[SUBROUTINE init_densmat]: <left_bm> allocation failure'
        end if
        allocate(right_bm(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &    '[SUBROUTINE init_densmat]: <right_bm> allocation failure'
        end if
        allocate(tmp_b(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &    '[SUBROUTINE init_densmat]: <tmp_b> allocation failure'
        end if

c       print '(/,a)',' TESTING NON-EQUILIBRIUM DENSITY MATRIX >>>'
        if (aims_input) then 
         print '(/,2x,a)','evaluating a density matrix'
	else  
         print '(/,2x,a)','evaluating a non-equilibrium density matrix'
	end if 
c       set-up chemical potentials of the leads
	mu_left   = av_mu - 0.5d0*diff_mu
	mu_right  = av_mu + 0.5d0*diff_mu

        call cpu_time(stime)
        do ispin = 1, nspin
         if (nspin == 2) then
          if (ispin == 1) then
           print '(5x,a)', '** alpha channel **'
          else
           print '(5x,a)', '** beta  channel **'
          end if
         end if
c        compute m_pp'
         print '(2x,a,$)', '-- computing Mpp matrices ...'
 	 call mppmat(mu_left,mu_right,ispin)
         print *, 'done'

         if (fintask) then
          print '(2x,a,$)', '-- building up the density matrix ...'
         else
          print '(2x,a,$)', '-- computing trace of the density matrix ...'
	 end if
c        print '(/,2x,a)', '== final step: build up density matrix =='	
c        save local copy of b-matrix 
         forall(n=1:nsaos,m=1:nsaos)
	  tmp_b(n,m) = heigvec(n,m,ispin)
	 end forall       

c        print *, ' -- B * M  --'
	 call cmplx_ab('n','n',nsaos,tmp_b,left_mpp,left_bm)
	 call cmplx_ab('n','n',nsaos,tmp_b,right_mpp,right_bm)
c        print '(/,2x,a,/)', ' done with B * M '

         if (fintask) then

c         first we care about occupation numbers: 
          occsum(ispin) = 0.0d0
          do p = 1, nsaos
	   tmpocc = 0.0d0
	   do n = 1, nsaos
            tmpocc=tmpocc + dble(left_bm(n,p)*dconjg(tmp_b(n,p)) +
     &                          right_bm(n,p)*dconjg(tmp_b(n,p)) )
	   end do
           if ( nspin.eq.1 ) then  
	    occn(p,ispin) = 2.0d0*dabs(tmpocc)
	    occsum(ispin) = occsum(ispin) + occn(p,ispin)
	   else 
	    occn(p,ispin) = dabs(tmpocc)
            occsum(ispin) = occsum(ispin) + occn(p,ispin)
	   end if 	   
	  end do

c         calculation is converged: a full density matrix is required 
c         print *, ' -- B * M * B+ --'
	  call cmplx_ab('n','c',nsaos,left_bm,tmp_b,left_dmat)
	  call cmplx_ab('n','c',nsaos,right_bm,tmp_b,right_dmat)
c         print '(/,2x,a,/)', ' done with B * M * B+ '

c         save result to the output array
          do m = 1, nsaos
	   do n = 1, m
c           due to a hermitian density matrix and real-valued 
c           basis functions, one needs only a real part of dens.matrix
	    densmtrx(n,m,ispin) = dble(left_dmat(n,m) + right_dmat(n,m))
	    if (n.ne.m) densmtrx(m,n,ispin) = densmtrx(n,m,ispin)
	   end do
	  end do      
          tmptrace = 0.0d0
          do n = 1, nsaos
	   tmptrace = tmptrace + densmtrx(n,n,ispin)
	  end do
         else 
c         we compute just number of electrons, that is trace(densmat)
          ctrace  = czero
          do m = 1, nsaos
	   do n = 1, nsaos
	    ctrace = ctrace
     &	           + (left_bm(n,m)+right_bm(n,m)) * dconjg(tmp_b(n,m)) 
           end do
	  end do
c         take just a real piece :
          tmptrace = dble(ctrace)
c         put dummy value to the output array densmtrx
          forall (n=1:nsaos,m=1:nsaos) densmtrx(n,m,ispin) = 0.0d0
         end if ! fintask
	         
c        put result to the output variable
         if (nspin.eq.1) then
           nel = 2.0d0*tmptrace	 ! factor 2 is due to spin    
c          we need that for updating of the occupation numbers
	   tmp_nel(ispin) = nel 
         else 
	   tmp_nel(ispin) = tmptrace 
         end if
         print *, 'done'

	end do ! ispin

c       take care about output electron number for nspin=2
	if (nspin.eq.2) then
	  nel = tmp_nel(1) + tmp_nel(2) 
	end if

        call cpu_time(ftime)
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts*60.0d0
        print '(2x,a,i3,a,f5.2,a)',
     &        'time spent:', mnts, ' min ', secs, ' sec'

        deallocate(left_mpp,right_mpp,left_dmat,right_dmat,
     & 	                        left_bm,right_bm,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', '  [SUBROUTINE init_densmat]: ',
     &            ' impossible to deallocate temporary arrays'
         print *, ' anyway, proceed further ...'
        end if

c       here we updating occupation numbers >>>
        if (fintask) then
         do ispin = 1, nspin
          do p = 1, nsaos
	   occn(p,ispin) = 
     &	      occn(p,ispin) * tmp_nel(ispin)/occsum(ispin)
	  end do
          tmpsum = 0.0d0
          do p = 1, nsaos
c	   write(40+ispin,*) p, occn(p,ispin)
           tmpsum = tmpsum + occn(p,ispin)
    	  end do
c          print *
c          print *, ' sum of occupation numbers        : ', tmpsum
c          print *
         end do	! ispin  
        end if
c       print '(/,a)',' <<< DONE WITH NON-EQUILIBRIUM DENSITY MATRIX'

       end subroutine compute_densmat

       
       subroutine fppmat(fpp,mu,ispin)
c      *******************************************
c      computing f_pp'(mu)-matrix : see math-notes
c      *******************************************
        implicit none
c       output matrix
        complex(8), intent(out) :: fpp(:,:)
c       chemical potential of the reservoir
        double precision, intent(in) :: mu
c       spin index       
        integer, intent(in) :: ispin 
	integer  p, q

	do p = 1, nsaos
	 do q = 1, p
          fpp(p,q) = fppfunc(p,q,mu,ispin) 	 
	  if (p.ne.q) fpp(q,p) = dconjg(fpp(p,q))
 	 end do
	end do
       end subroutine fppmat
       
       complex(8) function fppfunc(p,q,mu,ispin)
        implicit none
        integer, intent(in)          :: p, q
        double precision, intent(in) :: mu
        integer, intent(in)          :: ispin

        double precision eps_p, eps_q, eta_p, eta_q
	complex(8)       zp, zq, tmp_sum, zdiff
	double precision atan_p, atan_q, ln_pq, tmp_frac
	
	zp = zpoles(p,ispin)  ;  zq = zpoles(q,ispin)  
	eps_p = mu - dble(zp) ;	 eps_q = mu - dble(zq) 
	eta_p = -dimag(zp)    ;  eta_q = -dimag(zq)

        tmp_frac=(eps_p*eps_p+eta_p*eta_p)/(eps_q*eps_q+eta_q*eta_q)
	ln_pq  = dlog(tmp_frac)
	atan_p = datan(eps_p/eta_p)
	atan_q = datan(eps_q/eta_q)
c       ione --> imaginary unit=(0.0d0,1.0d0); declared in <globalvars.f>
	tmp_sum = 0.5d0*ln_pq - ione*(pi + atan_p + atan_q) 

        zdiff = zp - dconjg(zq)
        fppfunc = tmp_sum/(zdiff*2.0d0*pi)
	
       end function fppfunc

       subroutine mppmat(mu_l,mu_r,ispin)
c      **************************************************************  
c      computing m_pp'(mu) = [b^{-1} *gamma* b+^{-1}]_pp' * f_pp'(mu)
c      for reference: see math-notes
c      **************************************************************  
        implicit none
c       chemical potential of left/right reservoirs
        double precision, intent(in) :: mu_r, mu_l
c       spin index       
        integer, intent(in) :: ispin 

        integer ierr, n,m

        complex(8), allocatable :: b_gl(:,:), b_gr(:,:), 
     & 	                           right_bgb(:,:), left_bgb(:,:),
     &                             tmp_invb(:,:),
     & 	                           right_fpp(:,:), left_fpp(:,:)

c       complex(8) ctmp
c       logical :: conjtest = .true.

c       allocation of temporary arrays
        allocate(b_gl(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &    '[SUBROUTINE mppmat]: <bg_l> allocation failure'
        end if
        allocate(b_gr(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &    '[SUBROUTINE mppmat]: <bg_r> allocation failure'
        end if

        allocate(tmp_invb(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &    '[SUBROUTINE mppmat]: <tmp_invb> allocation failure'
        end if

        allocate(left_bgb(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &    '[SUBROUTINE mppmat]: <left_bgb> allocation failure'
        end if
        allocate(right_bgb(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &    '[SUBROUTINE mppmat]: <right_bgb> allocation failure'
        end if

        allocate(right_fpp(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &    '[SUBROUTINE mppmat]: <right_fpp> allocation failure'
        end if
        allocate(left_fpp(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        print * 
        stop
     &    '[SUBROUTINE mppmat]: <left_fpp> allocation failure'
        end if


c       print '(/,2x,a)', '== computing m_pp-matrices =='

c       step 1: form f_pp'-matrices
c       print *, ' -- Fpp call --'
        call fppmat(left_fpp,mu_l,ispin)
        call fppmat(right_fpp,mu_r,ispin)

c       step 2: initialize b^{-1}*gamma
c       print *, ' -- B^{-1} Gamma --'
        forall(n=1:nsaos,m=1:nsaos)
	 tmp_invb(n,m) = invheigvec(n,m,ispin)
	end forall       
	do m = 1, nsaos
	 do n = 1, nsaos
          b_gl(n,m) = tmp_invb(n,m) * gamma_left(m)
          b_gr(n,m) = tmp_invb(n,m) * gamma_right(m)
	 end do
	end do

c       step 3: compute b^{-1} * gamma * b+^{-1}  	 
c       print *, ' -- B^{-1} Gamma B+^{-1} --'
	call cmplx_ab('n','c',nsaos,b_gl,tmp_invb,left_bgb)
	call cmplx_ab('n','c',nsaos,b_gr,tmp_invb,right_bgb)

c       step 4: initialize m_pp'
c       print *, ' -- Mpp --'
        do m = 1, nsaos
	 do n = 1, m
	  left_mpp(n,m)  = left_bgb(n,m)  * left_fpp(n,m)
c         hermicity of m_pp' is kept up to 1.e-15
c         for even better precision, one can make a handwork >>>   
          if (n.ne.m) then 
	   left_mpp(m,n) = dconjg(left_mpp(n,m))
     	  else
	   left_mpp(n,n) = dble(left_mpp(n,n)) + ione * 0.0d0
	  end if
	  
	  right_mpp(n,m) = right_bgb(n,m) * right_fpp(n,m)
          if (n.ne.m) then 
	   right_mpp(m,n) = dconjg(right_mpp(n,m))
     	  else
	   right_mpp(n,n) = dble(right_mpp(n,n)) + ione * 0.0d0
	  end if

	 end do
	end do

c        checking whether m_pp' matrix is hermitian or not >>>
c        do n = 1, nsaos
c	 do m = 1, n
c          checking
c           ctmp = left_mpp(n,m) - dconjg(left_mpp(m,n))
c	   if (cmplxabs(ctmp)>myprec) then
c	    conjtest = .false.
c	    write(70+ispin,*) n, m, ctmp
c	   end if	
c	 end do    	
c        end do
c	if (conjtest) then
c	 print *, 'TEST FOR Mpp is OK'
c        else
c	 print *, 'TEST FOR Mpp FAILS!'
c	end if	
c	stop
	
        deallocate(b_gl,b_gr,left_bgb,right_bgb,tmp_invb,
     & 	                  right_fpp,left_fpp,stat=ierr)

c       print '(2x,a)', '== done with m_pp-matrices =='
	
       end subroutine mppmat

       subroutine save_occn(occnfile)
c       **********************************************
c       saves occupation numbers to file external file       
c       **********************************************
        character(*), intent(in) :: occnfile
        integer extfile, p, ierr

        extfile = 44 
        open(extfile,file=occnfile,status='replace',
     &	                 action='write',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a,a)',
     &         ' can not open file "', trim(occnfile), '" for writing'
         print *,
     &    'please, check your directory content or access rights'
         print *
         stop ': transport module is terminated now'
        end if
    
        if (nspin.eq.2) then
         do p = 1, nsaos
c         write(extfile,2000) occn(p,1), occn(p,2) 
          write(extfile,*) occn(p,1), occn(p,2) 
         end do
        else 
c        nspin == 1
         do p = 1, nsaos
c         write(extfile,2002) occn(p,1), '0.0D+00'
          write(extfile,*) occn(p,1), '  0.0D+00'
         end do
        end if
        close(extfile)

c2000    format(1x,d20.14,6x,d20.14)
c2002    format(1x,d20.14,6x,a)
       end subroutine save_occn


       subroutine qprint(tmpdmat)
c       ***************************************
c       prints out loewdin charges/magn.moments
c       for all atoms in the cluster
c       ***************************************

c       input: a non-eq.dens.matrix in the orthogonal basis
        double precision, intent(in) :: tmpdmat(:,:,:)

        double precision, allocatable :: lwcharge(:,:)

        double precision dq, qtot(2) 
	integer iat, iatype, iorb, norb, ispin, nu, ierr

        allocate(lwcharge(num_atoms,nspin),stat=ierr)
        if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE qprint]: <lwcharge> allocation failure'
        end if
     
        if (nspin.eq.1) then
	 print '(/,1x,a)', '=== LOEWDIN POPULATION ANALYSIS ==='
        else
         print '(/,13x,a)', '=== LOEWDIN POPULATION ANALYSIS ==='
	end if


        qtot(1) = 0.0d0
        qtot(2) = 0.0d0
        do ispin = 1, nspin
	 forall(iat=1:num_atoms) lwcharge(iat,ispin) = 0.0d0
	end do

	do ispin = 1, nspin
         nu = 0 
         do iat = 1, num_atoms
	  iatype = atom(iat)%atype
	  norb = n_basis_func(iatype)
	  do iorb = 1, norb
	   nu = nu + 1
	   dq = tmpdmat(nu,nu,ispin)
	   lwcharge(iat,ispin) = lwcharge(iat,ispin) + dq
	   qtot(ispin) = qtot(ispin) + dq
	  end do ! iorb 
 	 end do ! iat
	end do ! ispin  

        if (nspin.eq.1) then

          print '(4x,a)', '-------------------------'
          print '(4x,a)', '      atom      charge   '
          print '(4x,a)', '-------------------------'
          do iat = 1, num_atoms
           print '(4x,i4,3x,a2,2x,f12.6)', 
     &	          iat,atom(iat)%symbol(1:2), 2.0d0*lwcharge(iat,1)
          end do
          print '(4x,a)', '-------------------------'
          print '(9x,a,2x,f12.6,/)', 'SUM:', 2.0d0*qtot(1)

	else

          print '(1x,a)', 
     &	   '-----------------------------------------------------------'
          print '(1x,a)', 
     &	   '    atom      a.charge    b.charge     charge   magn.moment'
          print '(1x,a)', 
     &	   '-----------------------------------------------------------'
          do iat = 1, num_atoms
           print '(1x,i4,3x,a2,1x,2f12.6,1x,2f12.6)', 
     &	     iat,atom(iat)%symbol(1:2),lwcharge(iat,1),lwcharge(iat,2),
     &       lwcharge(iat,1)+lwcharge(iat,2), 
     &       lwcharge(iat,1)-lwcharge(iat,2)
          end do
          print '(1x,a)', 
     &	   '-----------------------------------------------------------'
          print '(6x,a,1x,2f12.6,1x,2f12.6,/)', 'SUM:', 
     &     	  qtot(1), qtot(2), qtot(1)+qtot(2), qtot(1)-qtot(2)
	
        end if

        deallocate(lwcharge,stat=ierr)
        if (ierr.ne.0) then
          print '(/,a,a)', ' [SUBROUTINE qprint]: ',
     &            ' impossible to deallocate array <lwcharge>'
          print *, ' anyway, proceed further ...'
        end if

       end subroutine qprint

      end module neq_density_matrix
