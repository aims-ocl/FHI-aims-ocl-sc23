c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           august 2008
c     last revision:  jan 2012
c###########################################################

      module resigma
c      ********************************************      
c      searches for a real piece of the self-energy  
c      ********************************************
       use globalvars
       use read_externals
       use neq_density_matrix
       use exthamiltonian
       use tools
       
       implicit none
       contains

       subroutine tune_sigma(rguess)
c       ******************************************	
c       adjust a real piece of the self-energy so 
c       that for the given fermi energy Tr(DS) = N 
c       ******************************************	
        implicit none
  
c       input/output 
        double precision rguess
	
        integer ::       ierr, i, n, nn, iterno = 0
        double precision ef0, diff_mu, tmp_conv, dn_over_n, 
     &	                 rguess_save, best_rguess,
     &                   rguess0, nelectr0, rguess1, nelectr1, best_nelectr

c       variables for self-consistent cycle   
        double precision tmpf0, tmpf1, tmpf2, n0, n1, n2, dn0, dn1  
	logical          convflag, limit_exceeded

c       temporary set of molecular orbitals (just alpha + beta channels)
c       which are not properly ordered 
        type(spectrum), allocatable :: tmpset(:)
        integer,        allocatable :: indarr(:)

        allocate(neqdmat(nsaos,nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE tune_sigma]: <neqdmat> allocation failure'
        end if
   
        allocate(fullsetmo(2*nsaos),stat=ierr)
        if (ierr.ne.0) then
        stop
     &   '[SUBROUTINE tune_sigma]: <fullsetmo> allocation failure'
        end if

c       save input value of 'rguess':
        rguess_save = rguess         
	best_rguess = rguess

c       prepair a full set of molecular orbitals
        print '(/,1x,a)', 
     &	      '*****  SEARCHING FOR A REAL PART OF THE SELF-ENERGY  *****'
        print '(/,1x,a,$)', 'merge mo-sets of two spin channels ...'
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
     &    '[SUBROUTINE tune_sigma]: <tmpset> allocation failure'
         end if
         allocate(indarr(2*nsaos),stat=ierr)
         if (ierr.ne.0) then
         stop
     &    '[SUBROUTINE tune_sigma]: <indarr> allocation failure'
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

         deallocate(tmpset,indarr,stat=ierr)
         if (ierr.ne.0) then
          print '(/,a,a)', '  [SUBROUTINE tune_sigma]: ',
     &             ' impossible to deallocate temporary arrays'
          print *, ' anyway, proceed further ...'
         end if
       
	end if ! nspin
        print *, 'done'
	if (.not.ecp) then
	  print '(/,1x,a,i5)',  'Overall number of electrons Ne = ', nelectr
	else
	  print '(/,1x,a,i5)',  'Number of valence electrons Ne = ', nelectr
	end if 
	print '(1x,a,e12.6)', 'Convergence criterium: |dN/Ne| = ', nelcnv
		
c       *** define fermi level here ***	

c       >>> caution! an fhi-aims user should be protected:
c                    she/he is not supposed to use $efermi keyword and
c                    have a freedom to modify a fermi level guess 'ef0'.
c                    default choice is always ef0 = (ehomo+elumo)/2
        if (efermi_guess.and.(.not.aims_input)) then
c        ... use for ef0 a value from <tcontrol>
	 ef0 = efermi
	else
c        ... take for ef0 = (ehomo+elumo)/2: 
c            this is a default (recommended!) choice 
c            if keyword 'efermi' is not found in <tcontrol>
         ef0 = 0.5d0*(fullsetmo(nelectr)%en + 
     &	                   fullsetmo(nelectr+1)%en)	
	end if	

	print '(1x,a,f16.12,a)', 'Fermi energy = ',ef0,' H'
        diff_mu = bias/hartree 
	print '(1x,a,f16.12,a)', 'Bias voltage = ',diff_mu,' eV'

c       first guess for re.sigma
        print '(/,1x,a)', 
     &	  '-------------------------------------------------------'
        iterno = iterno + 1
        print '(1x,a,i3,/)', 'ITERATION :',iterno
        if (expert) then
         rguess0 = rsigma(1) / isigma(1)
        else
         rguess0 = rguess
        end if 
        call update_sigma(rguess0)
        print '(1x,a,f16.12,a)',
     &	      'Re.Sigma(1st-layer) = ',rguess0*isigma(1),' H'
        print '(1x,a,f16.12,a)',
     &	      'Re.Sigma(2nd-layer) = ',rguess0*isigma(2),' H'
        print '(1x,a,f16.12,a)',
     &	      'Re.Sigma(3rd-layer) = ',rguess0*isigma(3),' H'
        call hfull
        call hspectrum
        call compute_densmat(ef0,diff_mu,neqdmat,nelectr0,.false.)
c       print '(1x,a,f16.12,a)','Re.Sigma(1st-layer) = ',rsigma(1),' H'
        print '(/,1x,a,f16.12,a,$)', 'EFermi = ', ef0, ' H'
	print '(6x,a,f17.12)', 'Ne = ', nelectr0         
        print '(1x,a)', 
     &	  '-------------------------------------------------------'

        IF (dabs(nelectr0-nelectr)/dabs(nelectr+nelcnv) > nelcnv ) THEN
c       second guess for re.sigma
        print '(/,1x,a)', 
     &	  '-------------------------------------------------------'
        iterno = iterno + 1
        print '(1x,a,i3,/)', 'ITERATION :',iterno
        if (expert) then
         rguess1 = rfactor *  rsigma(1) / isigma(1) 
        else
         rguess1 = 2.0d0 * rguess
        end if 
        call update_sigma(rguess1)
        print '(1x,a,f16.12,a)',
     &	      'Re.Sigma(1st-layer) = ',rguess1*isigma(1),' H'
        print '(1x,a,f16.12,a)',
     &	      'Re.Sigma(2nd-layer) = ',rguess1*isigma(2),' H'
        print '(1x,a,f16.12,a)',
     &	      'Re.Sigma(3rd-layer) = ',rguess1*isigma(3),' H'
        call hfull
        call hspectrum
        call compute_densmat(ef0,diff_mu,neqdmat,nelectr1,.false.)
c       print '(1x,a,f16.12,a)','Re.Sigma(1st-layer) = ',rsigma(1),' H'
        print '(/,1x,a,f16.12,a,$)', 'EFermi = ', ef0, ' H'
	print '(6x,a,f17.12)', 'Ne = ', nelectr1         
        print '(1x,a)', 
     &	  '-------------------------------------------------------'

ccc     ### self-consistent cycle for re.sigma ###
ccc     ... not too much work has to be done here ...
        convflag = .false.
	limit_exceeded = .false.
	tmp_conv = 10.0d0
	best_nelectr = 0.0d0
	do while ((.not.convflag).and.(.not.limit_exceeded))
c        put left/right boundaries for search
         if (nelectr0 < nelectr1) then
          tmpf0 = rguess0  ;  tmpf1 = rguess1 
	  n0    = nelectr0 ;  n1    = nelectr1   
 	 else
          tmpf0 = rguess1  ;  tmpf1 = rguess0 
	  n0    = nelectr1 ;  n1    = nelectr0   
	 end if

c        now, a better guess for efermi is >>>
         tmpf2 = tmpf0 + (nelectr-n0) * (tmpf1-tmpf0) / (n1-n0) 

c        searching for electron number with better guess >>>	 
         print '(/,1x,a)', 
     &	   '-------------------------------------------------------'
         iterno = iterno + 1
	 print '(1x,a,i3,/)', 'ITERATION :',iterno
	 call update_sigma(tmpf2)
         print '(1x,a,f16.12,a)',
     &	       'Re.Sigma(1st-layer) = ',tmpf2*isigma(1),' H'
         print '(1x,a,f16.12,a)',
     &	       'Re.Sigma(2nd-layer) = ',tmpf2*isigma(2),' H'
         print '(1x,a,f16.12,a)',
     &	       'Re.Sigma(3rd-layer) = ',tmpf2*isigma(3),' H'
         call hfull
         call hspectrum
         call compute_densmat(ef0,diff_mu,neqdmat,n2,.false.)
c        print '(1x,a,f16.12,a)','Re.Sigma(1st-layer) = ',rsigma(1),' H'
         print '(/,1x,a,f16.12,a,$)', 'EFermi = ', ef0, ' H'
	 print '(6x,a,f17.12)', 'Ne = ', n2         
         print '(1x,a)', 
     &	   '-------------------------------------------------------'

c        new guess is taken as one boundary for next iteration
         rguess0 = tmpf2 ;  nelectr0 = n2
c        while the second one is the most closest 
c        among two previous boundaries
         dn0 = dabs(nelectr-n0)
         dn1 = dabs(nelectr-n1)
         if (dn0 < dn1) then
           rguess1 = tmpf0 ;  nelectr1 = n0
	 else
           rguess1 = tmpf1 ;  nelectr1 = n1
	 end if

         dn_over_n = dabs(n2-nelectr)/dabs(nelectr+nelcnv)
c        save current best estimates ...
         if (dn_over_n .lt. tmp_conv) then 
	  tmp_conv = dn_over_n
	  best_rguess = rguess0
	  best_nelectr = n2 
	 end if 
	 
c        checking whether we achieved convergence >>>
	 if ( dn_over_n < nelcnv ) then
	  convflag = .true.
c        check whether we have exceeded a limit of iterations >>>
         else if (iterno.ge.iter_limit) then	
	  limit_exceeded = .true. 
	 end if 
	 
        end do
ccc     ### end of the self-consistent cycle ###
   
        if (convflag) then
         print '(/,1x,a)','CYCLE CONVERGED!'
        else 
c        print out a warning message >>>
         print '(/,1x,a,i3,a)', 
     &	       'WARNING: a cycle has not converged within', iter_limit, ' iterations !'
         print *
         print '(1x,a,e12.6)', 
     &	       '         Smallest relative error in electron number is |dN/Ne|= ', tmp_conv 
         print '(1x,a,e12.6)', 
     &	       '         while a default convergence criterium was set to :     ', nelcnv 
         print *
         print '(1x,a)', 
     &	       '         Best found parameters (see below) for the model self-energy will be'
         print '(1x,a)', 
     &	       '         written to file <tcontrol>.'
c        replace output with the best found values >>>
         rguess0  = best_rguess
	 nelectr0 = best_nelectr
	end if
       
        END IF  ! first guess was OK

        print '(/,1x,a)', 
     &	   '======================================================='
        print '(1x,a,f16.12,a)',
     &	      'Re.Sigma(1st-layer) = ',rguess0*isigma(1),' H'
        print '(1x,a,f16.12,a)',
     &	      'Re.Sigma(2nd-layer) = ',rguess0*isigma(2),' H'
        print '(1x,a,f16.12,a)',
     &	      'Re.Sigma(3rd-layer) = ',rguess0*isigma(3),' H'
        print '(/,1x,a,f16.12,a)', 
     &	      '            E_Fermi = ', ef0, ' H'
        print '(1x,a,f16.11)', 
     &	      '                 Ne = ', nelectr0         
        print '(1x,a,3x,e13.7)', 
     &	      '            |dN/Ne| = ', dabs( (nelectr0 - nelectr)/dble(nelectr) )         
        print '(1x,a)', 
     &	   '======================================================='
     
        print '(/,1x,a)', 
     &        'parameters for the self-energy are written to <tcontrol>'

c       updating tcontrol-file >>>

c       update a group "$adjust_rsigma" in tcontrol-file
        if (.not.do_dmat_cycle) then
         call updatetcontrol(trim(tcntrl_file_name),'$adjust_rsigma',-1,0.0d0)
        else
         call updatetcontrol(trim(tcntrl_file_name),'$adjust_rsigma',1,0.0d0)
        end if 
        call updatetcontrol(trim(tcntrl_file_name),'rfactor',0,rfactor)
        call updatetcontrol(trim(tcntrl_file_name),'rguess',0,rguess_save)
        call updatetcontrol(trim(tcntrl_file_name),'nelcnv',0,nelcnv)
        call updatetcontrol(trim(tcntrl_file_name),'iterlimit',iter_limit,0.0d0)

c       update self-energy parameters
        forall (i = 1:deep) rsigma(i) = rguess0*isigma(i)
	
	call updatetcontrol(trim(tcntrl_file_name),'$s1r',0,rsigma(1))
	call updatetcontrol(trim(tcntrl_file_name),'$s1i',0,-isigma(1))

	call updatetcontrol(trim(tcntrl_file_name),'$s2r',0,rsigma(2))
	call updatetcontrol(trim(tcntrl_file_name),'$s2i',0,-isigma(2))

	call updatetcontrol(trim(tcntrl_file_name),'$s3r',0,rsigma(3))
	call updatetcontrol(trim(tcntrl_file_name),'$s3i',0,-isigma(3))

c       put the value of fermi-energy to <tcontrol>
        if (.not.fixed_efermi) call updatetcontrol(trim(tcntrl_file_name),'$efermi',0,ef0)

c       ... and update output for rguess >>> 
        rguess = rguess0
     
        print '(/,1x,a)',
     &	      '*****  DONE WITH A REAL PART OF THE SELF-ENERGY  *****'

        deallocate(fullsetmo,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', '  [SUBROUTINE tune_sigma]: ',
     &            ' impossible to deallocate temp. array <fullsetmo>'
         print *, ' anyway, proceed further ...'
        end if

        deallocate(neqdmat,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE tune_sigma]: ',
     &           'impossible to deallocate <neqdmat>'
         print *, 'nevertheless, proceed further ...'
        end if

c        if (aims_input) then 
c         print '(/,/,a,/,a,/,a,/)', 
c     &         ' to proceed with transport calculation (and, optionally, LDOS)' ,	 
c     &         ' switch on flags "$transmission" and/or "$ldos" to "on" ',
c     &         ' in your <tcontrol> file '
c        end if
     
       end subroutine tune_sigma
       
       subroutine update_sigma(fctr)
c        updates self-energy matrix       
         double precision, intent(in) :: fctr

         integer n 
         double precision  tmprs, tmpis, hex
         double precision, parameter :: epsilon = 1.0d-10

	 do n = 1, nsaos
	  tmpis = dimag(sigmann(n))
	  tmprs = tmpis * fctr
          if (.not.sp_leads) then
c          usual case: no exchange field in self-energy	  
	   sigmann(n) = cmplx(tmprs,tmpis,8)
	  else
c           case of spin-polarized electrodes	  
	    if (dabs(tmpis)<epsilon) then
c              actually, here we put zeros into arrays	    
	       sigmann(n) = cmplx(tmprs,tmpis,8)
	       bsigmann(n) = cmplx(tmprs,tmpis,8)
	    else 
c              here, we update alpha- & beta-pieces accordingly 	    
	       hex = 0.5d0*(bsigmann(n)-sigmann(n))
	       sigmann(n)  = cmplx(tmprs-hex,tmpis,8)
	       bsigmann(n) = cmplx(tmprs+hex,tmpis,8)
	    end if  
	  end if 
	 end do

       end subroutine update_sigma

      end module resigma
