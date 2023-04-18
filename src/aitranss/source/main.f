c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           may   2008
c     last revision:  march 2012
c###########################################################

      program main_module
c     **********************************
c     main program:
c     call step-by-step required modules
c     **********************************
      
       use globalvars
       use info

       use read_externals
       use make_overlap
       use density_matrix0

       use selfenergy
       use hamiltonian
       use int_out_core_states
       use exthamiltonian

       use neq_density_matrix
       use resigma
       use transmission 
       use pop_analysis

       use hubbardu
       use spectral_function

       use molsubsystem
       use ldos_rpoint
   
       implicit none

       character ddate*8, ttime*10, currdate*10, currtime*12
       
       real(4)  totstime, totftime, totsecs
       integer  totmnts, tothours

       logical       outinfo
       character*128 syscall

       call date_and_time(ddate,ttime)
       currdate = ddate(1:4)//'-'//ddate(5:6)//'-'//ddate(7:8)
       currtime = ttime(1:2)//':'//ttime(3:4)//':'//ttime(5:10)
       print '(/,1x,a,a,a)', currdate, ' ', currtime         

       call cpu_time(totstime)    
       call output_version_stamp
       call info_message

       call readtcontrol(trim(tcntrl_file_name))

       if (do_dmat_cycle) then
        print '(/,a,i4,a)', ' ************ DENSITY CYCLE ITERATION : ',
     &                 iter,' ************'  
       end if    

       print '(/,a)', ' READ EXTERNAL FILES >>>'

       if (aims_input) then
        call grepatoms(trim(crd_file_name),trim(crd_tmpfile_name))
	call readcoord(trim(crd_tmpfile_name))
c       in case of aims-input, remove temporary coord-file >>>
        if (.not.testing) then
         syscall = 'rm '//trim(crd_tmpfile_name) ; call system(trim(syscall))
	end if ! testing 
       else
        call readcoord(trim(crd_file_name))
       end if

       if (localexf) then
        call readhsource(trim(hsrc_file_name))
       end if
       
       if (aims_input) then
        call read_aims_basis(trim(basis_file_name),trim(tmpbasis_file_name))
        call readbasis(trim(tmpbasis_file_name))
       else
        call readbasis(trim(basis_file_name))
       end if

       call import_mos

c      reading overlap matrix is required, e.g. in case of fhi-aims
       if (read_omat) then
        call readomat(inomat_file_name)        
       end if

c      import and initialize a self-energy matrix >>>
       if (import_self_energy) then
         call read_self_energy(trim(self_energy_file),outinfo)
         if (outinfo) then 
c         fill a self-energy matrix
	  call fill_sigma_matrix
	 else
c         if reading was unsuccessful, switch off a flag 'import_self_energy': 
c         self-energy will be automatically initialized later 
	  import_self_energy = .false.
	 end if  ! info
       end if  ! import_self_energy   

       print '(/,a)', ' <<< DONE WITH EXTERNAL FILES'
   
c      in case of full-electron treatment (fhi-aims),
c      check electron number >>> 
       if (aims_input) call check_electron_number

c      initialize a self-energy matrix if it has not been read >>>
       if (.not.import_self_energy) then
         call buildsigma(dz_tolerance)
       end if

c      if we haven't read overlap matrix, we have to 
c      evaluate/initialize it, which is e.g. a case of TURBOMOLE 
       if (.not.read_omat) then
        call init_overlap
       end if

c      next call is for the equilibrium density matrix ;
c      allows to check consistency of molecular orbitals ;
c      keywords $dens_matrix (and subsequent ones, see externals.f)
c      should be present in <tcontrol> file
c      ! required in case of lda+u and no_leads flags
       call init_densmat0

       call build_smtrxs
       
       call h0
c       if (.not.symmetrize) then
c        call h0
c       else
c        call h0_symmetrize
c       end if    

c      if required, integrate out core states
       if (ecp) then
c       get value of the bare fermi energy 
        efermi = get_fermi_level(nelectr,1)
c       call updatetcontrol(trim(tcntrl_file_name),'#$efermi0',0,efermi)
c       update hamiltonian, molecular orbitals, overlap integrals, etc.
        call wrapper_valence_states
       end if    

c      lda+u corrections to the xc-potential
c      and the modified hamiltonian are evaluated here
       if (ldau) call make_lda_u

       if (do_pop.or.do_lmpop) then 
        call lwcycle
       else 
c       do the rest >>>

c       in case of 'no_leads' flag, switch to the 
c       evaluation of the equilibrium density-matrix
        if (no_leads) then

c        in case of lda+u, update and save the equilibrium dens.matrix 
         call update_eq_dmat
c        if required, output a spectral function 
         if (spectr_func) then 
	  call hfull
	  call hspectrum
	  call corr_ldos_cycle
	 end if
 	       
        else if (landauer) then

c        here comes only Landauer transmission function
         call hfull
	 call hspectrum
         efermi = get_fermi_level(nelectr,1)

c        performs calculation of T(E) and (optional) of LDOS >>>
         if (do_landauer) call tcycle 
         if (do_ldos)     call ldoscycle
	 
        else
c        here come a usual calculation, with leads, 
c        options for the self-energy, transmission, ldos, dens.matrix, etc.
	
c        optional job is done here if $tune_rigma is "on"
c        >>> this flag can not be active if a self-consistent
c            cycle for the density matrix is required !
         if (tune_rsigma) then
           call tune_sigma(rguess)
         else
c          otherwise, use parameters for self-energy from <tcontrol>
           call hfull
           call hspectrum
	 end if  

c        main job is done here ...      
c        >>>  we will proceed with further calculation, even if 
c             just in the step above the self-energy parameters 
c             have been tuned  ($tune_rsigma flag)  
 
         if (do_dmat_cycle) then
          call make_dmat
c         in case of lda+u, print out a spectral function
c         of the correlated subspace to external file
          if (ldau.and.spectr_func) call corr_ldos_cycle
         
         else if (mol_subsystem) then
c         performs analysis of the "molecular" subsystem 
c         embedded into extended system: computes 
c         spectrum, occupation numbers, local dos ...
          call make_dmat
          call execute_molecular_subsystem

         else if (do_ldos3d.or.do_ldos3d_ewindow) then
c         compute a space resolved spectral function
          call make_dmat 
	  call make_rldos 
	  
         else 
c         here do_trans_calc==.true. or do_ldos==.true. 
c         so we are doing energy cycle          
          call ecycle

         end if ! do_dmat_cycle

        end if ! no_leads
	
       end if ! do_pop
       
       call deallocall

       call cpu_time(totftime)
       totmnts  = int((totftime-totstime)/60.0d0,4)
       totsecs  = totftime - totstime - totmnts*60.0d0
       tothours = int(totmnts/60.0d0,4)
       
       if (tothours.eq.0) then
        print '(/,1x,a,i2,a,f5.2,a,/)',
     &        '** TOTAL CPU TIME: ', 
     &        totmnts, ' min ', totsecs, ' sec **'
       else    
        totmnts =  totmnts - tothours*60.0
        print '(/,1x,a,i3,a,i2,a,f5.2,a,/)',
     &      '** TOTAL CPU TIME:', 
     &      tothours, ' h ', totmnts, ' min ', totsecs, ' sec **'
       end if

       print '(/,a,/)', ' ** aitranss : all done **'

       call date_and_time(ddate,ttime)
       currdate = ddate(1:4)//'-'//ddate(5:6)//'-'//ddate(7:8)
       currtime = ttime(1:2)//':'//ttime(3:4)//':'//ttime(5:10)
       print '(1x,a,a,a)', currdate, ' ', currtime         

       if (do_dmat_cycle) then
        print '(/,a,i4,a,/)',
     &        ' ******* END OF DENSITY CYCLE ITERATION : ',iter,
     &        ' *******'
       end if  

      end program main_module 

