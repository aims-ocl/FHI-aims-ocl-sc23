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

      module pop_analysis
c      *************************************      
c      performs atom- & eigenstate-projected  
c      loewdin population analysis 
c      *************************************
       use globalvars
       implicit none
           
       contains
       
      subroutine orbloewdin(p,ispin,lwout)
c       *********************************
c       splits to atoms a loewding charge 
c       coming from the eigenstate p  
c       ********************************* 
        implicit none
        integer,          intent(in)  :: p, ispin
        double precision, intent(out) :: lwout(:)
       
        integer isort, n, iat, iatype, norb, iorb
        double precision c_n
    
	forall(isort=1:num_atom_sorts) lwout(isort) = 0.0d0
        n = 0
        do iat = 1, num_atoms
         isort  = atom(iat)%asort
         iatype = atom(iat)%atype
         norb = n_basis_func(iatype)
         do iorb = 1, norb
          n = n + 1
          c_n = mo_orth(n,p,ispin)
          lwout(isort) = lwout(isort) + c_n*c_n
         end do
        end do
              
      end subroutine orbloewdin

      subroutine lmloewdin(p,ispin,lwout,lmlwout)
c       ***************************************************************
c       splits to atoms a loewding charge coming from the eigenstate p; 
c       decomposes charge to lm-contributions         
c       ***************************************************************
        implicit none
        integer,          intent(in)  :: p, ispin
        double precision, intent(out) :: lwout(:), lmlwout(:,:)
       
        integer          isort, n, iat, iatype, norb, iorb, lm, lmmax
        double precision c_n
        type(cgto)       tmporb
    
        lmmax = (lmax0+1)*(lmax0+1)
        forall(isort=1:num_atom_sorts,lm=1:lmmax)
     &        lmlwout(isort,lm) = 0.0d0
    	forall(isort=1:num_atom_sorts) lwout(isort) = 0.0d0
        n = 0
        do iat = 1, num_atoms
         isort  = atom(iat)%asort
         iatype = atom(iat)%atype
         norb = n_basis_func(iatype)
         do iorb = 1, norb
c         ! take lm index of the orbital !
          tmporb = aos(iatype,iorb)
          lm = tmporb%lm
          n = n + 1
          c_n = mo_orth(n,p,ispin)
          lwout(isort) = lwout(isort) + c_n*c_n
          lmlwout(isort,lm) = lmlwout(isort,lm) + c_n*c_n
         end do
        end do
              
      end subroutine lmloewdin

       
      subroutine lwcycle
c      *********************************************      
c      header routine: 
c      -- takes only orbitals in given energy window
c      -- saves output to external files
c      *********************************************
       implicit none
       character(80),    allocatable :: outfilename(:)
       integer,          allocatable :: outfile(:)
       double precision, allocatable :: lwpop(:), lmlwpop(:,:)

       integer           p, ispin, lm, lmmax, ierr, isort 
       double precision  ep
       logical           flag
       character(128)    lmfmt

c      output file names
       allocate(outfilename(num_atom_sorts),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE lwcycle]: <outfilename> allocation failure'
       end if
       do isort = 1, num_atom_sorts
        outfilename(isort) = 'loewdin.'//trim(ref_sort(isort))//'.dat'
       end do

c      output files
       allocate(outfile(num_atom_sorts),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE lwcycle]: <outfile> allocation failure'
       end if

c      loewdin populations
       allocate(lwpop(num_atom_sorts),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop
     &    '[SUBROUTINE lwcycle]: <lwpop> allocation failure'
       end if
       
       lmmax = (lmax0+1)*(lmax0+1)
       if (do_lmpop) then
        allocate(lmlwpop(num_atom_sorts,lmmax),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop
     &    '[SUBROUTINE lwcycle]: <lwpop> allocation failure'
        end if
       end if

       if (do_lmpop) then
        print '(/,a,/)', ' PERFORMING lm-decomposed POPULATION ANALYSIS >>>'
       else
        print '(/,a,/)', ' PERFORMING POPULATION ANALYSIS >>>'
       end if	
     
       lmfmt='(1x,f14.10,f16.8,2x,E20.10,2x,E16.6,2x,3(E16.6),2x,5(E16.6))'
c      we will print out only decomposition up to d-functions
       lmmax = 9
       do isort = 1, num_atom_sorts

        outfile(isort) = 20 + isort
        open(outfile(isort),file=outfilename(isort),
     &               status='unknown',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a,a,/,a)',
     &         ' can not open file <',
     &         trim(outfilename(isort)),'> for writing',
     &         ' please, check your output directory access rights'
         stop  ' transport module is terminated now'
        else
         write(outfile(isort),fmt='(a,a)')
     &        '#loewdin population : atoms = ', trim(atom_sort(isort))
         if (efermi_guess) then
          write(outfile(isort),fmt='(a,3x,a,7x,a,9x,a)')
     &      '#', 'E [Hartree]', 'E-EF [eV]', 'population'
         else
           write(outfile(isort),fmt='(a,3x,a,7x,a,12x,a)')
     &      '#', 'E [Hartree]', 'E [eV]', 'population'
         end if
        end if 
       end do

       do ispin = 1, nspin
        if (nspin.eq.2) then 
         do isort = 1, num_atom_sorts
	   if (ispin.eq.1) then
            write(outfile(isort),fmt='(a)'), '#alpha-channel'
	   else
            write(outfile(isort),fmt='(a)'), '#beta-channel'
	   end if
         end do  
        end if

c       cycle over molecular orbitals
        do p = 1, nsaos
	  ep   = mo_en(p,ispin)
	  flag = (ep>=ener).and.(ep<=eend)
          if (flag) then
           print '(2x,a,f14.10,a,$)', 'E = ',ep,' H    '
           if (do_pop) then 
            call orbloewdin(p,ispin,lwpop)
	    do isort = 1, num_atom_sorts
             if (ispin.eq.2) lwpop(isort) = -lwpop(isort)
             write(outfile(isort),fmt='(1x,f14.10,f16.8,2x,E20.10)')
     &             ep, (ep-efermi)*hartree, lwpop(isort)
	    end do ! isort
           else
c           lm-decomposed population	   
            call lmloewdin(p,ispin,lwpop,lmlwpop)
	    do isort = 1, num_atom_sorts
             if (ispin.eq.2) then
	      lwpop(isort) = -lwpop(isort)
	      do lm = 1, lmmax
	       lmlwpop(isort,lm) = -lmlwpop(isort,lm)
	      end do
	     end if  ! ispin==2
             write(outfile(isort),fmt=trim(lmfmt))  
     &	                   ep, (ep-efermi)*hartree, 
     &	                   lwpop(isort), (lmlwpop(isort,lm),lm=1,lmmax)
            end do ! isort
           end if ! do_pop vs. do_lmpop
	   print *, 'done'
	  end if ! flag
	end do ! p-cycle
       end do ! ispin

       do isort = 1, num_atom_sorts
        write(outfile(isort),fmt='(a)'), '#end'
        close(outfile(isort))
       end do 
       
       print '(/,a,/)', ' <<< DONE WITH POPULATION ANALYSIS'

       print '(a)', ' to proceed with transport calculation, please'
       print '(a)', ' use "$population off" in your <tcontrol> file'

       deallocate(outfilename,outfile,lwpop,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE lwcycle]: ',
     &           'impossible to deallocate temp. arrays'
        print *, 'nevertheless, proceed further ...'
       end if

       if (do_lmpop) then
        deallocate(lmlwpop,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE lwcycle]: ',
     &           'impossible to deallocate <lmlwpop>'
         print *, 'nevertheless, proceed further ...'
        end if
       end if

      end subroutine lwcycle
      
      end module pop_analysis
      
