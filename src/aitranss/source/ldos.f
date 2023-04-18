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

      module ldos
c     ********************************      
c     computes local density of states 
c     ********************************      
       
       use globalvars
       use tools
       implicit none

c      local array to store the green's function for given energy
       complex(8), private, allocatable :: gf(:)
       
c      local density of states projected on different atoms       
       complex(8), private, allocatable :: 
c                        -- alpha/default channel --
     &                   locdos(:), lmdos(:,:), 
c                        -- beta channel --       
     &                   b_locdos(:), b_lmdos(:,:)
       
      contains

      subroutine split_ldos(gfe,ispin,ldos_out)
c      *********************************************************
c      -- evaluates local dos 
c      -- splits contributions according to different atom sorts
c      *********************************************************
       complex(8), intent (in)  :: gfe(:)
       integer,    intent (in)  :: ispin
       complex(8), intent (out) :: ldos_out(:)
    
       integer     isort, iatype, iat, norb, iorb, n, p 
       complex(8)  tmpgr
       logical     reservoir
       
       forall(isort=1:num_atom_sorts) ldos_out(isort) = czero
       n = 0
       do iat = 1, num_atoms

         isort  = atom(iat)%asort
	 iatype = atom(iat)%atype
	 norb = n_basis_func(iatype)

         reservoir = ( (atom(iat)%llead).or.(atom(iat)%rlead) )
         if (.not.reservoir) then
          do iorb = 1, norb
  	    n = n + 1 
 	    do p = 1, nsaos
	     tmpgr  = heigvec(n,p,ispin) * gfe(p) * invheigvec(p,n,ispin)
             ldos_out(isort) = ldos_out(isort) + 	
     &	      ione*(tmpgr - dconjg(tmpgr))/(2.0d0*pi)
	    end do
	  end do
	 else
	  n = n + norb
	 end if 

       end do ! iat
           
      end subroutine split_ldos

      subroutine split_lmdos(gfe,ispin,ldos_out,lmdos_out)
c      ***********************************************
c      -- evaluates local dos & lm-projected local dos 
c      ***********************************************
       complex(8), intent (in)  :: gfe(:)
       integer,    intent (in)  :: ispin
       complex(8), intent (out) :: ldos_out(:), lmdos_out(:,:)
    
       integer     isort, iatype, iat, norb, iorb, n, p, lm, lmmax
       complex(8)  tmpgr
       type(cgto)  tmporb
       logical     reservoir
       
       lmmax = (lmax0+1)*(lmax0+1)       
       forall(isort=1:num_atom_sorts,lm=1:lmmax) 
     &                                lmdos_out(isort,lm) = czero
       forall(isort=1:num_atom_sorts) ldos_out(isort) = czero

       n = 0
       do iat = 1, num_atoms
        isort  = atom(iat)%asort
	iatype = atom(iat)%atype
	norb = n_basis_func(iatype)

        reservoir = ( (atom(iat)%llead).or.(atom(iat)%rlead) )
        if (.not.reservoir) then

         do iorb = 1, norb
 	  n = n + 1 
c         ! take lm index of the orbital !
          tmporb = aos(iatype,iorb) 
	  lm = tmporb%lm          
	  do p = 1, nsaos
	   tmpgr  = heigvec(n,p,ispin) * gfe(p) * invheigvec(p,n,ispin)
c          update ldos contribution
           ldos_out(isort) = ldos_out(isort) + 
     &	      ione*(tmpgr - dconjg(tmpgr))/(2.0d0*pi)
c          update lm-projected contributions
           lmdos_out(isort,lm) = lmdos_out(isort,lm) + 
     &	      ione*(tmpgr - dconjg(tmpgr))/(2.0d0*pi)
	  end do ! p
	 end do ! iorb
	else 
	 n = n + norb
	end if 
	
       end do ! iat
           
      end subroutine split_lmdos

      subroutine ldoscycle
c      *****************************************
c      performs energy cycle : evaluates ldos(e)
c      *****************************************

       double precision tmpldos, b_tmpldos, energy,
     &                  ener_save, estep_save, eend_save,
     &                  tmplmdos((lmax0+1)*(lmax0+1)), 
     &                  b_tmplmdos((lmax0+1)*(lmax0+1)) 
       integer  ispin, ierr, n, isort, outfile, 
     &          lmmax, tmp_lmax0, lm
       character(80), allocatable :: outfilename(:)
       character(128) lmfmt
           
       real(4)  stime, ftime, secs
       integer  mnts

c      allocation of global array with GF:
       allocate(gf(nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE ldoscycle]: <gf> allocation failure'
       end if

c      allocation of arrays with local dos   
       allocate(locdos(num_atom_sorts),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE ldoscycle]: <locdos> allocation failure'
       end if
c      additional array for beta channel
       if (nspin.eq.2) then
         allocate(b_locdos(num_atom_sorts),stat=ierr)
         if (ierr.ne.0) then
          print *
          stop
     &    '[SUBROUTINE ldoscycle]: <b_locdos> allocation failure'
         end if
       end if	 
       
c      next, we define some variables controlling an output 
c      of the lm-projected local density of states
       if (do_lmdos) then
        if (lmax0 > 2) then 
	  tmp_lmax0 = 2 ! only s,p & d contributions are counted
	 else 
	  tmp_lmax0 = lmax0
	 end if
	 
         select case (tmp_lmax0)
          case (0)
c          marginal case of s-orbitals only: nothing to decompose
           do_lmdos = .false.      
	  case (1)
c          s- & p-functions 
	   if (nspin.eq.1) then
	    lmfmt = '(1x,f14.10,f16.8,2x,2(E20.10),'//
     &	            '2x,E16.6,2x,3(E16.6))'
           else
	    lmfmt = '(1x,f14.10,f16.8,2x,3(E20.10),'//
     &	            '2x,2(E16.6),2x,6(E16.6))'
	   end if
	  case (2)
c          s-,p- & d-functions 
	   if (nspin.eq.1) then
	    lmfmt = '(1x,f14.10,f16.8,2x,2(E20.10)'//
     &	            ',2x,E16.6,2x,3(E16.6),2x,5(E16.6))'
           else
	    lmfmt = '(1x,f14.10,f16.8,2x,3(E20.10),'//
     &	            '2x,2(E16.6),2x,6(E16.6),2x,10(E16.6))'
	   end if
          case default ; 
           print *
           stop '[SUBROUTINE ldoscycle]: tmp_lmax0 is out of range!'
         end select
	 
       end if ! do_lmdos

c      in case of lm-dos, dynamical allocation of additional arrays is required
       lmmax = (lmax0+1)*(lmax0+1)
       if (do_lmdos) then
         allocate(lmdos(num_atom_sorts,lmmax),stat=ierr)
         if (ierr.ne.0) then
          print *
          stop
     &     '[SUBROUTINE ldoscycle]: <lmdos> allocation failure'
         end if
c        additional array for the beta channel
         if (nspin.eq.2) then
          allocate(b_lmdos(num_atom_sorts,lmmax),stat=ierr)
          if (ierr.ne.0) then
           print *
           stop
     &    '[SUBROUTINE ldoscycle]: <b_lmdos> allocation failure'
          end if
         end if
       end if	 

c      update lmmax: needed for output information
       if (lmax0 > 2) then
        lmmax = 9  ! just s, p, and d-functions
       end if

c      output file names
       allocate(outfilename(num_atom_sorts),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE ldoscycle]: <outfilename> allocation failure'
       end if
       do isort = 1, num_atom_sorts
        outfilename(isort) = 'ldos.'//trim(ref_sort(isort))//'.dat'    
       end do 
       
       call cpu_time(stime)
 
       if (.not.do_lmdos) then
        print '(/,a,/)', ' CALCULATING LOCAL DENSITY OF STATES >>>'      
       else 
        print '(/,a,/)', ' CALCULATING LM-PROJECTED LOCAL DENSITY OF STATES >>>'      
       end if

c      write headers to output files
       do isort = 1, num_atom_sorts
       
        outfile = 20 + isort
        open(outfile,file=outfilename(isort),
     &	             status='unknown',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a,a,/,a)',
     &         ' can not open file <',
     &         trim(outfilename(isort)),'> for writing',
     &         ' please, check your output directory access rights'  
         stop  ' transport module is terminated now'
        else
         write(outfile,fmt='(a,a)') 
     &	              '#local dos : atoms = ', trim(atom_sort(isort))
         if (nspin ==  1) then
          write(outfile,fmt='(a)') '#non-spin-polarized calculation'
          write(outfile,'(a,f10.6,a,f17.12,a)')
     &                  '#bias =',bias,' V    efermi =', efermi, ' H'
           write(outfile,fmt='(a,3x,a,7x,a,4x,a,4x,a)')
     &      '#', 'E [Hartree]', 'E-EF [eV]',
     &      'LDOS per spin [1/eV]', 'LDOS [1/eV]'
         else
          write(outfile,fmt='(a)') '#spin-polarized calculation'
          write(outfile,'(a,f10.6,a,f17.12,a)')
     &                  '#bias =',bias,' V    efermi =', efermi, ' H'
           write(outfile,fmt='(a,3x,a,7x,a,6x,a,3x,a,6x,a)')
     &     '#','E [Hartree]','E-EF [eV]',
     &     'LDOS.alpha [1/eV]','LDOS.beta [1/eV]','LDOS [1/eV]'
         end if
        end if           
      
       end do
       
c      take care about energy mesh, if input is in eV
       if (evunits) then
        ener_save = ener ; estep_save = estep ; eend_save = eend
        ener  = efermi + ener/hartree
        estep = estep/hartree
        eend =  efermi + eend/hartree
       end if

c      cycle over energy points      
       n = 0  ; energy = ener    
       do while (energy <= eend)  
       
        if (evunits) then
	 print '(2x,a,f14.10,a,$)', 'E - EF = ',(energy-efermi)*hartree,' eV    '
        else
         print '(2x,a,f14.10,a,$)', 'E = ',energy,' H    '
        end if

c       cycle over spin channels
        do ispin = 1, nspin
c        evaluating local dos
         call greensf(energy,ispin,gf)
         if (ispin.eq.1) then
 	   if (.not.do_lmdos) then
	     call split_ldos(gf,1,locdos)
	   else
	     call split_lmdos(gf,1,locdos,lmdos)
	   end if     
         else
	  if (.not.do_lmdos) then
            call split_ldos(gf,2,b_locdos)
	  else
	    call split_lmdos(gf,2,b_locdos,b_lmdos)
	  end if
	 end if  
	end do ! ispin

c       write data in external output-files

        if (.not.do_lmdos) then
c        just ldos contributions are written to output files
         do isort = 1, num_atom_sorts
	  outfile = 20 + isort
	  if (nspin == 1) then
 	   tmpldos = dble(locdos(isort))/hartree
 	   write(outfile,fmt='(1x,f14.10,f16.8,2x,2(E20.10))')
     &           energy, (energy-efermi)*hartree, 
     &           tmpldos, tmpldos*2.0d0
          else
	   tmpldos   = dble(locdos(isort))/hartree
	   b_tmpldos = dble(b_locdos(isort))/hartree
	   write(outfile,fmt='(1x,f14.10,f16.8,2x,3(E20.10))')
     &           energy, (energy-efermi)*hartree, 
     &           tmpldos, b_tmpldos, tmpldos+b_tmpldos      	   
       	  end if
	 end do ! isort

        else
c       here we deal with lm-projected dos
         do isort = 1, num_atom_sorts
          outfile = 20 + isort
	  if (nspin == 1) then
c          one spin channel
  	   tmpldos = dble(locdos(isort))/hartree
 	   do lm = 1, lmmax
	    tmplmdos(lm) = dble(lmdos(isort,lm))/hartree
	   end do
 	   write(outfile,trim(lmfmt)) 
     &	         energy,(energy-efermi)*hartree,tmpldos,tmpldos*2.0d0, 
     &           (tmplmdos(lm),lm=1,lmmax)
          else
c          two spin channels	  
	   tmpldos   = dble(locdos(isort))/hartree
	   b_tmpldos = dble(b_locdos(isort))/hartree
 	   do lm = 1, lmmax
	    tmplmdos(lm) = dble(lmdos(isort,lm))/hartree
	    b_tmplmdos(lm) = dble(b_lmdos(isort,lm))/hartree
	   end do
	   write(outfile,trim(lmfmt))
     &           energy,(energy-efermi)*hartree,
     &           tmpldos, b_tmpldos, tmpldos+b_tmpldos,
     &           (tmplmdos(lm),b_tmplmdos(lm), lm=1,lmmax) 
       	  end if
	 end do ! isort
	
	end if

        print *, ' done'

        n = n + 1 
	energy = ener + n * estep
       end do ! e-points
 
       do isort = 1, num_atom_sorts
	 outfile = 20 + isort
         write(outfile,fmt='(a)') '#end'
         close(outfile)
       end do

       print '(/,a)', ' <<< DONE WITH LOCAL DOS'
       call cpu_time(ftime)
       mnts = int((ftime-stime)/60.0d0,4)
       secs = ftime - stime - mnts*60.0d0
       print '(1x,a,i3,a,f5.2,a)',
     &       'time spent:', mnts, ' min ', secs, ' sec'

c      restore parameters for energy mesh
       if (evunits) then
        ener = ener_save ; estep = estep_save ; eend = eend_save
       end if

       if (.not.do_lmdos) then
         deallocate(gf,locdos,outfilename,stat=ierr)
       else 
         deallocate(gf,locdos,lmdos,outfilename,stat=ierr)
       end if	 
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE ldoscycle]: ',
     &           'impossible to deallocate temp. arrays'
        print *, 'nevertheless, proceed further ...'
       end if

c      take care about the beta channel if exists ...
       if (nspin.eq.2) then
        if (.not.do_lmdos) then 
          deallocate(b_locdos,stat=ierr)
	else 
          deallocate(b_locdos,b_lmdos,stat=ierr)
	end if 
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE ldoscycle]: ',
     &           'impossible to deallocate b-channel temp. array'
         print *, 'nevertheless, proceed further ...'
        end if
       end if

      end subroutine ldoscycle

      end module ldos 

