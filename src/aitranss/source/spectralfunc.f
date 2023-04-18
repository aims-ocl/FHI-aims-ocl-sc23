c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           may 2009
c     last revision:  jan 2012
c###########################################################

      module spectral_function
c     *****************************************      
c     spectral function for correlated subspace
c     in case of lda+u type of calculation 
c     *****************************************      
       
       use globalvars
       use tools
       use hubbardu
       implicit none

c      local array to store the green's function for given energy
       complex(8), private, allocatable :: gf(:)

c      spectral function of correlated d-ion site(s) is 
c      projecteed on d-orbitals only 
       integer, private, parameter :: dsize = 5
       
c      local dos at correlated site(s) projected on d-orbitals       
       complex(8), private :: 
c                        -- alpha/default channel --
     &                   locdos(dsize),
c                        -- beta channel --       
     &                   b_locdos(dsize)

c      temporary arrays
       complex(8), private, allocatable :: osb(:,:,:), bso(:,:,:)
       
      contains

      subroutine corr_aw(gfe,ispin,corr_ldos_out)
c      ****************************************************
c      -- evaluates spectral function at correlated site(s)
c      -- splits it into 5 d-orbital symmetries
c      ****************************************************
       complex(8), intent (in)  :: gfe(:)
       integer,    intent (in)  :: ispin
       complex(8), intent (out) :: corr_ldos_out(:)

       integer    iasort, iatype, lm, ilm, iat, norb, iorb,
     &            dcount, is, n, ia, p
       logical    flag
       complex(8) tmpgr

       forall(lm=1:dsize) corr_ldos_out(lm) = czero

       n = 0 ; ia = 0 
       do iat = 1, num_atoms
        iatype = atom(iat)%atype
        iasort = atom(iat)%asort

        flag = .false.
        do is = 1, num_ratoms
         flag = ( flag.or.(iasort == ratom_types(is)) )
        end do

        if (flag) then
c        given atom <iat> belongs to a "reservoir"
         n = n + n_basis_func(iatype)
        else
c        given atom <iat> is correlated site
c        here we pick up only the first d-type CGTOs
c        with large exponents but not the diffuse d-orbitals
         norb = n_basis_func(iatype)
         dcount = 0
         do iorb = 1, norb
          n  = n  + 1
          lm = aos(iatype,iorb)%lm
          if ((lm.ge.5).and.(lm.le.9).and.(dcount.lt.5)) then
c         if ((lm.ge.5).and.(lm.le.9)) then
           dcount = dcount + 1
           ia = ia + 1
ccc        compute lm-contribution to A(w)
           ilm = lm - 4
           do p = 1, nsaos
            tmpgr  = osb(ia,p,ispin) * gfe(p) * bso(p,ia,ispin)
            corr_ldos_out(ilm) = corr_ldos_out(ilm) +
     &                           ione*(tmpgr - dconjg(tmpgr))/(2.0d0*pi)
           end do
          end if  ! lm-choice
         end do ! orb

        end if ! flag

       end do ! iat

c      testing: print out 'n'
c      print '(1x,a,i4,a,i4,$)', '[corr_aw]: n =', n, '  ia =', ia

      end subroutine corr_aw

      subroutine corr_ldos_cycle
c      *****************************************
c      performs energy cycle : evaluates ldos(e)
c      *****************************************

       double precision tmpldos(dsize), b_tmpldos(dsize), energy,
     &                  ener_save, estep_save, eend_save  
       integer          ispin, ierr, n, outfile, lm
       character(64)    outfilename, lmfmt

       real(4)  stime, ftime, secs
       integer  mnts

c      allocation of array with GF:
       allocate(gf(nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE corr_ldos_cycle]: <gf> allocation failure'
       end if

       allocate(osb(ion_size,nsaos,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE corr_ldos_cycle]: <osb> allocation failure'
       end if

       allocate(bso(nsaos,ion_size,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE corr_ldos_cycle]: <osb> allocation failure'
       end if

       call cpu_time(stime)
       print '(/,a,/)', ' <HUBBARD-U>: OUTPUT A SPECTRAL FUNCTION OF THE d-ION(S) >>>'      

       print '(2x,a,$)', 'projector operators ... '      
c      evaluate projector operators >>>
       do ispin = 1, nspin

c        1. compute matrix osb = omat12inv * smat12 * heigvec
         call xcmplx_ab('n','n',ion_size,nsaos,nsaos,cmplx_osmat,heigvec(:,:,ispin),osb(:,:,ispin))

c        2. compute matrix bso = invheigvec * smat12 * omat12inv
         call xcmplx_ab('n','t',nsaos,nsaos,ion_size,invheigvec(:,:,ispin),cmplx_osmat,bso(:,:,ispin))

       end do ! ispin
       print '(a,/)', ' are initialized '      
   
c      output file name
       outfilename = 'corr.ldos.dat'    

       outfile = 20
       open(outfile,file=outfilename, status='unknown',iostat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a,a,/,a)',
     &         ' can not open file <',
     &         trim(outfilename),'> for writing',
     &         ' please, check your output directory access rights'  
         stop  ' transport module is terminated now'
       else
         write(outfile,fmt='(a)')  '#d-like local dos : correlated site(s)'
         if (nspin ==  1) then
          write(outfile,fmt='(a)') '#non-spin-polarized calculation'
         else
          write(outfile,fmt='(a)') '#spin-polarized calculation'
         end if ! ispin
         write(outfile,'(a,f16.10,a)') '#u-j ', u_j, ' eV'
         write(outfile,'(a,f10.6,a,f17.12,a)')
     &                  '#bias =',bias,' V    efermi =', efermi, ' H'
         write(outfile,fmt='(a,3x,a,7x,a)') '#', 'E [Hartree]', 'E-EF [eV]'
       end if           

c      choose a format of output line 
       if (nspin.eq.1) then
         lmfmt = '(1x,f14.10,f16.8,2x,5(E16.6))'
       else
         lmfmt = '(1x,f14.10,f16.8,2x,10(E16.6))'
       end if

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
c        evaluating local dos
c        print '(2x,a,i1,$)', ' ispin = ',ispin
         call greensf(energy,ispin,gf)
c	 print '(a,$)', ' ... [greensf] done '

c        compute spectral function
         if (ispin.eq.1) then        
	   call corr_aw(gf,1,locdos)
	 else
	   call corr_aw(gf,2,b_locdos)
	 end if 
c        print '(a,$)', ' ... [corr_aw] done '

	end do ! ispin

c       write data to external output-files
        if (ispin.eq.1) then
 	  tmpldos = dble(locdos)/hartree
          write(outfile,trim(lmfmt)) 
     &	         energy,(energy-efermi)*hartree,(tmpldos(lm),lm=1,dsize)
	else
 	  tmpldos   = dble(locdos)/hartree
 	  b_tmpldos = dble(b_locdos)/hartree
          write(outfile,trim(lmfmt)) 
     &	         energy,(energy-efermi)*hartree,(tmpldos(lm),b_tmpldos(lm),lm=1,dsize)
	end if        
        print *, ' done'

        n = n + 1 
	energy = ener + n * estep
       end do ! e-points

       write(outfile,'(a)') '#end'
       close(outfile)

       print '(/,a)', ' <HUBBARD-U>: <<< DONE WITH OUTPUT OF THE SPECTRAL FUNCTION'

       call cpu_time(ftime)
       mnts = int((ftime-stime)/60.0d0,4)
       secs = ftime - stime - mnts*60.0d0
       print '(1x,a,i3,a,f5.2,a)',
     &       'time spent:', mnts, ' min ', secs, ' sec'

c      restore parameters for energy mesh
       if (evunits) then
        ener = ener_save ; estep = estep_save ; eend = eend_save
       end if

       deallocate(gf,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE corr_ldos_cycle]: ',
     &           'impossible to deallocate <gf>'
        print *, 'nevertheless, proceed further ...'
       end if

       deallocate(osb,bso,stat=ierr)
c      print '(a,i4)', 'stat =', ierr
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE corr_ldos_cycle]: ',
     &         'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
       end if 

      end subroutine corr_ldos_cycle

      end module spectral_function 

