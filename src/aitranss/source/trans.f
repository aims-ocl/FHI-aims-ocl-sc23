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
c     last revision:  jul 2012
c###########################################################

      module transmission
c     ******************************************      
c     -- computes a transmission function T(E,V) 
c        performs spin & energy cycle
c     -- if required, performs LDOS calculation  
c     ******************************************    
       
       use globalvars
       use neq_density_matrix
       use tools
c      >>> inserted by RK and AB:
       use eigenvalue_solver
c      <<< done with insert by RK and AB
       use cond_chann_wf
       use ldos
       implicit none

c      diagonal representaion of the green's function for a given spin
       complex(8), private, allocatable :: gf(:)

c      represenation of gamma-matrices in the basis of eigenvectors:
       complex(8), private, allocatable :: bgb_left(:,:,:), bgb_right(:,:,:)

c     >>> inserted by RK and AB:
c     temporary arrays for the eigenchannels evaluation:
      complex(8), private, allocatable       :: g12b_left(:,:,:)
      double precision, private, allocatable :: channels(:,:)
c     <<< done with insert by RK and AB

      contains

      subroutine ecycle  
c     ***************************************************
c     higher level routine: 
c        -- searches for the fermi energy (efermi)   
c        -- executes "ldos" & "transmission" cycle calls 
c     ***************************************************     
       implicit none
       integer          ierr
       double precision tmp_nelectr

c      non-equilibrium density matrix       
       if (.not.allocated(neqdmat)) then
        allocate(neqdmat(nsaos,nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop
     &   '[SUBROUTINE ecycle]: <neqdmat> allocation failure'
        end if
       end if
       	
c      searching for the fermi energy
c      call searchef(efermi,neqdmat)       

c      searching for the fermi-energy; evaluate dens.matrix
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

c      prints out atomic loewdin charges/magn.moments
       if (qoutput) call qprint(neqdmat)

       deallocate(neqdmat,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE ecycle]: ',
     &           'impossible to deallocate <neqdmat>'
        print *, 'nevertheless, proceed further ...'
       end if

c      performs ldos calculation if required >>>
       if (do_ldos) call ldoscycle
c      performs transmission calculation if required >>>
       if (do_trans_calc) call tcycle

      end subroutine ecycle

      subroutine tcycle
c      ***************************************************
c      -- transforms gamma-matrices 
c         to hamiltonian's eigenvectors basis
c      -- performs energy cycle : evaluates bias dependent 
c                                 transmission t(e,v)
c      ***************************************************
       implicit none

       double precision :: te, ate, bte, energy, sgn, esmall = 1.0d-3
       double precision    ener_save, estep_save, eend_save, tmp_sum
       integer :: ispin, ierr, n, k, outfile = 77, chfile = 50
       logical :: output = .true.
    
       real(4)  stime, ftime, secs
       integer  mnts
       
c      allocation of global array with GF:
        allocate(gf(nsaos),stat=ierr)
        if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE tcycle]: <gf> allocation failure'
       end if

c      allocation of arrays with transformed gamma-matrices:
       allocate(bgb_left(nsaos,nsaos,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE tcycle]: <bgb_left> allocation failure'
       end if
       allocate(bgb_right(nsaos,nsaos,nspin),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE tcycle]: <bgb_right> allocation failure'
       end if

c      >>> inserted by RK and AB:
       if (do_cond_channels) then
         allocate(g12b_left(nsaos,nsaos,nspin),stat=ierr)
         if (ierr.ne.0) then
          print *
          stop '[SUBROUTINE tcycle]: <g12b_left> allocation failure'
         end if
         allocate(channels(nsaos,nspin),stat=ierr)
         if (ierr.ne.0) then
          print *
          stop '[SUBROUTINE tcycle]: <nsaos> allocation failure'
         end if
       end if
c      <<< done with insert by RK and AB

       if (do_scatt_wave_func.and.(.not.allocated(utt)) ) then
         allocate(utt(nsaos,nsaos,nspin),stat=ierr)
         if (ierr.ne.0) then
          print *
          stop
     &   '[SUBROUTINE tcycle]: <utt> allocation failure'
         end if
       end if

c      transform gamma-matrices to basis of hamiltonian's eigenvectors
       print '(/,a)', ' TRANSFORMING GAMMA MATRICES >>>'      

       call cpu_time(stime)
       print '(/,2x,a,$)', 'left electrode ...'     
       do ispin = 1, nspin
        call transform_gamma(gamma_left,bgb_left(:,:,ispin),ispin,'L')  
c       >>> modified by RK and AB:
        if (do_cond_channels) call transform12_gamma(gamma_left,g12b_left(:,:,ispin),ispin)
       end do
       print *, 'done'

       print '(2x,a,$)', 'right electrode ...'
       do ispin = 1, nspin
        call transform_gamma(gamma_right,bgb_right(:,:,ispin),ispin,'R')
       end do
       print *, 'done'
     
       call cpu_time(ftime)
       mnts = int((ftime-stime)/60.0d0,4)
       secs = ftime - stime - mnts*60.0d0
       print '(2x,a,i3,a,f5.2,a)',
     &       'time spent:', mnts, ' min ', secs, ' sec'
       print '(/,a)', ' <<< DONE WITH GAMMA MATRICES'
       
c      transmission evaluation starts here ...

       open(outfile,file=trim(output_file_name),status='unknown',iostat=ierr)
       if (ierr.ne.0) then
        output = .false.
        print '(/,a,a,a,/,a)', 
     &         ' can not open output file <',
     &         trim(output_file_name),'> for writing',
     &         ' output data can only be found on console'       
       else
	if (nspin ==  1) then
	 write(outfile,fmt='(a)') '#non-spin-polarized calculation'
         write(outfile,'(a,f10.6,a,f17.12,a)') 
     &	               '#bias =',bias,' V    efermi =', efermi, ' H'
         write(outfile,fmt='(a,3x,a,7x,a,7x,a)') 
     &                 '#', 'E [Hartree]', 'E-EF [eV]', 'T(E) per spin'
	else
	 write(outfile,fmt='(a)') '#spin-polarized calculation'
         write(outfile,'(a,f10.6,a,f17.12,a)') 
     &	               '#bias =',bias,' V    efermi =', efermi, ' H'
	 write(outfile,fmt='(a,3x,a,7x,a,9x,a,9x,a)') 
     &      '#','E [Hartree]','E-EF [eV]','T(E).alpha','T(E).beta'
	end if
       end if

       if (do_cond_channels) then
        if (nspin == 1) then
	 write(chfile+1,fmt='(a)') '#non-spin-polarized calculation'
         write(chfile+1,'(a,f10.6,a,f17.12,a)') '#bias =',bias,' V    efermi =', efermi, ' H'
         write(chfile+1,fmt='(a,3x,a,7x,a,7x,a)')  '#', 'E [Hartree]', 'E-EF [eV]', 'T(E) per spin'
        else
	 write(chfile+1,fmt='(a)') '#spin-polarized calculation: alpha channel'
	 write(chfile+2,fmt='(a)') '#spin-polarized calculation: beta channel'
         write(chfile+1,'(a,f10.6,a,f17.12,a)') '#bias =',bias,' V    efermi =', efermi, ' H'
         write(chfile+2,'(a,f10.6,a,f17.12,a)') '#bias =',bias,' V    efermi =', efermi, ' H'
         write(chfile+1,fmt='(a,3x,a,7x,a,7x,a)')  '#', 'E [Hartree]', 'E-EF [eV]', 'T(E).alpha'
         write(chfile+2,fmt='(a,3x,a,7x,a,7x,a)')  '#', 'E [Hartree]', 'E-EF [eV]', 'T(E).beta'
        end if
       end if ! do_cond_channels    
       
c      take care about energy mesh, if input is in eV
       if (evunits) then
        ener_save=ener ; estep_save=estep ; eend_save=eend
        ener  = efermi + ener/hartree
        estep = estep/hartree
        eend  = efermi + eend/hartree
       end if       

       if (conductance) then
               print '(/,a,/,/,a,/)', 
     &	 ' CONDUCTANCE CALCULATION:',
     &   ' G = T(EFermi) * [e^2/h] per spin channel'
        n = 0
c       instead, we evaluate transmission at EF
	ener = efermi
c       define dummy mesh >>>	
	energy = ener
	estep = 2.0d0*esmall
	eend  = ener + esmall
       else if (do_scatt_wave_func) then
        print '(/,a,/,a,/)', 
     &	 ' TRANSMISSION AND SCATTERING WAFE FUNCTIONS ',
     &   ' WILL BE COMPUTED AT FIXED ENERGY POINT ONLY '
        n = 0
c       define dummy mesh >>> 
        energy = ener
	estep = 2.0d0*esmall
	eend  = ener + esmall
       else
        print '(/,a,/)', 
     &	 ' TRANSMISSION CALCULATION: CYCLE OVER E-POINTS'
        n = 0
        energy = ener
       end if

c      cycle over energy points      
       
       if (eend>=ener) then ; sgn = 1.0d0
       else ;  sgn = -1.0d0
       end if
       	
       do while ( sgn*(eend-energy)>=0 )
       
        if (evunits) then
         print '(a,f14.10,a,$)',' E - EF = ',(energy-efermi)*hartree,' eV   '
        else 
         print '(a,f14.10,a,$)',' E = ',energy,' H   '
        end if 

c       cycle over spin channels
        do ispin = 1, nspin

c        evaluating transmission
         call cpu_time(stime)
         te = trans(energy,ispin)
         call cpu_time(ftime)
c        >>> modified by RK and AB :
         if (do_cond_channels) then
           call evaluate_eigenchannels(energy,ispin,channels(:,ispin),tmp_sum)
         end if 
         mnts = int((ftime-stime)/60.0d0,4)
         secs = ftime - stime - mnts*60.0d0

         if (nspin == 1) then
          print '(a,E20.10)', ' T(E) =', te
         else
          if (ispin==1) then
           print '(a,E20.10,$)', ' alpha.T(E) =', te
           ate = te
          else
           print '(a,E20.10)', '   beta.T(E)  =', te
           bte = te
          end if
         end if

	end do ! ispin

c       write data in external output-file
        if (output) then
	  if (nspin == 1) then
	   write(outfile,fmt='(1x,f14.10,f16.8,2x,E20.10)')
     &           energy, (energy-efermi)*hartree, te
          else
	   write(outfile,fmt='(1x,f14.10,f16.8,2x,2(E20.10))')
     &           energy, (energy-efermi)*hartree, ate, bte
       	  end if
	end if

        if (do_cond_channels) then
         do ispin = 1, nspin
          write(chfile+ispin,fmt='(1x,f14.10,f16.8,2x,E20.10,4x,15(E18.10))')
     &          energy, (energy-efermi)*hartree, tmp_sum, (channels(nsaos+1-k,ispin),k=1,5)
         end do ! ispin
        end if ! do_cond_channels 

        n = n + 1 
	energy = ener + n * estep
       end do ! e-points
 
       if (output) then
         write(outfile,fmt='(a)') '#end'
         close(outfile)
       end if

       print '(/,a,a,a)', 
     &       ' transmission is written to a file "',trim(output_file_name),'"'

c      restore parameters for energy mesh:
       if (evunits) then
        ener = ener_save ; estep = estep_save ; eend = eend_save
       end if

c      if scattering wave functions are required, 
c      call a separate routine here ...
       if (do_scatt_wave_func) call compute_scatt_wave_func(5)

       deallocate(gf,bgb_left,bgb_right,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE tcycle]: ',
     &           'impossible to deallocate temp. arrays'
        print *, 'nevertheless, proceed further ...'
       end if

c      >>> inserted by RK and AB:
c      take care of conduction channels
       if (do_cond_channels) then
         deallocate(g12b_left,channels,stat=ierr)
         if (ierr.ne.0) then
           print '(/,a,a)', ' [SUBROUTINE tcycle]: ',
     &           'impossible to deallocate b-channel temp. arrays'
           print *, 'nevertheless, proceed further ...'
         end if
       end if ! do_cond_channels
c      <<< done with insert by RK and AB...

      end subroutine tcycle

c     >>> inserted by RK and AB:
      subroutine transform12_gamma(ingamma,outgamma,ispin)
c      ****************************
c      outgamma = ingamma^{1/2} * b
c      ****************************
       implicit none
       double precision, intent(in) :: ingamma(:)
       complex(8), intent(out)      :: outgamma(:,:)
       integer,    intent(in)       :: ispin

       integer n, m
       double precision tmp_gamma12

       do n = 1, nsaos
        tmp_gamma12 = dsqrt( dabs(ingamma(n)) )
        do m = 1, nsaos
         outgamma(n,m) = tmp_gamma12 * heigvec(n,m,ispin)
        end do
       end do

      end subroutine transform12_gamma
c     <<< done with insert by RK and AB

      subroutine transform_gamma(ingamma,outgamma,ispin,lead)
c      *****************************************************************
c      if lead='l'.or.'L' then : outgamma = b+ * ingamma * b  
c      if lead='r'.or.'R' then : outgamma = b^{-1} * ingamma * b+^{-1}  
c      with b=<heigvec> being a complex matrix of eigenvectors of h_full
c      *****************************************************************
       implicit none
       double precision, intent(in) :: ingamma(:)
       complex(8), intent(out)      :: outgamma(:,:)
       integer,    intent(in)       :: ispin  
       character,  intent(in)       :: lead  
       
c      temporary arrays
       complex(8), allocatable :: b(:,:), bg(:,:)  
       integer n,m, ierr
      
       allocate(b(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE transform_gamma]: <b> allocation failure'
       end if

       allocate(bg(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE transform_gamma]: <bg> allocation failure'
       end if

       if (lead=='l'.or.lead=='L') then
        forall (n=1:nsaos,m=1:nsaos)
         b(n,m) = dconjg(heigvec(m,n,ispin))        
        end forall
       else       
c       lead=='r'.or.lead=='R'        
        forall (n=1:nsaos,m=1:nsaos)
         b(n,m) = invheigvec(n,m,ispin)        
        end forall
       end if
       
       do m = 1, nsaos
        do n = 1, nsaos
         bg(n,m) = b(n,m) * ingamma(m)
	end do
       end do
       
       call cmplx_ab('n','c',nsaos,bg,b,outgamma)
       
       deallocate(b,bg,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE transform_gamma]: ',
     &           'impossible to deallocate temp. arrays'
        print *, 'anyway, proceed further ...'
       end if
       
      end subroutine transform_gamma

c     >>> inserted by RK and AB:
      subroutine evaluate_eigenchannels(epoint,ispin,tau,sum)
       double precision, intent (in)  :: epoint
       integer, intent (in)           :: ispin
       double precision, intent (out) :: tau(:), sum

c      some temporary arrays:
       complex(8), allocatable :: g_gamma_g(:,:), tmp_mat(:,:), tmp_tt(:,:)
       integer     n, m, k, ierr

       allocate(g_gamma_g(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUB. evaluate_eigenchannels]: <g_gamma_g> allocation failure'
       end if
       allocate(tmp_mat(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUB. evaluate_eigenchannels]: <tmp_mat> allocation failure'
       end if
       allocate(tmp_tt(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUB. evaluate_eigenchannels]: <tmp_tt> allocation failure'
       end if

c      evaluate GF
       call greensf(epoint,ispin,gf)

       do n = 1, nsaos
        do m = 1, nsaos
          g_gamma_g(n,m) =  gf(n) * bgb_right(n,m,ispin) * dconjg(gf(m))
        end do
       end do

       call cmplx_ab('n','n',nsaos,g12b_left(:,:,ispin),g_gamma_g,tmp_mat)
       call cmplx_ab('n','c',nsaos,tmp_mat,g12b_left(:,:,ispin),tmp_tt)

       if (do_scatt_wave_func) then
        call find_xspectrum(tmp_tt,tau,utt(:,:,ispin),nsaos)
       else
        call find_spectrum(tmp_tt,tau,nsaos)
       end if    

       sum = 0.0d0
       do k = 1, nsaos
        sum = sum + tau(k)
       end do
c      write(50+ispin,fmt='(1x,f14.10,f16.8,2x,E20.10,4x,15(E18.10))')
c    &           epoint, (epoint-efermi)*hartree, sum, (tau(nsaos+1-k),k=1,5)

       deallocate(g_gamma_g,tmp_mat,tmp_tt,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUB. evaluate_eigenchannels]: ',
     &    'impossible to deallocate temporary arrays'
        print *, 'nevertheless, proceed further ...'
       end if

      end subroutine evaluate_eigenchannels
c     >>> done with insert by RK and AB

      double precision function trans(epoint,ispin)
c      evaluate t(e) := trace[gamma_left * g(e) * gamma_right * g+(e)]      
c      basis of eigenvectors is used  
       implicit none
       
       double precision, intent (in) :: epoint
       integer, intent (in)          :: ispin
       
       complex(8), allocatable :: lmat(:,:),
     &                            rmat(:,:)
       complex(8) tmptrans
       integer :: n, m, ierr ! , monitor = 500
      
       allocate(lmat(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
	stop '[FUNCTION trans]: <lmat> allocation failure'
       end if

       allocate(rmat(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
        print *
	stop '[FUNCTION trans]: <rmat> allocation failure'
       end if
   
c      evaluate GF
       call greensf(epoint,ispin,gf)       
       
       do m = 1, nsaos
        do n = 1, nsaos
c         here bgb_left = b+ * gamma_left * b 	
          lmat(n,m) = bgb_left(n,m,ispin) * gf(m)
	end do
       end do

       do m = 1, nsaos
        do n = 1, nsaos
c         here bgb_right = b^{-1} * gamma_right * b+^{-1} 	
          rmat(n,m) = bgb_right(n,m,ispin) * dconjg(gf(m))
	end do
       end do
       
       tmptrans = czero
       do n = 1, nsaos
        do m = 1, nsaos
         tmptrans = tmptrans + lmat(n,m)*rmat(m,n)
	end do 
       end do

       deallocate(lmat,rmat,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [FUNCTION trans]: ', 
     &	  'impossible to deallocate temporary arrays'
	print *, 'nevertheless, proceed further ...' 
       end if
  
c       print *
c       print *, ' checking:  TRACE = ', tmptrans
c       print *
c      return T(E)
       trans = dble(tmptrans)
      
      end function trans

      end module transmission     

