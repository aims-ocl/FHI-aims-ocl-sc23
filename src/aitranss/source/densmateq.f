c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           march/april 2008
c     last revision:  jan 2012
c###########################################################

      module density_matrix0
c      ***********************************************      
c      computes an equilibrium density matrix D0
c      checks number of electrons: N = tr(D0*S)    
c  
c      -- update added related to the lda+u functional
c
c      ***********************************************
       use globalvars
       use read_externals
       use neq_density_matrix
       use tools
       implicit none

c      equilibrium density matrix (per spin); last index == ispin
       double precision, private, allocatable :: densmat0(:,:,:)
       
       real(4), private :: stime, ftime, secs
       integer, private :: mnts
           
       contains
       
       double precision function mtrace(d_mat,s_mat,n,ispin)
        implicit none
        double precision, intent (in) :: d_mat(:,:,:), s_mat(:,:)
        integer, intent (in) :: n, ispin
        integer i,j    
        double precision tmp_trace
	
	tmp_trace = 0.0d0
	do i = 1, n
	 do j = 1, n
	  tmp_trace = tmp_trace + d_mat(i,j,ispin)*s_mat(j,i)
	 end do
	end do
        mtrace = tmp_trace

       end function mtrace

       double precision function xmtrace(d_mat,n,ispin)
        implicit none
        double precision, intent (in) :: d_mat(:,:,:)
        integer, intent (in) :: n, ispin
        integer i    
        double precision tmp_trace
	
	tmp_trace = 0.0d0
	do i = 1, n
	  tmp_trace = tmp_trace + d_mat(i,i,ispin)
        end do
        xmtrace = tmp_trace

       end function xmtrace
       
       subroutine dens_matrix0(mo_cfs,d_mat0,ispin,nocc)
        implicit none
        double precision, intent (in)  :: mo_cfs(:,:,:)       
	double precision, intent (out) :: d_mat0(:,:,:)
        integer, intent (in) :: ispin, nocc
	integer  iorb, n, m
	double precision c_n, c_m

	forall (n=1:nsaos,m=1:nsaos)
	 d_mat0(n,m,ispin) = 0.0d0
	end forall

	do iorb = 1, nocc ! cycle over occupied orbitals
c        print '(a,i4)', '   molecular orbital :', iorb
	 do n = 1, nsaos 
          c_n = mo_cfs(n,iorb,ispin)
          do m = n, nsaos
           c_m = mo_cfs(m,iorb,ispin)
	   d_mat0(n,m,ispin) = d_mat0(n,m,ispin) + c_n*c_m
	   d_mat0(m,n,ispin) = d_mat0(n,m,ispin)
	  end do
         end do
	end do
	
       end subroutine dens_matrix0
       
       subroutine init_densmat0
c       this routine is called for testing purposes
        implicit none
        integer  ierr 
	double precision nel, nel_a, nel_b

        if (.not.calc_densmat0) then
c        print '(/,a)', ' density matrix calculation is switched off'
	 return
	end if

        if (ldau.and.no_leads) then
c         calculation of eq. density matrix is done later on
	  return
	end if  
	
	allocate(densmat0(nsaos,nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
        stop 
     &    '[SUBROUTINE init_densmat0]: <densmat0> allocation failure'
        end if
        
	print '(/,a,/)', ' CALCULATING EQUILIBRIUM DENSITY MATRIX >>>'

        if (nspin==1) then 
         call dens_matrix0(mo_coeff,densmat0,1,occ)
	 nel = mtrace(densmat0,smat,nsaos,1) 
	 print 500, ' number of electrons = ', 2.0d0*nel	
	else
c	 print *, '  === alpha-electrons ==='
	 call dens_matrix0(mo_coeff,densmat0,1,a_occ)  ! alpha electrons
	 nel_a = mtrace(densmat0,smat,nsaos,1) 
	 print 501, ' number of alpha-electrons = ', nel_a

c        print '(/,a)', '   === beta-electrons ==='
	 call dens_matrix0(mo_coeff,densmat0,2,b_occ)  ! beta  electrons
         nel_b = mtrace(densmat0,smat,nsaos,2)
	 print 501, ' number of beta-electrons  = ', nel_b	
	end if

        deallocate(densmat0,stat=ierr)
        if (ierr.ne.0) then
          print '(/,a,a)', '  [SUBROUTINE init_densmat0]: ',
     &             ' impossible to deallocate <densmat0>'
          print *, ' anyway, proceed further ...'
        end if

	print '(/,a)', ' <<< DONE WITH EQUILIBRIUM DENSITY MATRIX'

500    format(/,a,f12.6)
501    format(a,f12.6)
       end subroutine init_densmat0

       subroutine update_eq_dmat
c       ***********************************************
c       with a given number of electrons in the system
c       and spectrum of the hamiltonian, construct 
c       an equilibrium density matrix ;
c       updates group $dens_matrix in <tcontrol> ;
c       save <dmat> & <smat> matrices to external files
c       ***********************************************

        integer          n, nn, m, ispin, ierr
        double precision nel, nel_a, nel_b
c       temporary set of molecular orbitals (just alpha + beta channels)
c       which are not properly ordered
        type(spectrum), allocatable :: fullsetmo(:), tmpset(:)
        integer, allocatable        :: indarr(:)

c       temporary arrays
        double precision, allocatable :: tmp_dmat(:,:), tmp_sd(:,:)        
	
        allocate(fullsetmo(2*nsaos),stat=ierr)
        if (ierr.ne.0) then
        stop
     &   '[SUBROUTINE update_eq_dmat]: <fullsetmo> allocation failure'
        end if

	print '(/,a)', ' <HUBBARD-U>: EQUILIBRIUM DENSITY MATRIX IS TO BE MODIFIED >>>'
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
     &    '[SUBROUTINE update_eq_dmat]: <tmpset> allocation failure'
         end if
         allocate(indarr(2*nsaos),stat=ierr)
         if (ierr.ne.0) then
         stop
     &    '[SUBROUTINE update_eq_dmat]: <indarr> allocation failure'
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
         print *, 'done'

         deallocate(tmpset,indarr,stat=ierr)
         if (ierr.ne.0) then
          print '(/,a,a)', '  [SUBROUTINE update_eq_dmat]: ',
     &             ' impossible to deallocate temporary arrays'
          print *, ' anyway, proceed further ...'
         end if 
	end if ! nspin

c       define occupation numbers
        print '(/,2x,a,i5)',     '-- number of electrons: ', nelectr
        if (nspin.eq.1) then
	  occ = nelectr/2
	  print '(1x,a,i5,a,i5)', '-- occupation numbers:  na =', occ, '  nb =', occ 
c         update <tcontrol> file
          call updatetcontrol(trim(tcntrl_file_name),'$dens_matrix',occ,0.0d0)
        else 
c         nspin == 2
	  a_occ = 0 ; b_occ = 0 
	  do nn = 1, nelectr
	   ispin = fullsetmo(nn)%spin
	   if (ispin.eq.1) then ; a_occ = a_occ + 1 
	   else ;  b_occ = b_occ + 1
	   end if              
	  end do
	  print '(2x,a,i5,a,i5)', '-- occupation numbers:  na =', a_occ, '  nb =', b_occ 
c         update <tcontrol> file
          call updatetcontrol(trim(tcntrl_file_name),'$dens_matrix',a_occ,dble(b_occ))
        end if ! nspin.eq.1

c       further, we have to compute a modified equilibrium density matrix  ...
c       and save it to external files (which are to be read by TURBOMOLE or FHI-aims)

        allocate(densmat0(nsaos,nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
        stop
     &    '[SUBROUTINE update_eq_dmat]: <densmat0> allocation failure'
        end if
  
        if (nspin==1) then
         call dens_matrix0(mo_orth,densmat0,1,occ)
	 nel = xmtrace(densmat0,nsaos,1) 
         print '(/,2x,a,f12.6)', 'number of electrons = ', 2.0d0*nel
        else
c        alpha channel
         call dens_matrix0(mo_orth,densmat0,1,a_occ) 
	 nel_a = xmtrace(densmat0,nsaos,1) 
         print '(/,2x,a,f12.6)', 'number of alpha-electrons = ', nel_a
c        beta  channel	 
         call dens_matrix0(mo_orth,densmat0,2,b_occ) 
	 nel_b = xmtrace(densmat0,nsaos,2) 
         print '(2x,a,f12.6)', 'number of beta-electrons  = ', nel_b
        end if
        print '(/,a)', ' <HUBBARD-U>: <<< DONE WITH UPDATE OF EQUILIBRIUM DENSITY MATRIX'

c       testing >>>
c       in case of domain-wall solution, symmetrize the density matrix
c        if (z_symmetrize) then
c         call z_updatedmat(densmat0)
c        end if ! symmetrize
c        if (x_symmetrize) then
c         call x_updatedmat(densmat0)
c        end if ! symmetrize
c        if (y_symmetrize) then
c         call y_updatedmat(densmat0)
c        end if ! symmetrize
c       <<< testing

c       prints out atomic loewdin charges/magn.moments
        if (qoutput) call qprint(densmat0)

c       transform dens.matrix to the non-orthogonal basis
        allocate(tmp_dmat(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        stop
     &    '[SUBROUTINE update_eq_dmat]: <tmp_dmat> allocation failure'
        end if

        allocate(tmp_sd(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
        stop
     &    '[SUBROUTINE update_eq_dmat]: <tmp_sd> allocation failure'
        end if

        print '(/,1x,a,$)',
     &  'transform density-matrix to the non-orthogonal basis ...'
        call cpu_time(stime)
        do ispin = 1, nspin
c        put a local copy of dens.matrix to temporary array  >>>
         forall (n=1:nsaos,m=1:nsaos)
           tmp_dmat(n,m) = densmat0(n,m,ispin)
         end forall
c        tmp_dmat --> smat12inv * tmp_dmat * smat12inv
c        1st step: tmp_sd = sma12inv * tmp_dmat
         call dble_ab('n','n',nsaos,smat12inv,tmp_dmat,tmp_sd)
c        2nd step: tmp_dmat = tmp_sd * smat12
         call dble_ab('n','n',nsaos,tmp_sd,smat12inv,tmp_dmat)

c        save a local copy of dens.matrix >>>
         forall (n=1:nsaos,m=1:nsaos)
c          we symmetrize results to avoid numerical uncertainties	 
           densmat0(n,m,ispin) = 0.5d0*(tmp_dmat(n,m) + tmp_dmat(m,n))
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
c          call updatedmat(densmat0,1)
c       end if ! symmetrize

c       -- admix a new density with the previous one
c       -- save result to external file
        call mixdensmat(densmat0)

        deallocate(densmat0,tmp_dmat,stat=ierr)
        if (ierr.ne.0) then
          print '(/,a,a)', '  [SUBROUTINE update_eq_dmat]: ',
     &             ' impossible to deallocate temporary arrays '
          print *, ' anyway, proceed further ...'
        end if

        deallocate(fullsetmo,stat=ierr)
        if (ierr.ne.0) then
          print '(/,a,a)', '  [SUBROUTINE update_eq_dmat]: ',
     &             ' impossible to deallocate <fullsetmo>'
          print *, ' anyway, proceed further ...'
        end if 

c500    format(/,a,f12.6)
c501    format(a,f12.6)
 
       end subroutine update_eq_dmat
      
       double precision function get_fermi_level(num_electrons,icall)
c       ************************************
c       searches for the fermi level, based 
c       on electron number, as (HOMO+LUMO)/2      
c       ************************************
        implicit none
c       input: number of electrons >>> 
        integer num_electrons, icall

c       local vars >>>	
        double precision efguess
        integer          n, nn, ierr
c       temporary set of molecular orbitals (just alpha + beta channels)
c       which are not properly ordered 
        type(spectrum), allocatable :: tmpset(:)
        integer, allocatable :: indarr(:)
c       character chann
   
        allocate(fullsetmo(2*nsaos),stat=ierr)
        if (ierr.ne.0) then
        stop
     &   '[FUNCTION get_fermi_level]: <fullsetmo> allocation failure'
        end if

c        print '(/,1x,a)', 
c     &	   '>>> define Fermi level, as (E_HOMO+E_LUMO)/2 for "extended molecule" :'
c 	 print '(/,1x,a,i5)', 'overall number of electrons = ', nelectr
c        print '(/,1x,a,$)',  'merge mo-sets of two spin channels ...'
        
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
     &    '[FUNCTION get_fermi_level]: <tmpset> allocation failure'
         end if
         allocate(indarr(2*nsaos),stat=ierr)
         if (ierr.ne.0) then
         stop
     &    '[FUNCTION get_fermi_level]: <indarr> allocation failure'
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
          print '(/,a,a)', '  [FUNCTION get_fermi_level]: ',
     &             ' impossible to deallocate temporary arrays'
          print *, ' anyway, proceed further ...'
         end if
       
	end if ! nspin
c       print *, 'done'

        if ( num_electrons.lt.(2*nsaos) ) then
         efguess  = 0.5d0*(fullsetmo(num_electrons)%en + 
     &	                   fullsetmo(num_electrons+1)%en)	
        else
c        this is impossible case, but to be on the safe side ...
         efguess = 0.0d0	
	end if

        if (icall.eq.1) print '(/,1x,a,f16.12,a)', 'Fermi level: EFermi = ', efguess, ' H' 
        
        deallocate(fullsetmo,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', '  [FUNCTION get_fermi_level]: ',
     &            ' impossible to deallocate temp. array <fullsetmo>'
         print *, ' anyway, proceed further ...'
        end if
  
c       put a new value of the fermi-energy to <tcontrol>
        call updatetcontrol(trim(tcntrl_file_name),'$efermi',0,efguess)

        get_fermi_level = efguess   

       end function get_fermi_level
       
      end module density_matrix0
