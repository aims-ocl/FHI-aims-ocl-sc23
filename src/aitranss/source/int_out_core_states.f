c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           july 2013
c     last revision:  july 2013
c###########################################################

      module int_out_core_states
c     ***********************************************      
c     integrate out core states, update global arrays 
c     only valence electrons are left further on
c     ***********************************************      
       use globalvars
       use tools
       use rmatrix
       use density_matrix0
       use hamiltonian
     
       implicit none
 
       double precision, private, allocatable :: 
     &                   h_val(:,:,:), u_val(:,:,:), valence_states(:,:),
     &                   h_core(:,:,:), u_core(:,:,:), core_states(:,:),
     &                   h_core_val(:,:,:), tmp_omat(:,:)

       type(cgto), private, allocatable :: tmp_aos(:,:)

       real(4), private ::  stime, ftime, secs
       integer, private ::  mnts  

      contains

      subroutine wrapper_valence_states
c      ******************************************************
c      wrapper routine: call for necessary tasks to 
c      integrate out core states and leave valence ones only
c      ******************************************************
       implicit none
       integer ispin, ierr, tmpfile, tmpfile1, tmpfile2, 
     &         n, i, j, itype, iorb, jorb
       double precision tmp_efermi

       print '(/,a)', ' INTEGRATING OUT CORE STATES >>>'

       allocate(tmp_omat(nvalence,nvalence),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <tmp_omat> allocation failure'
       end if

       allocate(h_val(nvalence,nvalence,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <h_val> allocation failure'
       end if
       allocate(u_val(nvalence,nvalence,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <u_val> allocation failure'
       end if
       allocate(valence_states(nvalence,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <val_states> allocation failure'
       end if
      
       allocate(h_core(ncore_states,ncore_states,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <h_core> allocation failure'
       end if
       allocate(u_core(ncore_states,ncore_states,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <h_core> allocation failure'
       end if
       allocate(core_states(ncore_states,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <h_core> allocation failure'
       end if

       allocate(h_core_val(ncore_states,nvalence,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <h_core_val> allocation failure'
       end if
      
       print '(/,2x,a,$)', 'split out hamiltonian in blocks ... '
       do ispin = 1, nspin
        call split_hamiltonian(ispin)
       end do
       print '(a)', 'done'

       do ispin = 1, nspin
       
        if (nspin.eq.2) then
         if (ispin.eq.1) then 
           print '(/,2x,a)', '*** alpha channel ***'
         else 
           print '(/,2x,a)', '*** beta  channel ***'
         end if
        end if

        call cpu_time(stime)
        print '(/,2x,a,$)', 'valence states: solving eigenvalue problem ... '
        call real_spectrum(h_val(:,:,ispin),valence_states(:,ispin),u_val(:,:,ispin),nvalence)
        print '(a)', 'done'

        call cpu_time(ftime)
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts*60.0d0
        print '(2x,a,i3,a,f5.2,a)', 'time spent:', mnts, ' min ', secs, ' sec'


        call cpu_time(stime)
        print '(/,2x,a,$)', 'core states: solving eigenvalue problem ... '
        call real_spectrum(h_core(:,:,ispin),core_states(:,ispin),u_core(:,:,ispin),ncore_states)
        print '(a)', 'done'

        call cpu_time(ftime)
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts*60.0d0
        print '(2x,a,i3,a,f5.2,a)', 'time spent:', mnts, ' min ', secs, ' sec'

        if (testing) then
    
         tmpfile1 = 31
         tmpfile2 = 32

         if (nspin == 1) then

          open(tmpfile1,file='mos.h_val.tmp',status='unknown',iostat=ierr)
          if (ierr.ne.0) then
           stop ': can not open mos.h_val.tmp!'
          end if
          write(tmpfile1,'(a)') '$mos.h_val.tmp'
          open(tmpfile2,file='mos.h_core.tmp',status='unknown',iostat=ierr)
          if (ierr.ne.0) then
           stop ': can not open mos.h_core.tmp!'
          end if
          write(tmpfile2,'(a)') '$mos.h_core.tmp'
    
         else  ! spin-polarized case
          if (ispin == 1) then
           open(tmpfile1,file='alpha.h_val.tmp',status='unknown',iostat=ierr)
           if (ierr.ne.0) then
            stop ': can not open alpha.h_val.tmp!'
           end if
           write(tmpfile1,'(a)') '$alpha.h_val.tmp'
           open(tmpfile2,file='alpha.h_core.tmp',status='unknown',iostat=ierr)
           if (ierr.ne.0) then
            stop ': can not open alpha.h_core.tmp!'
           end if
           write(tmpfile2,'(a)') '$alpha.h_core.tmp'
          else
           open(tmpfile1,file='beta.h_val.tmp',status='unknown',iostat=ierr)
           if (ierr.ne.0) then
            stop ': can not open beta.h_val.tmp!'
           end if
           write(tmpfile1,'(a)') '$beta.h_val.tmp'
           open(tmpfile2,file='beta.h_core.tmp',status='unknown',iostat=ierr)
           if (ierr.ne.0) then
            stop ': can not open beta.h_core.tmp!'
           end if
           write(tmpfile2,'(a)') '$beta.h_core.tmp'
          end if
         end if

         write(tmpfile1,'(a,i5)') '#nvalence =', nvalence
         write(tmpfile2,'(a,i5)') '#ncore    =', ncore_states

         do n = 1, nvalence
           write(tmpfile1,'(i6,4x,e20.14)') n, valence_states(n,ispin)
         end do
        write(tmpfile1,'(a)') '$end'

         do n = 1, ncore_states
           write(tmpfile2,'(i6,4x,e20.14)') n, core_states(n,ispin)
         end do
         write(tmpfile2,'(a)') '$end'
        
         close(tmpfile1)
         close(tmpfile2)

        end if ! testing

c       compute self-energy due to core-states and 
c       update "valence-valence" block of the hamiltonian

        call cpu_time(stime)
        print '(/,2x,a,$)', 'updating Kohn-Sham hamiltonian  ... '
        call update_hamiltonian(ispin)
        print '(a)', 'done'

        call cpu_time(ftime)
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts*60.0d0
        print '(2x,a,i3,a,f5.2,a)', 'time spent:', mnts, ' min ', secs, ' sec'

       end do ! ispin

c      update global arrays and solve
c      eigenvalue for the modified hamiltonian (with core self-energy) 

       if (allocated(h0mat)) deallocate(h0mat,stat=ierr)
       allocate(h0mat(nvalence,nvalence,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <h0mat> allocation failure'
       end if
       h0mat = h_val
 
       if (allocated(mo_orth)) deallocate(mo_orth,stat=ierr)
       allocate(mo_orth(nvalence,nvalence,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <mo_orth> allocation failure'
       end if

       if (allocated(mo_en)) deallocate(mo_en,stat=ierr)
       allocate(mo_en(nvalence,nspin),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <mo_en> allocation failure'
       end if

       do ispin = 1, nspin

        if (nspin.eq.1) then
         print '(/,2x,a,$)', 'modified valence states: solving eigenvalue problem ... '
        else if (ispin.eq.1) then
         print '(/,2x,a)', '*** alpha channel *** '
         print '(2x,a,$)',   'modified valence states: solving eigenvalue problem ... '
        else
         print '(/,2x,a)', '*** beta  channel *** '
         print '(2x,a,$)',   'modified valence states: solving eigenvalue problem ... '
        end if
        
        call cpu_time(stime)
        call real_spectrum(h0mat(:,:,ispin),mo_en(:,ispin),mo_orth(:,:,ispin),nvalence)
        print '(a)', 'done' 
        
        call cpu_time(ftime)  
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts*60.0d0
        print '(2x,a,i3,a,f5.2,a)', 'time spent:', mnts, ' min ', secs, ' sec'

       end do 

ccc    !!! UPDATE global array's dimension nsaos >>>>
       nsaos0 = nsaos
       nsaos  = nvalence
       print '(/,2x,a,i6)', 'from now on, dimension of the hilbert space (nsaos) =', nsaos
       print '(2x,a,i6)', 'amount of valence electrons in the system =', nelectr-2*ncore_states
       tmp_efermi = get_fermi_level(nelectr-2*ncore_states,2)
       do ispin = 1, nspin
        do n = 1, nvalence
         h0mat(n,n,ispin) = h0mat(n,n,ispin) - tmp_efermi + efermi
         mo_en(n,ispin) = mo_en(n,ispin) - tmp_efermi + efermi
        end do
       end do    
c      tmp_efermi = get_fermi_level(nelectr-2*ncore_states,2)

ccc    !!! UPDATE global electron number >>>>
       nelectr0 = nelectr
       nelectr  = nelectr-2*ncore_states
       call updatetcontrol(trim(tcntrl_file_name),'$valence_electrons',nelectr,0.0d0)
              
       if (testing) then
        
         do ispin = 1, nspin
        
         tmpfile = 31
         if (nspin == 1) then

          open(tmpfile,file='mos.valence.tmp',status='unknown',iostat=ierr)
          if (ierr.ne.0) then
           stop ': can not open mos.valence.tmp!'
          end if
          write(tmpfile,'(a)') '$mos.valence.tmp'
    
         else  ! spin-polarized case

          if (ispin == 1) then
           open(tmpfile,file='alpha.valence.tmp',status='unknown',iostat=ierr)
           if (ierr.ne.0) then
            stop ': can not open alpha.valence.tmp!'
           end if
           write(tmpfile,'(a)') '$alpha.valence.tmp'
          else
           open(tmpfile,file='beta.valence.tmp',status='unknown',iostat=ierr)
           if (ierr.ne.0) then
            stop ': can not open beta.valence.tmp!'
           end if
           write(tmpfile,'(a)') '$beta.valence.tmp'
          end if

         end if ! nspin

         write(tmpfile,'(a,i5)') '#nvalence =', nvalence

         do n = 1, nvalence
           write(tmpfile,'(i6,4x,e20.14)') n, mo_en(n,ispin)
         end do
         write(tmpfile,'(a)') '$end'

         close(tmpfile)

        end do ! nspin
       end if ! testing

c      deallocalte temporary arrays 
       deallocate(h_val,h_core,h_core_val,u_val,valence_states,u_core,core_states,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,/,a)',
     &   ' [SUBR. wrapper_valence_states]: deallocation of temp. arrays failed ',
     &   ' ... will proceed further anyway'
       end if

c      >>> update aos array

       print '(/,2x,a,$)', 'updating atomic orbitals array ... '

       allocate(tmp_aos(num_atom_types,num_val_gto),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <tmp_aos> allocation failure'
       end if

       do itype = 1, num_atom_types
         jorb = 0
         do iorb = 1, n_basis_func(itype)
           if (aos(itype,iorb)%valence) then
             jorb = jorb + 1 
             tmp_aos(itype,jorb)%ngto = aos(itype,iorb)%ngto
             tmp_aos(itype,jorb)%icoeff = aos(itype,iorb)%icoeff
             tmp_aos(itype,jorb)%xi = aos(itype,iorb)%xi
             tmp_aos(itype,jorb)%lm = aos(itype,iorb)%lm
             tmp_aos(itype,jorb)%valence = .true.
           end if
         end do
       end do

c      update gloval variables
       num_at_gto = num_val_gto
       forall(itype=1:num_atom_types)
        n_basis_func(itype) = n_val_basis_func(itype)
       end forall    

       deallocate(aos,stat=ierr)
       allocate(aos(num_atom_types,num_val_gto),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <aos> allocation failure'
       end if
       
       do itype = 1, num_atom_types
         do iorb = 1, n_val_basis_func(itype)
             aos(itype,iorb)%ngto = tmp_aos(itype,iorb)%ngto
             aos(itype,iorb)%icoeff = tmp_aos(itype,iorb)%icoeff
             aos(itype,iorb)%xi = tmp_aos(itype,iorb)%xi
             aos(itype,iorb)%lm = tmp_aos(itype,iorb)%lm
             aos(itype,iorb)%valence = .true.
         end do
       end do

c      deallocalte temporary arrays 
       deallocate(tmp_aos,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,/,a)',
     &   ' [SUBR. wrapper_valence_states]: deallocation of <tmp_aos> failed ',
     &   ' ... will proceed further anyway'
       end if
  
       print '(a)', 'done'
c      <<< update of the aos array

       print '(/,a)', ' <<< DONE WITH CORE STATES'

c      update overlap matrix >>>>

       deallocate(smat,stat=ierr)
       allocate(smat(nvalence,nvalence),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <smat> allocation failure'
       end if

       forall(i=1:nvalence,j=1:nvalence)
        smat(i,j) = tmp_omat(i,j)
       end forall    

       deallocate(tmp_omat,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,/,a)',
     &   ' [SUBR. wrapper_valence_states]: deallocation of <tmp_omat> array failed ',
     &   ' ... will proceed further anyway'
       end if

       deallocate(smat12,stat=ierr)
       allocate(smat12(nvalence,nvalence),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <smat> allocation failure'
       end if

       deallocate(smat12inv,stat=ierr)
       allocate(smat12inv(nvalence,nvalence),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE wrapper_valence_states]: <smat> allocation failure'
       end if

       call sqrt_smat(smat,smat12,smat12inv,nvalence,4)    

c      <<< done with update of the overlap matrix

      end subroutine wrapper_valence_states

      subroutine update_hamiltonian(ispin)
c      ****************************************************
c      computes self-energy correction due to core orbitals
c      updates hamiltonian of the valence states
c      ****************************************************
       implicit none
       integer, intent(in) :: ispin
       integer  ierr
       
       double precision, allocatable :: core_self_energy(:,:), uv(:,:), euv(:,:)
c      uv  = u_val^+ * h_core_val
c      euv = (efermi-core_states)^{-1} * u_val^+ * h_core_val 
c      core_self_energy = (uv)^+ * euv 

       double precision, parameter :: ieta = 1.0d-8
       double precision, parameter :: xeps = 1.0d-8
       double precision  de
       integer           m,n
c      logical           flag

       allocate(core_self_energy(nvalence,nvalence),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE update_hamiltonian]: <core_self_energy> allocation failure'
       end if

       allocate(uv(ncore_states,nvalence),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE update_hamiltonian]: <uv> allocation failure'
       end if

       allocate(euv(ncore_states,nvalence),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE update_hamiltonian]: <euv> allocation failure'
       end if

       uv = 0.0d0
       call xdble_ab('t','n',ncore_states,ncore_states,nvalence,
     &               u_core(:,:,ispin),h_core_val(:,:,ispin),uv)

       euv = 0.0d0
       do m = 1, ncore_states 
        de = efermi-core_states(m,ispin)
        do n = 1, nvalence
c        euv(m,n) = (de/( de*de + ieta*ieta)) * uv(m,n)
         euv(m,n) = (1.0d0/de) * uv(m,n)
        end do
       end do

       core_self_energy = 0.0d0
       call xdble_ab('t','n',nvalence,ncore_states,nvalence,uv,euv,core_self_energy)

       do m = 1, nvalence
        do n = 1, nvalence
         h_val(m,n,ispin) = h_val(m,n,ispin) + core_self_energy(m,n)
        end do
       end do

c       flag = .true.
c       if (testing) then
c        do m = 1, nvalence
c         do n = 1, m
c          if (dabs(h_val(m,n,ispin)-h_val(n,m,ispin)) > xeps) then
c           write(50+ispin,*) m, n, dabs(h_val(m,n,ispin)-h_val(n,m,ispin))
c           flag = .false.
c          end if
c         end do
c        end do
c       end if
c       write(50+ispin,*) 'updated hamiltonian, flag = ', flag 

       deallocate(core_self_energy,uv,euv,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,/,a)',
     &   ' [SUBR. wrapper_valence_states]: deallocation of temp. arrays failed ',
     &   ' ... will proceed further anyway'
       end if

      end subroutine update_hamiltonian

      subroutine split_hamiltonian(ispin)
c      ********************************************
c      splits hamiltonian to blocks asocciated with 
c      valence & core orbitals:
c              | h_val  v+     |
c      hmat =  | v      h_val  |
c      ********************************************
       implicit none
       integer, intent(in) :: ispin
       integer i, j, i1, j1, iat, jat, iatype, jatype, 
     &         iorb, jorb, norb, morb
c      integer m, n  
       
       double precision, parameter :: xeps = 1.0d-8
c      logical flag

c      picks up "valence-valence" block of h0 and redirect it to h_val
       i = 0 ; i1  = 0
       do iat = 1, num_atoms
        iatype = atom(iat)%atype
        norb = n_basis_func(iatype)
        do iorb = 1, norb 
         i = i + 1
         if (aos(iatype,iorb)%valence) then
          i1 = i1 + 1 
          j = 0 ; j1  = 0
          do jat = 1, num_atoms
           jatype = atom(jat)%atype
           morb = n_basis_func(jatype)
           do jorb = 1, morb
            j = j + 1
            if (aos(jatype,jorb)%valence) then 
             j1 = j1 + 1 
             h_val(i1,j1,ispin) = h0mat(i,j,ispin)
c            save valence-valence block of the overlap matrix as well
             if (ispin.eq.1) tmp_omat(i1,j1) = smat(i,j)
            end if
           end do
          end do
         end if
        end do 
       end do
c       print *
c       print *, 'i1 = ', i1, ' j1 = ', j1
c       print *, 'i  = ', i , ' j = ', j

c       flag = .true.
c       if (testing) then
c        do m = 1, nvalence
c         do n = 1, m
c          if (dabs(h_val(m,n,ispin)-h_val(n,m,ispin)) > xeps) then
c           write(30+ispin,*) m, n, dabs(h_val(m,n,ispin)-h_val(n,m,ispin))
c           flag = .false.
c          end if
c         end do
c        end do
c       end if
c       write(30+ispin,*) 'valence block, flag = ', flag 

c      picks up "core-core" block of h0 and redirect it to h_val
       i = 0 ; i1  = 0
       do iat = 1, num_atoms
        iatype = atom(iat)%atype
        norb = n_basis_func(iatype)
        do iorb = 1, norb 
         i = i + 1
         if (.not.aos(iatype,iorb)%valence) then
          i1 = i1 + 1 
          j = 0 ; j1  = 0
          do jat = 1, num_atoms
           jatype = atom(jat)%atype
           morb = n_basis_func(jatype)
           do jorb = 1, morb
            j = j + 1
            if (.not.aos(jatype,jorb)%valence) then 
             j1 = j1 + 1 
             h_core(i1,j1,ispin) = h0mat(i,j,ispin)
            end if
           end do
          end do
         end if
        end do 
       end do
c       print *
c       print *, 'i1 = ', i1, ' j1 = ', j1
c       print *, 'i  = ', i , ' j = ', j

c       flag = .true.
c       if (testing) then
c        do m = 1, ncore_states
c         do n = 1, m
c          if (dabs(h_core(m,n,ispin)-h_core(n,m,ispin)) > xeps) then
c           write(40+ispin,*) m, n, dabs(h_core(m,n,ispin)-h_core(n,m,ispin))
c           flag = .false.
c          end if
c         end do
c        end do
c       end if
c       write(40+ispin,*) 'core block, flag = ', flag 

c      picks up "core-valence" block of h0 and redirect it to h_val
       i = 0 ; i1  = 0
       do iat = 1, num_atoms
        iatype = atom(iat)%atype
        norb = n_basis_func(iatype)
        do iorb = 1, norb 
         i = i + 1
         if (.not.aos(iatype,iorb)%valence) then
          i1 = i1 + 1 
          j = 0 ; j1  = 0
          do jat = 1, num_atoms
           jatype = atom(jat)%atype
           morb = n_basis_func(jatype)
           do jorb = 1, morb
            j = j + 1
            if (aos(jatype,jorb)%valence) then 
             j1 = j1 + 1 
             h_core_val(i1,j1,ispin) = h0mat(i,j,ispin)
            end if
           end do
          end do
         end if
        end do 
       end do
c       print *
c       print *, 'i1 = ', i1, ' j1 = ', j1
c       print *, 'i  = ', i , ' j = ', j

      end subroutine split_hamiltonian

      end module int_out_core_states

