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
c     last revision:  oct 2012
c###########################################################

      module hamiltonian
c     ***************************************************      
c     computes the hamiltonian matrix in orthogonal basis
c     ***************************************************      
       use globalvars
       use tools
       use read_externals
       use rmatrix
       use domainwall

       implicit none
 
       double precision, private, parameter :: myprec  = 1.0d-9
       double precision, private, parameter :: mysmall = 1.0d-6
       double precision, private, parameter :: tiny    = 1.0d-13

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      parameters to call lapack subroutine 'dsyevr'
       character, private, parameter ::
     &              task ='V', ! eigenvalues and eigenvectors
     &              range='A', ! all eigenvalues will be found
     &              uplo ='U'  ! upper triangle is stored

       double precision, private, parameter :: abstol = 0.0d0

       double precision, private :: vl=0.0d0, vu=1.0d0
c                          ! not actually used
       integer, private :: il=1, iu=10
c                          ! not used as well

       double precision, private, allocatable :: tmp_work(:)
       integer, private, allocatable          :: isupp(:), tmp_iwork(:)
       integer, private :: numeigval, lwork, liwork
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
       real(4), private ::  stime, ftime, secs
       integer, private ::  mnts  

      contains
      
      subroutine build_smtrxs
c      ********************************************
c      splits jobs according to read/calculate-flag 
c      for overlap matrices
c      ********************************************
       implicit none
       
       integer ierr
       logical outstat

       allocate(smat12(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE build_smtrxs]: <smat12> allocation failure'
       end if
       allocate(smat12inv(nsaos,nsaos),stat=ierr)
       if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE build_smtrxs]: <smat12inv> allocation failure'
       end if

       if (importsmat.and.storesmat) then
c       then we try to read matrices from external files
        print '(/,a)', ' REQUEST FOR OVERLAP MATRICES >>>'
        call readsmat(trim(smat12_file_name),
     &                trim(smat12inv_file_name),outstat)
        if (.not.outstat) then
         call sqrt_smat(smat,smat12,smat12inv,nsaos,1)
        else
         print '(/,a)', ' <<< DONE WITH OVERLAP MATRICES'
        end if
       else    
c       otherwise, we calculate them     
        call sqrt_smat(smat,smat12,smat12inv,nsaos,1)
       end if

      end subroutine build_smtrxs
       
      subroutine sqrt_smat(smt,smt12,smt12inv,msize,icall)
c      ***********************************************************************
c      finds out the square-root of s(overlap)-matrix and its inverse
c      -- most advanced routine to solve eigenvalue problem DSYEVR is used --
c 
c      -- icall = 1: overlap matrix for the extended molecule
c      -- icall = 2: overlap matrix for the 3d-ion(s)
c      -- icall = 3: overlap matrix for the molecule
c      -- icall = 4: overlap matrix for the valence states
c      ***********************************************************************

       implicit none

       double precision, intent(in)  :: smt(:,:)
       double precision, intent(out) :: smt12(:,:), smt12inv(:,:)
       integer, intent(in)           :: msize, icall

       double precision, allocatable :: tmps(:,:), seigval(:), 
     &                                  utrans(:,:)      
       double precision, allocatable :: tmp_mat(:,:), tmpe(:,:), 
     &                                  tmp_us(:,:), tmp_usinv(:,:)
       double precision tmp_eig, tmp_diff
       integer  ierr, n, m, k, tmpfile
       logical :: testunit = .true. , 
     &            stest = .true. ,
     &            etest = .true. 
       character(64) smtfilename
     
c      local copy of s-matrix
       allocate(tmps(msize,msize),stat=ierr) 
       if (ierr.ne.0) then
        print *
	stop
     &	 '[SUBROUTINE sqrt_smat]: <tmps> allocation failure'
       end if       

c      u-matrix with eigenvectors (in columns)
       allocate(utrans(msize,msize),stat=ierr) 
       if (ierr.ne.0) then
        print *
	stop
     &	 '[SUBROUTINE sqrt_smat]: <utrans> allocation failure'
       end if       

c      eigenvalues of s-matrix
       allocate(seigval(msize),stat=ierr) 
       if (ierr.ne.0) then
        print *
	stop
     &	 '[SUBROUTINE sqrt_smat]: <seigval> allocation failure'
       end if       

c      tmp unit-matrix
       allocate(tmpe(msize,msize),stat=ierr) 
       if (ierr.ne.0) then
        print *
	stop
     &	 '[SUBROUTINE sqrt_smat]: <tmpe> allocation failure'
       end if       

c      saving a copy of s-matrix
       forall (n=1:msize,m=1:msize) tmps(n,m) = smt(n,m)
  
       select case (icall)
        case(1)
         print '(/,a)', ' CALCULATING SQUARE-ROOT OF THE S-MATRIX >>>'
         print '(/,a,$)', '  solving eigenvalue problem ...'
        case(2)
         print '(/,a)', ' <HUBBARD-U>: CALCULATING SQUARE-ROOT OF THE O-MATRIX >>>'
         print '(/,a,$)', '  solving eigenvalue problem ...'
        case(3)
         print '(/,a)', ' <MOLECULE''s SUBSPACE>: CALCULATING SQUARE-ROOT OF THE O-MATRIX >>>'
         print '(/,a,$)', '  solving eigenvalue problem ...'
        case(4)
         print '(/,a)', ' CALCULATING SQUARE-ROOT OF THE REDUCED S-MATRIX >>>'
         print '(/,a,$)', '  solving eigenvalue problem ...'
       end select
    
       call cpu_time(stime)    
c      solve eigenvalue problem
       call realspectrum(tmps,seigval,utrans,msize)
       print *, 'done'

       call cpu_time(ftime)
       mnts = int((ftime-stime)/60.0d0,4)
       secs = ftime - stime - mnts*60.0d0

c      computing square-root ...
       if (testing) then

        tmpfile = 51
        select case (icall)
	 case(1) ; smtfilename = 'smat.tmp'
	 case(2) ; smtfilename = 'omat.tmp'
	 case(3) ; smtfilename = 'omat_molecule.tmp'
	 case(4) ; smtfilename = 'smat.valence.tmp'
	end select

        open(tmpfile,file=trim(smtfilename),status='unknown',iostat=ierr)
        select case (icall)
	 case(1) 
          write(tmpfile,fmt='(a)') '$smat.tmp'
	  write(tmpfile,fmt='(a)') '#eigenvalues of s-matrix'
	 case(2) 
          write(tmpfile,fmt='(a)') '$omat.tmp'
	  write(tmpfile,fmt='(a)') '#eigenvalues of o-matrix'
	 case(3) 
          write(tmpfile,fmt='(a)') '$omat.molecule.tmp'
	  write(tmpfile,fmt='(a)') '#eigenvalues of o-matrix'
	 case(4) 
          write(tmpfile,fmt='(a)') '$smat.valence.tmp'
	  write(tmpfile,fmt='(a)') '#eigenvalues of the reduced s-matrix'
        end select

       end if

       do n = 1, msize
        tmp_eig = seigval(n)
        if (tmp_eig < 0.0d0) then
	 print '(/,a,i5,a)', ' uups! ... eigenvalue', n, 
     & 	                     ' is NEGATIVE ???'
	 stop ': transport module is terminated now'
        elseif (tmp_eig < tiny) then
	 print '(/,a,i5,a)', ' uups! ... eigenvalue', n, 
     &                       ' is extremely small !'
	 stop ': transport module is terminated now'
	else
         if (testing) then
	  write(tmpfile,fmt='(a,i4,a,d12.6)')  
     & 	          '  n = ', n, ' : ', tmp_eig  
	 end if
	 seigval(n) = sqrt(tmp_eig)
	end if
       end do

       if (testing) then 
        write(tmpfile,fmt='(a)') '$end'
        close(tmpfile)
       end if    

       if ((icall.eq.1).or.(icall.eq.4)) then
        print '(2x,a,i3,a,f5.2,a)', 
     &        'time spent:', mnts, ' min ', secs, ' sec'
       end if
       
c      checking for orthogonality of eigenvectors ...
       if (testing) then
        print '(/,a)', '  checking orthogonality of eigenvectors ...'
        allocate(tmp_mat(msize,msize),stat=ierr)
        tmp_mat = 0.0d0
        
c       tmp_mat = utrans * (utrans)^T
        call dble_ab('n','t',msize,utrans,utrans,tmp_mat)     

        do n = 1, msize
         do m = 1, msize
	  if (n==m) then
 	   tmp_diff = tmp_mat(n,m) - 1.0d0
	  else
	   tmp_diff = tmp_mat(n,m)
	  end if 
	  if (abs(tmp_diff) > myprec) then
	   print '(2x,a,i4,a,i4,a)', 
     &	         'element (',n,',',m,') causes a problem !'
           testunit = .false.
	  end if
	 end do
        end do

        if (.not.testunit) then
         stop ': eigenvectors are not orthogonal !'     
        else
	 print *, ' check is OK!'   
        end if
        deallocate(tmp_mat,stat=ierr)
       end if ! testing
            
c      building up square-root of s-matrix and its inverse
c       allocate(smat12(msize,msize),stat=ierr)
c       if (ierr.ne.0) then
c        print *
c        stop
c     &   '[SUBROUTINE sqrt_smat]: <smat12> allocation failure'
c       end if
c
c       allocate(smat12inv(msize,msize),stat=ierr)
c       if (ierr.ne.0) then
c        print *
c        stop
c     &   '[SUBROUTINE sqrt_smat]: <smat12inv> allocation failure'
c       end if

       select case (icall)
        case (1)
         print '(/,2x,a,$)', 'building up square-root of s-matrix ...'
        case (2)
         if (testing) print *
         print '(2x,a,$)', 'building up square-root of o-matrix ...'
        case (3)
         if (testing) print *
         print '(2x,a,$)', 'building up square-root of o-matrix ...'
        case (4)
         print '(/,2x,a,$)', 'building up square-root of the reduced s-matrix ...'
       end select
     
       allocate(tmp_us(msize,msize),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE sqrt_smat]: <tmp_us> allocation failure'
       end if

       allocate(tmp_usinv(msize,msize),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE sqrt_smat]: <tmp_us> allocation failure'
       end if

       call cpu_time(stime)     

c      smt12 = utrans * seigval * (utrans)^T
c      but ...
c      first we build up matrix tmp_us = utrans * seigval
c                        and tmp_usinv = utrans / seigval
       tmp_us = 0.0d0    
       do k = 1, msize
        do n = 1, msize
c        ! be careful:  array 'seigval' contains 
c        ! square-roots from original eigenvalues!
         tmp_us(n,k)    = utrans(n,k)*seigval(k)
         tmp_usinv(n,k) = utrans(n,k)/seigval(k)
	end do
       end do

c      next we compute smt12 = tmp_us * (utrans)^T
       smt12 = 0.0d0
       call dble_ab('n','t',msize,tmp_us,utrans,smt12)     
       smt12inv = 0.0d0
       call dble_ab('n','t',msize,tmp_usinv,utrans,smt12inv)     

       deallocate(tmp_us,tmp_usinv,stat=ierr)
       print *, 'done'

       call cpu_time(ftime)
       mnts = int((ftime-stime)/60.0d0,4)
       secs = ftime - stime - mnts*60.0d0
       if ((icall.eq.1).or.(icall.eq.4)) then
        print '(2x,a,i3,a,f5.2,a)', 
     &        'time spent:', mnts, ' min ', secs, ' sec'
       end if
       
ccc    checking for consitency
       if (testing) then
        select case (icall)
	 case (1)  
	   print '(/,2x,a,$)', 'checking: smat12*smat12 = smat ? ...'    
         case (2)
	   print '(/,2x,a,$)', 'checking: omat12*omat12 = omat ? ...'    
         case (3)
	   print '(/,2x,a,$)', 'checking: omat12*omat12 = omat ? ...'    
	 case (4)  
	   print '(/,2x,a,$)', 'checking: smat12*smat12 = smat ? ...'    
	end select
        tmps = 0.0d0
        tmpe = 0.0d0
	
c	tmps = smt12 * smt12
        call dble_ab('n','n',msize,smt12,smt12,tmps)     
c	tmpe = smat12 * smat12inv
        call dble_ab('n','n',msize,smt12,smt12inv,tmpe)     

        do n = 1, msize
	 do m = 1, msize
c         checking smt12*smt12 = smt	 
	  tmp_diff = smt(n,m) - tmps(n,m)
	  if (abs(tmp_diff) > myprec) then
	   write(90,*) n, m, tmp_diff
           stest = .false.
	  end if
c         checking smt12*smt12inv = 1	 
	  if (n==m) then
	   tmp_diff = tmpe(n,m) - 1.0d0
	  else
	   tmp_diff = tmpe(n,m)
	  end if 
	  if (abs(tmp_diff) > myprec) then
	   write(91,*) n, m, tmp_diff
           etest = .false.
	  end if
	 end do
	end do

        if (.not.stest) then
	 print *
         stop ': square-root of s-matrix is invalid !'     
        else if (.not.etest) then
	 print *
         stop ': inv(smat12) is invalid !'     
        else
	 print *, 'OK!'
	 print '(2x,a)', 'inv(smat12) is OK as well! '   
        end if
       
       end if ! testing

c      storing overlap matrix, if allowed ...
       if (storesmat.and.(icall.eq.1)) then 
        call savesmat(trim(smat12_file_name), trim(smat12inv_file_name))
c       updating <tcontrol> file
        call updatetcontrol(trim(tcntrl_file_name),'$overlapint',0,0.0d0)   
       end if ! storesmat

       deallocate(tmps,tmpe,utrans,seigval,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE sqrt_smat]: ',
     &           'impossible to deallocate temporary arrays'
        print *, 'nevertheless, proceed further ...'
       end if

       select case (icall)
        case(1)   
          print '(/,a)', ' <<< DONE WITH THE SQUARE-ROOT OF S-MATRIX'
        case(2)
          print '(/,a)', ' <HUBBARD-U>: DONE WITH THE SQUARE-ROOT OF O-MATRIX <<<'
        case(3)
          print '(/,a)', ' <MOLECULE''s SUBSPACE>: DONE WITH THE SQUARE-ROOT OF O-MATRIX <<<'
        case(4)   
          print '(/,a)', ' <<< DONE WITH THE SQUARE-ROOT OF THE REDUCED S-MATRIX'
       end select 
      
      end subroutine sqrt_smat
 
      subroutine savesmat(ofile1,ofile2)
c     ***************************************************
c     saves square-root (and its inverse) of overlap  
c     matrix to external files
c     sqrt(smat) --> ofile1 ; inv(sqrt(smat)) --> ofile2
c     ***************************************************
       implicit none
            
       character(*), intent (in) :: ofile1, ofile2
       integer, parameter        :: monitor = 100
       integer extfile1, extfile2, n,k,j, mincol, maxcol, ierr, ierr1
       character(128)    syscall
     
       extfile1 = 29       
       open(extfile1,file=ofile1,status='replace',
     &               action='write',form='unformatted',iostat=ierr)
c      print *, ' ierr = ', ierr
       if (ierr.ne.0) then
        syscall = 'rm '//ofile1 ;  call system(trim(syscall))
        open(extfile1,file=ofile1,status='new',
     &                action='write',form='unformatted',iostat=ierr1)
c       print *, ' ierr1 = ', ierr1
        if (ierr1.ne.0) then
         print '(/,a,a,a,/)',
     &            ' can not open file "', trim(ofile1), '" for writing'
	 stop ': encounters a serious problem here!'
	end if
       end if

       extfile2 = 28
       open(extfile2,file=ofile2,status='replace',
     &               action='write',form='unformatted',iostat=ierr)
c       print *, ' ierr = ', ierr
       if (ierr.ne.0) then
        syscall = 'rm '//ofile2 ;  call system(trim(syscall))
        open(extfile2,file=ofile2,status='new',
     &               action='write',form='unformatted',iostat=ierr1)
c       print *, ' ierr1 = ', ierr1
        if (ierr1.ne.0) then
         print '(/,a,a,a,/)',
     &            ' can not open file "', trim(ofile2), '" for writing'
	 stop ': encounters a serious problem here!'
	end if
       end if

c       write(extfile1,fmt='(a,i5)') '$smat12      nsoas= ',nsaos
c       write(extfile1,fmt='(a)') 
c    &  '#upper triangular of sqrt(overlap):      format(4d22.16))'
       write(extfile1) nsaos

c       write(extfile2,fmt='(a,i5)') '$smat12inv   nsoas= ',nsaos
c       write(extfile2,fmt='(a)') 
c     & '#upper triangular of inv(sqrt(overlap)): format(4d22.16))'
       write(extfile2) nsaos

       print '(/,2x,a,$)', 'saving overlap matrix : | '    
       k = 0    
       do n = 1, nsaos 
c        write(extfile1,fmt='(a,i5)') '#column ', n
c        write(extfile2,fmt='(a,i5)') '#column ', n
        k = k + n
      	if ( mod(k,monitor)==0 ) print '(a,$)', '='  ! monitor
c        only upper triangular part of matrices is stored
 	 maxcol = 0
         do while (maxcol<n)
          mincol = maxcol+1
          maxcol = min(maxcol+4,n)
c          write(extfile1,fmt='(4d22.16)')
c     &                  (smat12(j,n),j=mincol,maxcol)
c          write(extfile2,fmt='(4d22.16)')
c     &                  (smat12inv(j,n),j=mincol,maxcol)
          write(extfile1) (smat12(j,n),j=mincol,maxcol)
          write(extfile2) (smat12inv(j,n),j=mincol,maxcol)
         end do
       end do
       print *, '| done'
c      write(extfile1,fmt='(a)') '$end'
c      write(extfile2,fmt='(a)') '$end'
       close(extfile1) ; close(extfile2)
 
      end subroutine savesmat 


      subroutine get_deltah(dh)
c      ************************************************************** 
c      computes matrix <dh> of exchange field sources: 
c      when added/substracted to/from the spin-dependent hamiltonian, 
c      that helps to impose the required magnetic solution, e.g.,
c      a "domain wall" inside the molecule with d-atoms
c      ************************************************************** 
       double precision, intent (out) :: dh(:,:)

       double precision, allocatable :: dj(:) ! , smat12dj(:,:)
       double precision  tmp_exf
       integer iat, iatype, iorb, norb, nu, ierr !, n, m

c      print '(/,a,$)', '  compute matrix dh ...'
       allocate(dj(nsaos),stat=ierr) 
       if (ierr.ne.0) then
        print *
	stop
     &	 '[SUBROUTINE get_deltah]: <dj> allocation failure'
       end if       

c       allocate(smat12dj(nsaos,nsaos),stat=ierr) 
c       if (ierr.ne.0) then
c        print *
c	stop
c     &	 '[SUBROUTINE get_deltah]: <smat12dj> allocation failure'
c       end if       

c      get diagonal matrix dj <-- source of ex.fields on atoms        
       nu = 0
       do iat = 1, num_atoms
         iatype  = atom(iat)%atype
         tmp_exf = atom(iat)%exfield
         norb = n_basis_func(iatype)
         do iorb = 1, norb
          nu = nu + 1
          dj(nu) = tmp_exf
         end do ! iorb
       end do ! iat

cc      transform dj to orth.basis: dj <-- smat12 * dj * smat12       
c       do m = 1, nsaos
c        do n = 1, nsaos
c	 smat12dj(n,m) = smat12(n,m) * dj(m) 
c	end do
c       end do
c       call dble_ab('n','n',nsaos,smat12dj,smat12,dh)
cc      print *, 'done'     

       dh = 0.0d0
       forall(nu=1:nsaos) dh(nu,nu) = dj(nu)
 
c      deallocate(dj,smat12dj,stat=ierr)
       deallocate(dj,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE get_deltah]: ',
     &           'impossible to deallocate temp. arrays'
        print *, 'nevertheless, proceed further ...'
       end if
     
      end subroutine get_deltah

      subroutine h0
c      *******************************************
c      computes initial (h0) hamiltonian of the
c      extended molecule; no self-energy included
c      *******************************************
       implicit none
    
c      molecular orbitals in orthogonal basis for given ispin
       double precision, allocatable :: mo_ort(:,:)         
c      temporary arrays ... 
       double precision, allocatable :: tmp_ce(:,:),
     &                                  tmp_h0(:,:), tmp_en(:),
     &                                  utrans(:,:),
     &                                  mo_in(:), mo_out(:),
     &                                  dh(:,:)
       integer, allocatable :: tmpindex(:)
   
       integer ierr, ispin, n, nu, mu, nu1, p
       double precision tmp_diff, dev
       logical :: htest = .true.

       allocate(mo_ort(nsaos,nsaos),stat=ierr) 
       if (ierr.ne.0) then
        print *
	stop
     &	 '[SUBROUTINE h0]: <mo_ort> allocation failure'
       end if       
       
       if (do_pop.or.do_lmpop) then
        allocate(mo_orth(nsaos,nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop
     &    '[SUBROUTINE h0]: <mo_orth> allocation failure'
        end if
       end if

       allocate(h0mat(nsaos,nsaos,nspin),stat=ierr) 
       if (ierr.ne.0) then
        print *
	stop
     &	 '[SUBROUTINE h0]: <h0mat> allocation failure'
       end if       

       allocate(tmp_ce(nsaos,nsaos),stat=ierr) 
       if (ierr.ne.0) then
        print *
	stop
     &	 '[SUBROUTINE h0]: <tmp_ce> allocation failure'
       end if       

       allocate(tmp_h0(nsaos,nsaos),stat=ierr) 
       if (ierr.ne.0) then
        print *
	stop
     &	 '[SUBROUTINE h0]: <tmp_h0> allocation failure'
       end if       

c      if (symmetrize) then
       if ( (nspin.eq.2).and.(symmetrize)) then
         call update_mo
       end if

ccccccccccccccccccccccccccccccccccccccccccccccccc
c      !!!  do not uncommnet next block: 
c      !!!  it does not functon like it should :(
c
c       if (z_symmetrize) then
c        call symmetrize_mo('z')
c       end if
c
c       if (x_symmetrize) then
c        call symmetrize_mo('x')
c       end if
c
c       if (y_symmetrize) then
c        call symmetrize_mo('y')
c       end if
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

       print '(/,a,$)', ' COMPUTING HAMILTONIAN OF THE EXTENDED MOLECULE ...'

       call cpu_time(stime)     
       h0mat = 0.0d0 
       do ispin = 1, nspin

c       ortogonalize mo's
c       mo_ort = smat12 * mo_coeff
        mo_ort = 0.0d0
        call dble_ab('n','n',nsaos,smat12,mo_coeff(:,:,ispin),mo_ort)     

c       save a copy of <mo_ort> to the global array <mo_orth>
        if (do_pop.or.do_lmpop) then
         forall (n=1:nsaos,p=1:nsaos)
     &    mo_orth(n,p,ispin) = mo_ort(n,p)
        end if

c       h0 = mo_ort * mo_en * (mo_ort)^T
c       but ...
c       first, we build up matrix tmp_ce = mo_ort * mo_en 
        tmp_ce = 0.0d0
        do n = 1, nsaos
	 do nu = 1, nsaos
	  tmp_ce(nu,n) = mo_ort(nu,n) * mo_en(n,ispin)
	 end do
	end do
	
c       next, we compute h0 for given ispin
        tmp_h0 = 0.0d0
        call dble_ab('n','t',nsaos,tmp_ce,mo_ort,tmp_h0)

c       finally, we save result in matrix h0mat
        do nu = 1, nsaos
	 do nu1 = nu, nsaos
c	 do nu1 = 1, nsaos
c         h0mat(nu,nu1,ispin) = tmp_h0(nu,nu1)         
c         optional >>> 
c         symmetrization to avoid numerical uncertainty
          h0mat(nu,nu1,ispin) = 0.5d0*(tmp_h0(nu,nu1)+tmp_h0(nu1,nu))         
          h0mat(nu1,nu,ispin) = h0mat(nu,nu1,ispin)         
	 end do
	end do 

       end do ! ispin

       deallocate(tmp_ce,mo_ort,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE h0]: ',
     &           'impossible to deallocate temp. arrays'
        print *, 'nevertheless, proceed further ...'
       end if
       call cpu_time(ftime)     
       mnts = int((ftime-stime)/60.0d0,4)
       secs = ftime - stime - mnts*60.0d0

       print *, 'DONE'
       print '(1x,a,i3,a,f5.2,a)', 
     &       'time spent:', mnts, ' min ', secs, ' sec'

c       save kohn-sham hamiltonian to external file 
        if (save_hmat) call savehmat(outhmat_file_name)

c      internal check: 
c      compute eigenvalues of H0 and compare them with MO energies
       if (testing) then
        print '(/,a)', 
     &        ' internal consistency: checking eigenvalues of H0'

c       u-matrix with eigenvectors (in columns)
        allocate(utrans(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
          print *
 	  stop '[SUBROUTINE h0]: <utrans> allocation failure'
        end if       

c       eigenvalues of H0
        allocate(tmp_en(nsaos),stat=ierr) 
        if (ierr.ne.0) then
          print *
	  stop '[SUBROUTINE h0]: <tmp_en> allocation failure'
        end if        

        allocate(mo_in(nsaos),stat=ierr) 
        if (ierr.ne.0) then
          print *
	  stop '[SUBROUTINE h0]: <mo_in> allocation failure'
        end if       

        allocate(mo_out(nsaos),stat=ierr) 
        if (ierr.ne.0) then
          print *
	  stop '[SUBROUTINE h0]: <mo_out> allocation failure'
        end if       

        allocate(tmpindex(nsaos),stat=ierr) 
        if (ierr.ne.0) then
          print *
	  stop '[SUBROUTINE h0]: <tmpindex> allocation failure'
        end if       

        do ispin = 1, nspin

         print '(a,i1)', '  - ispin = ', ispin
         utrans = 0.0d0

c        make local copies of mos and mo-energies
         do nu = 1, nsaos 
          mo_in(nu) = mo_en(nu,ispin)
	  do nu1 = 1, nsaos
           tmp_h0(nu,nu1) = h0mat(nu,nu1,ispin)         
          end do	 
	 end do
c        ordering mo-energies
         mo_out = 0.0d0
	 tmpindex = 0
         call sortarray(mo_in,mo_out,nsaos,tmpindex)  
   
         print '(4x,a,$)', 'computing eigenvalues ...'
         call realspectrum(tmp_h0,tmp_en,utrans,nsaos)
         print *, 'done'

	 dev = 0.0d0
	 do n = 1, nsaos
	  tmp_diff = dabs(mo_out(n) - tmp_en(n))
	  if (tmp_diff > dev) dev = tmp_diff
	  if (tmp_diff > mysmall) then
	   print '(4x,a,i4,a)', 
     &	         'eigenvalue ', n,' causes a problem !'
           print '(2(D20.10))', mo_en(n,ispin), tmp_en(n)
           htest = .false.
	  end if 
 	 end do ! nsaos

	end do ! ispin

        if (.not.htest) then
         print '(a,E12.6,a)', ' max. deviation = ', dev, ' H'
         stop ': eigenvalues of H0 do not match MO energies!'    
        else
	 print '(1x,a)', 'check is OK!'   
         print '(a,E12.6,a)', ' max. deviation = ', dev, ' H'
        end if
   
        deallocate(utrans,tmp_en,mo_in,mo_out,tmpindex,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE h0]: ',
     &            'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
        end if
	
       end if ! testing

c      optional: activate local exchange fields
       if (localexf) then
        write(*,'(/,a)') 
     &	     ' >>> CAUTION: local exchange fields are active! <<<'

        allocate(dh(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
          print *
	  stop '[SUBROUTINE h0]: <dh> allocation failure'
        end if       

        call get_deltah(dh)
        do mu = 1, nsaos
         do nu = 1, nsaos
	  h0mat(nu,mu,1) = h0mat(nu,mu,1) - 0.5d0*dh(nu,mu)
	  h0mat(nu,mu,2) = h0mat(nu,mu,2) + 0.5d0*dh(nu,mu)
         end do ! nu
        end do ! mu

        deallocate(dh,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE h0]: ',
     &            'impossible to deallocate <dh>'
         print *, 'nevertheless, proceed further ...'
        end if
       end if ! localexf

       deallocate(tmp_h0,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', ' [SUBROUTINE h0]: ',
     &           'impossible to deallocate <tmp_h0>'
        print *, 'nevertheless, proceed further ...'
       end if

      end subroutine h0

      subroutine savehmat(hfile)
c     *****************************************      
c     saves hamiltonian matrix to external file
c     *****************************************      
       character(*), intent (in) :: hfile
       integer, parameter        :: monitor = 100
       integer extfile, n,k,j, mincol, maxcol, ispin, ierr, ierr1
       character(128) syscall

       extfile = 39
       open(extfile,file=hfile,status='replace',action='write',iostat=ierr)
c      print *, ' ierr = ', ierr
       if (ierr.ne.0) then
        syscall = 'rm '//hfile ;  call system(trim(syscall))
        open(extfile,file=hfile,status='new',action='write',iostat=ierr1)
c       print *, ' ierr1 = ', ierr1
        if (ierr1.ne.0) then
         print '(/,a,a,a,/)',
     &            ' can not open file "', trim(hfile), '" for writing'
         stop ': encounters a serious problem here!'
        end if
       end if

       write(extfile,fmt='(a,i5)') '$hamiltonian   format(4d20.14)   nsaos= ',nsaos
       write(extfile,fmt='(a,i5)') '#nspins:', nspin
c       write(extfile1,fmt='(a)')
c    &  '#upper triangular of sqrt(overlap):      format(4d22.16))'

       print '(/,1x,a,$)', 'saving hamiltonian matrix : | '
       do ispin = 1, nspin
        write(extfile,fmt='(a,i5)') '#ispin: ', ispin
        k = 0
        do n = 1, nsaos
         write(extfile,fmt='(a,i5)') '#column ', n
         k = k + n
         if ( mod(k,monitor)==0 ) print '(a,$)', '='  ! monitor
c         only upper triangular part of matrices is stored
          maxcol = 0
          do while (maxcol<n)
           mincol = maxcol+1
           maxcol = min(maxcol+4,n)
           write(extfile,fmt='(4d20.14)') (h0mat(j,n,ispin),j=mincol,maxcol)
          end do
        end do
       end do ! ispin
       print *, '| done'
       write(extfile,fmt='(a)') '$end'
       close(extfile)

      end subroutine savehmat
     
      subroutine readhmat(hmatfile,exthmat,msize)
c     *********************************************
c     reads a hamiltonian matrix from external file
c     *********************************************
      
        character(*), intent (in)     :: hmatfile
        double precision, intent(out) :: exthmat(:,:,:)
        integer, intent(in)           :: msize
		      
        integer, parameter :: monitor = 100
        integer extfile, n,k,j, mincol, ispin, tmp_nspin, tmp_ispin, stmp, 
     &    	tmp_msize, maxcol, ierr
        character(64) tmpstring, tmpstring1, tmpstring2
					            
        extfile = 28  ! external omat-file
        open(extfile,file=hmatfile,status='old',action='read',iostat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a,a)',
     &            ' can not open file "', trim(hmatfile), '" for reading'
         print *,
     &    'please, check your directory content or access rights'
         print *
         stop ': transport module is terminated now'
        end if

c        if (.not.allocated(exthmat)) then
c         allocate(exthmat(nsaos,nsaos),stat=ierr)
c         if (ierr.ne.0) then
c          stop '[SUBROUTINE readomat]: <smat> allocation failure'
c         end if
c        end if

        print '(/,2x,a,$)', 'importing molecular hamiltonian : | '
        read(extfile,*) tmpstring, tmpstring1, tmpstring2, tmp_msize
	if (msize.ne.tmp_msize) then 
         print '(/,2x,a)', 'confusion due to dimension mismatch!'
	 print '(2x,a,i4,a,i4,/)', 'msize= ', msize, 
     &         ' does not fit dimenstion of external hamiltonian tmp_msize =', tmp_msize
	 stop ': transport module will be terminated now!'
        end if

        read(extfile,*) tmpstring, tmp_nspin
	if (nspin.ne.tmp_nspin) then
         print '(/,2x,a)', 'amount of spin degrees of freedom do not match each other!'
	 print '(2x,a,i4,a,i4,/)', 'nspin = ', nspin, 
     &         ' does not equal tmp_nspin =',  tmp_nspin
	 stop ': transport module will be terminated now!'
        end if
		
        k = 0
	do ispin = 1, nspin
         read(extfile,*) tmpstring, tmp_ispin
c        print '(/,a,i2,/)', 'ispin = ', tmp_ispin
	 do n = 1, msize
          read(extfile,*) tmpstring, stmp
c         print '(a,i5)', trim(tmpstring), stmp
          k = k + n
          if ( mod(k,monitor)==0 ) print '(a,$)', '='  ! monitor
c          only upper triangular part is saved to file
           maxcol = 0
           do while (maxcol<n)
            mincol = maxcol+1
            maxcol = min(maxcol+4,n)
            read(extfile,fmt='(4d20.14)') (exthmat(j,n,ispin),j=mincol,maxcol)
c           take care about transposed elements :
            forall(j = mincol:maxcol) exthmat(n,j,ispin) = exthmat(j,n,ispin)
           end do
         end do
        end  do ! ispin
	
	print *, '| done'
        close(extfile)

      end  subroutine readhmat

      subroutine savemos(en_in,umat_in,jspin)
c     ***************************************************************
c     saves updated (symmetrized) molecular orbitals to external file
c     *************************************************************** 
       double precision, intent(in) :: en_in(:), umat_in(:,:)
       integer, intent(in) ::          jspin

       integer, parameter  :: monitor = 100
       integer extfile, extfile1, n,k,j, mincol, maxcol, ierr, ierr1
       character(128) syscall

       character(72) tmpmosfile, mosfile, outformat
       character(80) headerline 

       if (nspin.eq.1) then
         mosfile = trim(mos_file_name)
         tmpmosfile = trim(mos_file_name)//'.tmp'
       else  ! nspin == 2
         if (jspin.eq.1) then
          mosfile = trim(alpha_file_name)
          tmpmosfile = trim(alpha_file_name)//'.tmp'
         else
          mosfile = trim(beta_file_name)
          tmpmosfile = trim(beta_file_name)//'.tmp'
         end if       
       end if

       extfile = 64
       open(extfile,file=trim(tmpmosfile),status='replace',action='write',iostat=ierr)
c      print *, ' ierr = ', ierr
       if (ierr.ne.0) then
        syscall = 'rm '//trim(tmpmosfile) ;  call system(trim(syscall))
        open(extfile,file=trim(tmpmosfile),status='new',action='write',iostat=ierr1)
c       print *, ' ierr1 = ', ierr1
        if (ierr1.ne.0) then
         print '(/,a,a,a,/)',
     &            ' can not open file "', trim(tmpmosfile), '" for writing'
         stop ': encounters a serious problem here!'
        end if
       end if

       extfile1 = 65  ! external mos-file
       open(extfile1,file=trim(mosfile),status='old',action='read',iostat=ierr)
       if (ierr.ne.0) then
          print '(/,a,a,a)', ' can not open file "', trim(mosfile), '" for reading'
          print *, 'please, check your directory content or access rights'
          print *
          stop ': transport module is terminated now'
       end if

c      reading mos-file header >>>
       read(extfile1,fmt='(a80)') headerline
       write(extfile,fmt='(a80)') headerline

       if (nsaos<10) then
        outformat = '(i6,2x,a,6x,a,d20.14,3x,a,i1)'
       else if (nsaos<100) then
        outformat = '(i6,2x,a,6x,a,d20.14,3x,a,i2)'
       else if (nsaos<1000) then
        outformat = '(i6,2x,a,6x,a,d20.14,3x,a,i3)'
       else if (nsaos<10000) then
        outformat = '(i6,2x,a,6x,a,d20.14,3x,a,i4)'
       else if (nsaos<100000) then
        outformat = '(i6,2x,a,6x,a,d20.14,3x,a,i5)'
       else 
        outformat = '(i6,2x,a,6x,a,d20.14,3x,a,i6)'
       end if

       print '(2x,a,$)', 'updating mos/alpha/beta file(s) : | '
       k = 0
       do n = 1, nsaos
         k = k + n
         if ( mod(k,monitor)==0 ) print '(a,$)', '='  ! monitor
!        output a KS orbital energy >>>
         write(extfile,fmt=trim(outformat)) 
     &         n, 'a', 'eigenvalue=', en_in(n), 'nsaos=', nsaos
!        output an eigenvector >>>
         maxcol = 0
         do while (maxcol<nsaos)
           mincol = maxcol+1
           maxcol = min(maxcol+4,nsaos)
           write(extfile,fmt='(4d20.14)') (umat_in(j,n),j=mincol,maxcol)
         end do
       end do ! nsaos

       print *, '| done'
       print * 

       if (jspin.eq.2) then
        write(extfile,fmt='(a)') '$end'
        close(extfile)
       else 
        if (.not.mo_restart) then
         write(extfile,fmt='(a)') '$end'
         close(extfile)
        else 
         ! mos or alpha AND restart
         close(extfile)
         syscall = 'tail -10 '//trim(mosfile)//' >> '//trim(tmpmosfile)
         call system(trim(syscall))
        end if
       end if

       close(extfile1)

      end subroutine savemos

      end module hamiltonian       
  
      