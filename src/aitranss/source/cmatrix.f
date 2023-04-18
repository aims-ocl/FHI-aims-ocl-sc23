c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           nov 2009
c     last revision:  jan 2012
c###########################################################

      module cmatrix
c     ***************************************      
c     set of routines to solve the eigenvalue
c     problem for complex-valued matrices 
c     ***************************************      

       use globalvars
       use tools
       implicit none

       double precision, private, parameter :: myprec = 1.0d-9  
       real(4), private ::  stime, ftime, secs
       integer, private ::  mnts  
 
      contains

      subroutine cmplxspectrum(cmat,cpoles,bmat,invbmat,msize,icall)
c      **********************************************************
c      takes a complex-valued matrix <cmat> of dimension "msize", 
c      computes its spectrum <cpoles>, its eigenvectors <bmat>, 
c      and the inverse eigenvector matrix <invbmat>
c      **********************************************************

       implicit none
c      input/ouput        
       complex(8), intent (in) :: cmat(:,:,:)
       complex(8), intent (out):: cpoles(:,:)
       complex(8), intent (out):: bmat(:,:,:), invbmat(:,:,:)
       integer,    intent (in) :: msize, icall
       
c      temporary arrays for lapack subroutines ...      
       complex(8), allocatable :: vl(:,:), vr(:,:), tmp_work(:),
     &                            tmp_cmat(:,:),zenergy(:),szenergy(:) 
       double precision, allocatable :: rwork(:)
       
       integer ispin, lwork, n, m, p, ptmp, info, ierr, tmpfile
       integer, allocatable :: indarray(:)

       logical  ostat 
       character(50) :: zmos_file, zalpha_file, zbeta_file

       allocate(tmp_cmat(msize,msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE cmplxspectrum]: <tmp_cmat> allocation failure'
       end if

       allocate(zenergy(msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE cmplxspectrum]: <zenergy> allocation failure'
       end if

       allocate(szenergy(msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE cmplxspectrum]: <szenergy> allocation failure'
       end if

       allocate(indarray(msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE cmplxspectrum]: <indarray> allocation failure'
       end if

       allocate(vl(msize,msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE cmplxspectrum]: <vl> allocation failure'
       end if
       forall(n=1:msize,m=1:msize) vl(n,m) = czero

       allocate(vr(msize,msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE cmplxspectrum]: <vr> allocation failure'
       end if
       forall(n=1:msize,m=1:msize) vl(n,m) = czero
           
c      size of tmp_work (work space) comes out from a query call of 'zgeev'
       lwork = 33*msize
       allocate(tmp_work(lwork),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE cmplxspectrum]: <tmp_work> allocation failure'
       end if

       allocate(rwork(2*msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE cmplxspectrum]: <rwork> allocation failure'
       end if

c       print '(/,a)', 
c     &       ' COMPUTING SPECTRUM OF EXTENDED HAMILTONIAN >>>'
       
c      cycle over spin channels
       do ispin = 1, nspin

        if (nspin == 2) then
	 if (ispin == 1) then 
	  print '(/,2x,a)', '-- alpha CHANNEL --' 
         else 
	  print '(/,2x,a)', '-- beta CHANNEL --'
	 end if	
	end if
           
c       save a local copy of cmat
        forall(n=1:msize,m=1:msize)
         tmp_cmat(n,m) = cmat(n,m,ispin)     
        end forall 

        print '(/,2x,a,$)', 'solving eigenvalue problem ...'
        call cpu_time(stime)
        	
	tmp_work = czero
c       searching for eigevalues <zenergy> & eigenvectors <vr> 
c       attention: left eigenvectors <vl> are not referenced !               

        call zgeev('N','V',msize,tmp_cmat,msize,zenergy,
     &  	   vl,msize,vr,msize,tmp_work,lwork,rwork,info)

c       query call >>>
c        call zgeev('N','V',msize,tmp_cmat,msize,zenergy,
c     &  	   vl,msize,vr,msize,tmp_work,-1,rwork,info)
c
c         print *
c         print *, ' optimal lwork = ', tmp_work(1)
c	stop

        call cpu_time(ftime)
        if (info.eq.0) then
         print *, 'done'
        else
         print *, 'failed' 
	 print '(2x,a,i3,/)', 'error info =', info
         stop ': transport module is terminated now'
        end if 
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts*60.0d0
        print '(2x,a,i3,a,f5.2,a)',
     &        'time spent:', mnts, ' min ', secs, ' sec'  
        
c       final step: sorting eigen-energies in ascending order
c                   reshuffling corresponding eigenvectors  

        call cmplxsort(zenergy,szenergy,msize,indarray)

c       save complex eigenvalues/eigenvectors to global arrays
        forall (n=1:msize,p=1:msize) 
	 bmat(n,p,ispin) = czero
	end forall       

        do p=1, msize
	 cpoles(p,ispin) = szenergy(p)
         ptmp = indarray(p)
         do n = 1, msize
          bmat(n,p,ispin) = vr(n,ptmp)
	 end do
	end do
     
        if (testing.or.zspectrum) then
         
	 select case (icall)
	  case (1)
c          full hamiltonian	  
	   zmos_file   = 'zmos.tmp'
	   zalpha_file = 'zalpha.tmp'
	   zbeta_file  = 'zbeta.tmp'
	  case (2)
c          reservoirs
	   zmos_file   = 'zmos.res.tmp'
	   zalpha_file = 'zalpha.res.tmp'
	   zbeta_file  = 'zbeta.res.tmp'
	  case (3)
c          molecule
	   zmos_file   = 'zmos.mol.tmp'
	   zalpha_file = 'zalpha.mol.tmp'
           zbeta_file  = 'zbeta.mol.tmp'
	  case default ;
	 end select
 
         tmpfile = 32
         if (nspin == 1) then
	  open(tmpfile,file=trim(zmos_file),status='unknown',iostat=ierr)
          if (ierr.ne.0) then
	   stop ': can not open zmos.tmp!'
          end if
	  write(tmpfile,'(a)') '$'//trim(zmos_file)
	 else  ! spin-polarized case
          if (ispin == 1) then	 
	   open(tmpfile,file=trim(zalpha_file),status='unknown',iostat=ierr)
           if (ierr.ne.0) then
	    stop ': can not open zalpha.tmp!'
           end if
	   write(tmpfile,'(a)') '$'//trim(zalpha_file)
	  else
	   open(tmpfile,file=trim(zbeta_file),status='unknown',iostat=ierr)
           if (ierr.ne.0) then
	    stop ': can not open zbeta.tmp!'
           end if
	   write(tmpfile,'(a)') '$'//trim(zbeta_file)

	  end if 
	 end if
	 write(tmpfile,'(a)') '#poles of the Greens function'	 
	 if (icall.eq.1) then	 
	  write(tmpfile,'(a,i5)') '#nsaos=', msize
         else
	  write(tmpfile,'(a,i5)') '#msize=', msize
         end if
      
         do n = 1, msize
           write(tmpfile,'(i6,6x,d20.14,6x,d20.14)') n,cpoles(n,ispin)
c          write(tmpfile,'(i6,6x,d14.8,6x,d14.8)') n,cpoles(n,ispin)
         end do
 	 write(tmpfile,'(a)') '$end'
	 close(tmpfile)

        end if 
c       ! testing .or. zspectrum

c       checking eigen-decomposition
        if (testing) then
	 call eigdec_check(cmat,bmat,cpoles,msize,ispin,info)
	 if (info.eq.1) then ! check failed
	  print *
	  stop ': further calculations may go wrong! will quit now.'
         end if 
	end if

c       inverting eigenvector matrix 
c       (needed further for the density matrix evaluation)
        call eigvecinv(bmat,invbmat,ispin,msize,ostat)

        if (.not.ostat) then
	 print *
         stop 
     &   ': was not able to compute the inverse of eigenvectors matrix'
        end if
     
       end do ! ispin
            
c       print '(/,a)', ' <<< DONE WITH SPECTRUM OF HAMILTONIAN'
       
       deallocate(tmp_cmat,vl,vr,tmp_work,rwork,
     &            zenergy,szenergy,indarray,stat=ierr)
       if (ierr.ne.0) then
       print '(/,a,/,a)',
     &  ' [SUBROUTINE cmplxspectrum]: deallocation of temp. arrays failed',
     &  ' ... will proceed further anyway'
     
       end if
               
      end subroutine cmplxspectrum


      subroutine eigvecinv(bmat,invbmat,ispin,msize,ostat)
c      ************************************************
c      computes the inverse of complex-valued matrix
c      <bmat> of size "msize" for given "ispin",
c      puts answer into <invbmat> 
c      ************************************************
       implicit none
c      input/output       
       complex(8), intent(in)  :: bmat(:,:,:)
       complex(8), intent(out) :: invbmat(:,:,:)
       integer, intent(in)     :: ispin, msize
       logical, intent(out)    :: ostat ! output status
       
       complex(8), allocatable :: tmpb(:,:), tmpbinv(:,:), tmp_mat(:,:)
       integer,    allocatable :: ipiv(:)
       complex(8), allocatable :: tmpwork(:)
       
       integer     n, m, lwork, ierr, info
       complex(8)  cdiff
       logical     testunit
       
       allocate(tmpb(msize,msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE eigvecinv]: <tmpb> allocation failure'
       end if

       allocate(tmpbinv(msize,msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE eigvecinv]: <tmpb> allocation failure'
       end if
    
       allocate(tmp_mat(msize,msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE eigvecinv]: <tmpb> allocation failure'
       end if
    
       allocate(ipiv(msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE eigvecinv]: <ipiv> allocation failure'
       end if

c      that comes out from query call of [zgetri]
       lwork = 64 * msize
       allocate(tmpwork(lwork),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE eigvecinv]: <tmpwork> allocation failure'
       end if

       print '(/,2x,a,$)', 'inverting eigenvalues matrix (B) ...'

c      save local copy of b-matrix 
       forall (n=1:msize,m=1:msize) 
	 tmpb(n,m) = bmat(n,m,ispin)
       end forall

c      inverting matrix
c      1st step: LU decomposition
       tmpbinv = tmpb
       call cpu_time(stime)
       call zgetrf(msize,msize,tmpbinv,msize,ipiv,info)
       if (info.ne.0) then
         print '(/,2x,a,a)', '[SUBROUTINE eigvecinv]: ',
     &                       'problems with LU decomposition'
         print '(2x,a,i4)',  'info =', info
         ostat = .false.
        return
       end if

c      2nd step: inversion :: tmpb --> inv(tmpb)
       call zgetri(msize,tmpbinv,msize,ipiv,tmpwork,lwork,info)
       call cpu_time(ftime)
c      query call
c      call zgetri(msize,tmpb,msize,ipiv,tmpwork,-1,info)
c      if (info.eq.0) then
c          print *
c          print *, ' optimal lwork = ', tmpwork(1)
c      end if
c      stop
       if (info.ne.0) then
         print '(/,2x,a,a)', '[SUBROUTINE eigvecinv]: ',
     &                       'calling <zgetri> fails'
         ostat = .false.
         return
       end if
       
       ostat = .true.
       print *, 'done'

       mnts = int((ftime-stime)/60.0d0,4)
       secs = ftime - stime - mnts*60.0d0
       print '(2x,a,i3,a,f5.2,a)',
     &        'time spent:', mnts, ' min ', secs, ' sec'  

       if (testing) then
        print '(/,2x,a,$)', '-- checking inv(B) ...'
        testunit = .true.       
c       tmp_mat = tmpb * inv(tmpb)
        tmp_mat = czero
        call cmplx_ab('n','n',msize,tmpb,tmpbinv,tmp_mat)
	
        do n = 1, msize
         do m = 1, msize
          if (n==m) then
           cdiff = tmp_mat(n,m) - cone
	  else
           cdiff = tmp_mat(n,m)
          end if
	  if (cmplxabs(cdiff) > myprec) then
c           print '(2x,a,i4,a,i4,a)',
c     &           'element (',n,',',m,') causes a problem !'
            write(90-ispin,*) n, m, cdiff, cmplxabs(cdiff)
           testunit = .false.
          end if
         end do
        end do
        
        if (.not.testunit) then
         print '(1x,a,/)', 'inv(B) is invalid !'
         stop ': further results may be wrong!'
        else
         print '(1x,a)', 'OK!'
        end if

       end if

c      save inv(tmpb) to global array      
       forall(n=1:msize,m=1:msize)
        invbmat(n,m,ispin) = tmpbinv(n,m) 
       end forall

       deallocate(tmpb,tmpbinv,tmpwork,ipiv,tmp_mat,stat=ierr)
       if (ierr.ne.0) then
        print '(/,a,a)', '  [SUBROUTINE eigvecinv]: ',
     &             ' impossible to deallocate temporary arrays'
        print *, ' anyway, proceed further ...'
       end if

      end subroutine eigvecinv

      subroutine unitcheck(qmat,msize,info)
c      *****************************************  
c      performs internal check: qmat * qmat+ = E
c      *****************************************  
       implicit none
        
       complex(8), intent(in) :: qmat(:,:)
       integer, intent(in)    :: msize   
       integer, intent(out)   :: info   
c                 ! exit status  = 0: OK
c                                = 1: Q-matrix fails        
           
       complex(8), allocatable :: tmp_mat(:,:)
       complex(8)  cdiff

       integer ierr, n, m
       logical testunit 
              
       print '(/,2x,a)', 
     &       '-- checking unitary of eigenvectors matrix ...'
       
       allocate(tmp_mat(msize,msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE eigdec_check]: <tmp_mat> allocation failure'
       end if

       testunit = .true.       
c      tmp_mat = qmat * qmat+
       tmp_mat = czero
       call cmplx_ab('n','c',msize,qmat,qmat,tmp_mat)
	
       do n = 1, msize
         do m = 1, msize
          if (n==m) then
           cdiff = tmp_mat(n,m) - cone
	  else
           cdiff = tmp_mat(n,m)
          end if
	  if (cmplxabs(cdiff) > myprec) then
c           print '(2x,a,i4,a,i4,a)',
c     &           'element (',n,',',m,') causes a problem !'
            write(98,*) n, m, cdiff, cmplxabs(cdiff)
           testunit = .false.
          end if
         end do
       end do
        
       if (.not.testunit) then
         print '(5x,a)', 'eigenvectors are not orthogonal !'
         info = 1
       else
         print '(5x,a)', 'check is OK!'
         info = 0
       end if

       deallocate(tmp_mat,stat=ierr)
       
      end subroutine unitcheck
      
      subroutine eigdec_check(amat,qmat,eigval,msize,ispin,info)
c      **********************************************
c      perform internal check :  A*Q = Q*diag(eigval)
c      **********************************************
       implicit none
        
       complex(8), intent(in) :: amat(:,:,:), qmat(:,:,:),
     &                           eigval(:,:)
       integer, intent(in)    :: msize, ispin
       integer, intent(out)   :: info   
c                 ! exit status  = 0: OK
c                                = 1: check fails
           
       complex(8), allocatable :: tmpaq(:,:), qeig(:,:)
       complex(8)  cdiff

       integer ierr, n, m
       logical testeig
       
       info = 0 

       allocate(tmpaq(msize,msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE eigdec_check]: <tmpaq> allocation failure'
       end if

       allocate(qeig(msize,msize),stat=ierr)
       if (ierr.ne.0) then
       print *
       stop
     &   '[SUBROUTINE eigdec_check]: <qeig> allocation failure'
       end if

       print '(/,2x,a)', 
     &       '-- checking eigenvectors & eigenvalues ...'
       
       testeig = .true.       

c      aq = amat * qmat
       call cmplx_ab('n','n',msize,amat(:,:,ispin),qmat(:,:,ispin),tmpaq)
	
c      qeig = qmat * diag(eigval)
       do m = 1, msize
	 do n = 1, msize
          qeig(n,m) = qmat(n,m,ispin) * eigval(m,ispin)	 
	 end do
       end do
     
       do n = 1, msize
         do m = 1, msize
          cdiff = tmpaq(n,m) - qeig(n,m)
	  if (cmplxabs(cdiff) > myprec) then
c           print '(2x,a,i4,a,i4,a)',
c     &           'element (',n,',',m,') causes a problem !'
            write(98-ispin,*) n, m, cdiff, cmplxabs(cdiff)
            testeig = .false.
          end if
         end do
       end do
        
       if (.not.testeig) then
         print '(5x,a)', 'eigen-decomposition is invalid !'
	 info = 1
       else
         print '(5x,a)', 'check is OK!'
	 if (info.ne.1) info = 0
       end if

       deallocate(tmpaq,qeig,stat=ierr)
       
      end subroutine eigdec_check

      end module cmatrix       

