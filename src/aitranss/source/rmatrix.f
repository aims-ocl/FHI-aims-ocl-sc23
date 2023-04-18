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

      module rmatrix
c     ****************************************************
c     eigenvalue solver for real-valued symmetric matrices 
c     ****************************************************
      implicit none

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

      contains

      subroutine realspectrum(amat,rpoles,umat,msize)
c      ****************************************************************
c      calls the most advanced routine 'dsyevr' from LAPACK
c      to solve eigenvalue problem for the real-valued symmetric matrix
c      ****************************************************************
       double precision, intent(in)  :: amat(:,:)
       double precision, intent(out) :: rpoles(:), umat(:,:)
       integer, intent(in) :: msize
       
       integer info, ierr
       
       lwork = 33*msize
       allocate(tmp_work(lwork),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE realspectrum]: <tmp_work> allocation failure'
       end if

       liwork = 10*msize
       allocate(tmp_iwork(liwork),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE realspectrum]: <tmp_iwork> allocation failure'
       end if

       allocate(isupp(2*msize),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE realspectrum]: <isupp> allocation failure'
       end if

c      solving eigenvalue problem
       call dsyevr(task,range,uplo,msize,amat,msize,vl,vu,il,iu,
     &             abstol,numeigval,rpoles,umat,msize,isupp,
     &             tmp_work,lwork,tmp_iwork,liwork,info)

c       call dsyevr(task,range,uplo,msize,amat,msize,vl,vu,il,iu,
c     &             abstol,numeigval,rpoles,umat,msize,isupp,
c     &             tmp_work,-1,tmp_iwork,-1,info)
c       print *
c       print *, 'optimal lwork  = ', tmp_work(1)
c       print *, 'optimal liwork = ', tmp_iwork(1)
c       stop

       if (info.eq.0) then ;
c       print *, 'done'
       else
          print '(/,a)', 'CAUTION: unsuccessful termination of DSYEVR !'
          print '(/,i3,/)', ' -- info =', info
          stop ': transport module is terminated now'
       end if

       deallocate(tmp_work,tmp_iwork,isupp,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE realspectrum]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
       end if

      end subroutine realspectrum	

      subroutine real_spectrum(amat,rpoles,umat,msize)
c      ****************************************************************
c      calls the most advanced routine 'dsyevr' from LAPACK
c      to solve eigenvalue problem for the real-valued symmetric matrix
c      ****************************************************************
       double precision  amat(:,:)
       double precision, intent(out) :: rpoles(:), umat(:,:)
       integer, intent(in) :: msize

       double precision, allocatable :: tmp_amat(:,:)
       integer info, ierr

       allocate(tmp_amat(msize,msize),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE real_spectrum]: <tmp_amat> allocation failure'
       end if
       tmp_amat = amat
       
       lwork = 33*msize
       allocate(tmp_work(lwork),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE real_spectrum]: <tmp_work> allocation failure'
       end if

       liwork = 10*msize
       allocate(tmp_iwork(liwork),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE real_spectrum]: <tmp_iwork> allocation failure'
       end if

       allocate(isupp(2*msize),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE real_spectrum]: <isupp> allocation failure'
       end if

c      solving eigenvalue problem
       call dsyevr(task,range,uplo,msize,tmp_amat,msize,vl,vu,il,iu,
     &             abstol,numeigval,rpoles,umat,msize,isupp,
     &             tmp_work,lwork,tmp_iwork,liwork,info)

c       call dsyevr(task,range,uplo,msize,amat,msize,vl,vu,il,iu,
c     &             abstol,numeigval,rpoles,umat,msize,isupp,
c     &             tmp_work,-1,tmp_iwork,-1,info)
c       print *
c       print *, 'optimal lwork  = ', tmp_work(1)
c       print *, 'optimal liwork = ', tmp_iwork(1)
c       stop

       if (info.eq.0) then ;
c       print *, 'done'
       else
          print '(/,a)', 'CAUTION: unsuccessful termination of DSYEVR !'
          print '(/,i3,/)', ' -- info =', info
          stop ': transport module is terminated now'
       end if

       deallocate(tmp_amat,tmp_work,tmp_iwork,isupp,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE real_spectrum]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
       end if

      end subroutine real_spectrum	


      end module rmatrix	
      