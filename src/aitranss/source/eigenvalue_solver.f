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
c     last revision:  jan 2012
c###########################################################

      module eigenvalue_solver
c     *************************************************************
c     for a hermitian matrix, a relevant (fastest) routine to solve 
c     eigevalue problem is called 'zheevr.f', see LAPACK online at
c     http://www.netlib.org/lapack/explore-html/d9/dd2/zheevr_8f.html
c     *************************************************************
      use globalvars
      implicit none

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      parameters to call lapack subroutine 'zheevr'
       character, private, parameter ::
c    &              task ='N', ! compute eigenvalues only
     &              range='A', ! all eigenvalues will be found
     &              uplo ='U'  ! upper triangle is stored

       double precision, private, parameter :: abstol = 1.0d-14

       double precision, private :: vl=0.0d0, vu=1.0d0
c                          ! not actually used
       integer, private :: il=1, iu=10
c                          ! not used as well

c      local temporary arrays
       complex(8), private, allocatable       :: u_matrix(:,:)

       complex(8), private, allocatable       :: work(:)
       double precision, private, allocatable :: rwork(:) 
       integer, private, allocatable          :: isupp(:), iwork(:)
       integer, private :: n_found_eigenvalues, lwork, lrwork, liwork
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      contains

      subroutine find_spectrum(h_matrix,eigenvalues,msize)
c     ! purpose: 
c     !   call the advanced routine 'zheevr' from LAPACK
c     !   to solve eigenvalue problem for the hermitian matrix
c     !
c     ! input: hamiltonian matrix & its size
       complex(8), intent(out)       :: h_matrix(:,:)
       integer, intent(out)          :: msize
c     !  ouput: array with real eigenvalues
       double precision, intent(out) :: eigenvalues(:)

c     !  other local variables       
       integer info, ierr
       
       allocate(u_matrix(msize,msize),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUBROUTINE realspectrum]: <u_matrix> allocation failure'
       end if
       
c      lwork = 33*msize
       lwork = 10*msize
       allocate(work(lwork),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUBROUTINE realspectrum]: <work> allocation failure'
       end if
       work = czero

c      lrwork = 33*msize
       lrwork = 24*msize
       allocate(rwork(lrwork),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUBROUTINE realspectrum]: <rwork> allocation failure'
       end if
       rwork = 0.0d0

       liwork = 10*msize
       allocate(iwork(liwork),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUBROUTINE realspectrum]: <iwork> allocation failure'
       end if
       iwork = 0

       allocate(isupp(2*msize),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUBROUTINE realspectrum]: <isupp> allocation failure'
       end if
       isupp = 0
        
       eigenvalues = 0.0d0
       u_matrix = czero

c      SUBROUTINE ZHEEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
c     $                   ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,
c     $                   RWORK, LRWORK, IWORK, LIWORK, INFO )
c		
       call zheevr('N',range,uplo,msize,h_matrix,msize,vl,vu,il,iu,
     &             abstol,n_found_eigenvalues,eigenvalues,u_matrix,msize,isupp,
     &             work,lwork,rwork,lrwork,iwork,liwork,info)
c       print *,  ' info = ', info

c       call zheevr('N',range,uplo,msize,h_matrix,msize,vl,vu,il,iu,
c     &             abstol,n_found_eigenvalues,eigenvalues,u_matrix,msize,isupp,
c     &             work,-1,rwork,-1,iwork,-1,info)
c        print *
c        print *, 'optimal lwork  = ', work(1)
c        print *, 'optimal lrwork = ', rwork(1)
c        print *, 'optimal liwork = ', iwork(1)
c        stop

       if (info.ne.0) then 
          print '(/,a)', 'CAUTION: unsuccessful termination of ZHEEVR !'
          print *, ' -- info = ', info
c         stop ': program will be terminated now'
       end if

       deallocate(u_matrix,work,rwork,iwork,isupp,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE realspectrum]: ',
     &            'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
       end if

      end subroutine find_spectrum	

      subroutine find_xspectrum(h_matrix,eigenvalues,umat,msize)
c     ! purpose: 
c     !   call the advanced routine 'zheevr' from LAPACK
c     !   to solve eigenvalue problem for the hermitian matrix
c     !
c     ! input: hamiltonian matrix & its size
       complex(8), intent(out)       :: h_matrix(:,:)
       integer, intent(out)          :: msize
c     !  ouput: arrays with eigenvectors and real eigenvalues
       double precision, intent(out) :: eigenvalues(:)
       complex(8), intent(out)       :: umat(:,:)

c     !  other local variables       
       integer info, ierr
       
c      lwork = 33*msize
       lwork = 10*msize
       allocate(work(lwork),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUBROUTINE realspectrum]: <work> allocation failure'
       end if
       work = czero

c      lrwork = 33*msize
       lrwork = 24*msize
       allocate(rwork(lrwork),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUBROUTINE realspectrum]: <rwork> allocation failure'
       end if
       rwork = 0.0d0

       liwork = 10*msize
       allocate(iwork(liwork),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUBROUTINE realspectrum]: <iwork> allocation failure'
       end if
       iwork = 0

       allocate(isupp(2*msize),stat=ierr)
       if (ierr.ne.0) then
        print *
        stop '[SUBROUTINE realspectrum]: <isupp> allocation failure'
       end if
       isupp = 0
        
       eigenvalues = 0.0d0
       umat = czero

c      SUBROUTINE ZHEEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
c     $                   ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,
c     $                   RWORK, LRWORK, IWORK, LIWORK, INFO )
c		
       call zheevr('V',range,uplo,msize,h_matrix,msize,vl,vu,il,iu,
     &             abstol,n_found_eigenvalues,eigenvalues,umat,msize,isupp,
     &             work,lwork,rwork,lrwork,iwork,liwork,info)
c       print *,  ' info = ', info

c       call zheevr('V',range,uplo,msize,h_matrix,msize,vl,vu,il,iu,
c     &             abstol,n_found_eigenvalues,eigenvalues,u_matrix,msize,isupp,
c     &             work,-1,rwork,-1,iwork,-1,info)
c        print *
c        print *, 'optimal lwork  = ', work(1)
c        print *, 'optimal lrwork = ', rwork(1)
c        print *, 'optimal liwork = ', iwork(1)
c        stop

       if (info.ne.0) then 
          print '(/,a)', 'CAUTION: unsuccessful termination of ZHEEVR !'
          print *, ' -- info = ', info
c         stop ': program will be terminated now'
       end if

       deallocate(work,rwork,iwork,isupp,stat=ierr)
       if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE realspectrum]: ',
     &            'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
       end if

      end subroutine find_xspectrum	

      end module eigenvalue_solver	
      
