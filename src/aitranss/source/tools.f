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

      module tools
c     ***************************************
c     contains a set of auxiliary subroutines 
c     called in different parts of the code
c     ***************************************

       use globalvars
       implicit none
      
       contains 

       subroutine greensf(epoint,ispin,gftmp)
c       **********************************************************
c       computes the green's function in the basis of eigenvectors
c       g(E) = B^{-1} * [E-h_full]^{-1} * B
c       **********************************************************
        implicit none

        double precision, intent (in) :: epoint
        integer,          intent (in) :: ispin
c       output array: greens function (GF) for E=epoint
        complex(8), intent (out) :: gftmp(:)

        integer p
c       print '(/,2x,a,$)', 'evaluating GF ...'
        do p = 1, nsaos
c                    E --> E + i*eta
          gftmp(p) = cone/(cmplx(epoint,0.0d0,8) - zpoles(p,ispin))
        end do
c       print *, 'done'

       end subroutine greensf

       subroutine sortarray(in_arr, out_arr, dima, index_arr)
c      *********************************************************       
c      sorting in ascending order, from smaller to larger values 
c      index_arr(:) -- shows the indices related to input array
c      *********************************************************
        implicit none
	
c       input/output	
        double precision, intent (in)  :: in_arr(:)
	integer, intent (in)           :: dima
        double precision, intent (out) :: out_arr(:)
	integer, intent (out)          :: index_arr(:)
	
c       local variables        
        integer           i, n, j_x, j_y
	double precision  xx, yy

        do i = 1, dima
	 out_arr(i) = in_arr(i)
	 index_arr(i) = i
	end do 

        do n = dima, 2, -1
	 do i = 2, n
	  xx = out_arr(i-1) 
          yy = out_arr(i)
	  j_x = index_arr(i-1)
	  j_y = index_arr(i) 
          if (xx .gt. yy) then
     	   out_arr(i-1) = yy
	   out_arr(i)   = xx
	   index_arr(i-1) = j_y
	   index_arr(i)   = j_x
	  end if
	 end do
	end do

       end subroutine sortarray    

       subroutine cmplxsort(in_arr, out_arr, dima, index_arr)
c      ******************************************************       
c      sorting a complex valued array 
c      in ascending order regarding a real part only 
c      index_arr(:) -- shows the indecies 
c                      related to the input array
c      ******************************************************
        implicit none
	
c       input/output	
        complex(8), intent (in)  :: in_arr(:)
	integer, intent (in)     :: dima
        complex(8), intent (out) :: out_arr(:)
	integer, intent (out)    :: index_arr(:)
	
c       local variables        
        integer           i, n, j_x, j_y
	complex(8)        xx, yy
	double precision  rexx, reyy

        do i = 1, dima
	 out_arr(i) = in_arr(i)
	 index_arr(i) = i
	end do 

        do n = dima, 2, -1
	 do i = 2, n
	  xx = out_arr(i-1) 
          yy = out_arr(i)
	  j_x = index_arr(i-1)
	  j_y = index_arr(i) 
	  rexx = dble(xx)
	  reyy = dble(yy)
          if (rexx .gt. reyy) then
     	   out_arr(i-1) = yy
	   out_arr(i)   = xx
	   index_arr(i-1) = j_y
	   index_arr(i)   = j_x
	  end if
	 end do
	end do

       end subroutine cmplxsort       

       subroutine sort_moset(in_arr,out_arr,dima,index_arr)
c      *********************************************************
c      sorting in ascending order the array of type 'spectrum'
c      >> from smaller to larger values <<
c      index_arr(:) -- shows the indecies related to input array
c      *********************************************************
        implicit none

c       input/output
        type(spectrum), intent (in)    :: in_arr(:)
        integer, intent (in)           :: dima
        type(spectrum), intent (out)   :: out_arr(:)
        integer, intent (out)          :: index_arr(:)

c       local variables
        integer           i, n, j_x, j_y
        double precision  xx, yy
        integer xspin, yspin, xindex, yindex

        do i = 1, dima
         out_arr(i)%en      = in_arr(i)%en
         out_arr(i)%spin    = in_arr(i)%spin
         out_arr(i)%moindex = in_arr(i)%moindex
         index_arr(i) = i
        end do

        do n = dima, 2, -1
         do i = 2, n
          xx = out_arr(i-1)%en
          yy = out_arr(i)%en

          xspin = out_arr(i-1)%spin
          yspin = out_arr(i)%spin
          xindex = out_arr(i-1)%moindex
          yindex = out_arr(i)%moindex

          j_x = index_arr(i-1)
          j_y = index_arr(i)

          if (xx .gt. yy) then
           out_arr(i-1)%en = yy
           out_arr(i)%en   = xx

           out_arr(i-1)%spin = yspin
           out_arr(i)%spin   = xspin
           out_arr(i-1)%moindex = yindex
           out_arr(i)%moindex   = xindex

           index_arr(i-1) = j_y
           index_arr(i)   = j_x
          end if
         end do
        end do
       end subroutine sort_moset

       subroutine crossprod(avec, bvec, cvec)
c      ********************************************       
c      computes a cross product: cvec = avec x bvec
c      ********************************************       
        implicit none        
	double precision, intent  (in) :: avec(3), bvec(3)
	double precision, intent (out) :: cvec(3)
        double precision  tmpa, tmpb, tmpc, tmpd 
c      ! to be on the safe side
        cvec = 0.0d0  
c      x-component
        tmpa = avec(2) ;  tmpb = avec(3)
        tmpc = bvec(2) ;  tmpd = bvec(3)
	cvec(1) = tmpa*tmpd - tmpb*tmpc 
c      y-component
        tmpa = avec(3) ;  tmpb = avec(1)
        tmpc = bvec(3) ;  tmpd = bvec(1)
	cvec(2) = tmpa*tmpd - tmpb*tmpc 
c      z-component        
        tmpa = avec(1) ;  tmpb = avec(2)
        tmpc = bvec(1) ;  tmpd = bvec(2)
	cvec(3) = tmpa*tmpd - tmpb*tmpc 
       end subroutine crossprod       


       double precision function cmplxabs(z)
c       abs(x+iy) = |x| + |y| 
        implicit none
	complex(8), intent(in) :: z
        cmplxabs = abs(dble(z)) + abs(dimag(z))
       end function cmplxabs

 
       subroutine dble_ab(opa,opb,ndim,amat,bmat,cmat)
c      *****************************************************       
c      performs multiplication of real-valued matrices
c      C = opa(A)*opb(B), where ndim is a matrix dimension,    
c      using blas routine  DGEMM;
c      notations for operations (op) are taken from blaslib
c      *****************************************************       
        implicit none 

        character,  intent(in)        :: opa, opb
	integer,    intent(in)        :: ndim
	double precision, intent(in)  :: amat(:,:), bmat(:,:)
	double precision, intent(out) :: cmat(:,:)

	call dgemm(opa,opb,ndim,ndim,ndim,1.0d0,
     &            amat,ndim,bmat,ndim,0.0d0,cmat,ndim)
    
       end subroutine dble_ab

       subroutine xdble_ab(opa,opb,m,k,n,amat,bmat,cmat)
c      *****************************************************       
c      performs multiplication of real-valued matrices
c      C = opa(A)*opb(B), where 
c      --- opa(A) is m x k matrix
c      --- opb(B) is k x n matrix
c      --- C      is m x n matrix
c      the blas routine DGEMM is used;
c      notations for operations (op) are taken from blaslib
c      *****************************************************       
        implicit none 

        character,  intent(in)        :: opa, opb
	integer,    intent(in)        :: m, k, n
	double precision, intent(in)  :: amat(:,:), bmat(:,:)
	double precision, intent(out) :: cmat(:,:)

        integer lda, ldb, ldc
	
c	get dimensions
        if ((opa.eq.'n').or.(opa.eq.'N')) then
	  lda = m  
	else if ((opa.eq.'t').or.(opa.eq.'T')) then
	  lda = k 
	else
	 stop ': opa(A) >> wrong input in [xdble_ab] !'
	end if

        if ((opb.eq.'n').or.(opb.eq.'N')) then
	  ldb = k
	else if ((opb.eq.'t').or.(opb.eq.'T')) then
	  ldb = n
	else
	 stop ': opa(B) >> wrong input in [xdble_ab] !'
	end if
	ldc = m

	call dgemm(opa,opb,m,n,k,1.0d0,
     &	           amat,lda,bmat,ldb,0.0d0,cmat,ldc)
    
       end subroutine xdble_ab
       
       subroutine cmplx_ab(opa,opb,ndim,amat,bmat,cmat)
c      *****************************************************       
c      performs multiplication of comlex-valued matrices
c      C = opa(A)*opb(B), where ndim is a matrix dimension,    
c      using blas routine ZGEMM;
c      notations for operations (op) are taken from blaslib
c      *****************************************************       
        implicit none 

        character,  intent(in)  :: opa, opb
	integer,    intent(in)  :: ndim
	complex(8), intent(in)  :: amat(:,:), bmat(:,:)
	complex(8), intent(out) :: cmat(:,:)

c       cone  = (1.0d0, 0.0d0)
c       czero = (0.0d0, 0.0d0)
	call zgemm(opa,opb,ndim,ndim,ndim,cone,
     &            amat,ndim,bmat,ndim,czero,cmat,ndim)
	
       end subroutine cmplx_ab	
       
       subroutine xcmplx_ab(opa,opb,m,k,n,amat,bmat,cmat)
c      *****************************************************       
c      performs multiplication of complex-valued matrices
c      C = opa(A)*opb(B), where 
c      --- opa(A) is m x k matrix
c      --- opb(B) is k x n matrix
c      --- C      is m x n matrix
c      the blas routine ZGEMM is used;
c      notations for operations (op) are taken from blaslib
c      *****************************************************       
        implicit none 

        character,  intent(in)  :: opa, opb
	integer,    intent(in)  :: m, k, n
	complex(8), intent(in)  :: amat(:,:), bmat(:,:)
	complex(8), intent(out) :: cmat(:,:)

        integer lda, ldb, ldc
	
c	get dimensions
        if ((opa.eq.'n').or.(opa.eq.'N')) then
	  lda = m  
	else if ((opa.eq.'t').or.(opa.eq.'T').or.
     &	         (opa.eq.'c').or.(opa.eq.'C')) then
	  lda = k 
	else
	 stop ': opa(A) >> wrong input in [xcmplx_ab] !'
	end if

        if ((opb.eq.'n').or.(opb.eq.'N')) then
	  ldb = k
	else if ((opb.eq.'t').or.(opb.eq.'T').or.
     &	         (opb.eq.'c').or.(opb.eq.'C')) then
	  ldb = n
	else
	 stop ': opa(B) >> wrong input in [xcmplx_ab] !'
	end if
	ldc = m

c       cone  = (1.0d0, 0.0d0)
c       czero = (0.0d0, 0.0d0)
	call zgemm(opa,opb,m,n,k,cone,
     &	           amat,lda,bmat,ldb,czero,cmat,ldc)
    
       end subroutine xcmplx_ab

       subroutine dcmplx_ab(opa,opb,ndim,amat,bmat,cmat)
c      *****************************************************       
c      performs multiplication of comlex-valued matricies
c      C = opa(A)*opb(B), where ndim is a matrix dimension,    
c      using blas routine DGEMM instead of ZGEMM;
c      notations for operations (op) are taken from blaslib
c      *****************************************************       
        implicit none 

        character,  intent(in)  :: opa, opb
	integer,    intent(in)  :: ndim
	complex(8), intent(in)  :: amat(:,:), bmat(:,:)
	complex(8), intent(out) :: cmat(:,:)

        double precision, allocatable :: re_a(:,:), im_a(:,:), 
     &	                                 re_b(:,:), im_b(:,:), 
     &                                   re_c(:,:), im_c(:,:),
     &                                   re1(:,:),  re2(:,:),
     &                                   im1(:,:),  im2(:,:),
     &                                   re(:,:),   im(:,:)
        integer ierr, n, m
	
        allocate(re_a(ndim,ndim),stat=ierr)   
        allocate(re_b(ndim,ndim),stat=ierr)   
        allocate(re_c(ndim,ndim),stat=ierr)   
        allocate(im_a(ndim,ndim),stat=ierr)   
        allocate(im_b(ndim,ndim),stat=ierr)   
        allocate(im_c(ndim,ndim),stat=ierr)   
      
        allocate(re1(ndim,ndim),stat=ierr)
        allocate(re2(ndim,ndim),stat=ierr)
        allocate(im1(ndim,ndim),stat=ierr)
        allocate(im2(ndim,ndim),stat=ierr)

        allocate(re(ndim,ndim),stat=ierr)
        allocate(im(ndim,ndim),stat=ierr)
   
        do n = 1, ndim
	 do m = 1, ndim

         if (opa=='n') then
	  re_a(n,m) =  dble(amat(n,m)) 
 	  im_a(n,m) = dimag(amat(n,m)) 
	 else if (opa=='c') then 
c         hermitian conjugation	 
	  re_a(n,m) =  dble(amat(m,n)) 
 	  im_a(n,m) = -dimag(amat(m,n)) 
	 else
	  print *
	  stop ': incorrect operation <opa> in cmplx_ab!'
	 end if

        if (opb=='n') then
	  re_b(n,m) =  dble(bmat(n,m)) 
 	  im_b(n,m) = dimag(bmat(n,m)) 
	 else if (opb=='c') then
c         hermitian conjugation	 
	  re_b(n,m) =  dble(bmat(m,n)) 
 	  im_b(n,m) = -dimag(bmat(m,n)) 
	 else
	  print *
	  stop ': incorrect operation <opb> in cmplx_ab!'
	 end if
	  
	 end do
	end do
   
c       real part
	call dgemm('n','n',ndim,ndim,ndim,1.0d0,
     &            re_a,ndim,re_b,ndim,0.0d0,re1,ndim)
	call dgemm('n','n',ndim,ndim,ndim,1.0d0,
     &            im_a,ndim,im_b,ndim,0.0d0,re2,ndim)
        re = re1 - re2

c       imaginary part
	call dgemm('n','n',ndim,ndim,ndim,1.0d0,
     &            im_a,ndim,re_b,ndim,0.0d0,im1,ndim)
	call dgemm('n','n',ndim,ndim,ndim,1.0d0,
     &            re_a,ndim,im_b,ndim,0.0d0,im2,ndim)
        im = im1 + im2

c       output: cmat = re + i* im
        forall(n=1:ndim,m=1:ndim)
         cmat(n,m) = cmplx(re(n,m),im(n,m),8)
	end forall 

        deallocate(re_a,re_b,re_c,im_a,im_b,im_c,
     &    	   re1,re2,im1,im2,re,im,stat=ierr)
	
       end subroutine dcmplx_ab


       subroutine deallocall
c       ************************************       
c       -- called from the main module
c       -- deallocate all global arrays used
c       ************************************

        implicit none
	integer :: jerr, icall = 1
    
        print '(/,a,$)', ' deallocating all arrays used ...'
c       # icall = 1
        deallocate(h0mat,smat12,smat12inv,
     &	           gamma_left,gamma_right,stat=jerr)
        if (jerr.ne.0) call errdealloc(icall,jerr) 
	icall = icall + 1

c       # icall = 2
        if (sp_leads) then
	  deallocate(sigmann,bsigmann,stat=jerr)
	else
	  deallocate(sigmann,stat=jerr)
	end if 
        if (jerr.ne.0) call errdealloc(icall,jerr)  
	icall = icall + 1

c       # icall = 3
        if (.not.import_self_energy) then
	  deallocate(at_self_energy,stat=jerr)
	end if 
        if (jerr.ne.0) call errdealloc(icall,jerr)  
	icall = icall + 1

c       # icall = 4
        if (allocated(hmat)) then 
          deallocate(hmat,stat=jerr)
          if (jerr.ne.0)  call errdealloc(icall,jerr) 
	end if
        icall = icall + 1

c       # icall = 5
        if (allocated(zpoles)) then 
          deallocate(zpoles,stat=jerr)
          if (jerr.ne.0)  call errdealloc(icall,jerr) 
	end if
        icall = icall + 1

c       # icall = 6
        if (allocated(heigvec)) then 
          deallocate(heigvec,stat=jerr)
          if (jerr.ne.0)  call errdealloc(icall,jerr) 
	end if
        icall = icall + 1

c       # icall = 7
        if (allocated(invheigvec)) then 
          deallocate(invheigvec,stat=jerr)
          if (jerr.ne.0)  call errdealloc(icall,jerr) 
	end if
        icall = icall + 1

c       # icall = 8
        deallocate(atom,stat=jerr)
        if (jerr.ne.0) call errdealloc(icall,jerr) 
	icall = icall + 1

c       # icall = 9
        deallocate(atom_type,atom_sort,ref_sort,stat=jerr)
        if (jerr.ne.0)  call errdealloc(icall,jerr) 
	icall = icall + 1

c       # icall = 10
        deallocate(n_basis_func,stat=jerr)
        if (jerr.ne.0)  call errdealloc(icall,jerr) 
	icall = icall + 1

c       # icall = 11
        deallocate(mo_en,stat=jerr)
        if (jerr.ne.0) call errdealloc(icall,jerr) 
	icall = icall + 1

c       # icall = 12        
        deallocate(mo_coeff,stat=jerr)
        if (jerr.ne.0) call errdealloc(icall,jerr) 
	icall = icall + 1

c       # icall = 13
        deallocate(smat,stat=jerr)
        if (jerr.ne.0) call errdealloc(icall,jerr) 
	icall = icall + 1

c       # icall = 14 
        if (do_pop.or.no_leads) then
         deallocate(mo_orth,stat=jerr)
         if (jerr.ne.0) call errdealloc(icall,jerr) 
        end if
        icall = icall + 1  
	
c       # icall = 15 
        if (ldau) then
         deallocate(osmat,stat=jerr)
         if (jerr.ne.0) call errdealloc(icall,jerr) 
        end if
        icall = icall + 1  
	
c       # icall = 16 
        if (ldau.and.spectr_func) then
         deallocate(cmplx_osmat,stat=jerr)
         if (jerr.ne.0) call errdealloc(icall,jerr)  
        end if
        icall = icall + 1  

c       # icall = 17
c       # takes care about non-eq. density
        if (allocated(neqdmat)) then
	 deallocate(neqdmat,stat=jerr)
         if (jerr.ne.0) call errdealloc(icall,jerr)  
	end if
	icall = icall + 1 

c       # icall = 18
c       # takes care about normalization factors 
        if (allocated(orbnorm)) then
	 deallocate(orbnorm,stat=jerr)
         if (jerr.ne.0) call errdealloc(icall,jerr)  
	end if
	icall = icall + 1 

c       # icall = 19
c       # takes care about conduction eigenchannels 
        if (allocated(utt)) then
	 deallocate(utt,stat=jerr)
         if (jerr.ne.0) call errdealloc(icall,jerr)  
	end if
	icall = icall + 1 
	
	print *, 'done'
	
       end subroutine deallocall

       subroutine errdealloc(icall,ierr)
c       *******************************
c       -- prints out error message --        
c       *******************************
        implicit none 
	integer, intent (in) :: icall, ierr
        print '(/,a,i1,a)', ' call ',icall,': failed with deallocation'
        print '(a,i2)',     ' stat =  ',ierr	
       end subroutine errdealloc


      end module tools



  