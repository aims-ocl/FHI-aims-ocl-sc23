c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           april 2009
c     last revision:  jan   2012
c###########################################################

      module domainwall
c     *************************************      
c     set of routines required to build up
c     a domain wall solution      
c     *************************************      
       use globalvars
       use tools
   
       implicit none
 
       real(4), private ::  stime, ftime, secs
       integer, private ::  mnts  

      contains
 
      subroutine build_z_ur(prs,nuindex,ur)
c      ****************************************************
c      builds up an orthogonal matrix ur (ur^{-1} = ur^T) 
c      which transforms a density matrix under reflection 
c      (z'=-z,x'=x,y'=y), so that dmat' = ur^T * dmat * ur 
c      ****************************************************
       integer, intent(in) :: prs(:), nuindex(:)
       double precision, intent(out) :: ur(:,:)
       
       integer iat, jat, attype, nu, mu, n, m, nn, norb, lm
       logical flag
       
       ur = 0.0d0
       do iat = 1, num_atoms
c       take the counterpart atom
	jat= prs(iat)
c       define initial indices referring 
c       to the upper coner of the (iat,jat)-block
        nu = nuindex(iat)
	mu = nuindex(jat)
	
c       fill up diagonal elements of the (iat,jat)-block of 
c       the matrix <ur> -- its size is norb x norb
	attype = atom(iat)%atype
	norb = n_basis_func(attype)
	do nn = 1, norb
c        define global indecies
	 n = nu + nn-1 
	 m = mu + nn-1
c        take lm-index of the current orbital
	 lm = aos(iat,nn)%lm 
	 flag = (lm.eq.4).or.                           ! p_z
     & 	        (lm.eq.6).or.(lm.eq.7).or.              ! d_xz d_yz
     &          (lm.eq.10).or.(lm.eq.13).or.(lm.eq.14)  ! f_z3, xyz, (xx-yy)z
         if (flag) then
	  ur(n,m) = -1.0d0
	 else
          ur(n,m) =  1.0d0
	 end if  
        end do
        
       end do
       
      end subroutine build_z_ur

      subroutine def_z_pairs(prs,refind)
c      ***********************************************************
c      -- searches for pairs of equivalent atoms 
c         under reflection: z --> -z
c      -- initializes the reference array <refind>:
c         for a given atom 'iat' it returns an index nu=refind(iat), 
c         such that next subgroup of indices (nu, nu+1, ...) 
c         refers to internal degrees of freedom of the atom iat
c      ***********************************************************
       
       integer, intent(out) :: prs(:), refind(:)
       
       integer          iat, jat, iatype, nu, norb       
       double precision x0,y0,z0, x1,y1,z1, dr
       logical          pfound
       double precision, parameter :: reps = 1.0e-3
       
       do iat = 1, num_atoms
       
        x0 = atom(iat)%pos(1)
        y0 = atom(iat)%pos(2)
        z0 = atom(iat)%pos(3)

        if (dabs(z0)<reps) then
c        counterpart is the same atom (since z=0) 	
	 prs(iat) = iat
        else
c        scanning all atoms, search for that one
c        with x1=x0, y1=y0, and z1=-z0
	 pfound = .false.
	 jat = 0
	 do while ((.not.pfound).and.(jat<num_atoms))
          jat = jat + 1
	  x1 = atom(jat)%pos(1)
	  y1 = atom(jat)%pos(2)
	  z1 = atom(jat)%pos(3)
	  dr = dabs(x1-x0)+dabs(y1-y0)+dabs(z1+z0)
	  if ( (dr < 3.0*reps).and.(jat.ne.iat)) then
	   pfound = .true.
	  end if
	 end do
         if (pfound.and.(atom(iat)%atype.eq.atom(jat)%atype)) then
c         counterpart for iat is jat >>>
	  prs(iat) = jat
	 else
c         no counterpart is found !
          write (*,'(/,a,i3,/)') 
     &	  '[SUBROUTINE def_z_pairs]: no pair atom is found for atom ',iat          
          stop 'transport module is terminated now'	 
	 end if
       
        end if
       end do

c      build up a reference index array <refind>
       nu = 1
       do iat = 1, num_atoms
        refind(iat) = nu
        iatype = atom(iat)%atype
	norb   = n_basis_func(iatype) 
        nu = nu + norb
       end do      
      
      end subroutine def_z_pairs

      subroutine z_updatedmat(dmatrix)
c      *******************************************************
c      in case of domain-wall solution, this routine takes 
c      care about symmetrization of the density matrix, namely 
c      z --> -z : left/spin-up   == right/spin-down 
c               & left/spin-down == right/spin-up
c      *******************************************************
c       input/output: density matrix (either non-eq or equilibrium one)  
        double precision, allocatable :: dmatrix(:,:,:)

        integer, allocatable :: atpairs(:), nuind(:)
	integer iat, jat, ispin, jspin, n, m, ierr        
        double precision, allocatable :: umat(:,:), udmat(:,:), 
     &	                                 dmatrix1(:,:,:), u1(:,:),
     &                                   dmat0(:,:), dmat1(:,:)
        double precision nmtmp
        double precision, parameter :: ueps = 1.0d-10
	logical herm	

        allocate(atpairs(num_atoms),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE z_updatedmat]: <atpairs> allocation failure'
        end if       
        forall (iat=1:num_atoms) atpairs(iat)=0

        allocate(nuind(num_atoms),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE z_updatedmat]: <nuind> allocation failure'
        end if       
        forall (iat=1:num_atoms) nuind(iat)=0

        write(*,'(/,1x,a,/)') 
     &       'DOMAIN-WALL SOLUTION: searching for atom pairs (z <--> -z)'
        call def_z_pairs(atpairs,nuind) 
        write(*,'(a)') '    atom     counterpart'
        write(*,'(a)') ' ------------------------'
      
        do iat = 1, num_atoms
         jat = atpairs(iat)
c         write(*,'(i5,a4,3x,a,i7,a4,i6)') 
c     &	   iat, atom(iat)%symbol(1:2),':', jat, atom(jat)%symbol(1:2),nuind(iat)
         write(*,'(i5,a4,3x,a,i7,a4)') 
     &	   iat, atom(iat)%symbol(1:2),':', jat, atom(jat)%symbol(1:2)
        end do
        write(*,'(a)') ' ------------------------'
    
        allocate(dmatrix1(nsaos,nsaos,nspin),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE z_updatedmat]: <dmatrix1> allocation failure'
        end if       
        dmatrix1 = 0.0d0
        
        allocate(dmat0(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE z_updatedmat]: <dmat0> allocation failure'
        end if       

        allocate(dmat1(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE z_updatedmat]: <dmat1> allocation failure'
        end if       

        allocate(umat(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE z_updatedmat]: <umat> allocation failure'
        end if       
        umat = 0.0d0

        allocate(u1(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
	  stop
     &	  '[SUBROUTINE z_updatedmat]: <u1> allocation failure'
        end if       

        allocate(udmat(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE z_updatedmat]: <udmat> allocation failure'
        end if       
        
	print '(/,1x,a)', 'SYMMETRIZING DENSITY MATRIX >>>'

        call cpu_time(stime)       
c       initialize transformation matrix >>>
        print '(/,2x,a,$)', 'building up symmetry transform matrix U ...'
        call build_z_ur(atpairs,nuind,umat)
        print *, 'done'

c       checking that:  umat * umat^T = 1  
        call dble_ab('n','t',nsaos,umat,umat,u1)
 	herm = .true.
 	do n = 1, nsaos 
	 do m = 1, nsaos
	  if (m==n) then ; nmtmp = 1.0d0
	  else           ; nmtmp = 0.0d0
	  end if
	  herm = (dabs(u1(n,m)-nmtmp)<ueps) 
	 end do
	end do

        if (herm) then
	 print *, ' checking: matrix U is orthogonal'
        else 
	 print *, ' WARNING: matrix U is NOT orthogonal'
	end if

c       here nspin == 2 ! 
        if (nspin.eq.1) 
     &   stop '[SUB. z_updatedmat]: nspin==1 ! ... confused and will quit now!'

        print '(2x,a,$)', 'updating density matrix ...'
        do ispin = 1, nspin
        
         if (ap_symmetrize) then
           if (ispin.eq.1) then ; jspin = 2
 	   else                 ; jspin = 1
 	   end if
 	 else if (p_symmetrize) then
           if (ispin.eq.1) then ; jspin = 1
 	   else                 ; jspin = 2
 	   end if
 	 else
 	  stop '[SUB. z_updatedmat]: i am confused, neither p_symmetrize nor ap_symmetrize are active!'
 	 end if

c        save a local copy of dmatrix for spin 'jspin' opposite to 'ispin'
         forall(n=1:nsaos,m=1:nsaos)
          dmat0(n,m) = dmatrix(n,m,jspin)
         end forall
     
c        updating a density matrix, to insure a symmetrical domain wall solution
c        -- dmat_(sigma)  <--  1/2 *(dmat_(sigma) + umat^T * dmat_(-sigma) * umat)  
c        -- where 'sigma' is spin index; <umat> is xy-plane reflection 

c        1. compute dmat1 = umat^T * dmat0 * umat        
         call dble_ab('t','n',nsaos,umat,dmat0,udmat)
	 call dble_ab('n','n',nsaos,udmat,umat,dmat1)

c        2. admix different spin channels
         do m = 1, nsaos
	  do n = 1, nsaos
            dmatrix1(n,m,ispin) = 0.5d0*(dmatrix(n,m,ispin) + dmat1(n,m))	  
	  end do
	 end do

        end do ! ispin
	
c       put <dmatrix1> to <neqdmat>
	do ispin = 1, nspin
         do m = 1, nsaos
	  do n = 1, nsaos
           dmatrix(n,m,ispin) = dmatrix1(n,m,ispin) 	  
	  end do
	 end do
	end do  ! ispin
        print *, 'done'
	
        call cpu_time(ftime)
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts*60.0d0
        print '(2x,a,i3,a,f5.2,a)',
     &        'time spent:', mnts, ' min ', secs, ' sec'

	print '(/,1x,a)', '<<< DONE WITH SYMMETRIZATION '

        deallocate(atpairs,nuind,umat,udmat,dmat0,dmat1,dmatrix1,u1,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE z_updatedmat]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
        end if

      end subroutine z_updatedmat

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine build_x_ur(prs,nuindex,ur)
c      ****************************************************
c      builds up an orthogonal matrix ur (ur^{-1} = ur^T) 
c      which transforms a density matrix under reflection 
c      (x'=-x,y'=y,z'=z), so that dmat' = ur^T * dmat * ur 
c      ****************************************************
       integer, intent(in) :: prs(:), nuindex(:)
       double precision, intent(out) :: ur(:,:)
       
       integer iat, jat, attype, nu, mu, n, m, nn, norb, lm
       logical flag
       
       do iat = 1, num_atoms
c       take the counterpart atom
	jat= prs(iat)
c       define initial indices referring 
c       to the upper coner of the (iat,jat)-block
        nu = nuindex(iat)
	mu = nuindex(jat)
	
c       fill up diagonal elements of the (iat,jat)-block of 
c       the matrix <ur> -- its size is norb x norb
	attype = atom(iat)%atype
	norb = n_basis_func(attype)
	do nn = 1, norb
c        define global indecies
	 n = nu + nn-1 
	 m = mu + nn-1
c        take lm-index of the current orbital
	 lm = aos(iat,nn)%lm 
	 flag = (lm.eq.2).or.                           ! p_x
     & 	        (lm.eq.6).or.(lm.eq.8).or.              ! d_xz, d_xy
     &          (lm.eq.11).or.(lm.eq.13).or.(lm.eq.15)  ! -xxx-xyy+4xzz, xyz, xxx-3xyy
         if (flag) then
	  ur(n,m) = -1.0d0
	 else
          ur(n,m) =  1.0d0
         end if  
        end do
        
       end do
       
      end subroutine build_x_ur

      subroutine def_x_pairs(prs,refind)
c      ***********************************************************
c      -- searches for pairs of equivalent atoms 
c         under reflection: x --> -x
c      -- initializes the reference array <refind>:
c         for a given atom 'iat' it returns an index nu=refind(iat), 
c         such that next subgroup of indices (nu, nu+1, ...) 
c         refers to internal degrees of freedom of the atom iat
c      ***********************************************************
       
       integer, intent(out) :: prs(:), refind(:)
       
       integer          iat, jat, iatype, nu, norb       
       double precision x0,y0,z0, x1,y1,z1, dr
       logical          pfound
       double precision, parameter :: reps = 1.0e-3
       
       do iat = 1, num_atoms
       
        x0 = atom(iat)%pos(1)
        y0 = atom(iat)%pos(2)
        z0 = atom(iat)%pos(3)

        if (dabs(x0)<reps) then
c        counterpart is the same atom (since x=0)
	 prs(iat) = iat
        else
c        scanning all atoms, search for that one
c        with x1=-x0, y1=y0, and z1=z0
	 pfound = .false.
	 jat = 0
	 do while ((.not.pfound).and.(jat<num_atoms))
          jat = jat + 1
	  x1 = atom(jat)%pos(1)
	  y1 = atom(jat)%pos(2)
	  z1 = atom(jat)%pos(3)
	  dr = dabs(x1+x0)+dabs(y1-y0)+dabs(z1-z0)
	  if ( (dr < 3.0*reps).and.(jat.ne.iat)) then
	   pfound = .true.
	  end if
	 end do
         if (pfound.and.(atom(iat)%atype.eq.atom(jat)%atype)) then
c         counterpart for iat is jat >>>
	  prs(iat) = jat
	 else
c         no counterpart is found !
          write (*,'(/,a,i3,/)') 
     &	  '[SUBROUTINE def_x_pairs]: no pair atom is found for atom ',iat          
          stop 'transport module is terminated now'	 
	 end if
       
        end if
       end do

c      build up a reference index array <refind>
       nu = 1
       do iat = 1, num_atoms
        refind(iat) = nu
        iatype = atom(iat)%atype
	norb   = n_basis_func(iatype) 
        nu = nu + norb
       end do      
      
      end subroutine def_x_pairs

      subroutine x_updatedmat(dmatrix)
c      **************************************************************
c      in case of atomic structure with yz-mirror symmetry, this
c      routine takes care about symmetrization of the density matrix: 
c      x --> -x : left/spin-up   == right/spin-up 
c               & left/spin-down == right/spin-down
c      **************************************************************
c       input/output: density matrix (either non-eq or equilibrium one)  
        double precision, allocatable :: dmatrix(:,:,:)

        integer, allocatable :: atpairs(:), nuind(:)
	integer iat, jat, ispin, n, m, ierr        
        double precision, allocatable :: umat(:,:), udmat(:,:), 
     &	                                 dmatrix1(:,:,:), u1(:,:),
     &                                   dmat1(:,:)
        double precision nmtmp
        double precision, parameter :: ueps = 1.0d-10
	logical herm	

        allocate(atpairs(num_atoms),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE x_updatedmat]: <atpairs> allocation failure'
        end if       
        forall (iat=1:num_atoms) atpairs(iat)=0

        allocate(nuind(num_atoms),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE x_updatedmat]: <nuind> allocation failure'
        end if       
        forall (iat=1:num_atoms) nuind(iat)=0

        write(*,'(/,1x,a,/)') 
     &       'YZ-MIRROR SYMMETRIZED SOLUTION: searching for atom pairs (x <--> -x)'
        call def_x_pairs(atpairs,nuind) 
        write(*,'(a)') '    atom     counterpart'
        write(*,'(a)') ' ------------------------'
      
        do iat = 1, num_atoms
         jat = atpairs(iat)
c         write(*,'(i5,a4,3x,a,i7,a4,i6)') 
c     &	   iat, atom(iat)%symbol(1:2),':', jat, atom(jat)%symbol(1:2),nuind(iat)
         write(*,'(i5,a4,3x,a,i7,a4)') 
     &	   iat, atom(iat)%symbol(1:2),':', jat, atom(jat)%symbol(1:2)
        end do
        write(*,'(a)') ' ------------------------'
    
        allocate(dmatrix1(nsaos,nsaos,nspin),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE x_updatedmat]: <dmatrix1> allocation failure'
        end if       
        dmatrix1 = 0.0d0
        
        allocate(dmat1(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE x_updatedmat]: <dmat1> allocation failure'
        end if       

        allocate(umat(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE x_updatedmat]: <umat> allocation failure'
        end if       
        umat = 0.0d0

        allocate(u1(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
	  stop
     &	  '[SUBROUTINE x_updatedmat]: <u1> allocation failure'
        end if       

        allocate(udmat(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE x_updatedmat]: <udmat> allocation failure'
        end if       

        call cpu_time(stime)       
c       initialize transformation matrix >>>
        print '(/,2x,a,$)', 'building up symmetry transform matrix U ...'
        call build_x_ur(atpairs,nuind,umat)
        print *, 'done'

c       checking that:  umat * umat^T = 1  
        call dble_ab('n','t',nsaos,umat,umat,u1)
 	herm = .true.
 	do n = 1, nsaos 
	 do m = 1, nsaos
	  if (m==n) then ; nmtmp = 1.0d0
	  else           ; nmtmp = 0.0d0
	  end if
	  herm = (dabs(u1(n,m)-nmtmp)<ueps) 
	 end do
	end do

        if (herm) then
	 print *, ' checking: matrix U is orthogonal'
        else 
	 print *, ' WARNING: matrix U is NOT orthogonal'
	end if

        print '(2x,a,$)', 'updating density matrix ...'
        do ispin = 1, nspin

c        updating a density matrix, to insure a symmetrical solution
c        -- dmat_(ispin)  <--  1/2 *(dmat_(ispin) + umat^T * dmat_(ispin) * umat)  
c        -- where 'ispin' is a spin index; <umat> is zy-plane reflection 

c        1. compute dmat1 = umat^T * dmatrix * umat        
         call dble_ab('t','n',nsaos,umat,dmatrix(:,:,ispin),udmat)
	 call dble_ab('n','n',nsaos,udmat,umat,dmat1)

c        2. admix different spin channels
         do m = 1, nsaos
	  do n = 1, nsaos
            dmatrix1(n,m,ispin) = 0.5d0*(dmatrix(n,m,ispin) + dmat1(n,m))	  
	  end do
	 end do

        end do ! ispin
	
c       put <dmatrix1> to <neqdmat>
	do ispin = 1, nspin
         do m = 1, nsaos
	  do n = 1, nsaos
           dmatrix(n,m,ispin) = dmatrix1(n,m,ispin) 	  
	  end do
	 end do
	end do  ! ispin
        print *, 'done'
	
        call cpu_time(ftime)
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts*60.0d0
        print '(2x,a,i3,a,f5.2,a)',
     &        'time spent:', mnts, ' min ', secs, ' sec'

	print '(/,1x,a)', '<<< DONE WITH SYMMETRIZATION '

        deallocate(atpairs,nuind,umat,udmat,dmat1,dmatrix1,u1,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE x_updatedmat]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
        end if

      end subroutine x_updatedmat

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine build_y_ur(prs,nuindex,ur)
c      ****************************************************
c      builds up an orthogonal matrix ur (ur^{-1} = ur^T) 
c      which transforms a density matrix under reflection 
c      (x'=x,y'=-y,z'=z), so that dmat' = ur^T * dmat * ur 
c      ****************************************************
       integer, intent(in) :: prs(:), nuindex(:)
       double precision, intent(out) :: ur(:,:)
       
       integer iat, jat, attype, nu, mu, n, m, nn, norb, lm
       logical flag
       
       do iat = 1, num_atoms
c       take the counterpart atom
	jat= prs(iat)
c       define initial indices referring 
c       to the upper coner of the (iat,jat)-block
        nu = nuindex(iat)
	mu = nuindex(jat)
	
c       fill up diagonal elements of the (iat,jat)-block of 
c       the matrix <ur> -- its size is norb x norb
	attype = atom(iat)%atype
	norb = n_basis_func(attype)
	do nn = 1, norb
c        define global indecies
	 n = nu + nn-1 
	 m = mu + nn-1
c        take lm-index of the current orbital
	 lm = aos(iat,nn)%lm 
	 flag = (lm.eq.3).or.                           ! p_y
     & 	        (lm.eq.7).or.(lm.eq.9).or.              ! d_yz, d_xy
     &          (lm.eq.12).or.(lm.eq.13).or.(lm.eq.16)  ! -yyy-xxy+4yzz, xyz, yyy-3xxy
         if (flag) then
	  ur(n,m) = -1.0d0
	 else
          ur(n,m) =  1.0d0
         end if  
        end do
        
       end do
       
      end subroutine build_y_ur

      subroutine def_y_pairs(prs,refind)
c      *************************************************************
c      -- searches for pairs of equivalent atoms 
c         under reflection: y --> -y
c      -- initializes the reference array <refind>:
c         for a given atom 'iat' it returns an index nu=refind(iat), 
c         such that next subgroup of indices (nu, nu+1, ...) 
c         refers to internal degrees of freedom of the atom iat
c      *************************************************************
       
       integer, intent(out) :: prs(:), refind(:)
       
       integer          iat, jat, iatype, nu, norb       
       double precision x0,y0,z0, x1,y1,z1, dr
       logical          pfound
       double precision, parameter :: reps = 1.0e-3
       
       do iat = 1, num_atoms
       
        x0 = atom(iat)%pos(1)
        y0 = atom(iat)%pos(2)
        z0 = atom(iat)%pos(3)

        if (dabs(y0)<reps) then
c        counterpart is the same atom (since x=0) 	
	 prs(iat) = iat
        else
c        scanning all atoms, search for that one
c        with z1=z0, y1=-y0, and x1=x0
	 pfound = .false.
	 jat = 0
	 do while ((.not.pfound).and.(jat<num_atoms))
          jat = jat + 1
	  x1 = atom(jat)%pos(1)
	  y1 = atom(jat)%pos(2)
	  z1 = atom(jat)%pos(3)
	  dr = dabs(x1-x0)+dabs(y1+y0)+dabs(z1-z0)
	  if ( (dr < 3.0*reps).and.(jat.ne.iat)) then
	   pfound = .true.
	  end if
	 end do
         if (pfound.and.(atom(iat)%atype.eq.atom(jat)%atype)) then
c         counterpart for iat is jat >>>
	  prs(iat) = jat
	 else
c         no counterpart is found !
          write (*,'(/,a,i3,/)') 
     &	  '[SUBROUTINE def_y_pairs]: no pair atom is found for atom ',iat          
          stop 'transport module is terminated now'	 
	 end if
       
        end if
       end do

c      build up a reference index array <refind>
       nu = 1
       do iat = 1, num_atoms
        refind(iat) = nu
        iatype = atom(iat)%atype
	norb   = n_basis_func(iatype) 
        nu = nu + norb
       end do      
      
      end subroutine def_y_pairs

      subroutine y_updatedmat(dmatrix)
c      **************************************************************
c      in case of atomic structure with zx-mirror symmetry, this
c      routine takes care about symmetrization of the density matrix: 
c      y --> -y : left/spin-up   == right/spin-up 
c               & left/spin-down == right/spin-down
c      **************************************************************
c       input/output: density matrix (either non-eq or equilibrium one)  
        double precision, allocatable :: dmatrix(:,:,:)

        integer, allocatable :: atpairs(:), nuind(:)
	integer iat, jat, ispin, n, m, ierr        
        double precision, allocatable :: umat(:,:), udmat(:,:), 
     &	                                 dmatrix1(:,:,:), u1(:,:),
     &                                   dmat1(:,:)
        double precision nmtmp
        double precision, parameter :: ueps = 1.0d-10
	logical herm	

        allocate(atpairs(num_atoms),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE y_updatedmat]: <atpairs> allocation failure'
        end if       
        forall (iat=1:num_atoms) atpairs(iat)=0

        allocate(nuind(num_atoms),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE y_updatedmat]: <nuind> allocation failure'
        end if       
        forall (iat=1:num_atoms) nuind(iat)=0

        write(*,'(/,1x,a,/)') 
     &       'ZX-MIRROR SYMMETRIZED SOLUTION: searching for atom pairs y <--> -y)'
        call def_y_pairs(atpairs,nuind) 
        write(*,'(a)') '    atom     counterpart'
        write(*,'(a)') ' ------------------------'
      
        do iat = 1, num_atoms
         jat = atpairs(iat)
c         write(*,'(i5,a4,3x,a,i7,a4,i6)') 
c     &	   iat, atom(iat)%symbol(1:2),':', jat, atom(jat)%symbol(1:2),nuind(iat)
         write(*,'(i5,a4,3x,a,i7,a4)') 
     &	   iat, atom(iat)%symbol(1:2),':', jat, atom(jat)%symbol(1:2)
        end do
        write(*,'(a)') ' ------------------------'
    
        allocate(dmatrix1(nsaos,nsaos,nspin),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE y_updatedmat]: <dmatrix1> allocation failure'
        end if       
        dmatrix1 = 0.0d0
        
        allocate(dmat1(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE y_updatedmat]: <dmat1> allocation failure'
        end if       

        allocate(umat(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE y_updatedmat]: <umat> allocation failure'
        end if       
        umat = 0.0d0

        allocate(u1(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
	  stop
     &	  '[SUBROUTINE y_updatedmat]: <u1> allocation failure'
        end if       

        allocate(udmat(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE y_updatedmat]: <udmat> allocation failure'
        end if       
        
	print '(/,1x,a)', 'SYMMETRIZING DENSITY MATRIX >>>'

        call cpu_time(stime)       
c       initialize transformation matrix >>>
        print '(/,2x,a,$)', 'building up symmetry transform matrix U ...'
        call build_y_ur(atpairs,nuind,umat)
        print *, 'done'

c       checking that:  umat * umat^T = 1  
        call dble_ab('n','t',nsaos,umat,umat,u1)
 	herm = .true.
 	do n = 1, nsaos 
	 do m = 1, nsaos
	  if (m==n) then ; nmtmp = 1.0d0
	  else           ; nmtmp = 0.0d0
	  end if
	  herm = (dabs(u1(n,m)-nmtmp)<ueps) 
	 end do
	end do

        if (herm) then
	 print *, ' checking: matrix U is orthogonal'
        else 
	 print *, ' WARNING: matrix U is NOT orthogonal'
	end if

        print '(2x,a,$)', 'updating density matrix ...'
        do ispin = 1, nspin

c        updating a density matrix, to insure a symmetrical solution
c        -- dmat_(ispin)  <--  1/2 *(dmat_(ispin) + umat^T * dmat_(ispin) * umat)  
c        -- where 'ispin' is a spin index; <umat> is zy-plane reflection 

c        1. compute dmat1 = umat^T * dmatrix * umat        
         call dble_ab('t','n',nsaos,umat,dmatrix(:,:,ispin),udmat)
	 call dble_ab('n','n',nsaos,udmat,umat,dmat1)

c        2. admix different spin channels
         do m = 1, nsaos
	  do n = 1, nsaos
            dmatrix1(n,m,ispin) = 0.5d0*(dmatrix(n,m,ispin) + dmat1(n,m))	  
	  end do
	 end do

        end do ! ispin
	
c       put <dmatrix1> to <neqdmat>
	do ispin = 1, nspin
         do m = 1, nsaos
	  do n = 1, nsaos
           dmatrix(n,m,ispin) = dmatrix1(n,m,ispin) 	  
	  end do
	 end do
	end do  ! ispin
        print *, 'done'
	
        call cpu_time(ftime)
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts*60.0d0
        print '(2x,a,i3,a,f5.2,a)',
     &        'time spent:', mnts, ' min ', secs, ' sec'

	print '(/,1x,a)', '<<< DONE WITH SYMMETRIZATION '

        deallocate(atpairs,nuind,umat,udmat,dmat1,dmatrix1,u1,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE y_updatedmat]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
        end if

      end subroutine y_updatedmat

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine update_mo
c      **********************************************************
c      in case of domain-wall solution, this routine takes 
c      care about symmetrization of the mo expansion coefficients 
c      z --> -z : left/spin-up   == right/spin-down 
c               & left/spin-down == right/spin-up
c      **********************************************************
c       input/output: MO expansion coefficients  
c       double precision, allocatable :: coeffs(:,:,:)

        integer, allocatable :: atpairs(:), nuind(:)
	integer iat, jat, n, m, ierr        
        double precision, allocatable :: umat(:,:), u1(:,:), tmp_su(:,:)
        
        double precision nmtmp
        double precision, parameter :: ueps = 1.0d-10
	logical herm	

        allocate(atpairs(num_atoms),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE update_mo]: <atpairs> allocation failure'
        end if       
        forall (iat=1:num_atoms) atpairs(iat)=0

        allocate(nuind(num_atoms),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE update_mo]: <nuind> allocation failure'
        end if       
        forall (iat=1:num_atoms) nuind(iat)=0

        write(*,'(/,1x,a,/)') 
     &          'DOMAIN-WALL SOLUTION: searching for atom pairs (z <--> -z)'
        call def_z_pairs(atpairs,nuind) 
        write(*,'(a)') '    atom     counterpart'
        write(*,'(a)') ' ------------------------'
      
        do iat = 1, num_atoms
         jat = atpairs(iat)
c         write(*,'(i5,a4,3x,a,i7,a4,i6)') 
c     &	   iat, atom(iat)%symbol(1:2),':', jat, atom(jat)%symbol(1:2),nuind(iat)
         write(*,'(i5,a4,3x,a,i7,a4)') 
     &	   iat, atom(iat)%symbol(1:2),':', jat, atom(jat)%symbol(1:2)
        end do
        write(*,'(a)') ' ------------------------'

        allocate(umat(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE update_mo]: <umat> allocation failure'
        end if       
        umat = 0.0d0

        allocate(u1(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
	  stop
     &	  '[SUBROUTINE update_mo]: <u1> allocation failure'
        end if       

        print '(/,1x,a)', 'SYMMETRIZING MOLECULAR ORBITALS >>>'

        call cpu_time(stime)       
c       initialize transformation matrix >>>
        print '(/,2x,a,$)', 'building up symmetry transform matrix U ...'
        call build_z_ur(atpairs,nuind,umat)
        print *, 'done'

c       checking that:  umat * umat^T = 1  
        call dble_ab('n','t',nsaos,umat,umat,u1)
 	herm = .true.
 	do n = 1, nsaos 
	 do m = 1, nsaos
	  if (m==n) then ; nmtmp = 1.0d0
	  else           ; nmtmp = 0.0d0
	  end if
	  herm = (dabs(u1(n,m)-nmtmp)<ueps) 
	 end do
	end do

        if (herm) then
	 print *, ' checking: matrix U is orthogonal'
        else 
	 print *, ' WARNING: matrix U is NOT orthogonal'
	end if

c       here nspin == 2 ! 
        if (nspin.eq.1) 
     &   stop '[SUB. update_mo]: nspin==1 ! ... confused and will quit now!'

        print '(2x,a,$)', 'updating mo coefficients ...'

        allocate(tmp_su(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE update_mo]: <tmp_su> allocation failure'
        end if

c       transform matrix (testing!!!)
        tmp_su = 0.0d0
        call dble_ab('n','n',nsaos,smat12inv,umat,tmp_su)
        umat = 0.0d0
        call dble_ab('n','n',nsaos,tmp_su,smat12,umat)

        deallocate(tmp_su,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE update_mo]: ',
     &           'impossible to deallocate temporary array <tmp_su>'
         print *, 'nevertheless, proceed further ...'
        end if

        call dble_ab('n','n',nsaos,umat,mo_coeff(:,:,1),mo_coeff(:,:,2))
        forall (n=1:nsaos) mo_en(n,2) = mo_en(n,1)
        
        print *, 'done'

        call cpu_time(ftime)
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts*60.0d0
        print '(2x,a,i3,a,f5.2,a)',
     &        'time spent:', mnts, ' min ', secs, ' sec'

	print '(/,1x,a)', '<<< DONE WITH SYMMETRIZATION '

        deallocate(atpairs,nuind,umat,u1,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE update_mo]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
        end if

      end subroutine update_mo

      subroutine symmetrize_mo(scall)
c      **********************************************************
c      in case of domain-wall solution, this routine takes 
c      care about symmetrization of the mo expansion coefficients 
c      **********************************************************
c       input/output: MO expansion coefficients  
c       double precision, allocatable :: coeffs(:,:,:)
        character, intent (in) :: scall

        integer, allocatable :: atpairs(:), nuind(:)
	integer iat, jat, n, m, iorb, p, ispin, jspin, ierr        
        double precision, allocatable :: umat(:,:), u1(:,:), 
     &                    tmp_su(:,:), tmp_mo(:,:,:), tmp_coeff(:,:,:)
        
        double precision nmtmp
        double precision, parameter :: ueps = 1.0d-10
	logical herm	

        allocate(atpairs(num_atoms),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE symmetrize_mo]: <atpairs> allocation failure'
        end if       
        forall (iat=1:num_atoms) atpairs(iat)=0

        allocate(nuind(num_atoms),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE symmetrize_mo]: <nuind> allocation failure'
        end if       
        forall (iat=1:num_atoms) nuind(iat)=0

        select case (scall)
         case('x')
          write(*,'(/,1x,a,/)') 
     &            'SYMMETRIZING MOs COEFFICIENTS: searching for atom pairs (x <--> -x)'
          call def_x_pairs(atpairs,nuind) 

         case('y')
          write(*,'(/,1x,a,/)') 
     &            'SYMMETRIZING MOs COEFFICIENTS: searching for atom pairs (y <--> -y)'
          call def_y_pairs(atpairs,nuind) 

         case('z')
          write(*,'(/,1x,a,/)') 
     &            'SYMMETRIZING MOs COEFFICIENTS: searching for atom pairs (z <--> -z)'
          call def_z_pairs(atpairs,nuind) 
    
        end select
      
        write(*,'(a)') '    atom     counterpart'
        write(*,'(a)') ' ------------------------'
      
        do iat = 1, num_atoms
         jat = atpairs(iat)
c         write(*,'(i5,a4,3x,a,i7,a4,i6)') 
c     &	   iat, atom(iat)%symbol(1:2),':', jat, atom(jat)%symbol(1:2),nuind(iat)
         write(*,'(i5,a4,3x,a,i7,a4)') 
     &	   iat, atom(iat)%symbol(1:2),':', jat, atom(jat)%symbol(1:2)
        end do
        write(*,'(a)') ' ------------------------'

        allocate(umat(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
 	 stop
     &	  '[SUBROUTINE symmetrize_mo]: <umat> allocation failure'
        end if       
        umat = 0.0d0

        allocate(u1(nsaos,nsaos),stat=ierr) 
        if (ierr.ne.0) then
         print *
	  stop
     &	  '[SUBROUTINE symmetrize_mo]: <u1> allocation failure'
        end if       

        print '(/,1x,a)', 'SYMMETRIZING MOLECULAR ORBITALS >>>'

        call cpu_time(stime)       
c       initialize transformation matrix >>>
        print '(/,2x,a,$)', 'building up symmetry transform matrix U ...'
        select case (scall)
         case('x')
          call build_x_ur(atpairs,nuind,umat)
         case('y')
          call build_y_ur(atpairs,nuind,umat)
         case('z')
          call build_z_ur(atpairs,nuind,umat)
        end select
        print *, 'done'

c       checking that:  umat * umat^T = 1  
        call dble_ab('n','t',nsaos,umat,umat,u1)
 	herm = .true.
 	do n = 1, nsaos 
	 do m = 1, nsaos
	  if (m==n) then ; nmtmp = 1.0d0
	  else           ; nmtmp = 0.0d0
	  end if
	  herm = (dabs(u1(n,m)-nmtmp)<ueps) 
	 end do
	end do

        if (herm) then
	 print *, ' checking: matrix U is orthogonal'
        else 
	 print *, ' WARNING: matrix U is NOT orthogonal'
	end if

c       here nspin == 2 ! 
        if ((nspin.eq.1).and.(z_symmetrize).and.(ap_symmetrize)) 
     &   stop 
     &    '[SUB. symmetrize_mo]: nspin==1 but AP symmetrization is required ! ... i am confused'

        print '(2x,a,$)', 'updating mo coefficients ...'

        allocate(tmp_su(nsaos,nsaos),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE symmetrize_mo]: <tmp_su> allocation failure'
        end if

c       transform matrix (testing!!!)
        tmp_su = 0.0d0
        call dble_ab('n','n',nsaos,smat12inv,umat,tmp_su)
        umat = 0.0d0
        call dble_ab('n','n',nsaos,tmp_su,smat12,umat)

        deallocate(tmp_su,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE symmetrize_mo]: ',
     &           'impossible to deallocate temporary array <tmp_su>'
         print *, 'nevertheless, proceed further ...'
        end if

        allocate(tmp_mo(nsaos,nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
         print *
         stop '[SUBROUTINE symmetrize_mo]: <tmp_mo> allocation failure'
        end if
        tmp_mo = 0.0d0
  
        select case (scall)
         case('x')
          do ispin = 1, nspin
           call dble_ab('n','n',nsaos,umat,mo_coeff(:,:,ispin),tmp_mo(:,:,ispin))
           do p = 1, nsaos
            do iorb = 1, nsaos
             mo_coeff(iorb,p,ispin) = 0.5d0*(mo_coeff(iorb,p,ispin)+tmp_mo(iorb,p,ispin))
            end do
           end do    
          end do ! ispin

         case('y')
          do ispin = 1, nspin
           call dble_ab('n','n',nsaos,umat,mo_coeff(:,:,ispin),tmp_mo(:,:,ispin))
           do p = 1, nsaos
            do iorb = 1, nsaos
             mo_coeff(iorb,p,ispin) = 0.5d0*(mo_coeff(iorb,p,ispin)+tmp_mo(iorb,p,ispin))
            end do
           end do    
          end do ! ispin

         case('z')
          if (nspin.eq.1) then

            call dble_ab('n','n',nsaos,umat,mo_coeff(:,:,1),tmp_mo(:,:,1))
            do p = 1, nsaos
             do iorb = 1, nsaos
              mo_coeff(iorb,p,ispin) = 0.5d0*(mo_coeff(iorb,p,ispin)+tmp_mo(iorb,p,ispin))
             end do
            end do    

          else ! nspin == 2

           allocate(tmp_coeff(nsaos,nsaos,2),stat=ierr)
           if (ierr.ne.0) then
           print *
            stop '[SUBROUTINE symmetrize_mo]: <tmp_coeff> allocation failure'
           end if
           tmp_coeff = 0.0d0

           do ispin = 1, nspin
            if (ap_symmetrize) then
             if (ispin.eq.1) then ; jspin = 2
  	     else                 ; jspin = 1
 	     end if
 	    else if (p_symmetrize) then
             if (ispin.eq.1) then ; jspin = 1
 	     else                 ; jspin = 2
 	     end if
 	    else
 	     stop '[SUB. symmetrize_mo]: i am confused, neither p_symmetrize nor ap_symmetrize are active!'
 	    end if
            call dble_ab('n','n',nsaos,umat,mo_coeff(:,:,ispin),tmp_mo(:,:,jspin))
           end do ! ispin
    
           do p = 1, nsaos
             do iorb = 1, nsaos
              tmp_coeff(iorb,p,1) = 0.5d0*(mo_coeff(iorb,p,1)+tmp_mo(iorb,p,2))
              tmp_coeff(iorb,p,2) = 0.5d0*(mo_coeff(iorb,p,2)+tmp_mo(iorb,p,1))
             end do
           end do    

           forall(p=1:nsaos,iorb=1:nsaos) 
            mo_coeff(iorb,p,1) = tmp_coeff(iorb,p,1)
            mo_coeff(iorb,p,2) = tmp_coeff(iorb,p,2)
           end forall

           deallocate(tmp_coeff,stat=ierr)
           if (ierr.ne.0) then
            print '(/,a,a)', ' [SUBROUTINE symmetrize_mo]: ',
     &                       'impossible to deallocate temporary array <tmp_coeff>'
            print *, 'nevertheless, proceed further ...'
           end if

         end if ! nspin
    
        end select

        print *, 'done'

        call cpu_time(ftime)
        mnts = int((ftime-stime)/60.0d0,4)
        secs = ftime - stime - mnts*60.0d0
        print '(2x,a,i3,a,f5.2,a)',
     &        'time spent:', mnts, ' min ', secs, ' sec'

	print '(/,1x,a)', '<<< DONE WITH SYMMETRIZATION '

        deallocate(atpairs,nuind,umat,u1,tmp_mo,stat=ierr)
        if (ierr.ne.0) then
         print '(/,a,a)', ' [SUBROUTINE symmetrize_mo]: ',
     &           'impossible to deallocate temporary arrays'
         print *, 'nevertheless, proceed further ...'
        end if

      end subroutine symmetrize_mo

      end module domainwall

