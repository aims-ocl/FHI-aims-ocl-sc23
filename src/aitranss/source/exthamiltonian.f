c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           may 2008
c     last revision:  jan 2012
c###########################################################

      module exthamiltonian
c     ************************************      
c     computes complex eigen-energies and 
c     eigenvectors of extended hamiltonian 
c     ************************************      

c      use globalvars
c      use tools
       use cmatrix
       implicit none

c      double precision, private, parameter :: myprec = 1.0d-9    
c      real(4), private ::  stime, ftime, secs
c      integer, private ::  mnts  
 
      contains
       
      subroutine hfull
c      *********************************************************
c      builds up an "extended hamiltonian": h_full = h0 + sigma
c      *********************************************************
       implicit none
       integer n, m, ispin, ierr 
       logical allochmat

       allochmat = allocated(hmat)
       if (.not.allochmat) then
        allocate(hmat(nsaos,nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
        print *
        stop
     &    '[SUBROUTINE hfull]: <hmat> allocation failure'
        end if
       end if	

       print '(/,a,$)', ' INITIALIZING H+Sigma ...'

       if (no_leads) then
        forall (ispin=1:nspin,m=1:nsaos,n=1:nsaos)
     $   hmat(n,m,ispin) = cmplx(h0mat(n,m,ispin),0.0d0,8)
        forall (ispin=1:nspin,n=1:nsaos)
     $   hmat(n,n,ispin) = cmplx(h0mat(n,n,ispin),-eta,8)
       else
        forall (ispin=1:nspin,m=1:nsaos,n=1:nsaos)
     $   hmat(n,m,ispin) = cmplx(h0mat(n,m,ispin),0.0d0,8)
        do ispin = 1, nspin
         do n = 1, nsaos
          if (.not.sp_leads) then 
c           non-spin-polarized electrodes
	    hmat(n,n,ispin) = hmat(n,n,ispin) + sigmann(n)
	  else 
c           spin-polarized electrodes
	    if (ispin.eq.1) then 
	      hmat(n,n,1) = hmat(n,n,1) + sigmann(n)
	    else
	      hmat(n,n,2) = hmat(n,n,2) + bsigmann(n)
	    end if     
	  end if   ! sp_leads
         end do  ! nsaos
        end do  ! ispin
       
       end if ! no_leads
       
       print *, 'DONE'

      end subroutine hfull


      subroutine hspectrum
c      *****************************************
c      computes spectrum of the operator h+sigma
c      *****************************************
       implicit none
       integer  ierr
       
       if (.not.allocated(zpoles)) then
        allocate(zpoles(nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE hspectrum]: <zpoles> allocation failure'
        end if
       end if
       
       if (.not.allocated(heigvec)) then
        allocate(heigvec(nsaos,nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
        print *
        stop
     &    '[SUBROUTINE hspectrum]: <heigvec> allocation failure'
        end if
       end if	

       if (.not.allocated(invheigvec)) then
        allocate(invheigvec(nsaos,nsaos,nspin),stat=ierr)
        if (ierr.ne.0) then
        print *
        stop
     &   '[SUBROUTINE hspectrum]: <heigvec> allocation failure'
        end if
       end if	

       print '(/,a)', 
     &       ' COMPUTING SPECTRUM OF THE OPERATOR H+Sigma >>>'

       call cmplxspectrum(hmat,zpoles,heigvec,invheigvec,nsaos,1)

       print '(/,a)', ' <<< DONE WITH SPECTRUM OF H+Sigma'
      
      end subroutine hspectrum

      end module exthamiltonian       

