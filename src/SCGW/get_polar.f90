!****s
!* FHI-aims/evaluate_polarisability_freq
!  NAME
!   evaluate_polarisability_freq
!  SYNOPSIS

      subroutine get_polar &
         ( green_fn_time, tau, &
           ntau, wtau, polar_time)

!  PURPOSE
!  Subroutine evaluate_polarisability_freq  evaluates the non-interacting 
!  polarisability, represented within auxiliary basis.  Here imaginary 
!  frequency domain is used. The polarizability is evaluated using the formula:
!
!  Chi^0_mn (t) =  ________________________________
! 
!  It is used the trick of the partial summation as for the Self Energy
!
!  USES

      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use physics
      use hartree_fock
      use localorb_io, only: use_unit

      implicit none

!  ARGUMENTS
      integer    ::  ntau 
!      real*8     ::  ovlp_3fn(n_basis_pairs,n_loc_prodbas)  
      real*8     ::  tau (-ntau:ntau)
      real*8     ::  polar_time(n_basbas, n_loc_prodbas, 0:ntau)
      real*8     ::  green_fn_time (n_basis,n_basis,-ntau:ntau,n_spin)
      real*8     ::  wtau(-ntau:ntau)

!
! - ntau               number of points of the time grid
! - ovlp_3fn             is the overlap integrals between 2 NAO basis and a product basis fn
! - green_fn_time        the Greeen's Function in the NAO basis and time domain
! - wtau                 -- integration weights of the time grid
! - nomega             number of points of the frequency grids
! - tau                time points
! - polar_time         the polarizability expressed in the product basis and time domain

!INTERNAL
      real*8, dimension (:,:,:) , allocatable  ::  partial_sum_loc_t_plus
      real*8, dimension (:,:) , allocatable  ::  partial_sum_t_min
      real*8, dimension (:,:,:) , allocatable  ::  partial_sum_loc_t_min
      real*8, dimension (:,:,:) , allocatable  ::  aux_ovlp_3fn ! ovlp on the current processor
!      real*8, dimension (:,:,:) , allocatable  ::  tmp_ovlp_3fn ! this is the global ovlp3fn
      integer n_first_min

!     counters

      integer :: i_tau

      integer :: j_state
      integer :: k_state

      integer :: i_basbas
      integer :: j_basbas

      integer :: i_basis
      integer :: j_basis
      integer :: k_basis
      integer :: m_basis
      logical output
 
      integer :: i_index
      integer :: i_spin
      integer :: i_task

 
      real*8 chi1, chi2,t1,t2

!     begin work
!write(use_unit,*) 'ntau', ntau
!write(use_unit,*) 'size au', size(tau)
!write(use_unit,*) 'size wtau', size(wtau)
      if(myid.eq.0)then
        write(use_unit,*)"  --- Evaluating the Polarizability "
      endif
!      polar_time(:,:,:)=0.d0

      if(.not. allocated (partial_sum_loc_t_plus) )then
         allocate (partial_sum_loc_t_plus(n_basis, n_basis, n_loc_prodbas))
      endif
      if(.not. allocated (partial_sum_loc_t_min) )then
         allocate (partial_sum_loc_t_min(n_basis, n_basis, n_loc_prodbas))
      endif
      if(.not. allocated (partial_sum_t_min) )then
         allocate (partial_sum_t_min( n_basis, n_basbas))
      endif
      if(.not. allocated( aux_ovlp_3fn ))then
         allocate (aux_ovlp_3fn (n_basis, n_basis, n_loc_prodbas))
      endif

!stop

! first construct the matrix O^mu_ij = \int P^mu(r) \phi_i(r) \phi_j(r) dr
      aux_ovlp_3fn (:,:,:) = 0.d0
!      tmp_ovlp_3fn (:,:,:) = 0.d0
     
     !unpack the matrix
      do i_basis =1 , n_basis, 1
        do j_basis =1, n_basis, 1
          if(basis_nghr(i_basis,j_basis).gt.0)then
            aux_ovlp_3fn (i_basis, j_basis,:) = &
            ovlp_3fn(basis_nghr(i_basis,j_basis),:)
          else
            aux_ovlp_3fn (i_basis, j_basis,:) = 0.d0
          endif
        enddo
      enddo
!stop
      polar_time (:,:,:) = 0.d0

       do i_spin = 1, n_spin
        do i_tau = 1, ntau, 1
!     write(use_unit,*)'i_spin,itau',i_spin ,i_tau
          partial_sum_loc_t_plus(:,:,:) = 0.d0
          partial_sum_loc_t_min(:,:,:) = 0.d0
  
          do i_basbas = 1, n_loc_prodbas, 1
            !product of Overlap_3fn x Green's function (s overlap matrix already in G!)
             call dgemm ('T','N', n_basis,n_basis, n_basis, 1.d0, &
               aux_ovlp_3fn(:,:, i_basbas), n_basis, &
               green_fn_time(:,:,i_tau,i_spin), n_basis, 1.d0, &
               partial_sum_loc_t_plus (:,:, i_basbas), n_basis )
     
             call dgemm ('T','N', n_basis,n_basis, n_basis, 1.d0, &
               aux_ovlp_3fn(:,:, i_basbas), n_basis, &
               green_fn_time(:,:,-i_tau,i_spin), n_basis, 1.d0, &
               partial_sum_loc_t_min (:,:, i_basbas), n_basis )
          enddo
          do j_basis =1, n_basis ,1 
  
             partial_sum_t_min(:,:) = 0.d0
             do i_task = 1, n_tasks, 1
              do i_basbas = 1, n_loc_prodbas
                 i_index = map_prodbas(i_basbas,i_task)
                 if(i_index.gt.0 .and. myid.eq.i_task-1)then
                   partial_sum_t_min (:,i_index) = &
                   partial_sum_loc_t_min (j_basis,:,i_basbas)
                endif
               enddo
             enddo 
  
             call sync_matrix(&
               partial_sum_t_min (1:n_basis,1:n_basbas), &
               n_basis, n_basbas)
  
              if(n_loc_prodbas .gt.0)then
                call dgemm ('T','N', n_basbas, n_loc_prodbas, n_basis,&
                  1.d0, partial_sum_t_min ,  n_basis,&
                  partial_sum_loc_t_plus (:,j_basis,: ), n_basis, 1.d0, &
                  polar_time (:,:,i_tau), n_basbas )
              endif
  
          enddo ! i_basis
        enddo !i_tau
      enddo
 
       do i_basbas = 1, n_basbas, 1 
        do j_basbas = 1, n_loc_prodbas,1 
           chi1 = polar_time (i_basbas,j_basbas,1) 
           chi2 = polar_time (i_basbas,j_basbas,2)
           t1 = tau(1)  
           t2 = tau(2)
           polar_time (i_basbas,j_basbas,0) = chi1-(chi2-chi1)/(t2-t1)*t1
        enddo
       enddo

        output = .true.
        if (myid.eq.0 .and. output)then
          open (22,file="polar_green.dat")   
            i_basbas = 1
            j_basbas = 1
            do i_tau = 0, ntau, 1
!i             do i_basbas =1, n_loc_prodbas, 5
              write(22,*) tau(i_tau), polar_time (i_basbas,j_basbas, i_tau), &
                polar_time (i_basbas,1, i_tau)
!             enddo
            enddo
          close(22)
        endif

      if( allocated (partial_sum_t_min) )then
         deallocate (partial_sum_t_min)
      endif
      if( allocated (partial_sum_loc_t_plus) )then
         deallocate (partial_sum_loc_t_plus)
      endif
      if( allocated (partial_sum_loc_t_min) )then
         deallocate (partial_sum_loc_t_min)
      endif
      if(allocated (aux_ovlp_3fn) )then
         deallocate (aux_ovlp_3fn)
      endif
      return
      end subroutine get_polar
!---------------------------------------------------------------------
!**********
!****s
!* FHI-aims/evaluate_polarisability_freq
!  NAME
!   evaluate_polarisability_freq
!  SYNOPSIS

      subroutine get_polar_old &
         ( green_fn_time, tau, &
           ntau, wtau, polar_time)

!  PURPOSE
!  Subroutine evaluate_polarisability_freq  evaluates the non-interacting 
!  polarisability, represented within auxiliary basis.  Here imaginary 
!  frequency domain is used. The polarizability is evaluated using the formula:
!
!  Chi^0_mn (t) =  ________________________________
! 
!  It is used the trick of the partial summation as for the Self Energy
!
!  USES

      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use physics
      use hartree_fock
      use localorb_io, only: use_unit

      implicit none

!  ARGUMENTS
      integer    ::  ntau 
!      real*8     ::  ovlp_3fn(n_basis_pairs,n_loc_prodbas)  
      real*8     ::  tau (-ntau:ntau)
      real*8     ::  polar_time(n_basbas, n_loc_prodbas, 0:ntau)
      real*8     ::  green_fn_time (n_basis,n_basis,-ntau:ntau,n_spin)
      real*8     ::  wtau(-ntau:ntau)

!
! - ntau               number of points of the time grid
! - ovlp_3fn             is the overlap integrals between 2 NAO basis and a product basis fn
! - green_fn_time        the Greeen's Function in the NAO basis and time domain
! - wtau                 -- integration weights of the time grid
! - nomega             number of points of the frequency grids
! - tau                time points
! - polar_time         the polarizability expressed in the product basis and time domain

!INTERNAL
!      real*8, dimension (:,:,:) , allocatable  ::  partial_sum_t_plus
      real*8, dimension (:,:,:) , allocatable  ::  partial_sum_loc_t_plus
      real*8, dimension (:,:,:) , allocatable  ::  partial_sum_t_min
!      real*8, dimension (:,:,:) , allocatable  ::  partial_sum_loc_t_min
      real*8, dimension (:,:,:) , allocatable  ::  aux_ovlp_3fn ! ovlp on the current processor
      real*8, dimension (:,:,:) , allocatable  ::  tmp_ovlp_3fn ! this is the global ovlp3fn
      integer n_first_min

!     counters

      integer :: i_tau

      integer :: j_state
      integer :: k_state

      integer :: i_basbas
      integer :: j_basbas

      integer :: i_basis
      integer :: j_basis
      integer :: k_basis
      integer :: m_basis
      logical output
 
      integer :: i_index
      integer :: i_spin
      integer :: i_task

      real*8 aux_part_sum (n_basis, n_basis, n_loc_prodbas)
 
      real*8 chi1, chi2,t1,t2

!     begin work

      if(myid.eq.0)then
        write(use_unit,*)"  --- Evaluating the Polarizability "
      endif
      polar_time(:,:,:)=0.d0

      if(.not. allocated (partial_sum_loc_t_plus) )then
         allocate (partial_sum_loc_t_plus(n_basis, n_basis, n_loc_prodbas))
      endif
!      if(.not. allocated (partial_sum_loc_t_min) )then
!         allocate (partial_sum_loc_t_min(n_basis, n_basis, n_loc_prodbas))
!      endif
!      if(.not. allocated (partial_sum_t_plus) )then
!         allocate (partial_sum_t_plus(n_basis, n_basis, n_basbas))
!      endif
      if(.not. allocated (partial_sum_t_min) )then
         allocate (partial_sum_t_min(n_basis, n_basis, n_basbas))
      endif
      if(.not. allocated( tmp_ovlp_3fn ))then
         allocate (tmp_ovlp_3fn (n_basis, n_basis, n_basbas))
      endif
      if(.not. allocated( aux_ovlp_3fn ))then
         allocate (aux_ovlp_3fn (n_basis, n_basis, n_loc_prodbas))
      endif


! first construct the matrix O^mu_ij = \int P^mu(r) \phi_i(r) \phi_j(r) dr
      aux_ovlp_3fn (:,:,:) = 0.d0
      tmp_ovlp_3fn (:,:,:) = 0.d0
     
     !unpack the matrix
      do i_basis =1 , n_basis, 1
        do j_basis =1, n_basis, 1
          if(basis_nghr(i_basis,j_basis).gt.0)then
            aux_ovlp_3fn (i_basis, j_basis,:) = &
            ovlp_3fn(basis_nghr(i_basis,j_basis),:)
          else
            aux_ovlp_3fn (i_basis, j_basis,:) = 0.d0
          endif
        enddo
      enddo

!upfolding 
      do i_task = 1, n_tasks
         do i_basbas = 1, n_loc_prodbas
           i_index = map_prodbas(i_basbas,i_task)
           if(i_index.gt.0 .and. myid.eq.i_task-1) then
             tmp_ovlp_3fn(:,:,i_index) = &
             aux_ovlp_3fn(:,:,i_basbas)
           endif
         enddo
      enddo

! synchronizing the upfolded matrix
      do i_basis = 1, n_basis, 1
        call sync_matrix(&
          tmp_ovlp_3fn (i_basis,1:n_basis,1:n_basbas), &
          n_basis, n_basbas)
      enddo! i_basis

         
      polar_time (:,:,:) = 0.d0

     do i_spin = 1, n_spin
      do i_tau = 1, ntau, 1

        partial_sum_loc_t_plus(:,:,:) = 0.d0
        partial_sum_t_min(:,:,:) = 0.d0
        
        do i_basbas = 1, n_basbas, 1

         !partial sum on the current processor
              if(i_basbas.le.n_loc_prodbas)then
                 call dgemm ('T','N', n_basis,n_basis, n_basis, 1.d0, &
                   aux_ovlp_3fn(:,:, i_basbas), n_basis, &
                   green_fn_time(:,:,i_tau,i_spin), n_basis, 1.d0, &
                   partial_sum_loc_t_plus (:,:, i_basbas), n_basis )
              endif
   
          !global partial sum
              call dgemm ('T','N', n_basis,n_basis, n_basis, 1.d0, &
                tmp_ovlp_3fn(:,:, i_basbas), n_basis, &
                green_fn_time(:,:,-i_tau,i_spin), n_basis, 1.d0, &
                partial_sum_t_min (:,:, i_basbas), n_basis )
        enddo 

        do j_basis = 1, n_basis, 1
            call dgemm ('T','N', n_basbas, n_loc_prodbas, n_basis,&
              1.d0, partial_sum_t_min (j_basis, :, : ),  n_basis,&
              partial_sum_loc_t_plus (:,j_basis,: ), n_basis, 1.d0, &
              polar_time (:,:,i_tau), n_basbas )

!            if (i_tau.eq.1) then
!             !for the polarizability at t=0
!              call dgemm ('T','N', n_basbas, n_loc_prodbas, n_basis,&
!              1.d0, aux_part_sum (j_basis, :,:),  n_basis,&
!              partial_sum_loc_t_plus (:,j_basis,: ), n_basis, 1.d0, &
!              polar_time (:,:,0), n_basbas )
!           endif

        enddo ! i_basis
      enddo !i_tau

! t=0 treated separately
        partial_sum_loc_t_plus(:,:,:) = 0.d0
        partial_sum_t_min(:,:,:) = 0.d0

        do i_basbas = 1, n_basbas, 1
      !partial sum on the current processor
           if(i_basbas.le.n_loc_prodbas)then
              call dgemm ('T','N', n_basis,n_basis, n_basis, 1.d0, &
                aux_ovlp_3fn(:,:, i_basbas), n_basis, &
                green_fn_time(:,:,0,i_spin), n_basis, 1.d0, &
                partial_sum_loc_t_plus (:,:, i_basbas), n_basis )
           endif

       !global partial sum
           call dgemm ('T','N', n_basis,n_basis, n_basis, 1.d0, &
             tmp_ovlp_3fn(:,:, i_basbas), n_basis, &
             green_fn_time(:,:,1,i_spin), n_basis, 1.d0, &
             partial_sum_t_min (:,:, i_basbas), n_basis )
        enddo
      enddo !spin
  
        do j_basis = 1, n_basis, 1
            call dgemm ('T','N', n_basbas, n_loc_prodbas, n_basis,&
              1.d0, partial_sum_t_min (j_basis, :, : ),  n_basis,&
              partial_sum_loc_t_plus (:,j_basis,: ), n_basis, 1.d0, &
              polar_time (:,:,0), n_basbas )

        enddo ! i_basis

       !test
!       do i_basbas = 1, n_basbas, 1 
!        do j_basbas = 1, n_loc_prodbas,1 
!           chi1 = polar_time (i_basbas,j_basbas,1) 
!           chi2 = polar_time (i_basbas,j_basbas,2)
!           t1 = tau(1)  
!           t2 = tau(2)
!           polar_time (i_basbas,j_basbas,0) = chi1-(chi2-chi1)/(t2-t1)*t1
!        enddo
!       enddo


      output = .false.
      if (myid.eq.0 .and. output)then
        open (22,file="polar_green.dat")   
          i_basbas = 1
          j_basbas = 1
          do i_tau = 0, ntau, 1
!           do i_basbas =1, n_loc_prodbas, 5
            write(22,*) tau(i_tau), polar_time (i_basbas,j_basbas, i_tau), &
              polar_time (i_basbas,1, i_tau)
!           enddo
          enddo
        close(22)
      endif

!      if( allocated (partial_sum_t_plus) )then
!         deallocate (partial_sum_t_plus)
!      endif
      if( allocated (partial_sum_t_min) )then
         deallocate (partial_sum_t_min)
      endif
      if( allocated (partial_sum_loc_t_plus) )then
         deallocate (partial_sum_loc_t_plus)
      endif
!      if( allocated (partial_sum_loc_t_min) )then
!         deallocate (partial_sum_loc_t_min)
!      endif
      if(allocated (aux_ovlp_3fn) )then
         deallocate (aux_ovlp_3fn)
      endif
      return
      end subroutine get_polar_old
!---------------------------------------------------------------------
!**********
