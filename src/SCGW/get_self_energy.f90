
!**********
!****s
!* FHI-aims/get_self_energy
!  NAME
!   get_self_energy
!  SYNOPSIS

      subroutine get_self_energy &
         ( green_fn_time,screened_coul_int,&
           self_energy_time)


!  PURPOSE
!
! Evaluate the GW self energy in the time domain to avoid the convolution
! in frequency, G is expressed in the NAO basis, W is in the product basis
!
! Sigma_ij = Sum_lkmn  O^m_li O^n_kj G_lk (t) W_mn(t)
!
! lk refers to NAO basis, mn to the product basis
! O^m_li are overlap integrals between two NAO basis function and a product basis fn
!  and are given by  ovlp_3fn
! To improve scaling, first a partial summation is calculated
! 
!  A^m_ik (t) = Sum_l O^m_li G_lk (t)  
!  B^m_kj (t) = Sum_n O^n_kj W_mn(t)
!
! in this way the self energy is calculated as
!
! Sigma_ij (t) = Sum_mk A^m_ik (t)  B^m_kj (t)
!
!
!  USES

!      use physics
      use constants, only: pi
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use hartree_fock
      use scgw_grid
 
      implicit none

!  ARGUMENTS
      real*8     ::  green_fn_time (n_basis,n_basis,-ntau:ntau)
      real*8     ::  screened_coul_int (n_basbas, n_loc_prodbas, 0:ntau)
      real*8 ::  self_energy_time (n_basis,n_basis, -ntau:ntau)

! INPUTS
!
! o  ntau                 -- is the number of points in the time grid
! o  tau                  -- is the time grid
! o  ovlp_3fn             is the overlap integrals between 2 NAO basis and a product basis fn
! o  green_fn_time        the Greeen's Function in the NAO basis and time domain
! o  screened_coul_int    The screened Coulomb Interaction in the Product Basis and Time Domain, 
!                             which is symmetric in the time axis


!INTERNAL
      real*8, dimension (:,:) , allocatable  ::  partial_sum_G
      real*8, dimension (:,:,:) , allocatable  ::  partial_sum_W
      real*8, dimension (:,:,:) , allocatable  ::  aux_ovlp_3fn
      real*8, dimension (:,:) , allocatable  ::  tmp_ovlp_3fn
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
 
      integer :: i_index
      integer :: i_spin
      integer :: i_task

      logical :: output
!     begin work

!      if(myid.eq.0)then
!        write(use_unit,*)"  --- Evaluating the Self Energy "
!      endif

      if(.not. allocated (partial_sum_W) )then
         allocate (partial_sum_W(n_basis, n_basis, n_loc_prodbas))
      endif
      if(.not. allocated (partial_sum_G) )then
         allocate (partial_sum_G(n_basis, n_basis))
      endif
      if(.not. allocated( tmp_ovlp_3fn ))then
         allocate (tmp_ovlp_3fn (n_basis, n_basbas)) !global
      endif
      if(.not. allocated( aux_ovlp_3fn ))then
         allocate (aux_ovlp_3fn (n_basis, n_basis, n_loc_prodbas))!local
      endif

      self_energy_time (:,:,:) = 0.d0
      aux_ovlp_3fn (:,:,:) = 0.d0

! filling the aux_ovlp_3fn matrix, (in this way we can use dgemm)

      do j_basis =1 , n_basis, 1
        do k_basis =1, n_basis, 1
          if(.not.basis_nghr(j_basis,k_basis).eq.0)then
           aux_ovlp_3fn (j_basis, k_basis,:) = &
           ovlp_3fn(basis_nghr(j_basis,k_basis),:)
          endif
        enddo
      enddo

!upfolding

!synchronizing
!      do i_basis = 1, n_basis, 1
!        call sync_matrix(&
!          tmp_ovlp_3fn (i_basis,1:n_basis,1:n_basbas), &
!          n_basis, n_basbas)
!      enddo! i_basis
 
      do i_tau = 0, ntau,1
        partial_sum_W (:,:,:)  = 0.d0

        !evaluate partial_sum_W
        do j_basis = 1, n_basis, 1

           tmp_ovlp_3fn (:,:) = 0.d0
           do i_task = 1, n_tasks
             do i_basbas = 1, n_loc_prodbas
               i_index = map_prodbas(i_basbas,i_task)
               if(i_index.gt.0 .and. myid.eq.i_task-1) then
                 tmp_ovlp_3fn(:,i_index) = &
                  aux_ovlp_3fn(j_basis,:,i_basbas)
               endif
             enddo
           enddo
 
           call sync_matrix (tmp_ovlp_3fn, n_basis, n_basbas)           

           if(n_loc_prodbas.gt.0)then
             call dgemm ('N','N', n_basis,n_loc_prodbas,n_basbas,1.d0,&
                tmp_ovlp_3fn, n_basis, &
                screened_coul_int(:,:, i_tau), n_basbas,1.d0, &
                partial_sum_W (j_basis,:,:), n_basis )
           endif            
        enddo

        do i_basbas = 1, n_loc_prodbas, 1

           partial_sum_G (:,:)  = 0.d0
                !this is only for positive times
           call dgemm ('N','N', n_basis,n_basis,n_basis,1.d0,&
               aux_ovlp_3fn (:,:,i_basbas), n_basis, &
               green_fn_time (:,:, i_tau), n_basis,1.d0, &
               partial_sum_G (:,:),n_basis )

           call dgemm ('N','N', n_basis,n_basis,n_basis,1.d0,&
                partial_sum_G (:,:), n_basis, &
                partial_sum_W (:,:,i_basbas), n_basis,1.d0, &
                self_energy_time (:,:,i_tau),n_basis)
        enddo

        if (i_tau.gt.0)then !this is for negative time
          do i_basbas = 1, n_loc_prodbas, 1
             partial_sum_G (:,:)  = 0.d0
             call dgemm ('N','N', n_basis,n_basis,n_basis,1.d0,&
               aux_ovlp_3fn (:,:,i_basbas), n_basis, &
               green_fn_time (:,:, -i_tau), n_basis,1.d0, &
               partial_sum_G (:,:),n_basis )

             call dgemm ('N','N', n_basis,n_basis,n_basis,1.d0,&
                partial_sum_G (:,:), n_basis, &
                partial_sum_W (:,:,i_basbas), n_basis,1.d0, &
                self_energy_time (:,:,-i_tau),n_basis)            
          enddo
        endif

      enddo


      do i_tau = -ntau, ntau, 1
        call sync_matrix(&
          self_energy_time (1:n_basis,1:n_basis,i_tau), &
          n_basis, n_basis)
 
      enddo

      !this set the usual prefactor
      self_energy_time (:,:,:) = -self_energy_time (:,:,:)/2.d0/pi

      output = .false.
      if (output)then
       if(myid.eq.0)then      
         open(55,file= "self_energy_time.dat")
         do i_tau = -ntau, ntau, 1
           write(55,*)tau(i_tau), self_energy_time(1,1,i_tau)
         enddo 
         close(55)
       endif
      endif


      if( allocated (partial_sum_G) )then
         deallocate (partial_sum_G)
      endif
      if( allocated (partial_sum_W) )then
         deallocate (partial_sum_W)
      endif
      if(allocated (aux_ovlp_3fn) )then
         deallocate (aux_ovlp_3fn)
      endif
      if(allocated( tmp_ovlp_3fn ))then
         deallocate (tmp_ovlp_3fn)
      endif

      return
      end subroutine get_self_energy
!---------------------------------------------------------------------
!**********

!****s
!* FHI-aims/get_self_energy
!  NAME
!   get_self_energy
!  SYNOPSIS

      subroutine get_self_energy_v1 &
         ( green_fn_time,screened_coul_int,&
           self_energy_time)


!  PURPOSE
!
! Evaluate the GW self energy in the time domain to avoid the convolution
! in frequency, G is expressed in the NAO basis, W is in the product basis
!
! Sigma_ij = Sum_lkmn  O^m_li O^n_kj G_lk (t) W_mn(t)
!
! lk refers to NAO basis, mn to the product basis
! O^m_li are overlap integrals between two NAO basis function and a product basis fn
!  and are given by  ovlp_3fn
! To improve scaling, first a partial summation is calculated
! 
!  A^m_ik (t) = Sum_l O^m_li G_lk (t)  
!  B^m_kj (t) = Sum_n O^n_kj W_mn(t)
!
! in this way the self energy is calculated as
!
! Sigma_ij (t) = Sum_mk A^m_ik (t)  B^m_kj (t)
!
!
!  USES

!      use physics
      use constants, only: pi
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use hartree_fock
      use scgw_grid
 
      implicit none

!  ARGUMENTS
      real*8     ::  green_fn_time (n_basis,n_basis,-ntau:ntau)
      real*8     ::  screened_coul_int (n_basbas, n_loc_prodbas, 0:ntau)
      real*8 ::  self_energy_time (n_basis,n_basis, -ntau:ntau)

! INPUTS
!
! o  ntau                 -- is the number of points in the time grid
! o  tau                  -- is the time grid
! o  ovlp_3fn             is the overlap integrals between 2 NAO basis and a product basis fn
! o  green_fn_time        the Greeen's Function in the NAO basis and time domain
! o  screened_coul_int    The screened Coulomb Interaction in the Product Basis and Time Domain, 
!                             which is symmetric in the time axis


!INTERNAL
      real*8, dimension (:,:,:) , allocatable  ::  partial_sum_G
      real*8, dimension (:,:,:) , allocatable  ::  partial_sum_W
      real*8, dimension (:,:,:) , allocatable  ::  aux_ovlp_3fn
!      real*8, dimension (:,:,:) , allocatable  ::  tmp_ovlp_3fn
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
 
      integer :: i_index
      integer :: i_spin
      integer :: i_task

      logical :: output
!     begin work

      if(.not. allocated (partial_sum_W) )then
         allocate (partial_sum_W(n_basis, n_basis, n_loc_prodbas))
      endif
      if(.not. allocated (partial_sum_G) )then
         allocate (partial_sum_G(n_basis, n_basis, n_loc_prodbas))
      endif
      if(.not. allocated( aux_ovlp_3fn ))then
         allocate (aux_ovlp_3fn (n_basis, n_basis, n_loc_prodbas))!local
      endif

      self_energy_time (:,:,:) = 0.d0
      aux_ovlp_3fn (:,:,:) = 0.d0

! filling the aux_ovlp_3fn matrix, (in this way we can use dgemm)

      do i_basis =1 , n_basis, 1
        do j_basis =1, n_basis, 1
          if(.not.basis_nghr(i_basis,j_basis).eq.0)then
           aux_ovlp_3fn (i_basis, j_basis,:) = &
           ovlp_3fn(basis_nghr(i_basis,j_basis),:)
          endif
        enddo
      enddo

      do i_tau = 0, ntau,1
        partial_sum_W (:,:,:)  = 0.d0
        !evaluate partial_sum_W
        do j_basis = 1, n_basis, 1
!          call dgemm ('N','T', n_basis,n_loc_prodbas,n_loc_prodbas,1.d0,&
!               aux_ovlp_3fn (j_basis,:,:), n_basis, &
!               screened_coul_int(:,:, i_tau),&
!               n_basbas,1.d0, partial_sum_W (j_basis,:,:), n_basis )
!          call dgemm ('N','T', n_basis,n_loc_prodbas,n_loc_prodbas,1.d0,&
!               aux_ovlp_3fn (j_basis,:,:), n_basis, &
!               screened_coul_int(:,:, i_tau),&
!               n_basbas,1.d0, partial_sum_W (j_basis,:,:), n_basis )
  
          call dgemm ('N','T', n_basis,n_loc_prodbas,n_loc_prodbas,1.d0,&
               screened_coul_int(:,:, i_tau),n_basis,&
               aux_ovlp_3fn (j_basis,:,:), &
!               screened_coul_int(:,:, i_tau),&
               n_basbas,1.d0, partial_sum_W (j_basis,:,:), n_basis )
!           call dgemm ('N','N', n_basis,n_loc_prodbas,n_basbas,1.d0,&
!                tmp_ovlp_3fn (j_basis,:,:), n_basis, &
!                screened_coul_int(:,:, i_tau), n_basbas,1.d0, &
!                partial_sum_W (j_basis,:,:), n_basis )

        ! enddo
        ! do j_basis = 1, n_basis, 1
           call sync_matrix (partial_sum_W (j_basis,1,1),&
                       n_basis, n_loc_prodbas )
        enddo

        partial_sum_G (:,:,:)  = 0.d0
        do i_basbas = 1, n_loc_prodbas, 1

                !this is only for positive times
           call dgemm ('N','N', n_basis,n_basis,n_basis,1.d0,&
               aux_ovlp_3fn (:,:,i_basbas), n_basis, &
               green_fn_time (:,:, i_tau), n_basis,1.d0, &
               partial_sum_G (:,:,i_basbas),n_basis )

           call dgemm ('N','N', n_basis,n_basis,n_basis,1.d0,&
                partial_sum_G (:,:,i_basbas), n_basis, &
                partial_sum_W (:,:,i_basbas), n_basis,1.d0, &
                self_energy_time (:,:,i_tau),n_basis)
        enddo

        partial_sum_G (:,:,:)  = 0.d0
        if (i_tau.gt.0)then !this is for negative time
          do i_basbas = 1, n_loc_prodbas, 1
             call dgemm ('N','N', n_basis,n_basis,n_basis,1.d0,&
               aux_ovlp_3fn (:,:,i_basbas), n_basis, &
               green_fn_time (:,:, -i_tau), n_basis,1.d0, &
               partial_sum_G (:,:,i_basbas ),n_basis )

             call dgemm ('N','N', n_basis,n_basis,n_basis,1.d0,&
                partial_sum_G (:,:,i_basbas), n_basis, &
                partial_sum_W (:,:,i_basbas), n_basis,1.d0, &
                self_energy_time (:,:,-i_tau),n_basis)            
          enddo
        endif

      enddo


      do i_tau = -ntau, ntau, 1
        call sync_matrix(&
          self_energy_time (1:n_basis,1:n_basis,i_tau), &
          n_basis, n_basis)
 
      enddo

      !this set the usual prefactor
      self_energy_time (:,:,:) = -self_energy_time (:,:,:)/2.d0/pi

      output = .false.
      if (output)then
       if(myid.eq.0)then      
         open(55,file= "self_energy_time.dat")
         do i_tau = -ntau, ntau, 1
           write(55,*)tau(i_tau), self_energy_time(1,1,i_tau)
         enddo 
         close(55)
       endif
      endif


      if( allocated (partial_sum_G) )then
         deallocate (partial_sum_G)
      endif
      if( allocated (partial_sum_W) )then
         deallocate (partial_sum_W)
      endif
      if(allocated (aux_ovlp_3fn) )then
         deallocate (aux_ovlp_3fn)
      endif

      return
      end subroutine get_self_energy_v1
!---------------------------------------------------------------------
!**********
!****s
!* FHI-aims/get_self_energy
!  NAME
!   get_self_energy
!  SYNOPSIS

      subroutine get_self_energy_v0 &
         ( green_fn_time,screened_coul_int,&
           self_energy_time)


!  PURPOSE
!
! Evaluate the GW self energy in the time domain to avoid the convolution
! in frequency, G is expressed in the NAO basis, W is in the product basis
!
! Sigma_ij = Sum_lkmn  O^m_li O^n_kj G_lk (t) W_mn(t)
!
! lk refers to NAO basis, mn to the product basis
! O^m_li are overlap integrals between two NAO basis function and a product basis fn
!  and are given by  ovlp_3fn
! To improve scaling, first a partial summation is calculated
! 
!  A^m_ik (t) = Sum_l O^m_li G_lk (t)  
!  B^m_kj (t) = Sum_n O^n_kj W_mn(t)
!
! in this way the self energy is calculated as
!
! Sigma_ij (t) = Sum_mk A^m_ik (t)  B^m_kj (t)
!
!
!  USES

!      use physics
      use constants, only: pi
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use hartree_fock
      use scgw_grid
 
      implicit none

!  ARGUMENTS
      real*8     ::  green_fn_time (n_basis,n_basis,-ntau:ntau)
      real*8     ::  screened_coul_int (n_basbas, n_loc_prodbas, 0:ntau)
      real*8 ::  self_energy_time (n_basis,n_basis, -ntau:ntau)

! INPUTS
!
! o  ntau                 -- is the number of points in the time grid
! o  tau                  -- is the time grid
! o  ovlp_3fn             is the overlap integrals between 2 NAO basis and a product basis fn
! o  green_fn_time        the Greeen's Function in the NAO basis and time domain
! o  screened_coul_int    The screened Coulomb Interaction in the Product Basis and Time Domain, 
!                             which is symmetric in the time axis


!INTERNAL
      real*8, dimension (:,:,:) , allocatable  ::  partial_sum_G
      real*8, dimension (:,:,:) , allocatable  ::  partial_sum_W
      real*8, dimension (:,:,:) , allocatable  ::  aux_ovlp_3fn
      real*8, dimension (:,:,:) , allocatable  ::  tmp_ovlp_3fn
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
 
      integer :: i_index
      integer :: i_spin
      integer :: i_task

      logical :: output
!     begin work

!      if(myid.eq.0)then
!        write(use_unit,*)"  --- Evaluating the Self Energy "
!      endif

      if(.not. allocated (partial_sum_W) )then
         allocate (partial_sum_W(n_basis, n_basis, n_loc_prodbas))
      endif
      if(.not. allocated (partial_sum_G) )then
         allocate (partial_sum_G(n_basis, n_basis, n_loc_prodbas))
      endif
      if(.not. allocated( tmp_ovlp_3fn ))then
         allocate (tmp_ovlp_3fn (n_basis,n_basis, n_basbas)) !global
      endif
      if(.not. allocated( aux_ovlp_3fn ))then
         allocate (aux_ovlp_3fn (n_basis, n_basis, n_loc_prodbas))!local
      endif

      self_energy_time (:,:,:) = 0.d0
      aux_ovlp_3fn (:,:,:) = 0.d0
      tmp_ovlp_3fn (:,:,:) = 0.d0

! filling the aux_ovlp_3fn matrix, (in this way we can use dgemm)

      do i_basis =1 , n_basis, 1
        do j_basis =1, n_basis, 1
          if(.not.basis_nghr(i_basis,j_basis).eq.0)then
           aux_ovlp_3fn (i_basis, j_basis,:) = &
           ovlp_3fn(basis_nghr(i_basis,j_basis),:)
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

!synchronizing
      do i_basis = 1, n_basis, 1
        call sync_matrix(&
          tmp_ovlp_3fn (i_basis,1:n_basis,1:n_basbas), &
          n_basis, n_basbas)
      enddo! i_basis
 
      do i_tau = 0, ntau,1
        partial_sum_W (:,:,:)  = 0.d0

        !evaluate partial_sum_W
        do j_basis = 1, n_basis, 1
           call dgemm ('N','N', n_basis,n_loc_prodbas,n_basbas,1.d0,&
                tmp_ovlp_3fn (j_basis,:,:), n_basis, &
                screened_coul_int(:,:, i_tau), n_basbas,1.d0, &
                partial_sum_W (j_basis,:,:), n_basis )
        enddo

        partial_sum_G (:,:,:)  = 0.d0
        do i_basbas = 1, n_loc_prodbas, 1

                !this is only for positive times
           call dgemm ('N','N', n_basis,n_basis,n_basis,1.d0,&
               aux_ovlp_3fn (:,:,i_basbas), n_basis, &
               green_fn_time (:,:, i_tau), n_basis,1.d0, &
               partial_sum_G (:,:,i_basbas),n_basis )

           call dgemm ('N','N', n_basis,n_basis,n_basis,1.d0,&
                partial_sum_G (:,:,i_basbas), n_basis, &
                partial_sum_W (:,:,i_basbas), n_basis,1.d0, &
                self_energy_time (:,:,i_tau),n_basis)
        enddo

        partial_sum_G (:,:,:)  = 0.d0
        if (i_tau.gt.0)then !this is for negative time
          do i_basbas = 1, n_loc_prodbas, 1
             call dgemm ('N','N', n_basis,n_basis,n_basis,1.d0,&
               aux_ovlp_3fn (:,:,i_basbas), n_basis, &
               green_fn_time (:,:, -i_tau), n_basis,1.d0, &
               partial_sum_G (:,:,i_basbas ),n_basis )

             call dgemm ('N','N', n_basis,n_basis,n_basis,1.d0,&
                partial_sum_G (:,:,i_basbas), n_basis, &
                partial_sum_W (:,:,i_basbas), n_basis,1.d0, &
                self_energy_time (:,:,-i_tau),n_basis)            
          enddo
        endif

      enddo


      do i_tau = -ntau, ntau, 1
        call sync_matrix(&
          self_energy_time (1:n_basis,1:n_basis,i_tau), &
          n_basis, n_basis)
 
      enddo

      !this set the usual prefactor
      self_energy_time (:,:,:) = -self_energy_time (:,:,:)/2.d0/pi

      output = .false.
      if (output)then
       if(myid.eq.0)then      
         open(55,file= "self_energy_time.dat")
         do i_tau = -ntau, ntau, 1
           write(55,*)tau(i_tau), self_energy_time(1,1,i_tau)
         enddo 
         close(55)
       endif
      endif


      if( allocated (partial_sum_G) )then
         deallocate (partial_sum_G)
      endif
      if( allocated (partial_sum_W) )then
         deallocate (partial_sum_W)
      endif
      if(allocated (aux_ovlp_3fn) )then
         deallocate (aux_ovlp_3fn)
      endif
      if(allocated( tmp_ovlp_3fn ))then
         deallocate (tmp_ovlp_3fn)
      endif

      return
      end subroutine get_self_energy_v0
!---------------------------------------------------------------------
