!------------------------------------------------------
      subroutine get_exchange_self_energy_v1 &
         (ovlp_3fn, green_fn_time,&
          exchange_self_energy) 

!  PURPOSE
!  Evaluate the exact-exchange operator
!  from a given Green's function
!
!  USES
      use runtime_choices
      use species_data
      use constants
      use physics
      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: use_unit

      implicit none

!  ARGreal*8,UMENTS
      real*8     ::  ovlp_3fn             (n_basis_pairs,n_loc_prodbas)
      real*8     ::  green_fn_time        (n_basis,n_basis,n_spin)
      real*8     ::  exchange_self_energy (n_basis,n_basis,n_spin)

! INPUTS
!
! o  ovlp_3fn             is the overlap integrals between 2 NAO basis and a product basis fn
! o  green_fn_time        the Greeen's Function in the NAO basis and time domain
! o  exchange_self_energy the exact-exchange operator 

!INTERNAL
      real*8, dimension (:,:)   , allocatable  ::  partial_sum_G
      real*8, dimension (:,:,:) , allocatable  ::  aux_ovlp_3fn
      real*8  trace 
      real*8  exchange_en 
      logical output

!     counters
      integer :: i_basbas
      integer :: i_basis
      integer :: j_basis
      integer :: i_spin

!     begin work

      if(.not. allocated (partial_sum_G) )then
         allocate (partial_sum_G(n_basis, n_basis))
      endif
      if(.not. allocated( aux_ovlp_3fn ))then
         allocate (aux_ovlp_3fn (n_basis,n_basis, n_loc_prodbas))!local
      endif

      exchange_self_energy (:,:,:) = 0.d0
      
      ! unpack the aux_ovlp_3fn matrix
      aux_ovlp_3fn (:,:,:) = 0.d0
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

      do i_spin =1, n_spin , 1
        do i_basbas = 1, n_loc_prodbas, 1

           partial_sum_G (:,:) = 0.d0
  
           call dgemm ('N','N', n_basis,&
              n_basis,n_basis,1.d0,&
              green_fn_time(:,:,i_spin) , n_basis, &
              aux_ovlp_3fn (:,:,i_basbas), &
              n_basis,0.d0, &
              partial_sum_G , n_basis )
  
           call dgemm ( 'N','N', n_basis,&
              n_basis,n_basis,1.d0,&
              aux_ovlp_3fn (:,:,i_basbas), &
              n_basis, partial_sum_G ,&
              n_basis,1.d0, &
              exchange_self_energy (:,:,i_spin) , n_basis )
  
        enddo
   
        call sync_matrix(&
            exchange_self_energy(1:n_basis,1:n_basis,i_spin), &
            n_basis, n_basis)

      enddo

      exchange_self_energy (:,:,:) = &
           - exchange_self_energy (:,:,:) !*(2.d0/n_spin) 

      ! evaluate the exchange energy 
      if(.false.)then
        exchange_en = 0.d0
        do i_spin = 1, n_spin
          do i_basis =1 , n_basis, 1
            do j_basis =1, i_basis, 1
              if (i_basis .eq. j_basis)then
                 exchange_en = exchange_en + &
                   green_fn_time (i_basis,j_basis,i_spin)*  &
                   exchange_self_energy (i_basis,j_basis,i_spin)
              else 
                 exchange_en = exchange_en + &
                   2* green_fn_time (i_basis,j_basis,i_spin)*  &
                   exchange_self_energy (i_basis,j_basis,i_spin)
              endif
            enddo
          enddo
        enddo 
        exchange_en = exchange_en * (2.d0/n_spin) ! sum over spin
        exchange_en = exchange_en / 2.d0          ! exchange energy prefactor
  
        if(myid.eq.0) write(use_unit,*) "Exchange energy = ", &
              exchange_en *hartree, "eV" , exchange_en , "Ha"
      endif

      if( allocated (partial_sum_G) )then
         deallocate (partial_sum_G)
      endif
      if(allocated (aux_ovlp_3fn) )then
         deallocate (aux_ovlp_3fn)
      endif

      return
      end subroutine get_exchange_self_energy_v1

!------------------------------------------------
!****s

      subroutine get_exchange_self_energy &
         (ovlp_3fn, green_fn_time,&
          exchange_self_energy, hartree_pot, green_0_t_0, &
          full_hartree_pot )

!
!  USES
      use runtime_choices
      use species_data
      use constants
      use physics
      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi

      implicit none

!  ARGreal*8,UMENTS
      real*8     ::  ovlp_3fn (n_basis_pairs,n_loc_prodbas)  
      real*8     ::  green_fn_time (n_basis,n_basis)
      real*8     ::  exchange_self_energy (n_basis,n_basis)
      real*8     ::  hartree_pot (n_basis,n_basis)! this is the VARIATION of v_H from the DFT result
      real*8     ::  green_0_t_0 (n_basis,n_basis)
      real*8     ::  full_hartree_pot (n_basis,n_basis)

! INPUTS
!
! o  ovlp_3fn             is the overlap integrals between 2 NAO basis and a product basis fn
! o  green_fn_time        the Greeen's Function in the NAO basis and time domain
! o  screened_coul_int    The screened Coulomb Interaction in the Product Basis and Time Domain, 
!                             which is symmetric in the time axis


!INTERNAL
      real*8, dimension (:,:) , allocatable  ::  partial_sum_G
      real*8, dimension (:,:) , allocatable  ::  partial_sum_delta_G
      real*8, dimension (:,:) , allocatable  ::  delta_G
      real*8, dimension (:,:,:) , allocatable  ::  aux_ovlp_3fn
      real*8 trace 
      logical output

!     counters
      integer :: i_basbas
      integer :: i_basis
      integer :: j_basis
      integer :: i_spin

!     begin work

!      if(myid.eq.0)then
!        write(use_unit,*)"         | Evaluating the Exchange part"
!      endif

!      if(.not. allocated (partial_sum_loc) )then
!         allocate (partial_sum_loc(n_basis, n_basis, n_loc_prodbas, -ntau:ntau))
!      endif
!      if(.not. allocated (partial_sum_W) )then
!         allocate (partial_sum_W(n_basis, n_basis, n_loc_prodbas,0:ntau))
!      endif
      if(.not. allocated (partial_sum_G) )then
         allocate (partial_sum_G(n_basis, n_basis))
      endif
      if(.not. allocated (partial_sum_delta_G) )then
         allocate (partial_sum_delta_G(n_basis, n_basis))
      endif
      if(.not. allocated (delta_G) )then
         allocate (delta_G(n_basis, n_basis))
      endif
      if(.not. allocated( aux_ovlp_3fn ))then
         allocate (aux_ovlp_3fn (n_basis,n_basis, n_loc_prodbas))!local
      endif

      exchange_self_energy (:,:) = 0.d0
      delta_G (:,:)  = 0.d0
      
      do i_basis = 1, n_basis, 1
       do j_basis =1, n_basis, 1
         delta_G (i_basis, j_basis) =    &
         green_fn_time(i_basis, j_basis) &
       - green_0_t_0  (i_basis, j_basis)
       enddo
      enddo
! filling the aux_ovlp_3fn matrix, (in this way we can use dgemm)
      aux_ovlp_3fn (:,:,:) = 0.d0
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

      hartree_pot(:,:)  = 0.d0
      full_hartree_pot(:,:) = 0.d0

      do i_basbas = 1, n_loc_prodbas, 1
         partial_sum_G (:,:) = 0.d0
         partial_sum_delta_G (:,:) = 0.d0

         call dgemm ('N','N', n_basis,&
            n_basis,n_basis,1.d0,&
            delta_G , n_basis, &
            aux_ovlp_3fn (:,:,i_basbas), &
            n_basis,0.d0, &
            partial_sum_delta_G , n_basis )

         call dgemm ('N','N', n_basis,&
            n_basis,n_basis,1.d0,&
            green_fn_time , n_basis, &
            aux_ovlp_3fn (:,:,i_basbas), &
            n_basis,0.d0, &
            partial_sum_G , n_basis )
!
         call dgemm ( 'N','N', n_basis,&
            n_basis,n_basis,1.d0,&
            aux_ovlp_3fn (:,:,i_basbas), &
            n_basis, partial_sum_G ,&
            n_basis,1.d0, &
            exchange_self_energy , n_basis )

         !this is for the full hartree potential
         trace = 0
         do i_basis = 1, n_basis , 1
            trace =  trace + partial_sum_G (i_basis, i_basis)
         enddo
         full_hartree_pot(:,:) = full_hartree_pot(:,:)+&
           trace*aux_ovlp_3fn (:,:,i_basbas)

! this evaluates the hartree potential  in NAO basis, (we need 
! the same stuff)
         trace = 0
         do i_basis = 1, n_basis , 1
            trace =  trace + partial_sum_delta_G (i_basis, i_basis)
         enddo
         hartree_pot(:,:) = hartree_pot(:,:)+&
           trace*aux_ovlp_3fn (:,:,i_basbas)

      enddo
 
      call sync_matrix(&
          hartree_pot (1:n_basis,1:n_basis), &
          n_basis, n_basis)
      call sync_matrix(&
          full_hartree_pot (1:n_basis,1:n_basis), &
          n_basis, n_basis)
      call sync_matrix(&
          exchange_self_energy(1:n_basis,1:n_basis), &
          n_basis, n_basis)

!       exchange_self_energy (:,:) = exchange_self_energy (:,:) &
!            - hartree_pot (:,:)

!      if(myid.eq.0)then
!        output = .false.
!        if(output)then
!        open (22,file="exchange_self_energy.dat")   
!!        i_index =0
!        do i_basis = 1, n_basis, 1
!!          do j_basis = 1, n_basis, 1 
!!            i_index = i_index+1
!            write(22,*) full_hartree_pot(i_basis,i_basis), exchange_self_energy(i_basis,i_basis)
!!          enddo
!        enddo   
!        close(22)
!        endif
!      endif


      exchange_self_energy (:,:) = - exchange_self_energy (:,:)
!      hartree_pot(:,:)  = hartree_pot(:,:)!*2.0

      if( allocated (delta_G) )then
         deallocate (delta_G)
      endif
      if( allocated (partial_sum_G) )then
         deallocate (partial_sum_G)
      endif
      if(allocated (aux_ovlp_3fn) )then
         deallocate (aux_ovlp_3fn)
      endif
      if( allocated (partial_sum_delta_G) )then
         deallocate (partial_sum_delta_G)
      endif

      return
      end subroutine get_exchange_self_energy


