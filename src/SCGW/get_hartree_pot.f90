!****s
!* FHI-aims/get_self_energy
!  NAME
!   get_self_energy
!  SYNOPSIS

      subroutine get_hartree_pot &
         (ovlp_3fn, green_fn_time,&
          hartree_pot )

!  PURPOSE
!  evaluate the hartree potential (in the NAO basis) 
!  from a given (non-packed) density matrix
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
      real*8     ::  ovlp_3fn (n_basis_pairs,n_loc_prodbas)  
      real*8     ::  green_fn_time (n_basis,n_basis,n_spin)
      real*8     ::  hartree_pot (n_basis,n_basis)! this is the VARIATION of v_H from the DFT result

! INPUTS
!
! o  ovlp_3fn             the overlap integrals between 2 NAO basis and a product basis fn
! o  green_fn_time        the Greeen's Function in the NAO basis and time domain


!INTERNAL
      real*8, dimension (:,:)   , allocatable  ::  partial_sum_G
      real*8, dimension (:,:,:) , allocatable  ::  aux_ovlp_3fn
      real*8  trace 
      real*8  hartree_en 
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

      !evaluate the Hartree potential 
      hartree_pot(:,:)  = 0.d0
      do i_spin = 1, n_spin , 1
        do i_basbas = 1, n_loc_prodbas, 1
  
           partial_sum_G (:,:) = 0.d0
  
           call dgemm ('N','N', n_basis,&
              n_basis,n_basis,1.d0,&
              green_fn_time (:,:,i_spin) , n_basis, &
              aux_ovlp_3fn (:,:,i_basbas), &
              n_basis,0.d0, &
              partial_sum_G , n_basis )
  
           trace = 0
           do i_basis = 1, n_basis , 1
              trace =  trace + partial_sum_G (i_basis, i_basis)
           enddo
           hartree_pot(:,:) = hartree_pot(:,:)+&
             trace*aux_ovlp_3fn (:,:,i_basbas)
        enddo
      enddo
 
      hartree_pot(:,:) = hartree_pot(:,:)*(2.d0/n_spin)
      call sync_matrix(&
          hartree_pot (1:n_basis,1:n_basis), &
          n_basis, n_basis)

      ! evaluate the hartree energy 
      if(.false.)then
        hartree_en = 0.d0
        do i_spin = 1, n_spin
          do i_basis =1 , n_basis, 1
            do j_basis =1, i_basis, 1
              if (i_basis .eq. j_basis)then
                 hartree_en = hartree_en + &
                   green_fn_time (i_basis,j_basis,i_spin)*  &
                   hartree_pot (i_basis,j_basis)
              else 
                 hartree_en = hartree_en + &
                   2* green_fn_time (i_basis,j_basis,i_spin)*  &
                   hartree_pot (i_basis,j_basis)
              endif
            enddo
          enddo
        enddo 
        hartree_en = hartree_en * (2.d0/n_spin) ! sum over spin
        hartree_en = hartree_en / 2.d0          ! hartree energy prefactor
  
        if(myid.eq.0) write(use_unit,*) "Hartree energy  = ", &
              hartree_en *hartree, "eV" , hartree_en , "Ha" 
      endif

      if( allocated (partial_sum_G) )then
         deallocate (partial_sum_G)
      endif
      if(allocated (aux_ovlp_3fn) )then
         deallocate (aux_ovlp_3fn)
      endif

      return
      end subroutine get_hartree_pot 
