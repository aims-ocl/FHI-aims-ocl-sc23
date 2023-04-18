!****s*  FHI-aims/evaluate_exchange_energy
!  NAME
!    evaluate_exchange_energy
!  SYNOPSIS

     subroutine evaluate_exchange_energy &
         ( n_KS_states,n_homo, &
           occ_numbers, ovlp_3KS, &
           exchange_self_energy, &
           eex_energy &
          )

!  PURPOSE
!  Subroutine evaluate_exchange_energy evaluates the exchange part of
!  the self energy.
!
!   Sigma^_ii(tau) = - \sum_l^occ \sum_mn O_ilm O_iln  V _mn,
!
!  V_mn is the bare Coulomb interaction matrix.
! USES

      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: use_unit
      implicit none

!  ARGUMENTS

      integer :: n_KS_states
      integer :: n_homo(n_spin)
!      real*8  KS_eigenvector(n_basis,n_states,n_spin,n_k_points)
      real*8  occ_numbers(n_states,n_spin,n_k_points)
      real*8  ovlp_3KS(n_loc_prodbas, n_states, n_KS_states,n_spin)
      real*8  exchange_self_energy(n_KS_states,n_spin)
      real*8  eex_energy


! INPUTS
!  o n_KS_states -- number of KS states included in the calculations
!  o occ_numbers -- the occupation number 
!  o ovlp_3KS -- the transformed 3-function overlap integral involving
!                   two KS states and one auxiliary basis
!
! OUTPUT      
!  o n_homo -- the HOMO level
!  o exchange_self_energy -- exchange part of the self energy
!  o eex_energy -- exact HF exchange energy
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!  local variables

      real*8 ddot
!     counters

      integer :: i_state
      integer :: j_state
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_spin
      logical output


!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A)')"Starts to calculate the exchange energy ..."
      endif

!  determine the HOMO level
      do i_spin = 1, n_spin
       do i_state = 1, n_states
         if(occ_numbers(i_state,i_spin,1).gt.1.d-6) then
           n_homo(i_spin)=i_state
         endif
       enddo
      enddo

!    rearrange the index of ovlp_3KS to get a two-dimensioal, working matrix

      exchange_self_energy(:,:)=0.d0

      do i_spin = 1, n_spin
        do i_state = 1, n_KS_states, 1

          do j_state = 1, n_homo(i_spin)

           exchange_self_energy(i_state, i_spin) = &
            exchange_self_energy(i_state, i_spin) - &
                ddot( n_loc_prodbas, &
                      ovlp_3KS(:, j_state, i_state, i_spin), 1, &
                      ovlp_3KS(:, j_state, i_state, i_spin), 1 &
                    ) * &
                occ_numbers(j_state,i_spin,1)*dble(n_spin)/2.d0
!       if(i_state.eq.1) then
!        write(use_unit,*) j_state, exchange_energy(i_state,i_spin),
!     +           "ovlp",   ovlp_3KS(:,j_state,i_state,i_spin),
!     +           "occ",    occ_numbers(j_state,i_spin,1)
!       endif
          enddo

        enddo
      enddo
        call sync_matrix( &
             exchange_self_energy, n_KS_states, n_spin )
  
      eex_energy = 0.d0
      do i_spin = 1, n_spin
        do i_state = 1, n_homo(i_spin), 1

          eex_energy = eex_energy + &
            exchange_self_energy(i_state, i_spin)* &
            occ_numbers(i_state, i_spin, 1)
!           write(use_unit,'(2I4,f16.8)') i_spin, i_state, exchange_self_energy(i_state, i_spin)
!          if(myid.eq.0) then
!            write(use_unit,'(I4,3f18.8)')i_state, occ_numbers(i_state,i_spin,1), exchange_self_energy(i_state, i_spin), eex_energy
!          endif
        enddo
      enddo
      eex_energy = 0.5d0*eex_energy
!      if(myid.eq.0) then
!        write(use_unit,'(2X,A,f18.8,A)')&
!           "| Exact exchange energy:", eex_energy, "eV"
!      endif

!      output = .false.
!      if (output)then
!       open(44,file='exchange_self_energy_old.dat')
!       do i_state = 1, n_states, 1
!         write(44,*) i_state , exchange_self_energy(i_state, 1)
!       enddo
!       close(44)
!      endif

      return

      end subroutine evaluate_exchange_energy
!---------------------------------------------------------------------
!******
!===================================================================================================
!****s*  FHI-aims/evaluate_exchange_energy
!  NAME
!    evaluate_exchange_energy
!  SYNOPSIS

     subroutine evaluate_exchange_energy_2 &
         ( n_KS_states,n_homo, &
           occ_numbers, ovlp_3KS, &
           exchange_self_energy, &
           eex_energy &
          )

!  PURPOSE
!  Subroutine evaluate_exchange_energy evaluates the exchange part of
!  the self energy.
!
!   Sigma^_ii(tau) = - \sum_l^occ \sum_mn O_ilm O_iln  V _mn,
!
!  V_mn is the bare Coulomb interaction matrix.
! USES

      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: use_unit
      implicit none

!  ARGUMENTS

      integer :: n_KS_states
      integer :: n_homo(n_spin)
      real*8  occ_numbers(n_states,n_spin,n_k_points)
      real*8  ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
      real*8  exchange_self_energy(n_KS_states,n_spin)
      real*8  eex_energy


! INPUTS
!  o n_KS_states -- number of KS states included in the calculations
!  o occ_numbers -- the occupation number 
!  o ovlp_3KS -- the transformed 3-function overlap integral involving
!                   two KS states and one auxiliary basis
!
! OUTPUT      
!  o n_homo -- the HOMO level
!  o exchange_self_energy -- exchange part of the self energy
!  o eex_energy -- exact HF exchange energy
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




!  local variables

      real*8 ddot
!     counters

      integer :: i_state
      integer :: j_state
      integer :: i_spin
      integer :: i_state_loc
      integer :: j_state_loc


!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A)')"Starts to calculate the exchange energy ..."
      endif

!  determine the HOMO level
      do i_spin = 1, n_spin
       do i_state = 1, n_states
         if(occ_numbers(i_state,i_spin,1).gt.1.d-6) then
           n_homo(i_spin)=i_state
         endif
       enddo
      enddo

!    rearrange the index of ovlp_3KS to get a two-dimensioal, working matrix

      exchange_self_energy(:,:)=0.d0

      do i_spin = 1, n_spin
        do i_state = 1, n_KS_states, 1

          if(own_dim2_o3KS(i_state) /= myp2_o3KS) cycle
          i_state_loc = loc_dim2_o3KS(i_state)

          do j_state = 1, n_homo(i_spin)

            if(own_dim1_o3KS(j_state) /= myp1_o3KS) cycle
            j_state_loc = loc_dim1_o3KS(j_state)

            exchange_self_energy(i_state, i_spin) = &
              exchange_self_energy(i_state, i_spin) - &
                ddot( n_basbas, &
                      ovlp_3KS(:, j_state_loc, i_state_loc, i_spin), 1, &
                      ovlp_3KS(:, j_state_loc, i_state_loc, i_spin), 1 &
                    ) * &
                occ_numbers(j_state,i_spin,1)*dble(n_spin)/2.d0
          enddo

        enddo
      enddo

      call sync_matrix(exchange_self_energy, n_KS_states, n_spin )
  
      eex_energy = 0.d0
      do i_spin = 1, n_spin
        do i_state = 1, n_homo(i_spin), 1

          eex_energy = eex_energy + &
            exchange_self_energy(i_state, i_spin)* &
            occ_numbers(i_state, i_spin, 1)
        enddo
      enddo
      eex_energy = 0.5d0*eex_energy
!      if(myid.eq.0) then
!        write(use_unit,'(2X,A,f18.8,A)')&
!           "| Exact exchange energy:", eex_energy, "eV"
!      endif
      return

      end subroutine evaluate_exchange_energy_2
!---------------------------------------------------------------------
!******
