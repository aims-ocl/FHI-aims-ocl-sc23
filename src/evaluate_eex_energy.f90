!****s*  FHI-aims/evaluate_eex_energy
!  NAME
!    evaluate_eex_energy
!  SYNOPSIS

      subroutine evaluate_eex_energy &
           (n_homo_max,n_homo, &
            occ_numbers, ovlp_3KS, &
            eex_energy &
           )

!  PURPOSE
!  Subroutine evaluate_exchange_energy evaluates the exchange part of
!  the self energy for a given set of single-particle orbitals
!
!  Sigma^_ii(tau) = - \sum_l^occ \sum_mn O_ilm O_iln  V _mn,
!  V_mn is the bare Coulomb interaction matrix.
!
!  USES
      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: use_unit

      implicit none

!  ARGUMENTS  

      integer :: n_homo_max
      integer :: n_homo(n_spin)

      real*8  occ_numbers(n_states,n_spin,n_k_points)
      real*8  ovlp_3KS(n_loc_prodbas, n_states, n_homo_max,n_spin)

      real*8  eex_energy

!  INPUTS
!  o n_homo_max -- the larger HOMO level of the two spin components in case of
!         spin-polarized calculation
!  o n_homo -- HOMO level, or precisely, the number of occupied single particle state
!  o occ_number -- the occupation number 
!  o ovlp_3KS -- transformed 3-index Coulomb integral involving two KS/HF orbital
!       and one auxiliary basis
!
!  OUTPUTS
!   o eex_energy -- exact exchange energy 
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

      real*8  exchange_self_energy(n_homo_max,n_spin)
      real*8  ddot
!     counters

      integer :: i_state
      integer :: j_state
      integer :: i_spin


!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A)')"Starts to calculate the exchange energy ..."
      endif

!    rearrange the index of ovlp_3KS to get a two-dimensioal, working matrix

      exchange_self_energy(:,:)=0.d0

      do i_spin = 1, n_spin
        do i_state = 1, n_homo_max, 1

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
             exchange_self_energy, n_homo_max, n_spin )

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
!        write(use_unit,'(2X,A,f13.5)')
!     +      "  Exact exchange energy        :", eex_energy
!      endif

      return

      end subroutine evaluate_eex_energy
!---------------------------------------------------------------------
!******
!****s*  FHI-aims/evaluate_eex_energy
!  NAME
!    evaluate_eex_energy
!  SYNOPSIS

      subroutine evaluate_eex_energy_2 &
           (n_homo_max,n_homo, &
            occ_numbers, ovlp_3KS, &
            eex_energy &
           )

!  PURPOSE
!  Subroutine evaluate_exchange_energy evaluates the exchange part of
!  the self energy for a given set of single-particle orbitals
!
!  Sigma^_ii(tau) = - \sum_l^occ \sum_mn O_ilm O_iln  V _mn,
!  V_mn is the bare Coulomb interaction matrix.
!
!  USES
      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: use_unit

      implicit none

!  ARGUMENTS  

      integer :: n_homo_max
      integer :: n_homo(n_spin)

      real*8  occ_numbers(n_states,n_spin,n_k_points)
      real*8  ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)

      real*8  eex_energy

!  INPUTS
!  o n_homo_max -- the larger HOMO level of the two spin components in case of
!         spin-polarized calculation
!  o n_homo -- HOMO level, or precisely, the number of occupied single particle state
!  o occ_number -- the occupation number 
!  o ovlp_3KS -- transformed 3-index Coulomb integral involving two KS/HF orbital
!       and one auxiliary basis
!
!  OUTPUTS
!   o eex_energy -- exact exchange energy 
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

      real*8  exchange_self_energy(n_homo_max,n_spin)
      real*8  ddot
!     counters

      integer :: i_state, i_state_loc
      integer :: j_state, j_state_loc
      integer :: i_spin


!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A)')"Starts to calculate the exchange energy ..."
      endif

!    rearrange the index of ovlp_3KS to get a two-dimensioal, working matrix

      exchange_self_energy(:,:)=0.d0

      do i_spin = 1, n_spin
        do i_state = 1, n_homo_max, 1

         if(own_dim2_o3ks(i_state) /= myp2_o3ks) cycle
         i_state_loc = loc_dim2_o3ks(i_state)

          do j_state = 1, n_homo(i_spin)

           if(own_dim1_o3ks(j_state) /= myp1_o3ks) cycle
           j_state_loc = loc_dim1_o3ks(j_state)

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
        call sync_matrix( &
             exchange_self_energy, n_homo_max, n_spin )

      eex_energy = 0.d0
      do i_spin = 1, n_spin
        do i_state = 1, n_homo(i_spin), 1

          eex_energy = eex_energy + &
            exchange_self_energy(i_state, i_spin)* &
            occ_numbers(i_state, i_spin, 1)
        enddo
      enddo
      eex_energy = 0.5d0*eex_energy
      !if(myid.eq.0) then
      !  write(use_unit,'(2X,A,f13.5)') "  Exact exchange energy        :", eex_energy
      !endif

      return

      end subroutine evaluate_eex_energy_2
!---------------------------------------------------------------------
!******
