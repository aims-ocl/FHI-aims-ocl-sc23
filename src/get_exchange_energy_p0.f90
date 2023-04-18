!****s* FHI-aims/get_exchange_energy_p0
!  NAME
!   get_exchange_energy_p0
!  SYNOPSIS
      subroutine get_exchange_energy_p0 &
           (alpha, KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, k_p_weights, occ_numbers, &
           chemical_potential,fock_energy,energy_xc)

!  PURPOSE
!  Subroutine get_fock_energy evaluates the exact exchange energy
!  for Hatree-Fock (HF) and HF-based calculations.
!
!  USES

      use dimensions
      use hartree_fock_p0
      use mpi_utilities
      use runtime_choices
      use synchronize_mpi
      use scalapack_wrapper
      implicit none

!  imported variables
!  now defind in modle hartree_fock
!   alpha = 1: Hartree-Fock; 0.25: PBE0

!  ARGUMENTS

      real*8 alpha
      real*8, dimension(n_basis,n_states,n_spin,n_k_points_task) ::  &
           KS_eigenvector
      complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: &
           KS_eigenvector_complex
      real*8, dimension(n_states, n_spin, n_k_points) :: KS_eigenvalue
      real*8, dimension(n_k_points) :: k_p_weights
      real*8, dimension(n_states,n_spin,n_k_points) :: &
           occ_numbers
      real*8 :: chemical_potential

      real*8  fock_energy 
      real*8  energy_xc
!  INPUTS
!  o  alpha -- the mixing parameter for HF-based calculations, 1.0 for HF,
!         0.25 for PBE0, etc.
!  o  KS_eigenvector -- real array,
!         the eigenvector of the single-particle calculation
!  o  KS_eigenvector_complex -- complex array,
!         the eigenvector of the single-particle calculation
!  o  KS_eigenvalue -- eigenvalues
!  o  k_p_weights -- k-point weights
!  o  occ_numbers -- real array,
!         the occupation number of the electrons for each eigenstate and each spin
!  o chemical_potential -- Fermi energy
!
!  OUTPUTS
!  o  fock_energy -- real number, here the exact exchange energy
!  o  en_xc -- real number, the exchange-correlation energy, differs form fock_energy
!          in cases of hybrid functional calculations
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
    real*8 :: energy_fock
    real*8 :: energy_fock_SR
    real*8, allocatable :: tmp_r(:,:)
    real*8, allocatable :: tmp_r_SR(:,:)
    complex*16, allocatable :: tmp_c(:,:)
    complex*16, allocatable :: tmp_c_SR(:,:)
    real*8, dimension(n_basis) :: tmpvec_r
    complex*16, dimension(n_basis) :: tmpvec_c

!   counters

    integer :: i_state
    integer :: i_basis
    integer :: i_basis_1
    integer :: i_basis_2
    integer :: i_index
    integer :: i_spin
    integer :: i_k_point, i_k

!   begin work

    energy_fock=0.d0
    energy_fock_SR=0.d0

    if(use_scalapack) then

      if(real_eigenvectors)then

        allocate(tmp_r(mxld,mxcol))
        if (use_lc_wpbeh .and. alpha /= 0.d0) then
        	allocate(tmp_r_SR(mxld,mxcol))
        end if

        do i_spin = 1, n_spin

          call pdgemm('N','N', n_basis, n_states, n_basis, 1.d0, &
                      hf_exchange_matr_real(1,1,1,i_spin), 1, 1, sc_desc, &
                      eigenvec(1,1,i_spin),                1, 1, sc_desc, &
                      0.d0, tmp_r, 1, 1, sc_desc)
                      
          if (use_lc_wpbeh .and. alpha /= 0.d0) then
          	call pdgemm('N','N', n_basis, n_states, n_basis, 1.d0, &
                      hf_exchange_matr_real_SR(1,1,1,i_spin), 1, 1, sc_desc, &
                      eigenvec(1,1,i_spin),                1, 1, sc_desc, &
                      0.d0, tmp_r_SR, 1, 1, sc_desc)
          end if

          do i_state = 1, n_states
            if(l_col(i_state) == 0) cycle
            do i_basis = 1, n_basis
              if(l_row(i_basis)>0) then
                energy_fock = energy_fock + &
                 eigenvec(l_row(i_basis),l_col(i_state),i_spin) * tmp_r(l_row(i_basis),l_col(i_state)) * &
                 occ_numbers(i_state,i_spin,my_k_point) * k_p_weights(my_k_point)
                if (use_lc_wpbeh .and. alpha /= 0.d0) then
                  energy_fock_SR = energy_fock_SR + &
		             eigenvec(l_row(i_basis),l_col(i_state),i_spin) * tmp_r_SR(l_row(i_basis),l_col(i_state)) * &
		             occ_numbers(i_state,i_spin,my_k_point) * k_p_weights(my_k_point)
                end if
              endif
            enddo
          enddo

        enddo

        deallocate(tmp_r)
        if (use_lc_wpbeh .and. alpha /= 0.d0) then
        	deallocate(tmp_r_SR)
        end if

      else

        allocate(tmp_c(mxld,mxcol))
        if (use_lc_wpbeh .and. alpha /= 0.d0) then
        	allocate(tmp_c_SR(mxld,mxcol))
        end if

        do i_spin = 1, n_spin

          call pzgemm('N','N', n_basis, n_states, n_basis, (1.d0,0.d0), &
                      hf_exchange_matr_complex(1,1,1,i_spin), 1, 1, sc_desc, &
                      eigenvec_complex(1,1,i_spin),           1, 1, sc_desc, &
                      (0.d0,0.d0), tmp_c, 1, 1, sc_desc)
          
          if (use_lc_wpbeh .and. alpha /= 0.d0) then
          	call pzgemm('N','N', n_basis, n_states, n_basis, (1.d0,0.d0), &
                      hf_exchange_matr_complex_SR(1,1,1,i_spin), 1, 1, sc_desc, &
                      eigenvec_complex(1,1,i_spin),           1, 1, sc_desc, &
                      (0.d0,0.d0), tmp_c_SR, 1, 1, sc_desc)
          end if

          do i_state = 1, n_states
            if(l_col(i_state) == 0) cycle
            do i_basis = 1, n_basis
              if(l_row(i_basis)>0) then
                energy_fock = energy_fock + &
                 dble(conjg(eigenvec_complex(l_row(i_basis),l_col(i_state),i_spin)) * &
                      tmp_c(l_row(i_basis),l_col(i_state))) * &
                 occ_numbers(i_state,i_spin,my_k_point) * k_p_weights(my_k_point)
                if (use_lc_wpbeh .and. alpha /= 0.d0) then
                	energy_fock_SR = energy_fock_SR + &
		             dble(conjg(eigenvec_complex(l_row(i_basis),l_col(i_state),i_spin)) * &
		                  tmp_c_SR(l_row(i_basis),l_col(i_state))) * &
		             occ_numbers(i_state,i_spin,my_k_point) * k_p_weights(my_k_point)
                end if
              endif
            enddo
          enddo

        enddo

        deallocate(tmp_c)
        if (use_lc_wpbeh .and. alpha /= 0.d0) then
        	deallocate(tmp_c_SR)
        end if

      endif

    else

      if(real_eigenvectors) then        
         i_k = 0
         do i_k_point = 1, n_k_points,1
            if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid<= n_k_points ) then
               i_k = i_k + 1
               do i_spin = 1, n_spin
                  do i_state = 1, n_states, 1
                     call dgemv('t', n_basis, n_basis, 1.d0, hf_exchange_matr_real(1,1,i_k,i_spin), n_basis, &
                          KS_eigenvector(1,i_state,i_spin,i_k), 1, 0.d0, tmpvec_r, 1)
                     energy_fock = energy_fock + sum(tmpvec_r*KS_eigenvector(:,i_state,i_spin,i_k)) &
                          *occ_numbers(i_state,i_spin,i_k_point)*k_p_weights(i_k_point)
                     if (use_lc_wpbeh .and. alpha /= 0.d0) then
                        call dgemv('t', n_basis, n_basis, 1.d0, hf_exchange_matr_real_SR(1,1,i_k,i_spin), n_basis, &
                             KS_eigenvector(1,i_state,i_spin,i_k), 1, 0.d0, tmpvec_r, 1)
                        energy_fock_SR = energy_fock_SR + sum(tmpvec_r*KS_eigenvector(:,i_state,i_spin,i_k)) &
                             *occ_numbers(i_state,i_spin,i_k_point)*k_p_weights(i_k_point)
                     end if
                  enddo
               enddo
            endif
         enddo
      else ! eigenvectors complex
         i_k = 0
         do i_k_point = 1, n_k_points,1
            if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid<= n_k_points ) then
               i_k = i_k + 1
               do i_spin = 1, n_spin
                  do i_state = 1, n_states, 1
                     call zgemv('n', n_basis, n_basis, (1.d0,1.0d0), hf_exchange_matr_complex(1,1,i_k,i_spin), n_basis, &
                          KS_eigenvector_complex(1,i_state,i_spin,i_k), 1, (0.d0,0.0d0), tmpvec_c, 1)
                     energy_fock = energy_fock + dble(sum(tmpvec_c*conjg(KS_eigenvector_complex(:,i_state,i_spin,i_k)))) &
                          *occ_numbers(i_state,i_spin,i_k_point)*k_p_weights(i_k_point)
                     
                     if (use_lc_wpbeh .and. alpha /= 0.d0) then
                        call zgemv('n', n_basis, n_basis, (1.d0,1.0d0), hf_exchange_matr_complex_SR(1,1,i_k,i_spin), n_basis, &
                             KS_eigenvector_complex(1,i_state,i_spin,i_k), 1, (0.d0,0.d0), tmpvec_c, 1)
                        energy_fock_SR = energy_fock_SR + dble(sum(tmpvec_c*conjg(KS_eigenvector_complex(:,i_state,i_spin,i_k)))) &
                             *occ_numbers(i_state,i_spin,i_k_point)*k_p_weights(i_k_point)
                     end if
                  enddo
               enddo
            endif
         enddo
      endif ! end checking real/complex eigenvectors

    endif

    call sync_real_number(energy_fock)
    call sync_real_number(energy_fock_SR)

    if (use_lc_wpbeh) then
      fock_energy = -0.5d0*(energy_fock_SR+energy_fock)
      energy_xc = energy_xc + 0.5d0*(alpha*energy_fock_SR +(1.d0/lc_dielectric_constant)* energy_fock)
    else
      fock_energy = -0.5d0*energy_fock
      energy_xc = energy_xc + 0.5d0*alpha*energy_fock
    endif

    end subroutine get_exchange_energy_p0
!---------------------------------------------------------------------
!******
