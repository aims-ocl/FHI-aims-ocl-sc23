
!  NAME compute_g_times_w_single_kpoint
!   
!  SYNOPSIS

subroutine compute_g_times_w_single_kpoint &
( n_low_state, n_KS_states, n_homo_kq, &
  i_freq, n_freq, n_full_freq, omega_full, womega_full, omega,  &
  KS_eigenvalue_kq, KS_eigenvector_kq, KS_eigenvector_complex_kq, &
  KS_eigenvector_k_current, &
  chemical_potential_spin, &
  lvl_tricoeff_recip_k_current, lvl_tricoeff_recip_kq, &
  screened_coulomb, gw_selfe_band, delta_gw_selfe_band)

  !  PURPOSE
  !  Subroutine compute_g_times_w_single_kpoint.f90 computes the GW selfenergy at
  !  a given k point
  !
  !  USES

  use dimensions
  use prodbas
  use hartree_fock
  use hartree_fock_p0
  use mpi_tasks
  use synchronize_mpi
  use synchronize_mpi_basic
  use pbc_lists
  use runtime_choices
  use timing
  use constants
  use basis, only: basis_atom
  use geometry, only: species
  implicit none

  !  ARGUMENTS

  integer :: n_low_state
  integer :: n_KS_states
  integer :: n_freq
  integer :: n_full_freq
  integer :: i_freq
  integer :: n_homo_kq(n_spin,n_k_points)
  real*8 :: omega_full(n_full_freq)
  real*8 :: womega_full(n_full_freq)
  real*8 :: omega(n_freq)
  real*8, dimension(n_states,n_spin,n_k_points) :: KS_eigenvalue_kq
  real*8, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_kq
  complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex_kq
  complex*16, dimension(n_basis,n_states,n_spin) :: KS_eigenvector_k_current
  real*8 :: chemical_potential_spin(n_spin)
  complex*16, dimension(n_basbas*(n_basbas+1)/2,n_k_points_task) :: screened_coulomb
  complex*16, dimension(max_n_basbas_sp,n_basis,n_basis) :: lvl_tricoeff_recip_k_current
  complex*16, dimension(max_n_basbas_sp,n_basis,n_basis,n_k_points_task) :: lvl_tricoeff_recip_kq
  complex*16, dimension(n_freq,n_low_state:n_KS_states,n_spin) :: gw_selfe_band
  complex*16, dimension(n_low_state:n_KS_states,n_spin) :: delta_gw_selfe_band

  !  INPUTS
  !  o  n_low_state -- the lowest single-particle states whose GW self-energy will be calculated
  !  o  n_KS_states -- the highest single-particle states whose GW self-energy will be calculated
  !  o  n_freq -- number of frequency points for the self-energy
  !  o  n_full_freq -- number of frequency points for the screened Coulomb matrix
  !  o  i_freq -- the i-th frequency point for the screened Coulomb matrix
  !  o  n_homo_kq -- the k-point- and spin-depdenent HOMO levels 
  !  o  omega_full -- the frequency grid for W (in the imaginary axis)
  !  o  womega_full -- the integration weight for the frequency grid for W
  !  o  omega -- full frequency grid for the self-energy
  !  o  occ_numbers -- real array,
  !       the occupation number of the electrons for each eigenstate and each spin
  !  o  KS_eigenvalue_kq -- the eigenvalues from single-particle KS/HF calculations, 
  !            on a mixed k-q mesh in reciprocal space
  !  o  KS_eigenvector_kq --  the real eigenvectors from effective single-particle calculations
  !            on a mixed k-q mesh in reciprocal space
  !  o  KS_eigenvector_complex_kq -- the complex eigenvectors from effective single-particle calculations
  !            on a mixed k-q mesh in reciprocal space
  !  o  KS_eigenvector_k_current -- the (complex) eigenvectors of a preceding single-particle (KS/HF) 
  !     self-consistent calculation on a specific k for band plotting
  !  o  screened_coulomb -- the screened coulomb matrix at a set of irreducible (regular) k grid points
  !
  !  INPUTS/OUTPUTS
  !  o gw_selfe_band -- complex, periodic GW self-energy on a special set of k points for band plotting, 
  !            each calling to this routine adds one more frequency point in the convolution
  !
  !  o delta_gw_selfe_band -- a correction term to the gw_selfe_band at the frequency point omega_n
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

! FIXME: need to change the data structure of lvl_tricoeff_kq_2
! FIXME: need to distribute the memory storage of screened_coulomb

  real*8  omega_n
  real*8  womega_n
  complex*16 omega_minus_en
  complex*16 delta_gw_tmp

  complex*16, allocatable :: KS_eigenvector_kq_current(:,:,:) ! KS eigenvector at a given k-q point in the BZ
  complex*16, allocatable :: lvl_tricoeff_sum(:,:) ! sum of LVL triple coefficient at k and q points
  complex*16, allocatable :: lvl_tricoeff_sum_all(:,:,:) ! lvl_tricoeff_sum for all babas values
  complex*16, allocatable :: coeff_tmp(:,:) ! temporary triple coefficient, one basis index is transformed to
                                               ! (occupied) molecular orbital index 
  complex*16, allocatable :: coeff_tmp_2(:,:) ! temporary triple coefficient, two basis indices are transformed to molecular orbital indices.
  complex*16, allocatable :: lvl_tricoeff_kq(:,:,:,:) ! Full LVL triple coefficient, one basis index is transformed to
                                                    ! (occupied) molecular orbital index 
  complex*16, allocatable :: screened_coulomb_current(:,:) ! temporary screened coulomb matrix
  complex*16, allocatable :: screened_coulomb_tmp(:,:) ! temporary screened coulomb matrix
!  complex*16, allocatable :: screened_coulomb_irk(:,:,:) ! full screened coulomb matrix
  complex*16, allocatable :: aux_screened_coulomb(:) 

  integer :: info, mpierr
  character(*), parameter :: func = 'compute_g_times_w_single_kpoint.f90'
  character*150 :: info_str
! intrinsic function
  complex*16 zdotc

!  timing info
!   real*8   ::  time_gw_1, clock_time_gw_1
!   real*8   ::  time_gw_2, clock_time_gw_2
!   real*8   ::  tot_time_gw_2, tot_clock_time_gw_2
!   real*8   ::  time_gw_3, clock_time_gw_3
!   real*8   ::  tot_time_gw_3, tot_clock_time_gw_3
!   real*8   ::  time_gw_4, clock_time_gw_4

! counter 
  integer i_q_point
  integer i_irq_point
  integer i_q_point_local
  integer i_irq_point_local
  integer i_basis_1
  integer i_basis_2
  integer i_prodbas_1
  integer i_prodbas_2
  integer i_state
  integer j_state
  integer i_freq_1
  integer i_spin
  integer i_index
  integer id_root
  integer i_atom_1, i_species_1, bboff_1, n_spbb_1
  integer i_atom_2, i_species_2, bboff_2, n_spbb_2
  integer i_task
  integer id_send


!  call get_timestamps(time_gw_1, clock_time_gw_1)

  allocate(KS_eigenvector_kq_current(n_basis,n_states,n_spin),stat=info) 
  call check_allocation(info, 'KS_eigenvector_k', func)
  allocate(lvl_tricoeff_sum(n_basis,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_sum', func)
  allocate(lvl_tricoeff_sum_all(n_basbas,n_basis,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_sum_sll', func)
  allocate(coeff_tmp(n_basis,n_KS_states),stat=info) 
  call check_allocation(info, 'coeff_tmp', func)
  allocate(coeff_tmp_2(n_states,n_KS_states),stat=info) 
  call check_allocation(info, 'coeff_tmp_2', func)
  allocate(lvl_tricoeff_kq(n_basbas,n_states,n_KS_states,n_spin),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_kq', func)
  allocate(screened_coulomb_current(n_basbas,n_basbas),stat=info) 
  call check_allocation(info, 'screened_coulomb_current', func)
  allocate(screened_coulomb_tmp(n_basbas,n_states),stat=info) 
  call check_allocation(info, 'screened_coulomb_tmp', func)
  allocate(aux_screened_coulomb(n_states),stat=info) 
  call check_allocation(info, 'aux_screened_coulomb', func)

!  call get_timestamps(rtime, clock_rtime)
!  time_gw_1=rtime-time_gw_1
!  tot_time_gw_2 = 0.d0
!  tot_time_gw_3 = 0.d0

  omega_n = omega_full(i_freq)
  womega_n = womega_full(i_freq)
  do i_q_point = 1, n_k_points, 1
!     if(myid.eq.0) then
!       write(use_unit,*) "i_q_point :", i_q_point
!     endif

!     call get_timestamps(time_gw_2, clock_time_gw_2)

     i_q_point_local = (i_q_point-1)/n_tasks + 1

     i_irq_point = irk_point_mapping(i_q_point)
     i_irq_point_local = (i_irq_point-1)/n_tasks + 1
     i_task = mod(i_q_point, n_tasks)
!     id_send = mod(i_irq_point, n_tasks)

!     if(myid.eq.id_send) then
      if(myid.eq.i_task) then
        i_index = 0
        do i_prodbas_1 = 1, n_basbas, 1
          do i_prodbas_2 = 1, i_prodbas_1
             i_index = i_index + 1
             screened_coulomb_current(i_prodbas_2,i_prodbas_1) = screened_coulomb(i_index,i_q_point_local)
!             screened_coulomb_current(i_prodbas_1,i_prodbas_2) = screened_coulomb(i_index,i_q_point_local)
          enddo
        enddo
      endif

! puting screened Coulomb matrix in the right place
!     if(i_task .ne. id_send) then
!        if(myid .eq. id_send) then
!           call send_complex_vector(screened_coulomb_current,n_basbas*n_basbas,i_task) 
!        elseif(myid .eq. i_task) then
!           call receive_complex_vector(screened_coulomb_current,n_basbas*n_basbas,id_send) 
!        endif
!     endif

     KS_eigenvector_kq_current(:,:,:)=(0.d0,0.d0) 
     if(myid.eq.mod(i_q_point, n_tasks)) then
       if(real_eigenvectors) then
          KS_eigenvector_kq_current(:,:,:) = KS_eigenvector_kq(:,:,:,i_q_point_local) 
       else
          KS_eigenvector_kq_current(:,:,:) = KS_eigenvector_complex_kq(:,:,:,i_q_point_local) 
       endif
     endif
!     call sync_vector_complex(KS_eigenvector_kq_current, n_basis*n_states*n_spin)

!     call get_timestamps(rtime, clock_rtime)
!     tot_time_gw_2 = tot_time_gw_2 + rtime - time_gw_2

     if(myid .ne. mod(i_q_point, n_tasks)) cycle

     lvl_tricoeff_sum_all = (0.d0,0.d0)
     do i_basis_1 = 1, n_basis, 1
       do i_basis_2 = 1, n_basis, 1
!          i_atom_1 = Cbasis_to_atom(i_basis_1)
          i_atom_1 = basis_atom(i_basis_1)
          i_species_1 = species(i_atom_1)
          bboff_1 = atom2basbas_off(i_atom_1)
          n_spbb_1 = sp2n_basbas_sp(i_species_1)

!          i_atom_2 = Cbasis_to_atom(i_basis_2)
          i_atom_2 = basis_atom(i_basis_2)
          i_species_2 = species(i_atom_2)
          bboff_2 = atom2basbas_off(i_atom_2)
          n_spbb_2 = sp2n_basbas_sp(i_species_2)

! in lvl_tricoeff_recip_q the q vector is associated with the second basis i_basis_2,
! for convenince we change this in lvl_tricoeff_sum_all
          lvl_tricoeff_sum_all(bboff_1+1:bboff_1+n_spbb_1,i_basis_2,i_basis_1) = &
          lvl_tricoeff_sum_all(bboff_1+1:bboff_1+n_spbb_1,i_basis_2,i_basis_1) + &
             conjg(lvl_tricoeff_recip_kq(1:n_spbb_1, i_basis_1, i_basis_2,i_q_point_local)) + &
             lvl_tricoeff_recip2(1:n_spbb_1, i_basis_2, i_basis_1)

          lvl_tricoeff_sum_all(bboff_2+1:bboff_2+n_spbb_2,i_basis_2,i_basis_1) = &
          lvl_tricoeff_sum_all(bboff_2+1:bboff_2+n_spbb_2,i_basis_2,i_basis_1) + &
             lvl_tricoeff_recip_k_current(1:n_spbb_2, i_basis_2, i_basis_1) + &
             lvl_tricoeff_recip2(1:n_spbb_2, i_basis_1, i_basis_2)

       enddo
     enddo
       
     lvl_tricoeff_kq(:,:,:,:) = (0.d0,0.d0)


     do i_prodbas_1 = 1, n_basbas, 1
       do i_basis_1 = 1, n_basis, 1
        do i_basis_2 = 1, n_basis, 1
           lvl_tricoeff_sum(i_basis_2, i_basis_1) = lvl_tricoeff_sum_all(i_prodbas_1, i_basis_2, i_basis_1)
!  end loop over i_basis_2
        enddo
!  end loop over i_basis_1
       enddo

       do i_spin = 1, n_spin, 1
          call zgemm('N', 'N', n_basis, n_KS_states, n_basis, (1.d0,0.d0), &
                      lvl_tricoeff_sum, n_basis, &
                      KS_eigenvector_k_current(1,1,i_spin), n_basis, (0.d0,0.d0), &
                      coeff_tmp, n_basis)

          call zgemm('C', 'N', n_states, n_KS_states, n_basis, (1.d0,0.d0), &
                      KS_eigenvector_kq_current(1,1,i_spin), n_basis,  &
                      coeff_tmp, n_basis, (0.d0,0.d0), &
                      coeff_tmp_2, n_states)

!            coeff_tmp_2 = (0.d0,0.d0)
!            do i_state = 1, n_KS_states, 1
!               do j_state = 1, n_states, 1
!                 do i_basis_1 =1, n_basis, 1
!                    coeff_tmp_2(j_state, i_state) = coeff_tmp_2(j_state, i_state) + &
!                    conjg(KS_eigenvector_kq_current(i_basis_1, j_state, i_spin)) * &
!                    coeff_tmp(i_basis_1, i_state)
!                    if( (i_state .eq. 1 .and. j_state .eq. 1 .and. i_prodbas_1 .eq. 1) ) then
!                       write(use_unit,'(4I4,6f11.6)') i_prodbas_1, i_state, j_state, i_basis_1, &
!                         conjg(KS_eigenvector_kq_current(i_basis_1, j_state, i_spin)), &
!                         coeff_tmp(i_basis_1, i_state), coeff_tmp_2(j_state, i_state)
!                    endif
!                 enddo
!               enddo
!            enddo

          lvl_tricoeff_kq(i_prodbas_1,:,:,i_spin) = coeff_tmp_2(:,:)

! end loop over i_spin
       enddo
! end loop over i_prodbas_1
     enddo

! mathematics behind: \sum_{i,j} \sum_{k,q} [(f_{i,k}-f_{j,q}) * C_{i,j}^\mu (k,q) * C_{j,i}^\mu (q,k)] / 
!                     (e_ik - e_jq - i*omega)
!            write(use_unit,*) "gw band screened_coulomb_current"
!            do i_prodbas_1 = 1, n_basbas, 1
!               write(use_unit,'(I4,4f18.6)') i_prodbas_1, screened_coulomb_current(1,i_prodbas_1), &
!                       screened_coulomb_current(i_prodbas_1,1)
!            enddo
!
     do i_spin = 1, n_spin, 1
      
! the k-dependence of the HOMO and LUMO levels is introduced to account for fractional occupations.
         do i_state = n_low_state, n_KS_states, 1

            call zhemm('L','U', n_basbas, n_states, (1.d0,0.d0), &
                        screened_coulomb_current, n_basbas, &
                        lvl_tricoeff_kq(1,1,i_state,i_spin), n_basbas, (0.d0,0.d0), &
                        screened_coulomb_tmp(1,1),n_basbas)
!               screened_coulomb_tmp(:,:)=(0.d0,0.d0)
!               do j_state = 1, n_states, 1
!                 do i_prodbas_1 = 1, n_basbas, 1
!                     do i_prodbas_2 = 1, n_basbas, 1
!                        screened_coulomb_tmp(i_prodbas_1, j_state) = &
!                           screened_coulomb_tmp(i_prodbas_1, j_state) + &
!                           screened_coulomb_current(i_prodbas_2, i_prodbas_1) * &
!                           lvl_tricoeff_kq(i_prodbas_2,j_state,i_state,i_spin)
!!                        if(i_state .eq. 1 .and. j_state.eq.1 .and. i_prodbas_1 .eq. 1) then
!!                           write(use_unit, '(I4,6f16.8)') i_prodbas_2, screened_coulomb_current(i_prodbas_2, i_prodbas_1), &
!!                                lvl_tricoeff_kq(i_prodbas_2,j_state,i_state,i_spin), &
!!                                screened_coulomb_tmp(i_prodbas_1, j_state)
!!                        endif
!                     enddo
!                  enddo
!              enddo

            do j_state = 1, n_states, 1
                 aux_screened_coulomb(j_state) =  &
                   zdotc(n_basbas, lvl_tricoeff_kq(1,j_state,i_state,i_spin), 1, &
                         screened_coulomb_tmp(1,j_state), 1)
!               aux_screened_coulomb(j_state)=(0.d0,0.d0)
!               do i_prodbas_1 = 1, n_basbas, 1
!                 aux_screened_coulomb(j_state) = aux_screened_coulomb(j_state) + &
!                     conjg(lvl_tricoeff_kq(i_prodbas_1,j_state,i_state,i_spin)) * &
!                      screened_coulomb_tmp(i_prodbas_1,j_state)
!
!                   if( i_state .eq. 1 .and. j_state .eq.1 ) then
!                     write(use_unit,'(2I4,6f16.6)')i_q_point, i_prodbas_1, lvl_tricoeff_kq(i_prodbas_1,j_state,i_state,i_spin), &
!                               screened_coulomb_tmp(i_prodbas_1, j_state), &
!                               aux_screened_coulomb(j_state)
!                   endif
!
!               enddo
                    
               do i_freq_1 = 1, n_freq, 1
                  omega_minus_en = (0.d0,1.d0)*omega(i_freq_1) + chemical_potential_spin(i_spin) - &
                       KS_eigenvalue_kq(j_state,i_spin,i_q_point)

! also including the negative frequency point -omega_n here
                  gw_selfe_band(i_freq_1,i_state,i_spin) = &
                    gw_selfe_band(i_freq_1,i_state,i_spin) - &
                    2.d0*omega_minus_en/(omega_minus_en*omega_minus_en + omega_n*omega_n) * &
                    aux_screened_coulomb(j_state)*womega_n*k_weights(i_q_point)

!                    if(i_freq_1.eq.1 .and. (i_state .eq. 1 ) .and. j_state .eq. 1 ) then
!                       write(use_unit,'(A,3I4,4f16.8)') "gw band:", i_q_point, i_state, j_state, &
!                         gw_selfe_band(i_freq_1,i_state,i_spin), aux_screened_coulomb(j_state)
!                    endif
                 
              enddo

              delta_gw_tmp=(0.d0,0.d0)
              do i_freq_1 = 1, n_full_freq, 1
! evaluate the correction term
                 omega_minus_en = (0.d0,1.d0)*omega_n + chemical_potential_spin(i_spin) - &
                    KS_eigenvalue_kq(j_state,i_spin,i_q_point)

                 delta_gw_tmp = delta_gw_tmp + &
                   2.d0*omega_minus_en / &
                   (omega_minus_en*omega_minus_en + omega_full(i_freq_1)*omega_full(i_freq_1)) * &
                   womega_full(i_freq_1)
              enddo

              if (j_state .le. n_homo_kq(i_spin,i_q_point)) then
                  delta_gw_selfe_band(i_state,i_spin) =  &
                   delta_gw_selfe_band(i_state,i_spin) + &
                   (delta_gw_tmp/2.d0/pi - 0.5d0) * &
                    aux_screened_coulomb(j_state)*k_weights(i_q_point)
              else
                  delta_gw_selfe_band(i_state,i_spin) = &
                   delta_gw_selfe_band(i_state,i_spin) + &
                   (delta_gw_tmp/2.d0/pi + 0.5d0) * &
                   aux_screened_coulomb(j_state)*k_weights(i_q_point)
              endif
              
!              write(use_unit,*) i_freq, i_k_point, i_state, j_state, n_homo_kq(i_spin,i_q_point)
! end of loop over j_state
           enddo
! end of loop over i_state
         enddo
! end of loop over i_spin
      enddo
!      write(use_unit,'(A, I4,4f18.8)')"band",  i_q_point, gw_selfe_band(5,1,1), gw_selfe_band(9, 1, 1) 
!      call get_timestamps(time_gw_3, clock_time_gw_3)
!      tot_time_gw_3 = tot_time_gw_3 + time_gw_3 - rtime

! end loop i_q_point
  enddo

!  call get_timestamps(time_gw_4, clock_time_gw_4)
  deallocate(KS_eigenvector_kq_current)
  deallocate(lvl_tricoeff_sum)
  deallocate(lvl_tricoeff_sum_all)
  deallocate(coeff_tmp)
  deallocate(coeff_tmp_2)
  deallocate(lvl_tricoeff_kq)
  deallocate(screened_coulomb_tmp)
  deallocate(screened_coulomb_current)
  deallocate(aux_screened_coulomb)
!  call get_timestamps(rtime, clock_rtime)

!  time_gw_4 = rtime - time_gw_4
!  if(myid.eq.0) then
!    write(use_unit,'(A,4f18.4)') " TIME : ", time_gw_1, tot_time_gw_2, tot_time_gw_3, time_gw_4
!  endif

  return

end subroutine compute_g_times_w_single_kpoint
!---------------------------------------------------------------------
!******
