!****s*  FHI-aims/evaluate_single_excitation_correction_p0
!  NAME
!    evaluate_single_excitation_correction_p0
!  SYNOPSIS
      subroutine evaluate_single_excitation_correction_p0 &
         ( n_KS_states, &
           occ_numbers, &
           hartree_potential,rho,rho_gradient, &
           kinetic_density, &
           partition_tab, basis_l_max, &
           hamiltonian, &
           KS_eigenvalue,KS_eigenvector, &
           KS_eigenvector_complex, &
           en_se, en_rse &
          )


!  PURPOSE
!  Subroutine evaluate_single_excitation_correction_p0 evaluates the single excitation
!  contribution to the second order correlation energy. This is zero for a HF reference
!  state, but non-zero for a KS reference state.
!
!  There is also a renormalized single excitation contribution which alleviate the
!  divergence problem of 2nd-order single excitation term.
!
!  It is given by \sum_{ia} |<i|V_HFx - V_KSxc|a>|^2 / (E_i - E_a)
!  where V_HFx is the HF exchange operator,  V_KSxc is the KS
!  exchange-correlation potential operator, and the summation of i, a 
!  runs over the occupied states and unoccupied states respectively.
! USES

      use dimensions
      use prodbas
      use hartree_fock_p0
      use pbc_lists
      use constants
      use mpi_tasks
      use synchronize_mpi
      use lapack_wrapper
      use scalapack_wrapper
      implicit none

!  ARGUMENTS

      integer :: n_KS_states
      real*8  occ_numbers(n_states,n_spin,n_k_points)
      real*8, dimension(n_full_points)            :: hartree_potential
      real*8, dimension(n_spin, n_full_points)    :: rho
      real*8, dimension(3, n_spin, n_full_points) :: rho_gradient
      real*8, dimension(n_full_points)            :: partition_tab
      integer ::  basis_l_max (n_species)
      real*8  :: hamiltonian(n_hamiltonian_matrix_size,n_spin)
      real*8  :: overlap_matrix(n_hamiltonian_matrix_size)

      real*8  KS_eigenvalue(n_states,n_spin,n_k_points)
      real*8  KS_eigenvector(n_basis,n_states,n_spin,n_k_points_task)
      complex*16  KS_eigenvector_complex(n_basis,n_states,n_spin,n_k_points_task)
      real*8  en_se
      real*8  en_rse

      real*8, dimension(n_spin, n_full_points)    :: kinetic_density
! INPUTS
!  o n_KS_states -- number of KS states included in the calculations
!  o occ_numbers -- the occupation number 
!  o KS_eigenvector -- Single-particle KS eigenvector
!  o ovlp_3KS -- the transformed 3-function overlap integral involving
!                   two KS states and one auxiliary basis
!  o xc_matr -- the matrix form the KS exchange-correlation potential operator
!            -- within the atomic orbital basis
!
! OUTPUT      
!  o en_se -- single excitation contribution to the 2nd order correlation energy 
!  o en_rse -- a summation of a ladder-type single excitation contribution to the 
!              correlation energy
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

      real*8  xc_matr(n_basis,n_basis,n_spin)
      real*8  ovlp_3KS(n_loc_prodbas,n_states,n_KS_states,n_spin)

      real*8  en_xc, en_pot_xc
      real*8  en_vdw, en_pot_vdw
      real*8  tmp_numerator
! working hamiltonian and overlap matrix array
      real*8,    dimension(:,:),allocatable :: hamiltonian_w
      real*8,    dimension(:),  allocatable :: overlap_matrix_w
      complex*16,dimension(:,:),allocatable :: hamiltonian_w_complex
      complex*16,dimension(:),  allocatable :: overlap_matrix_w_complex
      real*8,    dimension(:,:),allocatable :: hf_hamiltonian
      complex*16,dimension(:,:),allocatable :: hf_hamiltonian_complex

      real*8, dimension(:,:), allocatable :: hf_ham_matr  
      real*8, dimension(:,:), allocatable :: aux_hf_ham_matr  
      complex*16, dimension(:,:), allocatable :: hf_ham_matr_complex
      complex*16, dimension(:,:), allocatable :: aux_hf_ham_matr_complex
      real*8, dimension(:,:), allocatable :: KS_xc_matr  
      real*8, dimension(:,:), allocatable :: aux_xc_matr
      real*8, dimension(:,:,:), allocatable :: global_hfex_real 
      complex*16, dimension(:,:,:), allocatable :: global_hfex_complex


!  working arrays used in rSE calculations
!  real
      real*8,    dimension(:,:),allocatable :: HF_ham_occ
      real*8,    dimension(:,:),allocatable :: HF_ham_unocc
      real*8,    dimension(:,:),allocatable :: HF_ham_ou
      real*8,    dimension(:,:),allocatable :: trans_matr_occ
      real*8,    dimension(:,:),allocatable :: trans_matr_unocc
      real*8,    dimension(:),allocatable   :: HF_ev_occ
      real*8,    dimension(:),allocatable   :: HF_ev_unocc
      real*8,    dimension(:,:),allocatable   :: tmp_HF_ham_ou
      real*8,    dimension(:,:,:),allocatable :: work_ham
      real*8,    dimension(:,:),allocatable   :: work_ovl

!  complex
      complex*16,    dimension(:),allocatable :: HF_ham_occ_complex
      complex*16,    dimension(:),allocatable :: HF_ham_unocc_complex
      complex*16,    dimension(:,:),allocatable :: HF_ham_ou_complex
      complex*16,    dimension(:,:),allocatable :: trans_matr_occ_complex
      complex*16,    dimension(:,:),allocatable :: trans_matr_unocc_complex
      complex*16,    dimension(:,:),allocatable   :: tmp_HF_ham_ou_complex

      integer :: n_homo_k(n_spin,n_k_points)
      integer :: n_occ, n_unocc
      integer :: n_nonsingular
      integer :: info
      integer :: flag_xc_old
      real*8 ddot
!     counters

      integer :: i_state, i_state_1
      integer :: j_state, j_state_1
      integer :: k_state
      integer :: i_spin
      integer :: i_k_point
      integer :: i_k
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_index
      integer :: i_index_1


!     begin work

      if(flag_xc.eq.0) then
        if(myid.eq.0) then
          write(use_unit,*)
          write(use_unit,'(2X,A)') &
          "HF reference state: single excitation contribution is zero"
        endif
        en_se=0.d0
        en_rse=0.d0
        return
      endif 

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A)') &
            "Starting to calculate the single excitation correction ..."
      endif

!  determine the HOMO level at each k_points
      n_homo_k(:,:)=0
      do i_spin = 1, n_spin, 1
       do i_k_point = 1, n_k_points, 1
        do i_state = 1, n_states, 1

         if(occ_numbers(i_state,i_spin,i_k_point).gt.1.d-6) then
           n_homo_k(i_spin,i_k_point)=i_state
         endif

        enddo
       enddo
      enddo

!    rearrange the index of ovlp_3KS to get a two-dimensioal, working matrix

      allocate(KS_xc_matr(n_states, n_states))
      allocate(aux_xc_matr(n_states, n_basis))

!     evaluate the Hartree Hamiltonian, hence set flag_xc = 0
      flag_xc_old = flag_xc
      flag_xc = 0
      call integrate_real_hamiltonian_matrix_p2 &
           ( hartree_potential,  rho, rho_gradient, kinetic_density, partition_tab, &
             basis_l_max, en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw)
      flag_xc = flag_xc_old

!     Construct the periodic Hartree-Fock hamiltonian
      if(real_eigenvectors) then
!          write(use_unit,*) "NOW REAL EIGENVECTOR!"
          allocate(hamiltonian_w(n_basis*(n_basis+1)/2,n_spin),stat=info)
          call check_allocation(info, 'hamiltonian_w                 ')

          allocate(overlap_matrix_w(n_basis*(n_basis+1)/2),stat=info)
          call check_allocation(info, 'overlap_matrix_w              ')

          allocate(hamiltonian_w_complex(1,1),stat=info)
          call check_allocation(info, 'hamiltonian_w_complex         ')

          allocate(overlap_matrix_w_complex(1),stat=info)
          call check_allocation(info, 'overlap_matrix_w_complex      ')

          allocate(hf_hamiltonian(n_states,n_states),stat=info)
          call check_allocation(info, 'hf_hamiltonian  ')

          allocate(hf_ham_matr(n_basis,n_basis),stat=info)
          call check_allocation(info, 'hf_ham_matr  ')

          allocate(aux_hf_ham_matr(n_basis,n_states),stat=info)
          call check_allocation(info, 'aux_hf_ham_matr  ')

      else 
!          write(use_unit,*) "NOW COMPLEX EIGENVECTOR!"
          allocate(hamiltonian_w_complex(n_basis*(n_basis+1)/2,n_spin),stat=info)
          call check_allocation(info, 'hamiltonian_w_complex         ')

          allocate(overlap_matrix_w_complex(n_basis*(n_basis+1)/2),stat=info)
          call check_allocation(info, 'overlap_matrix_w_complex      ')

          allocate(hamiltonian_w(1,1),stat=info)
          call check_allocation(info, 'hamiltonian_w_complex         ')

          allocate(overlap_matrix_w(1),stat=info)
          call check_allocation(info, 'overlap_matrix_w_complex      ')

          allocate(hf_hamiltonian_complex(n_states,n_states),stat=info)
          call check_allocation(info, 'hf_hamiltonian_complex  ')

          allocate(hf_ham_matr_complex(n_basis,n_basis),stat=info)
          call check_allocation(info, 'hf_ham_matr_complex ')

          allocate(aux_hf_ham_matr_complex(n_basis,n_states),stat=info)
          call check_allocation(info, 'aux_hf_ham_matr_complex  ')


      end if
      allocate(work_ham(n_centers_basis_I , n_centers_basis_I,n_spin ))
      call check_allocation(info, 'work_ham')
      allocate(work_ovl(n_centers_basis_I , n_centers_basis_I ))   
      call check_allocation(info, 'work_ovl')

      if(use_scalapack) then
         if(real_eigenvectors) then
            allocate(global_hfex_real(n_basis,n_basis,n_spin))
            call check_allocation(info, 'global_hfex_real')
            do i_spin=1,n_spin
               call get_scalapack_global_rmatrix(hf_exchange_matr_real(:,:,1,i_spin), global_hfex_real(:,:,i_spin))
            end do
         else
            allocate(global_hfex_complex(n_basis,n_basis,n_spin))
            call check_allocation(info, 'global_hfex_complex')
            do i_spin=1,n_spin
               call get_scalapack_global_zmatrix(hf_exchange_matr_complex(:,:,1,i_spin), global_hfex_complex(:,:,i_spin))
            end do
         end if
      end if

! start to evaluate SE and rSE contributions
      en_se = 0.d0
      en_rse = 0.d0
      i_k=0
      do i_k_point = 1, n_k_points,1

!         write(use_unit,*) "i_k_point: ", i_k_point
         if (myid.eq.k_point_loc(1,i_k_point)) then

           i_k = i_k + 1

           call construct_hamiltonian_and_ovl(hamiltonian, overlap_matrix, &
                  hamiltonian_w, overlap_matrix_w, &
                  hamiltonian_w_complex, overlap_matrix_w_complex, i_k_point, &
                  work_ham, work_ovl)

!           if(real_eigenvectors) then
!               write(use_unit,*) "real hamil: before ", hamiltonian_w(:,:)
!               write(use_unit,*) "exchange: before ", hf_exchange_matr_real(:,:,:,:)
!               call get_hf_hamiltonian_real_p0(hf_exchange_matr_real,hamiltonian_w,i_k)
!               write(use_unit,*) "real hamil: after ", hamiltonian_w(:,:)
!               write(use_unit,*) "exchange: after ", hf_exchange_matr_real(:,:,:,:)
!           else
!               write(use_unit,*) "hamiltonian_w: ",  hamiltonian_w
!               write(use_unit,*) "hf_exchange_matr_real: ",  hf_exchange_matr_real
!               call get_hf_hamiltonian_complex_p0(hf_exchange_matr_complex,hamiltonian_w_complex,i_k)
!               write(use_unit,*) " hamiltonian complex: ", hamiltonian_w_complex(:,:)
!           endif

           do i_spin = 1, n_spin

!    Construct the HF Hamiltonian in NAO basis
            if(real_eigenvectors) then
              i_index = 0
              if(use_scalapack) then
                 do i_basis_1 = 1, n_basis, 1
                    do i_basis_2 = 1, i_basis_1, 1
                       i_index = i_index+1
                       hf_ham_matr(i_basis_2,i_basis_1)=hamiltonian_w(i_index,i_spin) - &
                            global_hfex_real(i_basis_2,i_basis_1,i_spin)
                       hf_ham_matr(i_basis_1,i_basis_2)=hf_ham_matr(i_basis_2,i_basis_1)
                    enddo
                 enddo
              else
                 do i_basis_1 = 1, n_basis, 1
                    do i_basis_2 = 1, i_basis_1, 1
                       i_index = i_index+1
                       hf_ham_matr(i_basis_2,i_basis_1)=hamiltonian_w(i_index,i_spin) - &
                            hf_exchange_matr_real(i_basis_2,i_basis_1,i_k,i_spin)
                       hf_ham_matr(i_basis_1,i_basis_2)=hf_ham_matr(i_basis_2,i_basis_1)
                    enddo
                 enddo
              end if

             call dgemm('N', 'N', n_basis, n_states, &
                  n_basis, 1.0d0, &
                  hf_ham_matr, n_basis, &
                  KS_eigenvector(1,1,i_spin,i_k), n_basis, 0.d0, &
                  aux_hf_ham_matr, n_basis) 

!    Construct the HF Hamiltonian in KS basis
             call dgemm('T', 'N', n_states, n_states, &
                  n_basis, 1.0d0, &
                  KS_eigenvector(1,1,i_spin,i_k), n_basis,  &
                  aux_hf_ham_matr, n_basis, 0.d0, & 
                  hf_hamiltonian, n_states)

             else

              i_index = 0
              if(use_scalapack) then
                 do i_basis_1 = 1, n_basis, 1
                    do i_basis_2 = 1, i_basis_1, 1
                       i_index = i_index+1
                       hf_ham_matr_complex(i_basis_2,i_basis_1)=hamiltonian_w_complex(i_index,i_spin)- &
                            global_hfex_complex(i_basis_2,i_basis_1,i_spin)                       
                       hf_ham_matr_complex(i_basis_1,i_basis_2)=conjg(hf_ham_matr_complex(i_basis_2,i_basis_1))
                    enddo
                 enddo
              else
                 do i_basis_1 = 1, n_basis, 1
                    do i_basis_2 = 1, i_basis_1, 1
                       i_index = i_index+1
                       hf_ham_matr_complex(i_basis_2,i_basis_1)=hamiltonian_w_complex(i_index,i_spin)- &
                            hf_exchange_matr_complex(i_basis_2,i_basis_1,i_k,i_spin)
                       
                       hf_ham_matr_complex(i_basis_1,i_basis_2)=conjg(hf_ham_matr_complex(i_basis_2,i_basis_1))
                    enddo
                 enddo
              endif

!  XR: Matrix multiplication does not work here. No idea what happened.
!             aux_hf_ham_matr_complex(:,:)=(0.d0,0.d0)
!             call zgemm('N', 'N', n_basis, n_states, &
!                  n_basis, (1.d0,0.d0), &
!                  hf_ham_matr_complex, n_basis, &
!                  KS_eigenvector_complex(1,1,i_spin,i_k), &
!                  n_basis, (0.d0,0.d0), &
!                  aux_hf_ham_matr_complex, n_basis) 
              aux_hf_ham_matr_complex(:,:)=(0.d0,0.d0)
              do i_state = 1, n_states, 1
                 do i_basis_1= 1, n_basis, 1
                    do i_basis_2 =1, n_basis, 1
                       aux_hf_ham_matr_complex(i_basis_1,i_state) = &
                          aux_hf_ham_matr_complex(i_basis_1,i_state) + &
                          hf_ham_matr_complex(i_basis_1,i_basis_2) * &
                          KS_eigenvector_complex(i_basis_2,i_state,i_spin,i_k)
                    enddo
!                   write(use_unit,'(2I4,2f16.8)') i_basis_1, i_state, aux_hf_ham_matr_complex(i_basis_1,i_state)
                 enddo
              enddo

!    Construct the HF Hamiltonian in KS basis -- complex
!             call zgemm('C', 'N', n_states, n_states, &
!                  n_basis, (1.d0,0.d0), &
!                  KS_eigenvector_complex(1,1,i_spin,i_k), n_basis,  &
!                  aux_hf_ham_matr_complex, &
!                  n_basis, (0.d0,0.d0), & 
!                  hf_hamiltonian_complex, n_states)
             hf_hamiltonian_complex(:,:)=(0.d0,0.d0)
             do i_state = 1, n_states, 1
                 do j_state = 1, n_states, 1
                    do i_basis_2 =1, n_basis, 1
                       hf_hamiltonian_complex(j_state,i_state) = &
                           hf_hamiltonian_complex(j_state,i_state) + &
                           conjg(KS_eigenvector_complex(i_basis_2,j_state,i_spin,i_k)) *  &
                           aux_hf_ham_matr_complex(i_basis_2,i_state) 
                    enddo
                 enddo
              enddo

! end of if real_eigenvectors
            endif

            do j_state = n_homo_k(i_spin,i_k_point)+1, n_states, 1
              do i_state = 1, n_homo_k(i_spin,i_k_point), 1
               if(real_eigenvectors) then
                   en_se = en_se + &
                   hf_hamiltonian(i_state,j_state)*hf_hamiltonian(j_state,i_state) / &
                   (KS_eigenvalue(i_state,i_spin,i_k_point)-KS_eigenvalue(j_state,i_spin,i_k_point)) * &
                   occ_numbers(i_state,i_spin,i_k_point)*k_weights(i_k_point)
               else
                   en_se = en_se + &
                   dble(hf_hamiltonian_complex(i_state,j_state)*hf_hamiltonian_complex(j_state,i_state)) / &
                   (KS_eigenvalue(i_state,i_spin,i_k_point)-KS_eigenvalue(j_state,i_spin,i_k_point)) * &
                   occ_numbers(i_state,i_spin,i_k_point)*k_weights(i_k_point)
               endif
              enddo
            enddo

            n_occ=n_homo_k(i_spin,i_k_point)
            n_unocc=n_states-n_occ
            if(n_occ.le.0) cycle

            if(real_eigenvectors) then

              allocate(HF_ham_occ(n_occ,n_occ))
              allocate(HF_ham_unocc(n_unocc,n_unocc))
              allocate(HF_ham_ou(n_occ,n_unocc))
              allocate(trans_matr_occ(n_occ,n_occ))
              allocate(trans_matr_unocc(n_unocc,n_unocc))
              allocate(tmp_HF_ham_ou(n_occ,n_unocc))
              allocate(HF_ev_occ(n_occ))
              allocate(HF_ev_unocc(n_unocc))

              do i_state = 1, n_states, 1
               do j_state = 1, n_states, 1
                if((i_state .le. n_occ) .and. (j_state .le. n_occ)) then
                  HF_ham_occ(j_state,i_state) = hf_hamiltonian(j_state,i_state)
                elseif ((i_state .gt. n_occ) .and. (j_state .gt. n_occ)) then
                  HF_ham_unocc(j_state - n_occ, i_state - n_occ) = hf_hamiltonian(j_state,i_state)
                elseif ((j_state .le. n_occ) .and. (i_state .gt. n_occ)) then
                  HF_ham_ou(j_state,i_state - n_occ) = hf_hamiltonian(j_state,i_state)
                endif
               enddo
              enddo

! Diagonalize the occupied and unoccupied HF Hamiltonian (represented in KS orbital basis)
              call diagonalize_auxmat_lapack &
              ( n_occ, HF_ham_occ, safe_minimum, &
                -1.d9, n_nonsingular, HF_ev_occ, trans_matr_occ, &
                "RSE-HF-Occ" )

              call diagonalize_auxmat_lapack &
              ( n_unocc, HF_ham_unocc, safe_minimum, &
                -1.d9, n_nonsingular, HF_ev_unocc, trans_matr_unocc, &
                "RSE-HF-Unocc" )

              call dgemm('N','N', n_occ, n_unocc, n_unocc, 1.d0, &
                     HF_ham_ou, n_occ, trans_matr_unocc, n_unocc, 0.d0, &
                     tmp_HF_ham_ou, n_occ)

              call dgemm('N','N', n_occ, n_unocc, n_occ, 1.d0, &
                     trans_matr_occ, n_occ, tmp_HF_ham_ou, n_occ, 0.d0, &
                     HF_ham_ou, n_occ)

              do j_state = 1, n_unocc, 1
                do i_state = 1, n_occ, 1
                   en_rse = en_rse + &
                   HF_ham_ou(i_state,j_state)*HF_ham_ou(i_state,j_state) / &
                   (HF_ev_occ(i_state) - HF_ev_unocc(j_state)) * &
                   occ_numbers(i_state,i_spin,i_k_point)*k_weights(i_k_point)
                enddo
              enddo

              deallocate(HF_ham_occ)
              deallocate(HF_ham_unocc)
              deallocate(HF_ham_ou)
              deallocate(trans_matr_occ)
              deallocate(trans_matr_unocc)
              deallocate(tmp_HF_ham_ou)
              deallocate(HF_ev_occ)
              deallocate(HF_ev_unocc)
          else
             !complex case
             allocate(HF_ham_occ_complex(n_occ*(n_occ+1)/2))
             allocate(HF_ham_unocc_complex(n_unocc*(n_unocc+1)/2))
             allocate(HF_ham_ou_complex(n_occ,n_unocc))
             allocate(trans_matr_occ_complex(n_occ,n_occ))
             allocate(trans_matr_unocc_complex(n_unocc,n_unocc))
             allocate(tmp_HF_ham_ou_complex(n_occ,n_unocc))
             allocate(HF_ev_occ(n_occ))
             allocate(HF_ev_unocc(n_unocc))

              do i_state = 1, n_states, 1
               do j_state = 1, i_state, 1
                i_index = i_state*(i_state-1)/2 + j_state
                if(i_state .le. n_occ) then
                  HF_ham_occ_complex(i_index) = hf_hamiltonian_complex(j_state,i_state)
                elseif ( j_state .gt. n_occ ) then
                  i_state_1 = i_state - n_occ
                  j_state_1 = j_state - n_occ
                  i_index_1 = (i_state_1-1)*i_state_1/2 + j_state_1
                  HF_ham_unocc_complex(i_index_1) = hf_hamiltonian_complex(j_state,i_state)
                elseif ((j_state .le. n_occ) .and. (i_state .gt. n_occ)) then
                  HF_ham_ou_complex(j_state,i_state-n_occ) = hf_hamiltonian_complex(j_state,i_state)
                endif
               enddo
              enddo

               call diagonalize_rse_hamiltonian_complex &
              ( n_occ,  n_occ, HF_ham_occ_complex, &
               safe_minimum, HF_ev_occ,  trans_matr_occ_complex )

               call diagonalize_rse_hamiltonian_complex &
              ( n_unocc, n_unocc, HF_ham_unocc_complex, &
               safe_minimum, HF_ev_unocc,  trans_matr_unocc_complex )

              call zgemm('N','N', n_occ, n_unocc, n_unocc, (1.d0,0.d0), &
                     HF_ham_ou_complex, n_occ, trans_matr_unocc_complex, n_unocc, (0.d0,0.d0), &
                     tmp_HF_ham_ou_complex, n_occ)

              call zgemm('C','N', n_occ, n_unocc, n_occ, (1.d0,0.d0), &
                     trans_matr_occ_complex, n_occ, tmp_HF_ham_ou_complex, n_occ, (0.d0,0.d0), &
                     HF_ham_ou_complex, n_occ)

              do j_state = 1, n_unocc, 1
                do i_state = 1, n_occ, 1
                   en_rse = en_rse + &
                        dble(HF_ham_ou_complex(i_state,j_state)*conjg(HF_ham_ou_complex(i_state,j_state))) / &
                   (HF_ev_occ(i_state) - HF_ev_unocc(j_state)) * &
                   occ_numbers(i_state,i_spin,i_k_point)*k_weights(i_k_point)
!                   write(use_unit,*) HF_ham_ou_complex(i_state,j_state) 
!                   write(use_unit,*) HF_ev_occ(i_state), HF_ev_unocc(j_state)
                enddo
              enddo

              deallocate(HF_ham_occ_complex)
              deallocate(HF_ham_unocc_complex)
              deallocate(HF_ham_ou_complex)
              deallocate(trans_matr_occ_complex)
              deallocate(trans_matr_unocc_complex)
              deallocate(tmp_HF_ham_ou_complex)
              deallocate(HF_ev_occ)
              deallocate(HF_ev_unocc)

          endif
            
!  end of loop over i_spin
         enddo
 
!  end of  if (myid.eq. k_point_loc(1,i_k_point))
       endif
!  end of loop over i_k_point 
      enddo

      call sync_real_number(en_se)
      call sync_real_number(en_rse)

      if(myid.eq.0) then
        write(use_unit,*) 
        write(use_unit,'(2X,A,f13.8, A)') &
           "| Single excitation contribution to the correlation energy:", &
           en_se*hartree, " eV"
        write(use_unit,'(2X,2A,f13.8, A)') &
           "| Renormalized single excitation contribution (diag)", &
           " to the correlation energy:", &
           en_rse*hartree, "  eV"
        write(use_unit,*) 
      endif

      if(allocated(hamiltonian_w)) then
        deallocate(hamiltonian_w)
      endif
      if(allocated(overlap_matrix_w)) then
        deallocate(overlap_matrix_w)
      endif
      if(allocated(hamiltonian_w_complex)) then
        deallocate(hamiltonian_w_complex)
      endif
      if(allocated(overlap_matrix_w_complex)) then
        deallocate(overlap_matrix_w_complex)
      endif
      if(allocated(hf_hamiltonian)) then
        deallocate(hf_hamiltonian)
      endif
      if(allocated(hf_ham_matr)) then
        deallocate(hf_ham_matr)
      endif
      if(allocated(aux_hf_ham_matr)) then
        deallocate(aux_hf_ham_matr)
      endif
      if(allocated(hf_hamiltonian_complex)) then
        deallocate(hf_hamiltonian_complex)
      endif
      if(allocated(hf_ham_matr_complex)) then
        deallocate(hf_ham_matr_complex)
      endif
      if(allocated(aux_hf_ham_matr_complex)) then
        deallocate(aux_hf_ham_matr_complex)
      endif
      if(allocated(KS_xc_matr))then
        deallocate(KS_xc_matr)
      endif
      if(allocated(aux_xc_matr))then
        deallocate(aux_xc_matr)
      endif
      if(allocated(work_ham)) then
        deallocate(work_ham)
      endif
      if(allocated(work_ovl)) then
        deallocate(work_ovl)
      endif
      if(allocated(global_hfex_real)) then
         deallocate(global_hfex_real)
      end if      
      if(allocated(global_hfex_complex)) then
         deallocate(global_hfex_complex)
      end if

      return

      end subroutine evaluate_single_excitation_correction_p0
