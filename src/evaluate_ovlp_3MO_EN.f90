!****s* FHI-aims/evaluate_ovlp_3MO_EN
!  NAME
!   evaluate_ovlp_3MO_EN
!  SYNOPSIS

      subroutine evaluate_ovlp_3MO_EN &
               ( KS_eigenvector, ovlp_3KS, &
                 n_lumo_min, n_lumo )
!  PURPOSE
!  Subroutine evaluate_3KS_ovlp_matrix evaluates the three-KS-orbital overlap
!  integral, produced by matrix multiplication of the KS eigenvectors
!  and 3-basis integrals.
!  S^i_jk = \sum_lmn c_li c_mj c_nk O^l_mn
!
!  USES
      use dimensions
      use runtime_choices
      use prodbas
      use hartree_fock
      use mpi_tasks
      use mpi_utilities

      implicit none

!  ARGUMENTS

      integer :: n_lumo_min
      integer :: n_lumo(n_spin)
      real*8 KS_eigenvector(n_basis,n_states,n_spin)

      !real*8 ovlp_3KS( n_loc_prodbas,n_lumo_min:n_states, &
      !                 i_start_mp2:n_homo_max,n_spin )
      real*8 ovlp_3KS( n_loc_prodbas,i_start_mp2:n_states, &
                       i_start_mp2:n_states,n_spin )
!     local variables
      real*8, dimension(:,:), allocatable :: aux_ovlp_matr
      real*8, dimension(:,:), allocatable :: aux_ovlp_matr_2
      real*8, dimension(:), allocatable :: ovlp_matr

      character*20 filename

!     counters
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_basis_3
      integer :: i_index
      integer :: i_spin

!     error indicator
      integer :: mpierr
!  INPUTS
!   o n_lumo -- LUMO level
!   o KS_eigenvector -- Kohn-Sham eigenvectors
!   o n_lumo_min -- lower LUMO level
!  OUTPUT 
!   o ovlp_3KS -- overlap 3 centers Kohn-Sham integrals
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society.
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

!     begin work


      allocate(aux_ovlp_matr(n_basis,n_basis),stat=i_index)
      call check_allocation(i_index, 'aux_ovlp_matr                 ')

      !allocate(aux_ovlp_matr_2(n_basis,n_lumo_min:n_states),stat=i_index)
      allocate(aux_ovlp_matr_2(n_basis,i_start_mp2:n_states),stat=i_index)
      call check_allocation(i_index, 'aux_ovlp_matr_2               ')



      if(use_ovlp_swap) then

        if(use_mpi) then
         write(filename,'(A,I0)') 'OVLP', myid
        else
         write(filename,'(A)') 'OVLP'
        endif

        allocate (ovlp_matr(n_basis_pairs),stat=i_index)
        call check_allocation(i_index, 'ovlp_matr                     ')


        open (60,file=filename,form='unformatted', status='old', &
                 access ='direct',recl=8*n_basis_pairs)
      endif

      if(use_mpi) then
         call MPI_Bcast(KS_eigenvector, &
              n_basis*n_states*n_spin*n_k_points, &
              MPI_DOUBLE_PRECISION, &
              0, mpi_comm_global, mpierr)
      endif

      do i_basis_1 = 1, n_loc_prodbas, 1

         if(use_ovlp_swap) then
           read(60,rec=i_basis_1) ovlp_matr(:)
         endif

         aux_ovlp_matr(:,:) = 0.d0
         do i_basis_2 = 1, n_basis, 1
           do i_basis_3 = 1, i_basis_2, 1

             i_index = basis_nghr(i_basis_3, i_basis_2)
             if(i_index .gt. 0) then

               if(use_ovlp_swap) then
                 aux_ovlp_matr(i_basis_3, i_basis_2) = &
                    ovlp_matr(i_index)
               else
                 aux_ovlp_matr(i_basis_3, i_basis_2) = &
                    ovlp_3fn(i_index, i_basis_1 )
               endif

             endif
            enddo
         enddo

         do i_spin = 1, n_spin

           aux_ovlp_matr_2(:,:) = 0.d0
           call dsymm('L', 'U', n_basis, &
                      !n_states-n_lumo_min+1, 1.0d0, &
                      n_states-i_start_mp2+1, 1.0d0, &
                      aux_ovlp_matr, n_basis, &
                      KS_eigenvector(:,i_start_mp2:n_states,i_spin), &
                      n_basis, 0.d0, &
                      aux_ovlp_matr_2(:,i_start_mp2:n_states), &
                      n_basis &
                      )

           !call dgemm('T', 'N', n_states-n_lumo_min+1, &
           !            n_homo_max-i_start_mp2+1, n_basis, 1.0d0, &
           !            aux_ovlp_matr_2(:,n_lumo_min:n_states), &
           !            n_basis, &
           !            KS_eigenvector(:,i_start_mp2:n_homo_max,i_spin), &
           !            n_basis, 0.d0, &
           !            ovlp_3KS(i_basis_1,n_lumo_min:n_states, &
           !                     i_start_mp2:n_homo_max,i_spin), &
           !            n_states-n_lumo_min+1 &
           !           )
           call dgemm('T', 'N', n_states-i_start_mp2+1, &
                       n_states-i_start_mp2+1, n_basis, 1.0d0, &
                       aux_ovlp_matr_2(:,i_start_mp2:n_states), &
                       n_basis, &
                       KS_eigenvector(:,i_start_mp2:n_states,i_spin), &
                       n_basis, 0.d0, &
                       ovlp_3KS(i_basis_1,i_start_mp2:n_states, &
                                i_start_mp2:n_states,i_spin), &
                       n_states-i_start_mp2+1 &
                      )

       enddo
      enddo

      close(60,status='delete')

      if (use_ovlp_swap) then
        deallocate(ovlp_matr)
      endif
      deallocate(aux_ovlp_matr_2)
      deallocate(aux_ovlp_matr)

      end subroutine evaluate_ovlp_3MO_EN
!---------------------------------------------------------------------
!****** 
