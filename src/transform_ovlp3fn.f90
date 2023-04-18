!****s*  FHI-aims/transform_ovlp3fn
!  NAME
!    transform_ovlp3fn
!  SYNOPSIS
   subroutine transform_ovlp3fn &
         (n_KS_states, KS_eigenvector, ovlp_3fn, &
          ovlp_3KS &
          )

!  PURPOSE
!  Subroutine evaluate_3KS_ovlp_matrix evaluates the three-KS-orbital overlap
!  integral, produced by matrix multiplication of the KS eigenvectors
!  and 3-basis integrals.
!  S^l_jk = \sum_mn c_mj c_nk O^l_mn
!
! USES

      use dimensions
      use runtime_choices
      use prodbas
      use mpi_tasks
      use mpi_utilities
      implicit none

!  ARGUMENTS

      integer :: n_KS_states
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin)
      real*8, intent(IN) :: ovlp_3fn(n_basis_pairs, n_loc_prodbas)
      real*8  :: ovlp_3KS(n_loc_prodbas,n_states,n_KS_states,n_spin)

! INPUTS  
!  o n_KS_states -- the number of  KS states used in the transformation (and in later calculations) 
!   KS_eigenvector -- KS eigenvector
!   ovlp_3fn -- 3-function overlap integrals (orthonormalized)
!
! OUTPUT 
!  o ovlp_3KS -- transformed 3-function overlap integrals
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



      real*8, dimension(:,:), allocatable :: aux_ovlp_matr
      real*8, dimension(:,:), allocatable :: aux_ovlp_matr_2
      real*8, dimension(:), allocatable   :: ovlp_matr

      character*20 filename

!     counters
      integer :: i_state
      integer :: j_state

      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_basis_3
      integer :: i_index
      integer :: i_spin

!     error variable
      integer :: mpierr

!     begin work



      allocate(aux_ovlp_matr(n_basis,n_basis))
      allocate(aux_ovlp_matr_2(n_basis,n_states))


      if(use_ovlp_swap) then

        if(use_mpi) then
         write(filename,'(A,I0)') 'OVLP', myid
        else
         write(filename,'(A)') 'OVLP'
        endif

        allocate (ovlp_matr(n_basis_pairs))

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
                      n_states, 1.0d0, &
                      aux_ovlp_matr, n_basis, &
                      KS_eigenvector(:,:,i_spin), &
                      n_basis, 0.d0, &
                      aux_ovlp_matr_2(:,:), &
                      n_basis &
                      )

           call dgemm('T', 'N', n_states, &
                    n_KS_states, n_basis, 1.0d0, &
                    aux_ovlp_matr_2, &
                    n_basis, &
                    KS_eigenvector(:,1:n_KS_states,i_spin), &
                    n_basis, 0.d0, &
                    ovlp_3KS(i_basis_1,1:n_states, &
                             1:n_KS_states,i_spin), &
                    n_states &
                   )


         enddo
      enddo

      if (use_ovlp_swap) close(60,status='delete')

      if(allocated(ovlp_matr))then
        deallocate(ovlp_matr)
      endif
      deallocate(aux_ovlp_matr)
      deallocate(aux_ovlp_matr_2)

      end subroutine transform_ovlp3fn
!---------------------------------------------------------------------
!******
!===================================================================================================
!****s*  FHI-aims/transform_ovlp3fn
!  NAME
!    transform_ovlp3fn
!  SYNOPSIS
   subroutine transform_ovlp3fn_2 &
         (n_KS_states, KS_eigenvector, ovlp_3fn, &
          ovlp_3KS &
          )

!  PURPOSE
!  Subroutine evaluate_3KS_ovlp_matrix evaluates the three-KS-orbital overlap
!  integral, produced by matrix multiplication of the KS eigenvectors
!  and 3-basis integrals.
!  S^l_jk = \sum_mn c_mj c_nk O^l_mn
!
! USES

      use dimensions
      use runtime_choices
      use prodbas
      use mpi_tasks
      use mpi_utilities
      implicit none

!  ARGUMENTS

      integer :: n_KS_states
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin)
      real*8, intent(IN) :: ovlp_3fn(n_basis_pairs, n_loc_prodbas)
      real*8  :: ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)

! INPUTS  
!  o n_KS_states -- the number of  KS states used in the transformation (and in later calculations) 
!   KS_eigenvector -- KS eigenvector
!
! OUTPUT 
!  o ovlp_3KS -- transformed 3-function overlap integrals
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


      real*8, dimension(:,:), allocatable :: aux_ovlp_matr
      real*8, dimension(:), allocatable   :: ovlp_col
      real*8, dimension(:), allocatable   :: dist_ovlp_col
      real*8, dimension(:,:,:), allocatable :: aux1
      real*8, dimension(:,:), allocatable :: aux2, aux3
      real*8, dimension(:,:), allocatable :: aux_eigenvector_1(:,:,:)
      real*8, dimension(:,:), allocatable :: aux_eigenvector_2(:,:,:)

!     counters
      integer :: i_state
      integer :: j_state

      integer :: i_basis_2
      integer :: i_basis_3
      integer :: i_index
      integer :: i_spin
      integer :: i_state_loc
      integer :: j_state_loc

      integer :: ip_row, ip_col
      integer :: n_blk, n_col, i_basbas, i_glob


!     error variable
      integer :: mpierr

!     begin work



      allocate(aux_ovlp_matr(n_basis,n_basis))

      allocate(ovlp_col(n_basis_pairs))
      allocate(dist_ovlp_col(n_basis_pairs))

      allocate(aux1(ndim1_o3KS,n_basis,n_spin))
      allocate(aux2(ndim1_o3KS,n_basis))
      allocate(aux3(ndim1_o3KS,ndim2_o3KS))

      allocate(aux_eigenvector_1(n_basis,ndim1_o3KS,n_spin))
      allocate(aux_eigenvector_2(n_basis,ndim2_o3KS,n_spin))

      if(use_mpi) then
         call MPI_Bcast(KS_eigenvector, &
              n_basis*n_states*n_spin*n_k_points, &
              MPI_DOUBLE_PRECISION, &
              0, mpi_comm_global, mpierr)
      endif


      ! The last column of aux_eigenvector_1/2 may not be set by the code below
      ! so set these arrays to 0 before starting
      aux_eigenvector_1 = 0
      aux_eigenvector_2 = 0

      do i_spin = 1, n_spin

        do j_state = 1, n_states

          if(own_dim1_o3ks(j_state) /= myp1_o3KS) cycle
          j_state_loc = loc_dim1_o3ks(j_state)

          aux_eigenvector_1(:,j_state_loc,i_spin) = KS_eigenvector(:,j_state,i_spin)
        enddo

        do i_state = 1, n_KS_states

          if(own_dim2_o3ks(i_state) /= myp2_o3KS) cycle
          i_state_loc = loc_dim2_o3ks(i_state)

          aux_eigenvector_2(:,i_state_loc,i_spin) = KS_eigenvector(:,i_state,i_spin)
        enddo

      enddo


      ovlp_3KS = 0


      do n_blk=0,(n_basbas-1)/(nprow_aux_2d*npcol_aux_2d*nb_aux) ! blocks of columns on ovlp3_fn
        do ip_row=0,nprow_aux_2d-1 ! processor rows within one processor column
          do n_col=1,nb_aux ! columns of ovlp3_fn within one block

            ! get global column number of current column
            i_glob = global_id(ip_row,mypcol_aux_2d)
            i_basbas = n_col + i_glob*nb_aux + n_blk*nprow_aux_2d*npcol_aux_2d*nb_aux

            IF(i_basbas<=n_basbas) then ! if column exists

              ! Broadcast column within processor row

              if(ip_row==myprow_aux_2d) ovlp_col(:) = ovlp_3fn(:,n_col+n_blk*nb_aux)
              call MPI_Bcast(ovlp_col,n_basis_pairs,MPI_DOUBLE_PRECISION,ip_row, mpi_comm_rows_aux_2d, mpierr)

              ! form aux_ovlp_matr and multiply with aux_eigenvector_1 from left

              aux_ovlp_matr(:,:) = 0.d0
              do i_basis_2 = 1, n_basis, 1
                do i_basis_3 = 1, i_basis_2, 1

                  i_index = basis_nghr(i_basis_3, i_basis_2)
                  if(i_index .gt. 0) then
                    aux_ovlp_matr(i_basis_3, i_basis_2) = ovlp_col(i_index)
                    aux_ovlp_matr(i_basis_2, i_basis_3) = ovlp_col(i_index)
                  endif
                enddo
              enddo

              do i_spin = 1, n_spin
                call dgemm('T', 'N', ndim1_o3KS, n_basis, n_basis, 1.d0, &
                        aux_eigenvector_1(1,1,i_spin), ubound(aux_eigenvector_1,1), &
                        aux_ovlp_matr,ubound(aux_ovlp_matr,1),0.d0,aux1(1,1,i_spin),ubound(aux1,1))
              enddo
            ENDIF

            ! Broadcast "half ready" aux1 to all processor columns for multiplying with aux_eigenvector_2
            ! This is done by every processor column which still has work

            do ip_col = 0, npcol_aux_2d-1

              i_glob = global_id(ip_row,ip_col)
              i_basbas = n_col + i_glob*nb_aux + n_blk*nprow_aux_2d*npcol_aux_2d*nb_aux

              if(i_basbas <= n_basbas) then

                do i_spin = 1, n_spin
                  if(ip_col==mypcol_aux_2d) aux2(:,:) = aux1(:,:,i_spin)
                  call MPI_Bcast(aux2, ndim1_o3KS*n_basis, MPI_DOUBLE_PRECISION, ip_col, mpi_comm_cols_aux_2d, mpierr)
                  call dgemm('N', 'N', ndim1_o3KS, ndim2_o3KS, n_basis, 1.d0, aux2, ubound(aux2,1), &
                          aux_eigenvector_2(1,1,i_spin), ubound(aux_eigenvector_2,1), &
                          0.d0, aux3, ubound(aux3,1))
                  ovlp_3KS(i_basbas,:,:,i_spin) = aux3(:,:)
                enddo
              endif
            enddo

          enddo
        enddo
      enddo

      deallocate(ovlp_col)
      deallocate(dist_ovlp_col)
      deallocate(aux_ovlp_matr)
      deallocate(aux1)
      deallocate(aux2)
      deallocate(aux3)
      deallocate(aux_eigenvector_1)
      deallocate(aux_eigenvector_2)

      end subroutine transform_ovlp3fn_2
!---------------------------------------------------------------------
!******
