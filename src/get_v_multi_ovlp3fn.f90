!****s* FHI-aims/get_v_multi_ovlp3fn
!  NAME
!   get_v_multi_ovlp3fn
!  SYNOPSIS
      subroutine get_v_multi_ovlp3fn &
                 ( inv_coulomb_matr, &
                   ovlp_3fn &
                  )
!  PURPOSE
!  The routine multiply the 3-center overlap matrix with the
!  inverse square root of the coulomb matirx (inv_coulomb_matr)
!
!  USES

      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: localorb_info, use_unit, OL_norm

      implicit none

! ARGUMENTS
      real*8, dimension(n_basbas, n_loc_prodbas) :: inv_coulomb_matr
      real*8, dimension(n_basis_pairs, n_loc_prodbas) :: ovlp_3fn

! INPUTS
!  o ovlp_3fn -- real array, the three-center overlap integral (O integral) over
!           two NAO basis functions and one auxiliary basis function. This is the
!           central quantity of the entire formalism of post-Hartree-Fock calculation.
!           Later on, there is a transformation from ovlp_3fn to ovlp_3KS, the latter
!           being the integral over two single-particle orbital (KS or HF) and
!           one auxiliary basis. O == (ij|\nu)
!  o inv_coulomb_matr -- real array, this is actually the square root of the inverse
!           Coulomb matrix. V^(-1/2) == (\nu|\mu)^(-1/2)
! OUTPUTS
!  o ovlp_3fn -- real array, the three center integral times the the square root of
!            the inverse Coulomb interaction, O * V^{-1/2}
!
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
      real*8, allocatable, dimension(:,:) :: temp_ovlp_matr
      real*8, allocatable, dimension(:,:) :: aux_ovlp_matr

!  counter
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_basis_3
      integer :: i_index_1
      integer :: i_index_2
      integer :: i_prodbas_1
      integer :: n_compute
      integer :: i_task

!  start to work

      call localorb_info('Multiplying V^-0.5 x ovlp_3fn (1d-scalapack)', &
      &                  use_unit, '(2X,A)', OL_norm)

      allocate (temp_ovlp_matr(n_loc_prodbas, n_basis))
      allocate (aux_ovlp_matr(n_basbas, n_basis))

      i_index_1 = 0
      i_index_2 = 0
      do i_basis_1 =1, n_basis, 1

          temp_ovlp_matr (:,:) = 0
          n_compute = n_nghr_basis(i_basis_1)
          do i_basis_2 = 1, n_compute, 1
            i_index_1 = i_index_1 + 1

            temp_ovlp_matr(:,i_basis_2) = &
               ovlp_3fn(i_index_1,:)
          enddo

          aux_ovlp_matr (:,:) = 0
          ! mark by igor for the illegal case with n_loc_prodbas=0
          if (n_loc_prodbas .gt. 0) then
              call dgemm('N', 'N', n_basbas, n_compute, &
                          n_loc_prodbas, 1.0d0, &
                          inv_coulomb_matr, n_basbas, &
                          temp_ovlp_matr, n_loc_prodbas, 0.d0, &
                          aux_ovlp_matr, n_basbas &
                         )
          end if

          call sync_matrix(aux_ovlp_matr, n_basbas, n_compute)

          temp_ovlp_matr(:,:) =0.d0
          i_prodbas_1 = 0
          do i_task = 1, n_tasks, 1
            do i_basis_2 = 1, n_loc_prodbas, 1

               i_prodbas_1 = map_prodbas(i_basis_2,i_task)
               if(i_prodbas_1.gt.0 .and. myid.eq.i_task-1) then

                  temp_ovlp_matr(i_basis_2, 1:n_compute) = &
                  aux_ovlp_matr(i_prodbas_1, 1:n_compute)
               endif

            enddo
          enddo

          do i_basis_2 = 1, n_compute, 1

            i_index_2 = i_index_2 + 1
            ovlp_3fn(i_index_2,:) = &
              temp_ovlp_matr(:,i_basis_2)
          enddo
!  end of loop over i_basis_1
      enddo

!       do i_basis_3 =1, 1, 1
!         do i_basis_1 =1, n_basis, 1
!          i_index_1 = 0
!          do i_basis_2 = 1, i_basis_1, 1
!          i_index_1 = basis_nghr(i_basis_2,i_basis_1)
!           if(i_index_1 .gt. 0) then
!             write(use_unit, '(4X,I5,5X,I5,5X,I5,5X,F16.10)')
!     +            i_basis_1, i_basis_2, i_basis_3,
!     +            ovlp_3fn(i_index_1, i_basis_3)
!           endif
!         enddo
!        enddo
!       enddo


      deallocate (temp_ovlp_matr)
      deallocate (aux_ovlp_matr)

      end subroutine get_v_multi_ovlp3fn
!******
! --------------------------------------------------------------------------------------------------
!****s* FHI-aims/get_v_multi_ovlp3fn
!  NAME
!   get_v_multi_ovlp3fn
!  SYNOPSIS
      subroutine get_v_multi_ovlp3fn_2 &
                 ( inv_coulomb_matr, &
                   ovlp_3fn &
                  )
!  PURPOSE
!  The routine multiply the 3-center overlap matrix with the
!  inverse square root of the coulomb matrix (inv_coulomb_matr)
!
!  USES

      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: localorb_info, use_unit, OL_norm
      implicit none

! ARGUMENTS
      real*8, dimension(max_row_2d, max_col_2d) :: inv_coulomb_matr
      real*8, dimension(n_basis_pairs, n_loc_prodbas) :: ovlp_3fn

! INPUTS
!  o ovlp_3fn -- real array, the three-center overlap integral (O integral) over
!           two NAO basis functions and one auxiliary basis function. This is the
!           central quantity of the entire formalism of post-Hartree-Fock calculation.
!           Later on, the is a transformation from ovlp_3fn to ovlp_3KS, the latter
!           being the the integral over two single-particle orbital (KS or HF) and
!           one auxiliary basis. O == (ij|\nu)
!  o inv_coulomb_matr -- real array, this is actually the square root of the inverse
!           Coulomb matrix. V^(-1/2) == (\nu|\mu)^(-1/2)
! OUTPUTS
!  o ovlp_3fn -- real array, the three center integral times the the square root of
!            the inverse Coulomb interaction, O * V^{-1/2}
!
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
      real*8, allocatable, dimension(:,:) :: aux1, aux2, aux3
      integer, allocatable, dimension(:) :: send_cnt, send_loc, send_dsp, idx
      integer, allocatable, dimension(:) :: recv_cnt, recv_loc, recv_dsp

      ! Please note: If the auxilliary arrays should need too much memory,
      ! n_strip can safely be reduced. PDGEMM performance might suffer, however!!!

      integer, parameter :: n_strip = 2048 ! strip size for stripmining multiplication

      integer :: n_off, n_off1, n_len1
      integer :: j, np1, npc, npr, npg, ns, nr, n_col, n_rows, mpierr, desc_ovlp_3fn(9)

      call localorb_info('Multiplying ovlp_3fn x V^-0.5 (2d-scalapack)', &
      &                  use_unit, '(2X,A)', OL_norm)

      ! Set up data for mpi_alltoall

      allocate(send_cnt(0:n_tasks-1), send_loc(0:n_tasks-1), send_dsp(0:n_tasks-1))
      allocate(recv_cnt(0:n_tasks-1), recv_loc(0:n_tasks-1), recv_dsp(0:n_tasks-1))
      allocate(idx(0:n_tasks-1))

      ! Get counts how many blocks we send/receive

      recv_cnt(:) = 0
      send_cnt(:) = 0
      do npr = 0, nprow_aux_2d-1
        do j = 0, n_basbas-1
          np1 = MOD(j/nb_aux,n_tasks)         ! Processor in 1D distribution
          npc = MOD(j/nb_aux_2d,npcol_aux_2d) ! Column processor in 2D distribution
          npg = global_id(npr,npc)            ! Global ID in 2D distribution
          if(npc==mypcol_aux_2d .and. npr==myprow_aux_2d) recv_cnt(np1) = recv_cnt(np1) + 1
          if(np1==myid) send_cnt(npg) = send_cnt(npg) + 1
        enddo
      enddo

      ! Set up location where data of remote processor starts

      ns = 0
      nr = 0
      do j=0,n_tasks-1
         send_loc(j) = ns
         recv_loc(j) = nr
         ns = ns + send_cnt(j)
         nr = nr + recv_cnt(j)
      enddo

      ! Scale from rows to variable counts

      send_dsp = send_loc*n_strip
      recv_dsp = recv_loc*n_strip
      send_cnt = send_cnt*n_strip
      recv_cnt = recv_cnt*n_strip

      ! Allocate auxiliary arrays for communictaion/matrix multiply

      allocate(aux1(n_strip, n_loc_prodbas*nprow_aux_2d))
      allocate(aux2(n_strip, max_col_2d))
      allocate(aux3(n_strip, max_col_2d))

      aux1 = 0
      aux2 = 0
      aux3 = 0

      do n_off = 0, n_basis_pairs-1, n_strip*nprow_aux_2d

        ! Redistribute strip of ovlp_3fn from 1D to 2D distribution
        ! =========================================================

        ! Put current strip of ovlp_3fn into send buffers

        idx(:) = send_loc(:)
        do npr = 0, nprow_aux_2d-1

          n_off1 = n_off + npr*n_strip
          if(n_off1>=n_basis_pairs) exit
          n_len1 = MIN(n_strip, n_basis_pairs-n_off1) ! length of current strip

          do j = 1, n_loc_prodbas
            n_col = INDXL2G(j, nb_aux, myid, n_tasks)
            npc = MOD((n_col-1)/nb_aux_2d,npcol_aux_2d) ! Column processor in 2D distribution
            npg = global_id(npr,npc)                    ! Global ID in 2D distribution
            idx(npg) = idx(npg) + 1
            aux1(1:n_len1,idx(npg)) = ovlp_3fn(n_off1+1:n_off1+n_len1,j)
          enddo
        enddo

        ! Redistribute strip of ovlp_3fn using mpi_alltoallv

        call mpi_alltoallv(aux1,send_cnt,send_dsp,MPI_REAL8, &
                           aux2,recv_cnt,recv_dsp,MPI_REAL8, &
                           mpi_comm_global,mpierr)

        ! Sort received data to 2D distribution

        idx(:) = recv_loc(:)
        do j = 1, max_col_2d
          n_col = INDXL2G(j, nb_aux_2d, mypcol_aux_2d, npcol_aux_2d)
          np1 = MOD((n_col-1)/nb_aux,n_tasks) ! Processor in 1D distribution
          idx(np1) = idx(np1) + 1
          aux3(:,j) = aux2(:,idx(np1))
        enddo

        ! Calculate ovlp_3fn * inv_coulomb_matr for current strip
        ! =======================================================

        n_rows = min(n_strip*nprow_aux_2d, n_basis_pairs-n_off)

        call descinit(desc_ovlp_3fn, n_rows, n_basbas, n_strip, nb_aux_2d, &
                      0, 0, my_blacs_ctxt_aux_2d, n_strip, mpierr)

        call pdgemm('N','N',n_rows,n_basbas,n_basbas,1.d0, &
                    aux3, 1, 1, desc_ovlp_3fn,                           &
                    inv_coulomb_matr, 1, 1, aux_sc_desc_2d, 0.d0,        &
                    aux2, 1, 1, desc_ovlp_3fn)

        ! Redistribute result to 1D distribution
        ! ======================================

        ! put 2D distributed ovlp_3fn into send buffers

        aux3(:,:) = aux2(:,:)
        idx(:) = recv_loc(:)
        do j = 1, max_col_2d
          n_col = INDXL2G(j, nb_aux_2d, mypcol_aux_2d, npcol_aux_2d)
          np1 = MOD((n_col-1)/nb_aux,n_tasks) ! Processor in 1D distribution
          idx(np1) = idx(np1) + 1
          aux2(:,idx(np1)) = aux3(:,j)
        enddo

        ! Redistribute strip of ovlp_3fn using mpi_alltoallv
        ! recv_cnt/dsp and send_cnt/dsp have the reverse meaning here!

        call mpi_alltoallv(aux2,recv_cnt,recv_dsp,MPI_REAL8, &
                           aux1,send_cnt,send_dsp,MPI_REAL8, &
                           mpi_comm_global,mpierr)

        ! Sort received data to 1D distribution

        idx(:) = send_loc(:)
        do npr = 0, nprow_aux_2d-1

          n_off1 = n_off + npr*n_strip
          if(n_off1>=n_basis_pairs) exit
          n_len1 = MIN(n_strip, n_basis_pairs-n_off1) ! length of current strip

          do j = 1, n_loc_prodbas
            n_col = INDXL2G(j, nb_aux, myid, n_tasks)
            npc = MOD((n_col-1)/nb_aux_2d,npcol_aux_2d) ! Column processor in 2D distribution
            npg = global_id(npr,npc)                    ! Global ID in 2D distribution
            idx(npg) = idx(npg) + 1
            ovlp_3fn(n_off1+1:n_off1+n_len1,j) = aux1(1:n_len1,idx(npg))
          enddo

        enddo

      enddo

      deallocate (aux1)
      deallocate (aux2)
      deallocate (aux3)
      deallocate (send_cnt, send_loc, send_dsp, idx)
      deallocate (recv_cnt, recv_loc, recv_dsp)

contains
      INTEGER FUNCTION INDXL2G( INDXLOC, NB, IPROC, NPROCS )
      ! from Scalapack TOOLS directory with ISRCPROC omitted
      INTEGER INDXLOC, IPROC, NB, NPROCS
      INDXL2G = NPROCS*NB*((INDXLOC-1)/NB) + MOD(INDXLOC-1,NB) + IPROC*NB + 1
      END FUNCTION INDXL2G

      end subroutine get_v_multi_ovlp3fn_2
!******
