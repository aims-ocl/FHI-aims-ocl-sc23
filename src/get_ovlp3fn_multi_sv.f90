!****s* FHI-aims/get_ovlp3fn_multi_sv
!  NAME
!   get_ovlp3fn_multi_sv
!  SYNOPSIS
      subroutine get_ovlp3fn_multi_sv &
                 ( invs_times_sqrtv, &
                   ovlp_3fn &
                  )
!  PURPOSE
!  The routine multiply the 3-center overlap matrix with the
!  inverse square root of the coulomb matirx (invs_times_sqrtv)
!
!  USES

      use dimensions
      use mpi_tasks
      use localorb_io, only: localorb_info, use_unit, OL_norm
      use prodbas
      use synchronize_mpi

      implicit none

! ARGUMENTS
      real*8, dimension(n_basbas, n_loc_prodbas) :: invs_times_sqrtv
      real*8, dimension(n_basis_pairs, n_loc_prodbas) :: ovlp_3fn

! INPUTS
!  o ovlp_3fn -- real array, the three-center overlap integral (O integral) over
!           two NAO basis functions and one auxiliary basis function. This is the
!           central quantity of the entire formalism of post-Hartree-Fock calculation.
!           Later on, the is a transformation from ovlp_3fn to ovlp_3KS, the latter
!           being the the integral over two single-particle orbital (KS or HF) and
!           one auxiliary basis.
!  o invs_times_sqrtv -- real array, this is actually the inverted overlap matrix
!           times the squre root of the Coulomb matrix
! OUTPUTS
!  o ovlp_3fn -- real array, the three center integral times the the square root of
!            the inverse Coulomb interaction, O * S^(-1)*V^{1/2}
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
      real*8, allocatable, dimension(:,:) :: sync_ovlp_matr

!  counter
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_basis_3
      integer :: i_index_1
      integer :: i_index_2
      integer :: i_index
      integer :: i_prodbas_1
      integer :: n_compute
      integer :: i_task

!  start to work

      call localorb_info('Multiplying ovlp_3fn x (S^-1 V^0.5) (lapack)', &
      &                  use_unit, '(2X,A)', OL_norm)

      allocate (temp_ovlp_matr(n_basis, n_loc_prodbas))
      allocate (sync_ovlp_matr(n_basis, n_basbas))

      i_index_1 = 0
      i_index_2 = 0
      i_task = myid+1
      do i_basis_1 =1, n_basis, 1

          temp_ovlp_matr (:,:) = 0
          n_compute = n_nghr_basis(i_basis_1)

          do i_basis_2 = 1, n_compute, 1
            i_index_1 = i_index_1 + 1

            temp_ovlp_matr(i_basis_2,:) = &
               ovlp_3fn(i_index_1,:)
          enddo

          sync_ovlp_matr(:,:) = 0.d0
          do i_prodbas_1 = 1, n_loc_prodbas, 1
             i_index = map_prodbas(i_prodbas_1,i_task) 
             if(i_index .gt. 0 ) then
               sync_ovlp_matr(:,i_index) = temp_ovlp_matr(:,i_prodbas_1)
             endif
          enddo

          call sync_matrix(sync_ovlp_matr,n_basis,n_basbas)

          temp_ovlp_matr (:,:) = 0
          call dgemm('N', 'N', n_compute, n_loc_prodbas, &
                      n_basbas, 1.0d0, &
                      sync_ovlp_matr, n_basis, &
                      invs_times_sqrtv, n_basbas, 0.d0,&
                      temp_ovlp_matr, n_basis &
                     )

          do i_basis_2 = 1, n_compute, 1

            i_index_2 = i_index_2 + 1
            ovlp_3fn(i_index_2,:) = &
              temp_ovlp_matr(i_basis_2,:)
          enddo
!  end of loop over i_basis_1
      enddo

!      if(myid.eq.0) then
!       do i_basis_3 =1, n_loc_prodbas, 1
!         do i_basis_1 =1, 1, 1
!          i_index_1 = 0
!          do i_basis_2 = 1, i_basis_1, 1
!          i_index_1 = basis_nghr(i_basis_2,i_basis_1)
!           if(i_index_1 .gt. 0) then
!             write(use_unit, '(4X,I5,5X,I5,5X,I5,5X,F16.10)') &
!                 i_basis_1, i_basis_2, i_basis_3, &
!                 ovlp_3fn(i_index_1, i_basis_3)
!           endif
!         enddo
!        enddo
!       enddo
!      endif


      deallocate (temp_ovlp_matr)
      deallocate (sync_ovlp_matr)

      end subroutine get_ovlp3fn_multi_sv
!******
