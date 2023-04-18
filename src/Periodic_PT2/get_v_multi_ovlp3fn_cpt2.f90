!****s* FHI-aims/get_v_multi_ovlp3fn_complex
!  NAME
!   get_v_multi_ovlp3fn_complex
!  SYNOPSIS
      subroutine get_v_multi_ovlp3fn_complex &
                 ( coulomb_matr, &
                   ovlp_3fn, iop )
!  PURPOSE
!  The routine multiply the 3-center overlap matrix with the
!  square root of the complex coulomb matirx (coulomb_matr)
!
!  USES

      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi

      implicit none

! ARGUMENTS
      complex*16, dimension(n_basbas, n_basbas) :: coulomb_matr
      complex*16, dimension(n_basbas,n_basis,n_basis) :: ovlp_3fn
      integer :: iop

! INPUTS
!  o ovlp_3fn -- complex array, the three-center overlap integral (O integral) over
!           two NAO basis functions and one auxiliary basis function. This is the
!           central quantity of the entire formalism of post-Hartree-Fock calculation.
!           Later on, there is a transformation from ovlp_3fn to ovlp_3KS, the latter
!           being the the integral over two single-particle orbital (KS or HF) and
!           one auxiliary basis. O == (ij|\nu) * V_{LVL}^{-1} * V^{1/2} (IOP=1) 
!                             or   == V^{1/2} * V_{LVL}^{-1} * (\nu|ij) (IOP=2)
!  o coulomb_matr -- complex array, this is actually the square root of the complex
!           Coulomb matrix. V^(1/2) == (\nu|\mu)^(1/2)
!  o iop -- 1 :: O == (ij|\nu) * V_{LVL}^{-1} * V^{1/2}
!           2 :: O == V^{1/2} * V_{LVL}^{-1} * (\nu|ij)
!
! OUTPUTS
!  o ovlp_3fn -- complex array, the three center integral times the the square root of
!            the complex Coulomb interaction, which is the key to reduce the
!            scaling of periodic-PT2 from N^6 to N^5
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
      complex*16, allocatable, dimension(:,:) :: temp_ovlp_matr_a
      complex*16, allocatable, dimension(:,:) :: temp_ovlp_matr_b

!  counter
      integer :: i_basis_1

!  start to work

      allocate (temp_ovlp_matr_a(n_basbas, n_basis))
      !if (iop .eq. 1) then
      !    allocate (temp_ovlp_matr_b(n_basis, n_basbas))
      !endif

      do i_basis_1 =1, n_basis, 1


          temp_ovlp_matr_a(:,:) = ovlp_3fn(:,:,i_basis_1)

          if (iop .eq. 1) then
            call zgemm('T', 'N', n_basbas, n_basis, &
                        n_basbas, (1.0d0,0.0d0), &
                        coulomb_matr, n_basbas, &
                        temp_ovlp_matr_a, n_basbas, (0.d0,0.d0), &
                        ovlp_3fn(1,1,i_basis_1), n_basbas &
                       )
            !call zgemm('T', 'N', n_basis, n_basbas, &
            !            n_basbas, (1.0d0,0.0d0), &
            !            temp_ovlp_matr_a, n_basbas, &
            !            coulomb_matr, n_basbas, (0.d0,0.d0), &
            !            temp_ovlp_matr_b, n_basis &
            !           )

            !ovlp_3fn(:,:,i_basis_1) = transpose(temp_ovlp_matr_b)
          else
            call zgemm('N', 'N', n_basbas, n_basis, &
                        n_basbas, (1.0d0,0.0d0), &
                        coulomb_matr, n_basbas, &
                        temp_ovlp_matr_a, n_basbas, (0.d0,0.d0), &
                        ovlp_3fn(1,1,i_basis_1), n_basbas &
                       )
          endif
      enddo

      deallocate (temp_ovlp_matr_a)
      !if (allocated(temp_ovlp_matr_b)) deallocate (temp_ovlp_matr_b)
      end subroutine get_v_multi_ovlp3fn_complex
!******

      subroutine get_v_multi_ovlp3fn_complex_blacs &
                 ( n_low_state, n_homo_max_curr, n_homo_curr, coulomb_matr, &
                   ovlp_3fn, iop )
!  PURPOSE
!  The routine multiply the 3-center overlap matrix with the
!  square root of the complex coulomb matirx (coulomb_matr)
!
!  USES

      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use cpt2_blacs

      implicit none

! ARGUMENTS
      integer :: n_low_state, n_homo_max_curr, n_homo_curr, iop
      complex*16, intent(in), dimension(lbb_row:ubb_row,lbb_col:ubb_col) :: coulomb_matr
      complex*16, intent(inout), dimension(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max_curr), target :: ovlp_3fn

! INPUTS
!  o ovlp_3fn -- complex array, the three-center overlap integral (O integral) over
!           two NAO basis functions and one auxiliary basis function. This is the
!           central quantity of the entire formalism of post-Hartree-Fock calculation.
!           Later on, there is a transformation from ovlp_3fn to ovlp_3KS, the latter
!           being the the integral over two single-particle orbital (KS or HF) and
!           one auxiliary basis. O == (ij|\nu) * V_{LVL}^{-1} * V^{1/2} (IOP=1) 
!                             or   == V^{1/2} * V_{LVL}^{-1} * (\nu|ij) (IOP=2)
!  o coulomb_matr -- complex array, this is actually the square root of the complex
!           Coulomb matrix. V^(1/2) == (\nu|\mu)^(1/2)
!  o iop -- 1 :: O == (ij|\nu) * V_{LVL}^{-1} * V^{1/2}
!           2 :: O == V^{1/2} * V_{LVL}^{-1} * (\nu|ij)
!
! OUTPUTS
!  o ovlp_3fn -- complex array, the three center integral times the the square root of
!            the complex Coulomb interaction, which is the key to reduce the
!            scaling of periodic-PT2 from N^6 to N^5
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
      complex*16, dimension(lbb_row:ubb_row,lpb_col:upb_col), target :: temp_ovlp_matr_a

      complex*16, dimension(:,:), pointer :: ptr_ovlp_matr_a
      complex*16, dimension(:,:), pointer :: ptr_ovlp_matr_b

!  counter
      integer :: i_basis_1

!  start to work

      do i_basis_1 =n_low_state, n_homo_curr, 1


          temp_ovlp_matr_a(lbb_row:ubb_row, lpb_col:upb_col) = &
              ovlp_3fn(lbb_row:ubb_row, lpb_col:upb_col,i_basis_1)

          ptr_ovlp_matr_a => temp_ovlp_matr_a(lbb_row:ubb_row, lpb_col:upb_col)
          ptr_ovlp_matr_b => ovlp_3fn(lbb_row:ubb_row, lpb_col:upb_col,i_basis_1)

          if (iop .eq. 1) then
            call pzgemm('T', 'N', n_basbas, n_unocc_max, &
                        n_basbas, (1.0d0,0.0d0), &
                        coulomb_matr, 1,1, bb2desc, &
                        ptr_ovlp_matr_a, 1,1, pb2desc, (0.d0,0.d0), &
                        ptr_ovlp_matr_b, 1,1, pb2desc &
                       )
                

          else
            call pzgemm('N', 'N', n_basbas, n_unocc_max, &
                        n_basbas, (1.0d0,0.0d0), &
                        coulomb_matr, 1,1, bb2desc, &
                        ptr_ovlp_matr_a, 1,1, pb2desc, (0.d0,0.d0), &
                        ptr_ovlp_matr_b, 1,1, pb2desc &
                       )
          endif
      enddo

      end subroutine get_v_multi_ovlp3fn_complex_blacs

      subroutine get_v_multi_ovlp3fn_real_blacs &
                 ( n_low_state, n_homo_max_curr, n_homo_curr, coulomb_matr, &
                   ovlp_3fn, iop )
!  PURPOSE
!  The routine multiply the 3-center overlap matrix with the
!  square root of the complex coulomb matirx (coulomb_matr)
!
!  USES

      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use cpt2_blacs

      implicit none

! ARGUMENTS
      integer :: n_low_state, n_homo_max_curr, n_homo_curr, iop
      real*8, intent(in), dimension(lbb_row:ubb_row,lbb_col:ubb_col) :: coulomb_matr
      real*8, intent(inout), &
          dimension(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max_curr), target :: ovlp_3fn

! INPUTS
!  o ovlp_3fn -- real array, the three-center overlap integral (O integral) over
!           two NAO basis functions and one auxiliary basis function. This is the
!           central quantity of the entire formalism of post-Hartree-Fock calculation.
!           Later on, there is a transformation from ovlp_3fn to ovlp_3KS, the latter
!           being the the integral over two single-particle orbital (KS or HF) and
!           one auxiliary basis. O == (ij|\nu) * V_{LVL}^{-1} * V^{1/2} (IOP=1) 
!                             or   == V^{1/2} * V_{LVL}^{-1} * (\nu|ij) (IOP=2)
!  o coulomb_matr -- real array, this is actually the square root of the coulomb matrix
!           Coulomb matrix. V^(1/2) == (\nu|\mu)^(1/2)
!  o iop -- 1 :: O == (ij|\nu) * V_{LVL}^{-1} * V^{1/2}
!           2 :: O == V^{1/2} * V_{LVL}^{-1} * (\nu|ij)
!
! OUTPUTS
!  o ovlp_3fn -- real array, the three center integral times the the square root of
!            the real Coulomb interaction, which is the key to reduce the
!            scaling of periodic-PT2 from N^6 to N^5
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
      real*8, dimension(lbb_row:ubb_row,lpb_col:upb_col), target :: temp_ovlp_matr_a

      real*8, dimension(:,:), pointer :: ptr_ovlp_matr_a
      real*8, dimension(:,:), pointer :: ptr_ovlp_matr_b

!  counter
      integer :: i_basis_1

!  start to work

      do i_basis_1 =n_low_state, n_homo_curr, 1


          temp_ovlp_matr_a(lbb_row:ubb_row, lpb_col:upb_col) = &
              ovlp_3fn(lbb_row:ubb_row, lpb_col:upb_col,i_basis_1)

          ptr_ovlp_matr_a => temp_ovlp_matr_a(lbb_row:ubb_row, lpb_col:upb_col)
          ptr_ovlp_matr_b => ovlp_3fn(lbb_row:ubb_row, lpb_col:upb_col,i_basis_1)

          if (iop .eq. 1) then
            call pdgemm('T', 'N', n_basbas, n_unocc_max, &
                        n_basbas, (1.0d0,0.0d0), &
                        coulomb_matr, 1,1, bb2desc, &
                        ptr_ovlp_matr_a, 1,1, pb2desc, (0.d0,0.d0), &
                        ptr_ovlp_matr_b, 1,1, pb2desc &
                       )
                

          else
            call pdgemm('N', 'N', n_basbas, n_unocc_max, &
                        n_basbas, (1.0d0,0.0d0), &
                        coulomb_matr, 1,1, bb2desc, &
                        ptr_ovlp_matr_a, 1,1, pb2desc, (0.d0,0.d0), &
                        ptr_ovlp_matr_b, 1,1, pb2desc &
                       )
          endif
      enddo

      end subroutine get_v_multi_ovlp3fn_real_blacs
!******
