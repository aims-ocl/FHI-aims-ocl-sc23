!****s* FHI-aims/output_rrs_pbc_matrices
!  NAME
!    output_real_matrices
!  SYNOPSIS

subroutine output_rrs_pbc_matrices &
     ( overlap_matrix, hamiltonian )

!  PURPOSE
!  Subroutine output_matrices writes out all basis function properties, and the
!  previously calculated overlap and Hamiltonian matrices.
!
!  Output functionality for complex and real versions of matrices.
!  Split into two subroutines due to variable types, and to adjust
!  formatting.
!
!  USES

  use dimensions
  use basis
  use localorb_io
  implicit none

!  ARGUMENTS
!  imported variables

  real*8 overlap_matrix(n_basis*(n_basis+1)/2)
  real*8 hamiltonian(n_basis*(n_basis+1)/2,n_spin)

!  INPUTS
!   o overlap_matrix -- overlap matrix
!   o hamiltonian -- Hamiltonian matrix
!
!  OUTPUT
!    none
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

      real*8 output_element(n_basis)

      character*18,dimension(n_spin) :: file_name

      character*100 :: info_str

!  counters

      integer i_basis, i_basis_1, i_basis_2, i_index, i_spin, i_fn

!  begin work

      write (info_str,'(2X,A)') &
        "Writing basis function properties and matrices."
      call localorb_info(info_str,use_unit,'(A)')

      write (info_str,'(2X,I5,A)') n_basis, " basis functions used."
      call localorb_info(info_str,use_unit,'(A)')

!  write basis function properties

      open (50, file="basis-indices.out")

      write(info_str,*)
      call localorb_info(info_str,50,'(A)')

      write(info_str,'(A5,1X,A8,1X,A3,1X,A3,1X,A3,1X,A3)') &
        "fn.", "  type  ", "at.", "n", &
        "l", "m"
      call localorb_info(info_str,50,'(A)')

      do i_basis = 1, n_basis, 1
        i_fn = basis_fn(i_basis)

        write(info_str,'(I5,1X,A8,1X,I3,1X,I3,1X,I3,1X,I3)') &
          i_basis, basisfn_type(i_fn), &
          basis_atom(i_basis), basisfn_n(i_fn), basis_l(i_basis), &
          basis_m(i_basis)
        call localorb_info(info_str,50,'(A)')

      enddo

      close (50)

!  now, write the overlap matrix

      open (50, file="overlap-matrix.out")

      i_basis_2 = 1

      do while (i_basis_2.le.n_basis)

        write(info_str,'(8X,5(4X,I4,4X))') &
          (i_basis, i_basis = i_basis_2, i_basis_2+4, 1)
        call localorb_info(info_str,50,'(A)')

        do i_basis_1 = 1, n_basis, 1

            if (i_basis_2.le.(n_basis-4)) then

              do i_basis = i_basis_2, i_basis_2+4,1
                if (i_basis.ge.i_basis_1) then
                  i_index = i_basis_1 + (i_basis-1)*i_basis/2
                  output_element(i_basis) = overlap_matrix(i_index)
                else
                  output_element(i_basis) = 0.
                end if
              end do

              write(info_str,'(I5, 3X, 5(F10.6,2X))') i_basis_1, &
                (output_element(i_basis), &
                 i_basis = i_basis_2, i_basis_2+4, 1)
              call localorb_info(info_str,50,'(A)')

            else

              do i_basis = i_basis_2, n_basis, 1
                if (i_basis.ge.i_basis_1) then
                  i_index = i_basis_1 + (i_basis-1)*i_basis/2
                  output_element(i_basis) = overlap_matrix(i_index)
                else
                  output_element(i_basis) = 0.
                end if
              end do

              write(info_str,'(I5, 3X, 5(F10.6,2X))') i_basis_1, &
                (output_element(i_basis), &
                 i_basis = i_basis_2, n_basis, 1)
              call localorb_info(info_str,50,'(A)')

            end if

        enddo

        i_basis_2 = i_basis_2+5
        write(info_str,*)
        call localorb_info(info_str,50,'(A)')

      enddo

      close(50)

!  now, write the Hamiltonian matrix

      if (n_spin.eq.1) then
        write(file_name(1),'(A)') "hamiltonian.out"
      else if (n_spin.eq.2) then
        write(file_name(1),'(A)') "hamiltonian_up.out"
        write(file_name(2),'(A)') "hamiltonian_dn.out"
      end if

      do i_spin = 1, n_spin, 1

      open (50, file=file_name(i_spin))

      i_basis_2 = 1

      do while (i_basis_2.le.n_basis)

        write(info_str,'(8X,5(4X,I4,4X))') &
          (i_basis, i_basis = i_basis_2, i_basis_2+4, 1)
        call localorb_info(info_str,50,'(A)')

        do i_basis_1 = 1, n_basis, 1

            if (i_basis_2.le.(n_basis-4)) then

              do i_basis = i_basis_2, i_basis_2+4,1
                if (i_basis.ge.i_basis_1) then
                  i_index = i_basis_1 + (i_basis-1)*i_basis/2
                  output_element(i_basis) = hamiltonian(i_index,i_spin)
                else
                  output_element(i_basis) = 0.
                end if
              end do

              write(info_str,'(I5, 3X, 5(F12.6,2X))') i_basis_1, &
                (output_element(i_basis), &
                 i_basis = i_basis_2, i_basis_2+4, 1)
              call localorb_info(info_str,50,'(A)')

            else

              do i_basis = i_basis_2, n_basis,1
                if (i_basis.ge.i_basis_1) then
                  i_index = i_basis_1 + (i_basis-1)*i_basis/2
                  output_element(i_basis) = hamiltonian(i_index,i_spin)
                else
                  output_element(i_basis) = 0.
                end if
              end do

              write(info_str,'(I5, 3X, 5(F12.6,2X))') i_basis_1, &
                (output_element(i_basis), &
                 i_basis = i_basis_2, n_basis, 1)
              call localorb_info(info_str,50,'(A)')

            end if

        enddo

        i_basis_2 = i_basis_2+5
        write(info_str,*)
        call localorb_info(info_str,50,'(A)')

      enddo

      close(50)

      ! spin
      enddo

    end subroutine output_real_matrices
      
!----------------------------------------------------------------------
!******	
