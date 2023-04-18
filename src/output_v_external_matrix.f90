!****s* FHI-aims/output_v_external_matrix
!  NAME
!    output_v_external_matrix
!  SYNOPSIS

subroutine output_v_external_matrix &
     (  hamiltonian )

!  PURPOSE
!  Writes out the nuclear potential matrix
!  Proof of concept only! Works only for non-periodic case, non-packed matrices.
!  USES

  use dimensions
  use basis
  use localorb_io
  use runtime_choices, only : out_t_plus_v_matrix
  implicit none

!  ARGUMENTS
!  imported variables

  real*8 hamiltonian(n_basis*(n_basis+1)/2,n_spin)

!  INPUTS
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
!    Computer Physics Communications (2009), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2009).
!  SOURCE


!  local variables

      real*8 output_element(n_basis)

      character*40,dimension(n_spin) :: file_name

      character*100 :: info_str

!  counters

      integer i_basis, i_basis_1, i_basis_2, i_index, i_spin

!  begin work

      write (info_str,'(2X,A)') &
        "Writing electron-nuclear potential matrix."
      call localorb_info(info_str,use_unit,'(A)')

      write (info_str,'(2X,I5,A)') n_basis, " basis functions used."
      call localorb_info(info_str,use_unit,'(A)')

!  now, write the external potential matrix
      if (out_t_plus_v_matrix) then
        write(file_name(1),'(A)') trim('t_plus_v_matrix.out')
      else
        write(file_name(1),'(A)') trim('nuclear_potential_matrix.out')
      end if

      open (50, file=file_name(1))

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
                  output_element(i_basis) = hamiltonian(i_index,1)
                else
                  output_element(i_basis) = 0.
                end if
              end do

              write(info_str,'(I5, 3X, 5(F15.9,2X))') i_basis_1, &
                (output_element(i_basis), &
                 i_basis = i_basis_2, i_basis_2+4, 1)
              call localorb_info(info_str,50,'(A)')

            else

              do i_basis = i_basis_2, n_basis,1
                if (i_basis.ge.i_basis_1) then
                  i_index = i_basis_1 + (i_basis-1)*i_basis/2
                  output_element(i_basis) = hamiltonian(i_index,1)
                else
                  output_element(i_basis) = 0.
                end if
              end do

              write(info_str,'(I5, 3X, 5(F15.9,2X))') i_basis_1, &
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

    end subroutine output_v_external_matrix
      
!----------------------------------------------------------------------
!******	
