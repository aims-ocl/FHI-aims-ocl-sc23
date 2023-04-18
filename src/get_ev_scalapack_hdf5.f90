!****s* FHI-aims/get_angular_grid
!  NAME
!    get_angular_grid
!  SYNOPSIS

subroutine get_ev_scalapack_hdf5 &
     (KS_eigenvalue, eigenvec_complex, ovlp_complex, l_row, l_col, mxld, &
       mxcol,my_k_point,chemical_potential)

!  PURPOSE
!
!  Wrapper function for outputting the overlap_matrix and the KS_eigenvectors
!  in case of a scalapack run using hdf5. The eigenvectors don't have to be 
!  collected.
!  All data to calculate the MODOS is written to the file KS_EV.h5
!
!  USES
  use calculate_mommat_base
  use dimensions
  use runtime_choices
  use localorb_io
  use hdf5_tools, only: HID_T, HSIZE_T, open_hdf5, open_hdf5_dataset, &
      out_dimensions, out_basis_atom, out_KS_eigenvalue, out_KS_vec_scalapack, &
      out_ovlp_scalapack, close_hdf5
  use mpi_utilities
  
  implicit none 

!  ARGUMENTS
  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) :: KS_eigenvalue
  complex*16, dimension(mxld, mxcol, n_spin), INTENT(IN):: eigenvec_complex
  complex*16, dimension(mxld, mxcol), INTENT(IN):: ovlp_complex
                                                               
  integer, dimension(n_basis), INTENT(IN) :: l_row
  integer, dimension(n_basis), INTENT(IN) :: l_col
  integer, INTENT(IN) :: mxld
  integer, INTENT(IN) :: mxcol 
  integer, INTENT(IN) :: my_k_point
  real*8, INTENT(IN) :: chemical_potential

!  INPUTS
!   o KS_eigenvalue -- KS-Eigenvalue
!   o eigenvec_complex -- KS_eigenvectors in scalapack arrays
!   o ovlp_complex -- (stored) scalapack overlap matrix
!   o l_row/l_col -- local row/column positions
!   o mxld/mxcol -- size of scalapack arrays
!   o my_k_point -- local k_point
!   o chemical_potential -- \epsilon_F
!  OUTPUT
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

  INTEGER(HID_T) :: file_id                  ! File identifier
  INTEGER(HID_T) :: plist_id                 ! Property list identifier

  character*128 :: info_str
  integer :: info

  !  counters

  integer(HSIZE_T) :: dim_KS(5)
  integer(HSIZE_T) :: dim_ovl(4)
  integer :: i_col, i_row, i_spin
  !  begin work

    write(info_str,'(6X,A,1X,I4)') "Output of eigenvector - scalapack"
    call localorb_info ( info_str )
!    call get_state_minmax_k(KS_eigen, n_state_min_in, n_state_max_in)	
    dim_KS(1) =  n_basis
    dim_KS(2) =  n_states
    dim_KS(3) =  n_spin
    dim_KS(4) =  n_k_points
    dim_KS(5) =  2
    
    dim_ovl(1) = n_basis
    dim_ovl(2) = n_basis
    dim_ovl(3) = n_k_points
    dim_ovl(4) = 2

    call open_hdf5('KS_EV.h5', file_id, plist_id)
    call open_hdf5_dataset('KS_Eigenvector', file_id,5, dim_KS)
    call open_hdf5_dataset('overlap_matrix', file_id,4, dim_ovl)
    call out_dimensions(file_id,plist_id)    
    call out_basis_atom(file_id,plist_id)
    call out_KS_Eigenvalue((KS_eigenvalue- chemical_potential)*hartree,file_id,&
                            plist_id)
    call mpi_barrier(mpi_comm_global,info)
    do i_spin = 1, n_spin
	do i_col = 1, n_states
	  if(l_col(i_col)==0) cycle
	  do i_row = 1, n_basis
	      if(l_row(i_row)>0) then
		call out_KS_vec_scalapack(eigenvec_complex(l_row(i_row), &
                          l_col(i_col), i_spin), file_id, plist_id, &
			  i_row, i_col, i_spin, my_k_point, 'KS_Eigenvector')
	      end if
	  end do
	end do
    end do
    do i_col = 1, n_basis
      if(l_col(i_col)==0) cycle
	do i_row = 1, i_col
	  if(l_row(i_row)>0) then
               call out_ovlp_scalapack(ovlp_complex(l_row(i_row), &
                          l_col(i_col)), file_id, plist_id, &
			  i_row, i_col, my_k_point, 'overlap_matrix')
          endif
        enddo
    enddo
    call mpi_barrier(mpi_comm_global,info)
    call close_hdf5(file_id,plist_id)
    write(info_str,'(6X,A,1X,I4)') "Output finished"
 
end subroutine get_ev_scalapack_hdf5
!******	
