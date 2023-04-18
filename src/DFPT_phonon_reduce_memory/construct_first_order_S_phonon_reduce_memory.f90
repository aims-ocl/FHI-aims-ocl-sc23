!****s* FHI-aims/construct_first_order_S_phonon_reduce_memory
!  NAME
!   construct_first_order_S_phonon_reduce_memory
!  SYNOPSIS

subroutine construct_first_order_S_phonon_reduce_memory( & 
                                first_order_S_sparse, &
                                first_order_S_complex)


  !  PURPOSE
  !   Contruct first_order_S_complex on all k_point
  !   first_order_S_complex(i_basis,i_basis)=sum{R1,R2}[first_order_S(I_Cbasis,I_Cbasis)] 
  !   This is needed in periodic systems both for PM_none and PM_index
  !   shanghui 2012.10.06
  !   shanghui 2013.12.30 change to p0 version
  !   shanghui 2015.07.30 change to phonon_reduce_meory

  ! USES

  use pbc_lists
  use geometry
  use dimensions
  use mpi_tasks, only: myid, n_tasks
  implicit none

  !  ARGUMENTS

  complex*16 :: first_order_S_sparse(n_hamiltonian_matrix_size)
  complex*16 :: first_order_S_complex(n_basis, n_basis,n_k_points_task)

  !  INPUTS
  !    o first_order_S_sparse -- overlap in PM_index
  !  OUTPUT
  !    o first_order_S_complex -- overlap matrix in kpoint
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !
  ! SOURCE



  integer::     i_k_point, i_k_task


  i_k_task=0
  do i_k_point=1,n_k_points

   if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
      i_k_task = i_k_task + 1
   call construct_first_order_matrix_phonon_reduce_memory(first_order_S_sparse, &
       first_order_S_complex(1:n_basis,1:n_basis,i_k_task), &
       i_k_point)
   endif
 
 enddo  !i_k_point




end subroutine construct_first_order_S_phonon_reduce_memory
!******
