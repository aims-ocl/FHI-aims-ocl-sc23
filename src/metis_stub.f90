!  FHI-aims/metis_stub
!  NAME
!    metis_stub
!  PURPOSE
!    This file contains the stub routines for metis, TetGen and qhull wrapper calls
!    in FHI-aims. These routines are provided only for the purposes of name resolution in
!    the linking stage when the actual libraries are not present. Calling these
!    routines will result in termination of the program.
!  USES
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
!
subroutine metis_tetgen_wrapper(coords_x, coords_y, coords_z, grid_partition, &
     n_points_in_grid, n_grid_batches)

  use localorb_io

  implicit none

  integer :: n_points_in_grid
  real*8 :: coords_x(n_points_in_grid)
  real*8 :: coords_y(n_points_in_grid)
  real*8 :: coords_z(n_points_in_grid)
  integer :: grid_partition(n_points_in_grid)
  integer :: n_grid_batches
  
  character*100 :: info_str

  write (info_str,'(1X,A)') &
       "* Metis library called without proper linking of the library. Aborting."
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine metis_tetgen_wrapper

subroutine metis_qhull_wrapper(coords_x, coords_y, coords_z, grid_partition, &
     n_points_in_grid, n_grid_batches)

  use localorb_io

  implicit none

  integer :: n_points_in_grid
  real*8 :: coords_x(n_points_in_grid)
  real*8 :: coords_y(n_points_in_grid)
  real*8 :: coords_z(n_points_in_grid)
  integer :: grid_partition(n_points_in_grid)
  integer :: n_grid_batches
  
  character*100 :: info_str

  write (info_str,'(1X,A)') &
       "* Metis library called without proper linking of the library. Aborting."
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine metis_qhull_wrapper

subroutine metis_nearest_wrapper(xadj, adjncy, grid_partition, &
     n_points_in_grid, n_grid_batches)

  use localorb_io

  implicit none

  integer :: n_points_in_grid
  integer :: xadj(n_points_in_grid + 1)
  integer :: adjncy
  integer :: grid_partition(n_points_in_grid)
  integer :: n_grid_batches
  
  character*100 :: info_str

  write (info_str,'(1X,A)') &
       "* Metis library called without proper linking of the library. Aborting."
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine metis_nearest_wrapper

subroutine metis_qhull_batches_wrapper(coms_x, coms_y, coms_z, weights, batch_partition, &
     n_grid_batches, n_threads)

  use localorb_io

  implicit none

  integer :: n_grid_batches
  real*8 :: coms_x(n_grid_batches)
  real*8 :: coms_y(n_grid_batches)
  real*8 :: coms_z(n_grid_batches)
  integer :: batch_partition(n_grid_batches)
  integer :: weights(n_grid_batches)
  integer :: n_threads
  
  character*100 :: info_str

  write (info_str,'(1X,A)') &
       "* Metis library called without proper linking of the library. Aborting."
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine metis_qhull_batches_wrapper
