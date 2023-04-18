!****s* FHI-aims/auxiliary_matrix_multi
!  NAME
!   auxiliary_matrix_multi
!  SYNOPSIS

      subroutine auxiliary_matrix_multi &
           (n_dim, n_loc_dim, map_basis, a, b, c ) 

!  PURPOSE
!  evaluate the dipole polarisability from the dipole moments of
!  KS state pairs.
!
!
!  USES

      use dimensions, only: n_max_loc_prodbas
      use mpi_tasks
      use synchronize_mpi
      use prodbas

      implicit none

!  ARGUMENTS

      integer :: n_dim
      integer :: n_loc_dim
      integer :: map_basis(n_max_loc_prodbas,n_tasks)
      real*8  a(n_dim, n_loc_dim)
      real*8  b(n_dim, n_loc_dim)
      real*8  c(n_dim, n_loc_dim)

!  INPUTS
!  o n_dim -- the global dimension of the array 
!  o n_loc_dim -- the local dimension of the array 
!  o map_basis -- mapping the local index of the array to its global index
!  o a -- input matrix a, its columns are distibuted over processors
!         according to the mapping relation defined by map_basis
!  o b -- input matrix b, its columns are distibuted over processors
!         according to the mapping relation defined by map_basis
!
!  OUTPUTS
!  o c -- output matrix given by a*b, its columns are distibuted over processors
!         according to the mapping relation defined by map_basis
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

      integer :: i_index
      real*8, dimension(:,:), allocatable :: aux_loc_matr
      real*8, dimension(:,:), allocatable :: aux_matr

!     counters


      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_basis_3

      integer :: i_task

!     begin work

      allocate (aux_loc_matr(n_loc_dim,n_loc_dim),stat=i_index)
      call check_allocation(i_index, 'aux_loc_matr                  ')

      allocate (aux_matr(n_dim,n_dim),stat=i_index)
      call check_allocation(i_index, 'aux_matr                      ')

!      do i_task = 1, n_tasks, 1
! distribute matrix b into different processors, for parallel calc.
          aux_matr(:,:) = 0.d0
          i_task = myid+1
          do i_basis_1 = 1, n_loc_dim, 1
             i_basis_3 = map_basis(i_basis_1, i_task)
             if(i_basis_3.gt.0) then
                aux_matr(:,i_basis_3) = &
                  a(:,i_basis_1)
             endif
          enddo
          call sync_matrix(aux_matr,n_dim,n_dim)

          call dgemm('N', 'N', n_dim, n_loc_dim, n_dim, 1.d0, &
                 aux_matr, n_dim, b, n_dim, 0.d0, &
                 c, n_dim)

!          call sync_matrix(aux_matr,n_dim,n_loc_dim)
!
!          if(myid.eq.i_task-1) then 
!           c(:,:) = aux_matr(:,:)
!          endif

!   end of loop over i_tasks
!      enddo


      deallocate(aux_loc_matr)
      deallocate(aux_matr)
      return
      end subroutine auxiliary_matrix_multi
!---------------------------------------------------------------------
!**********
