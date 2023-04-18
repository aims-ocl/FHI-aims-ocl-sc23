!****s* FHI-aims/construct_first_order_H_p0
!  NAME
!   construct_first_order_H_p0
!  SYNOPSIS

subroutine construct_first_order_H_p0(first_order_H, &
      first_order_H_complex, k_point)


  !  PURPOSE
  !   Contruct first_order_H_complex on the given k_point
  !   first_order_H_complex(i_basis,i_basis)=sum{R1,R2}[first_order_H(I_Cbasis,I_Cbasis)] 
  !   This is needed in periodic systems both for PM_none and PM_index
  !   shanghui 2012.10.06
  ! USES

  use pbc_lists
  use geometry
  use dimensions
  use mpi_tasks, only: check_allocation
  implicit none

  !  ARGUMENTS

  real*8 :: first_order_H(3, n_atoms, n_Cbasis, n_Cbasis)
  complex*16:: first_order_H_complex(3,n_atoms, n_basis, n_basis)
  integer:: k_point

  !  INPUTS
  !    o first_order_H -- first_order_H in PM_none
  !    o k_point -- k-point wanted to be calculated 
  !  OUTPUT
  !    o first_order_H_complex -- first_order_H matrix in kpoint
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


  complex*16, dimension(:,:),allocatable :: work
  real*8, dimension(:,:),allocatable :: work_ma

  integer:: i_index, i_coord, i_atom, i_basis_2,  i_basis_1, i_basis


  first_order_H_complex = (0.0d0,0.0d0)

  allocate(work(n_centers_basis_I, n_basis),stat=i_index)
  call check_allocation(i_index, 'work                          ')
  allocate(work_ma(n_centers_basis_I, n_centers_basis_I),stat=i_index)
  call check_allocation(i_index, 'work_ma                          ')


  do i_coord = 1,3
  do i_atom = 1, n_atoms

  do i_basis_1 = 1,n_centers_basis_I,1 

     work_ma(1:n_centers_basis_I, i_basis_1) = first_order_H(i_coord, i_atom, 1:n_centers_basis_I, i_basis_1)

  end do


!------------sum_{R1,R2}-----------------------------------
   work=(0.0d0,0.0d0)
   work(1:n_centers_basis_I, 1:n_basis) = work_ma(1:n_centers_basis_I,1:n_basis)
      
  do i_basis_2 = n_basis+1,n_centers_basis_I,1 

     work(1:n_centers_basis_I, Cbasis_to_basis(i_basis_2)) =   &
     work(1:n_centers_basis_I, Cbasis_to_basis(i_basis_2)) &
            + work_ma(1:n_centers_basis_I, i_basis_2) &
            *  dconjg(k_phase(center_to_cell(Cbasis_to_center( i_basis_2 )),k_point))
          !my_(R2-R1)  * k_phase(center_to_cell(Cbasis_to_center( i_basis_2 )),k_point)
      
  end do

  do i_basis_1 = n_basis+1,n_centers_basis_I,1 

     work(Cbasis_to_basis(i_basis_1), 1:n_basis) =   &
     work(Cbasis_to_basis(i_basis_1), 1:n_basis) &
          + work(i_basis_1, 1:n_basis) &
          * k_phase(center_to_cell(Cbasis_to_center( i_basis_1 )),k_point)
         !my_(R2-R1) *  dconjg(k_phase(center_to_cell(Cbasis_to_center( i_basis_1 )),k_point))

  end do
!------------end sum_{R1,R2}-----------------------------------

  do i_basis_2 = 1,n_basis
     first_order_H_complex(i_coord, i_atom, 1:n_basis, i_basis_2) =  &
                      work( 1:n_basis, i_basis_2)
  enddo

!-------------sum_{R3}------------------------------------
!--------->shanghui move this part to evaluate_first_order_S.f90, 2012.10.16
!  do i_basis_2 = 1,n_basis
!
!     first_order_S_complex(i_coord, center_to_atom(i_center), 1:n_basis, i_basis_2) =  &
!     first_order_S_complex(i_coord, center_to_atom(i_center), 1:n_basis, i_basis_2)+   &
!                      work( 1:n_basis, i_basis_2)
!  end do
!------------end sum_{R3}---------------------------------


  enddo  !i_atom
  enddo  !i_coord


  deallocate(work)
  deallocate(work_ma)


end subroutine construct_first_order_H_p0
!******
