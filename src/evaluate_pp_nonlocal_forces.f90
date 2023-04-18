!****s* FHI-aims/evaluate_nonlocal_energy
!  NAME
!   evaluate_nonlocal_energy
!  SYNOPSIS

subroutine evaluate_pp_nonlocal_forces (KS_coefficents, occ_numbers)

!  PURPOSE
!     Subroutine evaluate_nonlocal_engery
!     calculates energy term contributed by a nonlocal embedding potential.
!
!  USES

  use dimensions
  use pseudodata
  use runtime_choices, only: real_eigenvectors
  use mpi_tasks, only: aims_stop
  use basis, only: basis_atom
  use synchronize_mpi
  use synchronize_mpi_basic

!  use constants

  implicit none

!  ARGUMENTS



!  INPUTS
!  o KS_coefficents -- coefficient matrix of KS_states
!  o occ_numbers -- occupation numbers of KS_states

!  
!  OUTPUT
!  o en_nonlocal -- value of energy contributed by nonlocal embedding potential
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


   real*8, dimension(n_states, n_spin, n_k_points)             :: occ_numbers
   real*8, dimension(n_basis, n_states,n_spin,n_k_points_task) :: KS_coefficents
 
  !     counters

   integer :: i_spin, i, k, l, i_pp_basis_fns, i_pp_basis, i_pp_atom, i_atom, i_coord
   integer :: i_basis_1, i_basis_2


  ! for parallalization on pp_basis
   integer , dimension (:,:), allocatable :: map_index_basis
   integer :: n_remain, n_loc_pp_basis, i_index, i_task, i_loc_pp_basis


  ! other local variables
   character(*), parameter :: func = 'evaluate_pp_nonlocal_forces'

  pp_nonlocal_forces = 0.d0
  pp_nonlocal_forces_2 = 0.d0

  if(real_eigenvectors)then


      n_remain = MOD(n_pp_basis, n_tasks)
      if (n_remain.eq.0) then
        n_loc_pp_basis = n_pp_basis / n_tasks
      else
        n_loc_pp_basis = n_pp_basis / n_tasks + 1
      endif
      allocate(map_index_basis(n_tasks, n_loc_pp_basis))


       map_index_basis = 0
       i_index = 0
       do i_task = 1, n_tasks

         if((n_remain.eq.0) .or. &
            i_task.le.n_remain ) then

            do i = 1, n_loc_pp_basis
             i_index = i_index + 1
             map_index_basis( i_task, i) = i_index
            enddo

         else

            do i = 1, n_loc_pp_basis -1
             i_index = i_index + 1
             map_index_basis(i_task, i) = i_index
            enddo

         endif
       enddo

       if (i_index .ne. n_pp_basis) then
         call aims_stop(" * Error in the distributing pseudopotantial basis over the cores. Stop! ", func)
         stop
       endif


!   if(myid==0) then

    do i_spin = 1, n_spin
    do i_loc_pp_basis = 1, n_loc_pp_basis, 1
       i_pp_basis = map_index_basis(myid+1, i_loc_pp_basis)
       if(i_pp_basis.gt.0) then

       i_pp_basis_fns = pp_basis_fn(i_pp_basis)
       if (pp_basisfn_l(i_pp_basis_fns).ne.pp_local_component(pp_basisfn_species(i_pp_basis_fns))) then
!           if(i_pp_atom .eq.pp_basis_atom(i_pp_basis)) then
       i_pp_atom = pp_basis_atom(i_pp_basis)


     do k = 1, n_max_basis_overlap(i_pp_basis)
           i_basis_1 = nonzero_overlap_entries(k,i_pp_basis)
           i_atom = basis_atom(i_basis_1)



        do l = 1, n_max_basis_overlap(i_pp_basis)

           i_basis_2 = nonzero_overlap_entries(l,i_pp_basis)



        do i = 1, n_states

              do i_coord = 1,3,1

!              pp_nonlocal_forces(i_coord,i_pp_atom) = pp_nonlocal_forces(i_coord,i_pp_atom) + & 
!                   occ_numbers(i, i_spin, 1) * KS_coefficents(i_basis_1,i, i_spin, 1)* &
!                   KS_coefficents(i_basis_2,i, i_spin, 1)*&
!                   (-d_basiswave_pp_overlap(i_coord,i_basis_1, i_pp_basis)*E_l_KB(i_pp_basis_fns)*&
!                   pp_basiswave_overlap(i_pp_basis,i_basis_2) + &
!                   basiswave_pp_overlap(i_basis_1, i_pp_basis)*E_l_KB(i_pp_basis_fns)*&
!                   d_pp_basiswave_overlap(i_coord,i_pp_basis,i_basis_2))

! we can make use of symmetry which makes calc of d_pp_basiswave_overlap redundent
              pp_nonlocal_forces(i_coord,i_pp_atom) = pp_nonlocal_forces(i_coord,i_pp_atom) + & 
                   occ_numbers(i, i_spin, 1) * KS_coefficents(i_basis_1,i, i_spin, 1)* &
                   KS_coefficents(i_basis_2,i, i_spin, 1)*&
                   (-d_basiswave_pp_overlap(i_coord,i_basis_1, i_pp_basis)*E_l_KB(i_pp_basis_fns)*&
                   basiswave_pp_overlap(i_basis_2,i_pp_basis) - &
                   basiswave_pp_overlap(i_basis_1, i_pp_basis)*E_l_KB(i_pp_basis_fns)*&
                   d_basiswave_pp_overlap(i_coord,i_basis_2,i_pp_basis))



!              pp_nonlocal_forces_2(i_coord,i_atom) = pp_nonlocal_forces_2(i_coord,i_atom) + & 
!                   occ_numbers(i, i_spin, 1) * KS_coefficents(i_basis_1,i, i_spin, 1)* &
!                   KS_coefficents(i_basis_2,i, i_spin, 1)*&
!                   (d_basiswave_pp_overlap(i_coord,i_basis_1, i_pp_basis)*E_l_KB(i_pp_basis_fns)*&
!                   pp_basiswave_overlap(i_pp_basis,i_basis_2) - &
!                   basiswave_pp_overlap(i_basis_1, i_pp_basis)*E_l_KB(i_pp_basis_fns)*&
!                   d_pp_basiswave_overlap(i_coord,i_pp_basis,i_basis_2))

! same here
              pp_nonlocal_forces_2(i_coord,i_atom) = pp_nonlocal_forces_2(i_coord,i_atom) + & 
                   occ_numbers(i, i_spin, 1) * KS_coefficents(i_basis_1,i, i_spin, 1)* &
                  KS_coefficents(i_basis_2,i, i_spin, 1)*&
                   (d_basiswave_pp_overlap(i_coord,i_basis_1, i_pp_basis)*E_l_KB(i_pp_basis_fns)*&
                   basiswave_pp_overlap(i_basis_2,i_pp_basis) + &
                   basiswave_pp_overlap(i_basis_1, i_pp_basis)*E_l_KB(i_pp_basis_fns)*&
                   d_basiswave_pp_overlap(i_coord,i_basis_2,i_pp_basis))



              enddo

        enddo
        enddo
     enddo
     endif
     endif

     enddo
     enddo

  else 
! FIXME: not yet implemented how to do this in the case of complex eigenvectors
       call aims_stop('This is not yet implemented for complex eigenvectors!', func)

  end if   

  deallocate(map_index_basis)

  call sync_vector(pp_nonlocal_forces,3*n_pp_atoms)
  call sync_vector(pp_nonlocal_forces_2,3*n_real_atoms)

  return

end subroutine evaluate_pp_nonlocal_forces

!---------------------------------------------------------------------
!******	
