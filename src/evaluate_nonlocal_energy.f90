!****s* FHI-aims/evaluate_nonlocal_energy
!  NAME
!   evaluate_nonlocal_energy
!  SYNOPSIS

subroutine evaluate_nonlocal_energy (KS_coefficents, occ_numbers)

!  PURPOSE
!     Subroutine evaluate_nonlocal_engery
!     calculates energy term contributed by a nonlocal embedding potential.
!
!  USES

  use dimensions
  use pseudodata
  use runtime_choices, only: real_eigenvectors
  use mpi_tasks, only: aims_stop
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

   integer :: i_spin, i, k, l


  ! for parallalization over basis
   integer , dimension (:,:), allocatable :: map_index_basis
   integer :: n_remain, n_loc_basis, i_index, i_task, i_loc_basis,n_basis_mytask


  ! other local variables
   character(*), parameter :: func = 'evaluate_nonlocal_energy'

  en_nonlocal = 0.d0
 

  if(real_eigenvectors)then


      n_remain = MOD(n_basis, n_tasks)
      if (n_remain.eq.0) then
        n_loc_basis = n_basis / n_tasks
      else
        n_loc_basis = n_basis / n_tasks + 1
      endif
      allocate(map_index_basis(n_tasks, n_loc_basis))


       map_index_basis = 0
       i_index = 0
       n_basis_mytask = n_loc_basis
       do i_task = 1, n_tasks

         if((n_remain.eq.0) .or. &
            i_task.le.n_remain ) then

            do i = 1, n_loc_basis
             i_index = i_index + 1
             map_index_basis( i_task, i) = i_index
            enddo

         else

            do i = 1, n_loc_basis -1
             i_index = i_index + 1
             map_index_basis(i_task, i) = i_index
            enddo
            if ((i_task.eq.(myid+1))) then
                 n_basis_mytask = n_loc_basis-1
            endif

         endif
       enddo

       if (i_index .ne. n_basis) then
         call aims_stop(" * Error in the distributing basis over the cores. Stop! ", func)
         stop
       endif


    do i_spin = 1, n_spin
    do i = 1, n_states    
    do i_loc_basis = 1, n_basis_mytask, 1
       k = map_index_basis(myid+1, i_loc_basis)
!               if(k.gt.0) then

        
         do l = 1, n_basis
            en_nonlocal = en_nonlocal + & 
                        occ_numbers(i, i_spin, 1) * KS_coefficents(k,i, i_spin, 1)* &
                        KS_coefficents(l,i, i_spin, 1)*nonlocal_matrix( k, l)

        
        enddo
!      endif
     enddo
     enddo
     enddo
  else 
! FIXME: not yet implemented how to do this in the case of complex eigenvectors
       call aims_stop('This is not yet implemented for complex eigenvectors!', func)

  end if   

  deallocate(map_index_basis)

  call sync_real_number(en_nonlocal)

!   return

end subroutine evaluate_nonlocal_energy

!---------------------------------------------------------------------
!******	
