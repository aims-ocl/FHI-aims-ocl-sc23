!------------------------------------------------------------------------------
!****s* restart/write_rrs_pbc_restart_info
!  NAME
!    write_rrs_pbc_restart_info
!  SYNOPSIS
  subroutine write_rrs_pbc_restart_info
!  PURPOSE
!    write restart information to file
!  USES
    use runtime_choices
    use mpi_tasks
    use mpi_utilities
    use physics

    implicit none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    integer :: i_basis, i_states, i_spin, i_k
    character*150 :: info_str


    if(use_scalapack .and. use_density_matrix)then
       return
    end if


    write(info_str,'(2X,2A)') 'Writing periodic restart information to file ', rrs_pbc_restart_write_file
    call localorb_info(info_str,use_unit,'(A)')
    open(file = rrs_pbc_restart_write_file, unit = 7, status = 'unknown', form = 'unformatted')
    
    ! checkpoint information 
    write(7) n_basis*n_states*n_spin*n_k_points_task
    
    do i_k = 1, n_k_points_task
       do i_basis = 1, n_basis
          do i_states = 1, n_states
             write(7) (KS_eigenvector_complex(i_basis,i_states,i_spin,i_k) ,i_spin = 1, n_spin)
          end do
       end do
    end do
    
    ! also output: 
    ! KS_eigenvalue(n_states,n_spin,n_k_points)
    ! occ_numbers(n_states,n_spin,n_k_points)
    do i_k = 1, n_k_points 
       do i_states = 1, n_states
          do i_spin = 1, n_spin 
             write(7) KS_eigenvalue(i_states,i_spin,i_k), occ_numbers(i_states,i_spin,i_k)
          end do
       end do
    end do
    close(unit = 7)
  end subroutine write_rrs_pbc_restart_info
