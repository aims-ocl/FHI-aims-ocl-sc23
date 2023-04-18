!****s* FHI-aims/multi_ovlp3KS_sqrtw
!  NAME
!   multi_ovlp3KS_sqrtw
!  SYNOPSIS
      subroutine multi_ovlp3KS_sqrtw &
                 ( screened_coulomb, &
                   ovlpw &
                  )
!  PURPOSE
!  The routine multiply the ovlp_3KS matrix with the
!  inverse square root of the screened coulomb matirx
!
!  USES

      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi

      implicit none

! ARGUMENTS
      real*8, dimension(n_basbas, n_loc_prodbas) :: screened_coulomb
      real*8, dimension(n_loc_prodbas, n_states) :: ovlpw

! INPUTS
!  o screened_coulomb -- square root of the screened Coulomb matrix (more precisely
!            the square root of dielectric function)
!  o ovlpw -- the ovlp_3KS matrix for a given i_state
!
! OUTPUTS
! o ovlpw -- the input ovlp_3KS multiplied with screened_coulomb
!
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
      real*8, allocatable, dimension(:,:) :: temp_ovlp_matr
      integer :: info

!  counter
      integer :: i_index
      integer :: i_prodbas_1

      character(*), parameter :: func = 'multi_ovlp3KS_sqrtw.f90'

!  start to work


       allocate (temp_ovlp_matr(n_basbas, n_states),stat=info)
       call check_allocation(info,'temp_ovlp_matrx',func)
      
       temp_ovlp_matr(:,:)=0.d0
       do i_prodbas_1 = 1, n_loc_prodbas, 1
          i_index = map_prodbas(i_prodbas_1,myid+1)
          temp_ovlp_matr(i_index,:) = ovlpw(i_prodbas_1,:)
       enddo

       call sync_matrix(temp_ovlp_matr,n_basbas,n_states)

       call dgemm('T', 'N', n_loc_prodbas, n_states, &
                     n_basbas, 1.0d0, &
                     screened_coulomb, n_basbas, &
                     temp_ovlp_matr, n_basbas, 0.d0,&
                     ovlpw, n_loc_prodbas &
                  )


      deallocate (temp_ovlp_matr)

      end subroutine multi_ovlp3KS_sqrtw
!******
