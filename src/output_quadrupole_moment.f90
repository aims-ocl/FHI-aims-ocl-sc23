!----------------------------------------------------------------------
!  Subroutine output_quadrupole_moment
!
!----------------------------------------------------------------------
      subroutine output_quadrupole_moment()

      use dimensions
      use grids
      use geometry
      use spline
      use free_atoms
      use mpi_tasks
      use localorb_io
      use constants, only: bohr
      use physics

      implicit none

!  local variables
      real*8 :: center_of_charge(3)
      real*8, dimension(3,3) :: quadrupole_moment
      real*8, dimension(3,3) :: quadrupole_ele_moment
      real*8, dimension(3,3) :: quadrupole_ion_moment
      real*8, dimension(3)   :: quadrupole_eigenvalues
      real*8, dimension(3)   :: electron_quadrupole_eigenvalues
      real*8, dimension(3)   :: ion_quadrupole_eigenvalues
      real*8, dimension(20)  :: work
      integer i_coord, info


!  begin work
      quadrupole_ele_moment = 0.0d0
      quadrupole_ion_moment = 0.0d0
      quadrupole_moment = 0.0d0


!  calculate quadrupole moment
      call evaluate_quadrupole_moment &
      (1d0,rho,partition_tab,quadrupole_ele_moment &
            ,quadrupole_ion_moment)

      quadrupole_moment = quadrupole_ion_moment - quadrupole_ele_moment

      call DSYEV ('N','L',3,quadrupole_moment,3, &
          quadrupole_eigenvalues,work, 20,info)

      call DSYEV ('N','L',3,quadrupole_ele_moment,3, &
            electron_quadrupole_eigenvalues,work,20,info)

      call DSYEV ('N','L',3,quadrupole_ion_moment,3, &
            ion_quadrupole_eigenvalues,work,20,info)

       if (myid .eq. 0) then
        write(use_unit,'(A)') &
        "------------------------------------------------------------"

        write (use_unit,'(2X,A)') "Quadrupole information:"
!        write(use_unit,'(2X,A,1X,3E30.15,A)')
!     +     "| Pure electronic quadrupole moment [eAng]:",
!     +      (quadrupole_ele_moment(i_coord)*bohr, i_coord=1,3,1)
!        write(use_unit,'(2X,A,1X,3E30.15,A)')
!     +     "| Pure ionic quadrupole moment moment [eAng]     :",
!     +     (quadrupole_ion_moment(i_coord)*bohr, i_coord=1,3,1)
        write (use_unit,'(2X,A,1X,3E30.15,A)') &
         "| electronic quadrupole moment [eAng^2]     :", &
         (electron_quadrupole_eigenvalues(i_coord)*bohr*bohr, &
          i_coord=1,3,1)
        write (use_unit,'(2X,A,1X,3E30.15,A)') &
         "| ionic quadrupole moment [eAng^2]          :", &
         (ion_quadrupole_eigenvalues(i_coord)*bohr*bohr, i_coord=1,3,1)
        write (use_unit,'(2X,A,1X,3E30.15,A)') &
         "| Total quadrupole moment [eAng^2]          :", &
         (quadrupole_eigenvalues(i_coord)*bohr*bohr, i_coord=1,3,1)
      end if
      return
      end
