!----------------------------------------------------------------------
!  Subroutine output_periodic_dipole_moment
!
!----------------------------------------------------------------------
!This is a clone of output_dipole_moment, targeted towards periodic systems.
!Although the dipole moment in all directions is computed, only the z-component makes
!any sense, unless vacuum_levels in x / y are specified too (which the code does not support yet).
!Also, this function makes only sense for uncharged systems, so no monopoles.

      subroutine output_periodic_dipole_moment(dipole_moment)

      use dimensions
      use runtime_choices, only: ipi_dip
      use grids
      use geometry
      use spline
      use free_atoms
      use mpi_tasks
      use localorb_io
      use constants, only: bohr
      use physics

      implicit none

      !arguments
      real*8 :: dipole_moment(3)
     
      ! dip -- dipole moment of the cell in e*Bohr

!  local variables
      real*8 :: center_of_charge(3)
      real*8 :: dipole_ele_moment(3)
      real*8 :: dipole_ion_moment(3)
      real*8 :: ion_charge(3)
      real*8 :: ele_charge(3)
      real*8 :: total_charge(3)
      integer i_coord
      character(len=100) :: info_str

!  begin work

       if (myid .eq. 0) then
         write(info_str,'(2X,A)') &
         "------------------------------------------------------------"
         call localorb_info(info_str)
        write (info_str,'(2X,A)') &
        " Computing dipole for complete unit cell:"
         call localorb_info(info_str)
       end if

      dipole_ele_moment = 0d0
      dipole_ion_moment = 0d0
      dipole_moment = 0d0


!  calculate total charge
      call evaluate_moment_p2 &
       (0d0,rho,partition_tab,ele_charge,ion_charge)
      total_charge = ion_charge - ele_charge

!     safeguard: should actually trigger earlier, but better safe than sorry
      if(abs(total_charge(1)).gt.1d-6) then
         write(info_str,*) 'Total charge of system is ', total_charge(1)
         call localorb_info(info_str)
         call localorb_info('Total charge is not zero.') 
         !call aims_stop &
         !('The dipole is not well defined in this case. Stopping.')
      endif


!  calculate center of charge density 
      call evaluate_moment_p2 &
       (1d0,rho,partition_tab,dipole_ele_moment,dipole_ion_moment)

      dipole_moment = dipole_ion_moment - dipole_ele_moment

       if (myid .eq. 0) then
        !write(use_unit,'(2X,A,1X,E30.15,A)') &
        !    "| Total electronic charge [e]         :", ele_charge(1)
        !write(use_unit,'(2X,A,1X,E30.15,A)') &
        !    "| Total ionic charge [e]              :", ion_charge(1)
        !write(use_unit,'(2X,A,1X,E30.15,A)') &
        !    "| Total charge [e]                    :", total_charge(1)
        !write(use_unit,'(2X,A,1X,3E30.15,A)') &
        !   "| Pure electronic dipole moment [eAng]:", &
        !    (dipole_ele_moment(i_coord)*bohr, i_coord=1,3,1)
        !write(use_unit,'(2X,A,1X,3E30.15,A)') &
        !   "| Pure ionic dipole moment [eAng]     :", &
        !   (dipole_ion_moment(i_coord)*bohr, i_coord=1,3,1)
        write(use_unit,'(2X,A,1X,1E30.15,A)') &
         "| Total dipole moment in z-direction [eAng]          :", &
         (dipole_moment(3)*bohr)
      end if
      if (use_pimd_wrapper .and. ipi_dip) then
         write (comm_string,'(2X,A,1X,3E30.15,A)') &
         "| Total dipole moment in z-direction [eAng]          :", &
         (dipole_moment(i_coord)*bohr, i_coord=1,3,1)
      endif
      return
      end subroutine output_periodic_dipole_moment


