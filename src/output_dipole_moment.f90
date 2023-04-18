   subroutine output_dipole_moment()

      use dimensions, only: use_pimd_wrapper
      use localorb_io, only: localorb_info, use_unit, comm_string
      use runtime_choices, only: ipi_dip
      use constants, only: bohr
      use physics, only: rho, partition_tab

      implicit none

!  local variables
      real*8 :: dipole_moment(3)
      real*8 :: abs_dipole_moment !absolute of the dipole moment
      real*8 :: dipole_ele_moment(3)
      real*8 :: dipole_ion_moment(3)
      real*8 :: ion_charge(3)
      real*8 :: ele_charge(3)
      real*8 :: total_charge(3)
      integer i_coord
      character*300 :: info_str
      character*300 :: tmp_str
      real*8, parameter :: Debye_To_eA = 0.20822678

!  begin work

      write(info_str, '(2X,A)') &
         "------------------------------------------------------------"
      call localorb_info(info_str,use_unit,'(A)')

      write(info_str, '(2X,A)') &
         "Computing monopole / dipole moments"
      call localorb_info(info_str,use_unit,'(A)')

      total_charge = 0.d0
      dipole_ele_moment = 0.d0
      dipole_ion_moment = 0.d0
      dipole_moment = 0.d0
      abs_dipole_moment = 0.d0

      ! calculate total charge
      call evaluate_moment_p1 &
              (0.d0,rho,partition_tab,ele_charge,ion_charge)

      total_charge = ion_charge - ele_charge

      ! calculate center of charge density
      call evaluate_moment_p1 &
              (1.d0,rho,partition_tab,dipole_ele_moment,dipole_ion_moment)

      dipole_moment = dipole_ion_moment - dipole_ele_moment

      abs_dipole_moment = sqrt(dipole_moment(1)**2+dipole_moment(2)**2+dipole_moment(3)**2)

      write(info_str,'(2X,A,X,E30.15)') &
         "| Total electronic charge [e]         :", ele_charge(1)
      call localorb_info(info_str,use_unit,'(A)')

      write(info_str,'(2X,A,X,E30.15)') &
         "| Total ionic charge [e]              :", ion_charge(1)
      call localorb_info(info_str,use_unit,'(A)')

      write(info_str,'(2X,A,X,E30.15)') &
         "| Total charge [e]                    :", total_charge(1)
      call localorb_info(info_str,use_unit,'(A)')

      if (abs(total_charge(1)) > 1.d-3) then
         write(info_str,'(1X,A,F15.8,A)') &
            "* Warning: The system has a total charge of ", total_charge(1), " ."
         call localorb_info(info_str,use_unit,'(A)')

         write(info_str,'(1X,A)') &
           "* The dipole moment depends on the choice of the origin."
         call localorb_info(info_str,use_unit,'(A)')
      end if

      write(info_str,'(2X,A,X,3E30.15)') &
         "| Total dipole moment [eAng]          :", &
         (dipole_moment(i_coord)*bohr, i_coord=1,3,1)
      call localorb_info(info_str,use_unit,'(A)')

      write(info_str,'(2X,A,X,E30.15,A,14X,A,10X,E30.15,A)') &
         "| Absolute dipole moment              :", &
         abs_dipole_moment*bohr, " eAng", "/", &
         abs_dipole_moment*bohr/Debye_To_eA, " Debye  ."
      call localorb_info(info_str,use_unit,'(A)')

      if (use_pimd_wrapper .and. ipi_dip) then
         write(tmp_str,'(2X,A,1X,3E30.15,A)') &
         "| Total dipole moment [eAng]          :", &
         (dipole_moment(i_coord)*bohr, i_coord=1,3,1)

         comm_string = trim(comm_string) // trim(tmp_str)
      end if

   end subroutine
