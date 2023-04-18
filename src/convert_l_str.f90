!---------------------------------------------------------------------
!  convert a string s,p,d,f,g, ... into the appropriate
!  angular momentum quantum number, and vice versa
!
!  provides
!  * subroutine convert_l_str: letter -> number
!  * function l_to_str : number -> letter
!  * function str_to_l : letter -> number
!
!---------------------------------------------------------------------

      subroutine convert_l_str(l_shell_str, l_shell)

      use localorb_io, only: use_unit
      implicit none

      character l_shell_str
      integer l_shell

      if (l_shell_str.eq."s") then

        l_shell = 0

      else if (l_shell_str.eq."p") then

        l_shell = 1

      else if (l_shell_str.eq."d") then

        l_shell = 2

      else if (l_shell_str.eq."f") then

        l_shell = 3

      else

        write(use_unit,*) "Unknown angular momentum state ", l_shell_str, "."
        stop

      end if

      return
      end

!---------------------------------------------------------------------

      integer function str_to_l(l_shell_str)

      use localorb_io, only: use_unit
      implicit none

      character l_shell_str

      if (l_shell_str.eq."s") then

        str_to_l = 0

      else if (l_shell_str.eq."p") then

        str_to_l = 1

      else if (l_shell_str.eq."d") then

        str_to_l = 2

      else if (l_shell_str.eq."f") then

        str_to_l = 3

      else if (l_shell_str.eq."g") then

        str_to_l = 4

      else if (l_shell_str.eq."h") then

        str_to_l = 5

      else if (l_shell_str.eq."i") then

        str_to_l = 6

      else

        write(use_unit,*) "Unknown angular momentum state ", l_shell_str, "."
        stop

      end if

      return
      end

!---------------------------------------------------------------------

      character function l_to_str(l_shell)

      use localorb_io, only: use_unit
      implicit none

      integer l_shell

      if (l_shell.eq.0) then

        l_to_str = 's'

      else if (l_shell.eq.1) then

        l_to_str = 'p'

      else if (l_shell.eq.2) then

        l_to_str = 'd'

      else if (l_shell.eq.3) then

        l_to_str = 'f'

      else if (l_shell.eq.4) then

        l_to_str = 'g'

      else if (l_shell.eq.5) then

        l_to_str = 'h'

      else if (l_shell.eq.6) then

        l_to_str = 'i'

      else

        write(use_unit,'(A,A,I4,A)') "* convert_l_str.f: No character ", &
          "to express angular momentum", l_shell, "."
        stop

      end if

      return
      end

