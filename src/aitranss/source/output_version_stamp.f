      subroutine output_version_stamp
       implicit none

       write(6,fmt='(/,1x,a)') 'invoking aitranss (version: 071813) ... '
       write(6,fmt='(1x,a)') 'compiled on 2013/07/30 at 10:43:36 on host intcool'

      end subroutine output_version_stamp
