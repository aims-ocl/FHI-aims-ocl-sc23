! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.
!

!>
!! Returns details about ELSI's versioning.
!!
subroutine elsi_version_info(version,datestamp,commit,hostname,datetime)

   implicit none

   character(len=8), intent(out) :: version
   character(len=8), intent(out) :: datestamp
   character(len=8), intent(out) :: commit
   character(len=40), intent(out) :: hostname
   character(len=20), intent(out) :: datetime

   version = "FHI-aims"
   datestamp = "N/A"
   commit = "N/A"
   hostname = "N/A"
   datetime = "N/A"

end subroutine

!>
!! Returns details about ELSI's installation.
!!
subroutine elsi_solver_info(have_pexsi)

   implicit none

   logical, intent(out) :: have_pexsi

   have_pexsi = .false.

end subroutine
