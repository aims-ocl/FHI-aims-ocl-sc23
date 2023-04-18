!****h* FHI-aims/adjust_scf_settings
!  NAME
!    adjust_scf_settings
!  SYNOPSIS
!
 subroutine adjust_scf_settings
!  PURPOSE
!    Subroutine will adjust s.c.f. settings on the fly based on simple estimated insights
!    about the system under study. The following parameters will be automatically adjusted
!    if they are not set explicitly in control.in:
!    - charge_mix_param
!    - occupation_width
!    Currently, the intention is to do this only once in a given run, namely in the first
!    s.c.f. iteration and based on the band gap of the system (gap or no gap) - systems without
!    a gap will not, as a tendency, need to be so careful about their mixer settings.
!    However, we can see how this experience evolves and if the adjustment should be 
!    considered in other s.c.f. steps, too.
! USES
   use localorb_io
   use runtime_choices
   use physics, only : estimate_low_gap
   use constants, only : hartree
   use elsi_wrapper, only : eh_scf,aims_elsi_set_mu_broaden_width

   implicit none

!  ARGUMENTS

!  INPUT

!  OUTPUT

!  local variables

!  AUTHOR
!    Volker Blum (2017)
!  HISTORY
!    Release version, FHI-aims (2017).
!  COPYRIGHT
!   Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement.
!  SOURCE

! local variables
    
   character*200 :: info_str

! begin work

    write(info_str,'(2X,A)') &
    " "
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A)') &
    "Checking to see if s.c.f. parameters should be adjusted."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

   ! Adjust default settings for the s.c.f. cycle based on whether or not the system is
   ! gapless.

   ! We can always change the mixing parameter (although it is not a priori clear if
   ! this will be helpful). We only do this for the Pulay mixer (default, used almost always)
   if ((.not.charge_mix_param_set) .and. (mixer.eq.MIX_PULAY)) then
      if (estimate_low_gap) then
         charge_mix_param(:) = 0.02
         prec_mix_param(1) = charge_mix_param(1)
         write(info_str,'(2X,A,A,F10.6,A)') &
         "Adjusted the Pulay mixing parameter (charge_mix_param) for an expected", &
         " low-gap system. New value: ", charge_mix_param(1), " ."
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
      else
         ! 
         charge_mix_param(:) = 0.2
         prec_mix_param(1) = charge_mix_param(1)
         write(info_str,'(2X,A,F10.6,A)') &
         "The system likely has a gap. Increased the default Pulay mixing parameter (charge_mix_param). Value: ", &
         charge_mix_param(1), " ."
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
      end if
   end if

   ! The occupation broadening is trickier. If no default was set in control.in, then we can adjust
   ! its value once, based on whether or not we think we are dealing with a metallic / open-shell 
   ! system. However, we cannot adjust it again during a relaxation run since the value of the 
   ! occupation broadening enters the free energy expression, which must remain consistent
   ! along the trajectory.
   if (.not.occupation_width_set) then
      if (estimate_low_gap) then
         occupation_width = 0.05/hartree
         write(info_str,'(2X,A,A,F15.8,A)') &
         "Adjusted the occupation broadening width for an expected", &
         " low-gap system. New value: ", occupation_width*hartree, " eV."
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
         ! Must also inform ELSI's scf cycle settings.
         call aims_elsi_set_mu_broaden_width(eh_scf,occupation_width)
      else
         write(info_str,'(2X,A,F10.6,A)') &
         "Kept the default occupation width. Value: ", &
         occupation_width*hartree, " eV."
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
      end if
      ! Do not change the occupation broadening ever again.
      occupation_width_set = .true.
   end if

 end subroutine adjust_scf_settings

