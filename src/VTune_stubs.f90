!****h* FHI-aims/F12/VTune_stubs
!  NAME
!    VTune
!  SYNOPSIS

module VTune

!  PURPOSE
!    stub version of the VTune module, has no purpose besides providing empty
!    dummy functions
!  USES
   use, intrinsic :: iso_c_binding, only: &
      c_char, &
      c_ptr

   implicit none
!  AUTHOR
!    Arvid Ihrig, FHI-aims team, Fritz-Haber-Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180, 2175-2196 (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2016).
!  SOURCE
   public

   contains

   subroutine vtune_resume()
   end subroutine vtune_resume

   subroutine vtune_pause()
   end subroutine vtune_pause

   subroutine vtune_detach()
   end subroutine vtune_detach

   type(c_ptr) function vtune_create_domain(name)
      character(kind=c_char) :: name(*)
   end function vtune_create_domain

   type(c_ptr) function vtune_create_string_handle(name)
      character(kind=c_char) :: name(*)
   end function vtune_create_string_handle

   subroutine vtune_frame_start(domain)
      type(c_ptr), value :: domain
   end subroutine vtune_frame_start

   subroutine vtune_frame_end(domain)
      type(c_ptr), value :: domain
   end subroutine vtune_frame_end

   subroutine vtune_task_start(domain, name)
      type(c_ptr), value :: domain
      type(c_ptr), value :: name
   end subroutine vtune_task_start

   subroutine vtune_task_end(domain)
      type(c_ptr), value :: domain
   end subroutine vtune_task_end

end module VTune
