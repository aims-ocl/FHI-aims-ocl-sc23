!****h* FHI-aims/F12/VTune
!  NAME
!    VTune
!  SYNOPSIS

module VTune

!  PURPOSE
!    provides a fortran interface to the ITT API used by Intel VTune, useful
!    for more fine-grained profiling
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
   private

   public :: &
      vtune_resume, &
      vtune_pause, &
      vtune_detach, &
      vtune_create_domain, &
      vtune_create_string_handle, &
      vtune_frame_start, &
      vtune_frame_end, &
      vtune_task_start, &
      vtune_task_end

   interface

      subroutine vtune_resume() &
         bind(c, name='fortran_itt_resume')
      end subroutine vtune_resume

      subroutine vtune_pause() &
         bind(c, name='fortran_itt_pause')
      end subroutine vtune_pause

      subroutine vtune_detach() &
         bind(c, name='fortran_itt_detach')
      end subroutine vtune_detach

      type(c_ptr) function vtune_create_domain(name) &
         bind(c, name='fortran_itt_domain_create')
         import c_char, c_ptr
         character(kind=c_char) :: name(*)
      end function vtune_create_domain

      type(c_ptr) function vtune_create_string_handle(name) &
         bind(c, name='fortran_itt_string_handle_create')
         import c_char, c_ptr
         character(kind=c_char) :: name(*)
      end function vtune_create_string_handle

      subroutine vtune_frame_start(domain) &
         bind(c, name='fortran_itt_frame_begin')
         import c_ptr
         type(c_ptr), value :: domain
      end subroutine vtune_frame_start

      subroutine vtune_frame_end(domain) &
         bind(c, name='fortran_itt_frame_end')
         import c_ptr
         type(c_ptr), value :: domain
      end subroutine vtune_frame_end

      subroutine vtune_task_start(domain, name) &
         bind(c, name='fortran_itt_task_begin')
         import c_ptr
         type(c_ptr), value :: domain
         type(c_ptr), value :: name
      end subroutine vtune_task_start

      subroutine vtune_task_end(domain) &
         bind(c, name='fortran_itt_task_end')
         import c_ptr
         type(c_ptr), value :: domain
      end subroutine vtune_task_end

   end interface

end module VTune
