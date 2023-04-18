!****s* FHI-aims/prepare_partition_tabs
!  NAME
!   prepare_partition_tabs
!  SYNOPSIS

      subroutine prepare_partition_tabs ( )

!  PURPOSE
!   Creates the pieces needed to create the partition tables (partition_tab and
!   hartree_partition_tab) used for integrals and the Hartree potential. Usually,
!   the necessary pieces are the densities of free atoms - but not necessarily for
!   the same xc functional that is used in the actual electronic-structure calculation.
!
!  USES
!
      use localorb_io
      use dimensions
      use runtime_choices
      use free_atoms
      use species_data,    only : species_pseudoized, no_basis
      use mpi_tasks, only: aims_stop
      implicit none

!  ARGUMENTS


! INPUTS 
!   none
! OUTPUTS
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
!    Release version, FHI-aims (2009).
!  SOURCE
!

!     local variables

      logical :: fill_derivative

      character*120 :: info_str

      integer :: i_species
! begin work


  write (info_str,'(2X,A)') "Preparing densities etc. for the partition functions (integrals / Hartree potential)."
  call localorb_info(info_str,use_unit,'(A)')

  ! hartree_partition_type is always evaluated first, separate evaluation for partition_type only if different.
  if ( (flag_hartree_partition_type .le. 5) .or. (flag_hartree_partition_type .eq. 7) & 
       .or. (flag_hartree_partition_type .eq. 8) .or. (flag_hartree_partition_type .eq. 9) ) then
    ! old-style partition tables : Use free-atom density for whatever xc functional
    ! is specified in control.in 
    !
    ! For some partition_types, partition_rho_spl is not even used, we fill it here only
    ! to avoid it being uninitialized ...

    hartree_partition_rho_spl = free_rho_spl

    if (multipole_interpolation_style.eq.1) then
      ! not the default, and not tested for a while. 
      ! partition_drho_dr_spl array is currently not used by default.
      hartree_partition_drho_dr_spl = free_drho_dr_spl

    end if

  else if (flag_hartree_partition_type.eq.6) then
    ! In this case, partitioning is done based on the LDA free-atom density. This is a purely
    ! technical use of the free-atom density, and should have no bearing whatsoever on any
    ! physical results - except for making the integrations in the code more resilient
    ! against unexpected corner cases.

    if (multipole_interpolation_style.eq.1) then
      fill_derivative = .true.
    else
      fill_derivative = .false.
    end if
    call get_lda_free_atom_densities &
    ( hartree_partition_rho_spl, hartree_partition_drho_dr_spl, fill_derivative )

  else
    ! Error trap - we should never reach this point unless someone implemented
    ! a new partition type and forgot to tell us ...

    write (info_str,'(1X,A)') "* Subroutine prepare_partition_tabs.f90 : Attention."
    call localorb_info(info_str,use_unit,'(A)')

    write (info_str,'(1X,A,I5,A)') "* Found hartree_partition_type ", flag_hartree_partition_type, & 
                                   ", which was unknown at the time of writing of this routine."
    call localorb_info(info_str,use_unit,'(A)')

    write (info_str,'(1X,A)') "* If a recent modification was made to the partition_types, please also "
    call localorb_info(info_str,use_unit,'(A)')

    write (info_str,'(1X,A)') "* update subroutine prepare_partition_tabs.f90 ."
    call localorb_info(info_str,use_unit,'(A)')

    write (info_str,'(1X,A)') "* Stopping the code for now."
    call localorb_info(info_str,use_unit,'(A)')

    call aims_stop("Unknown partition type","prepare_partition_tabs.f90") 

  end if

  if ( ( flag_hartree_partition_type .ne. partition_type ) .or. use_prodbas ) then
    ! same game for partition tab

    if ( (partition_type .le. 5) .or. (partition_type .eq. 7) .or. (partition_type .eq. 8) .or. (partition_type .eq. 9) ) then
      ! old-style partition tables : Use free-atom density for whatever xc functional
      ! is specified in control.in 
      !
      ! For some partition_types, partition_rho_spl is not even used, we fill it here only
      ! to avoid it being uninitialized ...

      partition_rho_spl = free_rho_spl

      ! not (yet) needed for integral partition tab, left here only as a reminder
      !
      ! if (multipole_interpolation_style.eq.1) then
      !  ! this criterion is in principle overkill, as the 
      !  ! partition_drho_dr_spl array is currently not used by default.
      !  partition_drho_dr_spl = free_drho_dr_spl
      ! end if

    else if (partition_type.eq.6) then
      ! In this case, partitioning is done based on the LDA free-atom density. This is a purely
      ! technical use of the free-atom density, and should have no bearing whatsoever on any
      ! physical results - except for making the integrations in the code more resilient
      ! against unexpected corner cases.

      ! derivative will not be touched by this call
      fill_derivative = .false.
      call get_lda_free_atom_densities &
      ( partition_rho_spl, hartree_partition_drho_dr_spl, fill_derivative )

    else
      ! Error trap - we should never reach this point unless someone implemented
      ! a new partition type and forgot to tell us ...

      write (info_str,'(1X,A)') "* Subroutine prepare_partition_tabs.f90 : Attention."
      call localorb_info(info_str,use_unit,'(A)')

      write (info_str,'(1X,A,I5,A)') "* Found partition_type ", partition_type, &
                                     ", which was unknown at the time of writing of this routine."
      call localorb_info(info_str,use_unit,'(A)')

      write (info_str,'(1X,A)') "* If a recent modification was made to the partition_types, please also "
      call localorb_info(info_str,use_unit,'(A)')

      write (info_str,'(1X,A)') "* update subroutine prepare_partition_tabs.f90 ."
      call localorb_info(info_str,use_unit,'(A)')

      write (info_str,'(1X,A)') "* Stopping the code for now."
      call localorb_info(info_str,use_unit,'(A)')

      call aims_stop("Unknown partition type","prepare_partition_tabs.f90") 

    end if

    ! end hartree partition tab special case
  end if



! FIXME: D.B. 22.11.
if (use_embedding_pp) then
!! from here on we need to set some free atom properties to zero if (species_pseudoized)
  do i_species =1, n_species

    if (species_pseudoized(i_species).or.no_basis(i_species)) then

        free_rho_spl(:,:,i_species) = 0.d0
        if (allocated(free_drho_dr_spl)) then
           free_drho_dr_spl(:,:,i_species) = 0.d0
        end if
    end if
  end do
end if

end subroutine prepare_partition_tabs
