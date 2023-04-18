!****s* FHI-aims/get_total_forces
!  NAME
!    get_total_forces
!  SYNOPSIS

subroutine get_total_forces(ionic_forces, pulay_forces, hellman_feynman_forces, multipole_forces, &
   gga_forces, gga_forces_on, total_forces, vdw_forces, d_exx_ene, EXX_forces_on, pseudocore_forces, &
   nlcc_forces, nlcc_forces_on,Gnonmf_forces,Gnonmf_forces_on )

!  PURPOSE
!    subroutine get_total_forces to sum up all the force contributions
!    that are:
!    o (1) hellman-feynman-forces (-> sum_up_whole_potential_p1.f90)
!    o (2) pulay-forces (-> update_density_and_forces_orbital, update_forces_densmat)
!    o (3) ionic-forces (-> get_free_superpos_energy_p1 in initialize_scf.f90)
!    o (4) multipole-forces (-> sum_up_whole_potential_p1.f90)
!    o (5) gga-forces (-> update_density_and_forces_orbital, update_forces_densmat)
!    o (6) meta-gga-forces (-> update_density_and_forces_orbital, update_forces_densmat)
!    o (7) external-forces as given in geometry.in file (stored in geometry)
!    o (8) EXX-forces (???)
!    o (9) van-der-Wall forces (various implementations throughout aims...)
!
!  USES

  use runtime_choices, only: force_new_functional, hybrid_coeff, &
                             use_symmetry_reduced_spg
  use dimensions
  use species_data
  use geometry
  use synchronize_mpi
  use constants
  use localorb_io
  use sym_base,only: symmetrize_forces_spg, map_atoms_all, destroy_symmetry_arrays
  use spglib_symmetry,only: write_symm_info, out_symm_mats, destroy_symmetry,&
                            destroy_symmats
!, destroy_sym_maps
  use lpb_solver_utilities, only: use_Gnonmf_forces
  implicit none

!  ARGUMENTS

  real*8, dimension(3, n_atoms), intent(in) :: ionic_forces
  real*8, dimension(3, n_atoms), intent(in) :: pulay_forces
  real*8, dimension(3, n_atoms), intent(in) :: hellman_feynman_forces
  real*8, dimension(3, n_atoms), intent(in) :: multipole_forces 
  real*8, dimension(3, n_atoms), intent(in) :: gga_forces 
  real*8, dimension(3, n_atoms), intent(in) :: vdw_forces
  real*8, dimension(3, n_atoms), intent(in) :: d_exx_ene
  real*8, dimension(3, n_atoms), intent(in) :: pseudocore_forces
  real*8, dimension(3, n_atoms), intent(in) :: nlcc_forces
  real*8, dimension(3, n_atoms), intent(in) :: Gnonmf_forces 
  
  logical :: gga_forces_on
  logical :: EXX_forces_on 
  logical :: nlcc_forces_on
  logical :: Gnonmf_forces_on
  
  real*8, dimension(3, n_atoms), intent(out) :: total_forces

!  INPUTS
!  o ionic_forces -- forces between atom nucleus
!  o pulay_forces -- Pulay forces
!  o hellman_feynman_forces -- Hellman-Feyman forces
!  o multipole_forces -- multipole forces
!  o gga_forces -- gga forces
!  o gga_forces_on -- are the gga forces calculated or not ?
!
!  OUTPUT
!   o total_forces -- total forces
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
!    Release version, FHI-aims (2008).
!  SOURCE


	

 ! local variables
 real*8 :: conversion
 real*8 :: pulay_force_sum(3)

 ! counter
 integer :: i_atom
 integer :: i_atom_2
 integer :: i_coords
 character*120 :: info_str
 character*31 :: pulay_string

!ctest
!write(info_str*) "hellman_feynman_forces: ", myid
!write(info_str*) hellman_feynman_forces
!ctest end

!call sync_forces(pulay_forces, hellman_feynman_forces, multipole_forces, gga_forces, gga_forces_on)

!ctest
!write(info_str*) "hellman_feynman_forces after sync: ", myid
!write(info_str*) hellman_feynman_forces
!write(info_str*) "ionic_forces after sync: ", myid
!write(info_str*) ionic_forces
!ctest end

 ! synchronize all forces here
 conversion = hartree / bohr
 write(info_str,'(A)') ''
 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
 write (info_str,'(2X, A)') "atomic forces [eV/Ang]:"
 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
 write (info_str,'(2X, A)') "-----------------------"
 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
 
 do i_atom = 1, n_atoms, 1
    write (info_str,'(2X, A, I4)') "atom # ", i_atom
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  ) 

    do i_coords = 1, 3, 1
       ! add force contributions
       total_forces(i_coords, i_atom) = ionic_forces(i_coords, i_atom) + hellman_feynman_forces(i_coords, i_atom) + &
            pulay_forces(i_coords, i_atom) + multipole_forces(i_coords, i_atom) + pseudocore_forces(i_coords, i_atom)
       if (gga_forces_on) then
          total_forces(i_coords, i_atom) = total_forces(i_coords, i_atom) + gga_forces(i_coords, i_atom)
       end if
       if (use_vdw_correction_hirshfeld .or. use_mbd_std .or. &
               use_mbd_dev .or. use_mbd_old .or. use_libmbd) then
           total_forces(i_coords, i_atom) = total_forces(i_coords, i_atom) + vdw_forces(i_coords, i_atom)
       end if
       if (EXX_forces_on) then
          total_forces(i_coords, i_atom) = total_forces(i_coords, i_atom) + hybrid_coeff*d_exx_ene(i_coords, i_atom)
       end if
       if(nlcc_forces_on) then
           total_forces(i_coords, i_atom) = total_forces(i_coords, i_atom) + nlcc_forces(i_coords, i_atom)
       end if
       if(external_forces_on) then
           total_forces(i_coords, i_atom) = total_forces(i_coords, i_atom) + external_forces(i_coords, i_atom)
       end if
       if (Gnonmf_forces_on) then
           total_forces(i_coords, i_atom) = total_forces(i_coords, i_atom) + Gnonmf_forces(i_coords, i_atom)
       end if
    end do

    if(n_periodic == 0 .and. (.not. force_new_functional)) then
!       write (info_str,'(3X, A, 3(1X, E13.6))') "ion-ion         :", (ionic_forces(i_coords, i_atom) * conversion, i_coords = 1, 3, 1)
       write (info_str,'(3X, A, 3(1X, E13.6))') "Hellman-Feynman :", & 
            ((hellman_feynman_forces(i_coords, i_atom) + ionic_forces(i_coords, i_atom)) &
            * conversion, i_coords = 1, 3, 1)
       call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       write (info_str,'(3X, A, 3(1X, E13.6))') "Pulay           :", & 
            (pulay_forces(i_coords, i_atom) * conversion, i_coords = 1, 3, 1)
       call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       write (info_str,'(3X, A, 3(1X, E13.6))') "Multipole       :",& 
            (multipole_forces(i_coords, i_atom)*conversion, i_coords = 1, 3, 1)
       call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       if (gga_forces_on) then
          write (info_str,'(3X, A, 3(1X, E13.6))') "GGA             :", & 
                (gga_forces(i_coords, i_atom) * conversion, i_coords = 1, 3, 1)
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       else if (use_gga) then
          write (info_str,'(3X, A)') "GGA             :    --- not yet available --- "
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       end if
       if (Gnonmf_forces_on) then
          write (info_str,'(3X, A, 3(1X, E13.6))') "MPB, Gnonmf     :", & 
                (Gnonmf_forces(i_coords, i_atom) * conversion, i_coords = 1, 3, 1)
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )       
       end if
       if(external_forces_on) then
          write (info_str,'(3X, A, 3(1X, E13.6))') "External force  :", & 
                (external_forces(i_coords, i_atom) * conversion, i_coords = 1, 3, 1)
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       end if
       write (info_str,'(3X, A)') "--------------------------------------------------"
       call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       if (use_gga.and.(.not.gga_forces_on)) then
          write (info_str,'(3X, A, I4, A)') "Total forces(" ,i_atom, "):    --- not yet available --- "
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       else
          write (info_str,'(3X, A, I4, A, 3(1X, E13.6))') "Total forces(" ,i_atom, "):", & 
               (total_forces(i_coords, i_atom) * conversion, i_coords = 1, 3, 1)
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       end if
    else !periodic case and lvl_fast cluster case
!       write (info_str,'(3X, A, 3(1X, E13.6))') "ion-ion         :", (ionic_forces(i_coords, i_atom) * conversion, i_coords = 1, 3, 1)
!       call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
!
! AJL/2014
! I have separated the output for all the force components.
! The overhead is tiny but the understanding is then greatly improved. And debugging is easier!
!
!       write (info_str,'(3X, A, 3(1X, E13.6))') "Hellmann-Feynman + Multipole  :", & 
       write (info_str,'(3X, A, 3(1X, E13.6))') "Hellmann-Feynman              :", &
            ((hellman_feynman_forces(i_coords, i_atom)) &
            * conversion, i_coords = 1, 3, 1)
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       write (info_str,'(3X, A, 3(1X, E13.6))') "Ionic forces                  :", & 
            ((ionic_forces(i_coords, i_atom)) &
            * conversion, i_coords = 1, 3, 1)
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       write (info_str,'(3X, A, 3(1X, E13.6))') "Multipole                     :", & 
            ((multipole_forces(i_coords, i_atom)) &
            * conversion, i_coords = 1, 3, 1)
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       if (EXX_forces_on) then
         write (info_str,'(3X, A, 3(1X, E13.6))') "Hartree-Fock exchange         :", &
             (hybrid_coeff * d_exx_ene(i_coords,i_atom) * conversion, i_coords = 1, 3, 1)
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       endif
       !PULAY FORCES
       do i_coords=1,3, 1
           pulay_force_sum(i_coords) = pulay_forces(i_coords, i_atom)
           if (use_gga) then
               pulay_force_sum(i_coords) = pulay_force_sum(i_coords)+gga_forces(i_coords,i_atom)
             if (use_Gnonmf_forces) then
               pulay_force_sum(i_coords) = pulay_force_sum(i_coords)+Gnonmf_forces(i_coords,i_atom)
               pulay_string="Pulay + GGA + MPB, Gnonmf     :"
             else
               pulay_string="Pulay + GGA                   :"
             end if
           else
             if (use_Gnonmf_forces) then
               pulay_force_sum(i_coords) = pulay_force_sum(i_coords)+Gnonmf_forces(i_coords,i_atom)
               pulay_string="Pulay + MPB, Gnonmf           :"
             else
               pulay_string="Pulay                         :"
             end if
           end if
       end do
       write (info_str,'(3X, A, 3(1X, E13.6))') pulay_string, &
           ( pulay_force_sum(i_coords)  * conversion, i_coords = 1, 3, 1)
       call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       !END PULAY
       if (use_vdw_correction_hirshfeld .or. use_mbd_std .or. use_mbd_dev .or. use_libmbd) then
         write (info_str,'(3X, A, 3(1X, E13.6))') "Van der Waals                 :", &
             ( vdw_forces(i_coords, i_atom)  * conversion, i_coords = 1, 3, 1)
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       end if
       if (use_mbd_old) then
         write (info_str,'(3X, A, 3(1X, E13.6))') "MBD@rsSCS                     :", &
             ( vdw_forces(i_coords, i_atom)  * conversion, i_coords = 1, 3, 1)
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       end if

       if(use_embedding_pp) then
         write (info_str,'(3X, A, 3(1X, E13.6))') "Pseudopotential               :", & 
            ((pseudocore_forces(i_coords, i_atom)) &
            * conversion, i_coords = 1, 3, 1)
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         if(nlcc_forces_on) then
            write (info_str,'(3X, A, 3(1X, E13.6))') "nlcc forces                   :", & 
               ((nlcc_forces(i_coords, i_atom)) &
               * conversion, i_coords = 1, 3, 1)
            call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         endif
       endif
       if(external_forces_on) then
          write (info_str,'(3X, A, 3(1X, E13.6))') "External Force                :", & 
                (external_forces(i_coords, i_atom) * conversion, i_coords = 1, 3, 1)
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       end if
       write (info_str,'(3X, A)') "----------------------------------------------------------------"
       call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
       write (info_str,'(3X, A, I4, A, 3(1X, E13.6))') "Total forces(" ,i_atom, ")            :", & 
            (total_forces(i_coords, i_atom) * conversion, i_coords = 1, 3, 1)
       call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

!       total_forces(1,:) =  total_forces(1,:) - sum(total_forces(1,:))
!       total_forces(2,:) =  total_forces(2,:) - sum(total_forces(2,:))
!       total_forces(3,:) =  total_forces(3,:) - sum(total_forces(3,:))
    end if ! periodic
 end do ! i_atom

 if (use_symmetry_reduced_spg) call symmetrize_forces_spg(total_forces)
 if (use_spglib .and. use_symmetric_forces .and..not.use_symmetry_reduced_spg)then
    call destroy_symmetry()
    call destroy_symmats()
!    call destroy_sym_maps()
    call destroy_symmetry_arrays()
    call write_symm_info()
    call out_symm_mats()
    call map_atoms_all()
    call symmetrize_forces_spg(total_forces)
 endif

 write(info_str,'(A)') ''
 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

end subroutine get_total_forces
!******	
