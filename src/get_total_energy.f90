!****s* FHI-aims/get_total_energy
!  NAME
!   get_total_energy
!  SYNOPSIS

      subroutine get_total_energy &
          ( ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
            en_ion_embed, en_density_embed, en_vdw, en_pot_vdw, &
            en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, &
            en_elec_free, &
            en_elec_delta, en_fock, &
            hartree_energy_free, hartree_delta_energy, &
            hartree_multipole_correction, &
            entropy_correction, penalty_energy, total_energy &
          )

! PURPOSE
! Subroutine get_total_energy.f
!
! Evaluates the total energy by summing up individual pieces calculated
! in prior subroutines.
!
! To understand the meaning of the terms, please read Section 4.6 in 
!   Blum et al., Computer Physics Communications 180 (2009) 2175-2196
! very carefully! The expression evaluated here is not simply the
! textbook expression for finite (molecular) systems, which would separate
! electrons and nuclei completely. Instead, the potential energy terms
! always refer to combined potentials from electrons and nuclei. 
! It is not possible to compute a convergent total energy for large and/or
! periodic systems without this transformation.
! Again, please follow the Equations in Section 4.6 .
!
! pieces of total energy are calculated in:
! * subroutine get_hartree_potential_and_density
!   ( -> get_hartree_potential_and_density.f )
! * subroutine sum_up_whole_potential
!    ( -> sum_up_whole_potential.f )
! * subroutine get_ev_sum
!     ( -> get_ev_sum.f )
! * subroutine get_entropy_correction
!     ( -> get_entropy_correction.f)
! here, just sum up these pieces and the ion-ion-interaction
!
!  USES

      use types, only: dp
      use dimensions
      use runtime_choices
      use geometry
      use species_data
      use runtime_choices
      use mpi_tasks
      use constants
      use localorb_io
      use pseudodata, only: en_nonlocal, en_xc_nlcc
      use lpb_solver_utilities, only: surface_and_volume_calc, surface_mpb, volume_mpb,&
	alphaplusgamma_mpb,beta_mpb, solve_lpbe_only,en_eps_rho,en_alpha_rho, &
	use_mpbe_free_energy, mpbe_no_ion, freeen_LPBE_nonelstat_solv_energy,&
	en_pot_Gnonmf, en_Gnonmf
      use mpb_solver_utilities, only: freeen_MPBE_nonelstat_solv_energy
      use mpe_interface, only: mpe_read_energy_contributions
      use plus_u, only: plus_u_energy_correction_term, plus_u_energy_correction
      use mbd_std_wrapper, only: mbd_self_consistent, mbd_first_step, &
          mbd_scf_converged
      implicit none

!  ARGUMENTS

      real*8, intent(in) :: ev_sum
      real*8, intent(in) :: ev_sum_shifted
      real*8, intent(in) :: en_xc
      real*8, intent(in) :: en_pot_xc
      real*8, intent(in) :: en_ion_ion
      real*8, intent(in) :: hartree_energy_free
      real*8, intent(in) :: hartree_delta_energy
      real*8, intent(in) :: hartree_multipole_correction
      real*8, intent(in) :: entropy_correction
      real*8, intent(IN) :: penalty_energy

      real*8, intent(in) :: en_ion_embed
      real*8, intent(in) :: en_density_embed
      real*8, intent(in) :: en_vdw
      real*8, intent(in) :: en_pot_vdw
      real*8, intent(in) :: en_ll_vdw, en_ll_vdw_err
      real*8, intent(in) :: en_lda_c, en_pbe_c
      real*8, intent(in) :: en_elec_free
      real*8, intent(in) :: en_elec_delta
      real*8, intent(in) :: en_fock 

      real*8, intent(out) :: total_energy

! INPUTS
! o ev_sum -- sum of eigenvalues
! o ev_sum_shifted -- shifted sum of eigenvalues
! o en_xc -- exchange and correlation energy
! o en_pot_xc -- potential part of the exchange and correlation energy
! o en_ion_ion -- ion ion energy
! o en_ion_embed -- embedded ion energy
! o en_density_embed -- energy of embedded charge density
! o en_vdw -- energy of van der Waals 
! o en_pot_vdw -- potential part of the vdW energy
! o en_elec_free -- free atoms energy
!   en_elec_free is the term written in the last line of Eq. (64) of 
!   Blum et al., Computer Physics Communications 180 (2009) 2175-2196
!   when evaluated for the superposition density of spherical free atoms.
! o en_elec_delta  -- Hartree energy from electron density only
! o en_fock -- Hartree Fock energy
! o hartree_energy_free -- free atoms Hartree energy
!   hartree_energy_free is the quantity given in Eq. (61) of the 
!   FHI-aims CPC paper, Blum et al., Computer Physics Communications 180 (2009) 2175-2196
!   when evaluated for the superposition density of spherical free atoms.
! o hartree_delta_energy -- Hartree energy minus free atoms Hartree energy
! o hartree_multipole_correction -- multipole correction to Hartree energy
! o entropy_correction -- entropy correction 
! o penalty_energy -- energy term from condition_penalty
!
!  OUTPUT
!   o total_energy -- total energy
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


	


!  local variables

      real*8 :: out_energy
      real*8 :: en_elecstat_tot
      real*8 :: tot_en_ll_vdw
      character*150 :: info_str
      real(dp) :: en_mpe_electrostatic
      real(dp) :: en_mpe_nonelectrostatic
      real(dp) :: en_mpe_total
      
!  begin work


      ! First, compute total energy only of the system itself - no embedding contribution
      total_energy = ev_sum + en_xc - &
                     en_pot_xc - &
                     0.5 * (hartree_energy_free + &
                     hartree_delta_energy) + &
                     en_vdw - en_pot_vdw + &
                     en_ll_vdw + en_lda_c - en_pbe_c


      if ( (n_periodic.eq.0) .and. (.not.force_new_functional) ) then
        total_energy = total_energy + en_ion_ion
      end if


      if ( ((use_embedding_potential) .and. (full_embedding)).or. use_embedding_pp ) then
        ! In this case, the external embedding potential is part of the sum of eigenvalues.
        ! For a total energy based on the density itself (discounting any external contributions),
        ! it needs to be subtracted out.
        total_energy = total_energy - en_density_embed - en_nonlocal
 
      end if

      if (solvent_method.eq.SOLVENT_MPB) then
	! We have to substract the terms added to the KS Hamiltonian due to 
	!dependence of dielectric function and ionic function on density
	!if dynamic_cavity and dynamic_ions = .false, these terms will be zero
	total_energy = total_energy - en_eps_rho - en_alpha_rho
	total_energy = total_energy - en_pot_Gnonmf
	!(the corresponding energy functional term is written out separately below)
      end if

      !MS Note: In case of MPE solvation, the _electrostatic_ interaction
      !  energy of the "reaction field" with the electron density is double
      !  counted in the sum of eigenvalues. This is already corrected for via
      !  also adding it to hartree_delta_energy.
      !  The _electrostatic_ interaction energy of the "reaction field" with
      !  nuclei is included in hartree_delta_energy.
      !  Non-electrostatic contributions are missing at this point and will
      !  be added in this routine.
      if (solvent_method.eq.SOLVENT_MPE) then
         call mpe_read_energy_contributions( &
            nonel=en_mpe_nonelectrostatic, &
            el=en_mpe_electrostatic, &
            tot=en_mpe_total )
         total_energy = total_energy + en_mpe_nonelectrostatic
      end if


      if ( use_plus_u ) then
        ! This is just a fix for the total energy.
        ! In the former DFT+U implementation a correction term to the total energy was missing.
        call plus_u_energy_correction_term
        total_energy = total_energy + plus_u_energy_correction
      end if

!  The total energy here is the Harris sum, where the Hartree term occurs both
!  in the sum of eigenvalues, and the double-counting Hartree term which is
!  subtracted.
!

         write(info_str, '(A)') ''
	 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X)') &
              "Total energy components:"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Sum of eigenvalues            :", ev_sum, " Ha", &
              ev_sum * hartree, " eV"
	 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| XC energy correction          :", en_xc, " Ha", &
              en_xc * hartree, " eV"
	 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| XC potential correction       :", - en_pot_xc, &
              " Ha", - en_pot_xc * hartree, " eV"
	 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         if ( (n_periodic.eq.0) .and. (.not.force_new_functional) ) then
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Free-atom Hartree contribution:", &
              - 0.5 * hartree_energy_free, &
              " Ha", - 0.5 * hartree_energy_free * hartree, " eV"
	 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         else
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Free-atom electrostatic energy:", &
              - 0.5 * hartree_energy_free, &
              " Ha", - 0.5 * hartree_energy_free * hartree, " eV"
	 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         end if
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Hartree energy correction     :", &
              - 0.5 * hartree_delta_energy, &
              " Ha", - 0.5 * hartree_delta_energy * hartree, " eV"
	 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         if ( (n_periodic.eq.0) .and. (.not.force_new_functional) ) then
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Ion-ion-interaction           :", en_ion_ion, &
              " Ha", en_ion_ion * hartree, " eV"
	 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         end if
         if ( use_plus_u ) then
            write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
               "| DFT+U energy correction       :", plus_u_energy_correction, " Ha", &
               plus_u_energy_correction * hartree, " eV"
      call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         end if
         if (use_embedding_potential.or.use_embedding_pp) then
           write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Ionic    embedding energy     :", &
              en_ion_embed , &
              " Ha", (en_ion_embed)*hartree, " eV"
	   call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
           write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Density  embedding energy     :", &
              en_density_embed , &
              " Ha", (en_density_embed)*hartree, " eV"
	   call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
           write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Nonlocal embedding energy     :", &
              en_nonlocal , &
              " Ha", (en_nonlocal)*hartree, " eV"
           call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
           if(use_nonlinear_core) then
             write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
                "| nonlinear core correction     :", &
                en_xc_nlcc , &
                " Ha", (en_xc_nlcc)*hartree, " eV"
             call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
           endif
           write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| External embedding energy     :", &
              en_ion_embed + en_density_embed + en_nonlocal + en_xc_nlcc, &
              " Ha", (en_ion_embed+en_density_embed + en_nonlocal + en_xc_nlcc)*hartree, " eV"
           call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

         end if
         if (use_vdw_correction .or. use_vdw_correction_hirshfeld.or.use_nlcorr_post&
              .or.use_vdw_post .or. use_libmbd) then   !SAG
           write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| vdW energy correction         :", &
              en_vdw, &
              " Ha", en_vdw * hartree, " eV"
           call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         end if
         if (use_vdw_correction_hirshfeld_sc) then
            write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
                 "| vdW energy correction         :", &
                 en_vdw, &
                 " Ha", en_vdw * hartree, " eV"
            call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
            write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
                 "| vdW potential correction      :", &
                 - en_pot_vdw, &
                 " Ha", - en_pot_vdw * hartree,   " eV"
            call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         end if 
         if (use_mbd_old .or. use_mbd_dev) then
             write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
                 "| MBD@rsSCS energy              :", &
                 en_vdw, " Ha", en_vdw*hartree, " eV"
             call localorb_info(info_str, use_unit, '(A)', OL_norm)
         end if
         if (use_mbd_std .and. (mbd_scf_converged .or. &
                mbd_self_consistent .and. .not. mbd_first_step)) then
             write (info_str, '(2X,A,1X,F20.8,A,1X,F20.8,A)') &
                 "| MBD@rsSCS energy              :", &
                 en_vdw, " Ha", en_vdw*hartree, " eV"
             call localorb_info(info_str, use_unit, '(A)', OL_norm)
             if (mbd_self_consistent) then
                write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
                    "| MBD@rsSCS potential           :", &
                    -en_pot_vdw, " Ha", -en_pot_vdw*hartree, " eV"
                call localorb_info(info_str, use_unit, '(A)', OL_norm)
             end if
         end if
         if (solvent_method.eq.SOLVENT_MPE) then
             write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
                 "| MPE non-electrostatic energy  :", &
                 en_mpe_nonelectrostatic , " Ha", &
                 en_mpe_nonelectrostatic*hartree , " eV"
             call localorb_info(info_str, use_unit, '(A)', OL_norm)
         end if


         if (use_ll_vdwdf) then
!
           tot_en_ll_vdw=en_ll_vdw+en_lda_c-en_pbe_c
           write (info_str,'(2X,A,2(1X,A,F9.5,A,F6.3,A,A))') &
              "| LL van der Waals energy corr. :", &
              "(",tot_en_ll_vdw, " +- ", en_ll_vdw_err, ")", "Ha",&
              "(",(tot_en_ll_vdw) * hartree, " +- ", en_ll_vdw_err*hartree, ")", " eV"
	   call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         end if
!
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Entropy correction            :", entropy_correction, &
              " Ha", entropy_correction * hartree, " eV"
	 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         if (penalty_energy /= 0.d0) then
            write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
            & "| Condition penalty             :", penalty_energy, &
            & " Ha", penalty_energy * hartree, " eV"
            call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         end if

!         write (info_str,'(2X,A,/,2X,A,1X,F20.8,A,1X,F20.8,A)') &
!              "| Energy correction for multipole", &
!              "| error in Hartree potential    :", &
!              0.5*hartree_multipole_correction, &
!              " Ha", 0.5*hartree_multipole_correction * hartree, " eV"
         write (info_str,'(2X,A,1X)') &
              "| ---------------------------"
	 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

! 	  write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)')&
! 	   "| Total energy w/o multipole corr.:",&
! 	    total_energy - 0.5* hartree_multipole_correction, " Ha",&
! 	   (total_energy - 0.5*hartree_multipole_correction)* hartree, " eV"

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Total energy                  :", &
              total_energy, " Ha", &
              (total_energy)* hartree, " eV"
	 call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

      if (occupation_type .eq. 2) then
!     methfessel-paxton
!         if (myid.eq.0) then
            write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,2A)') &
                 "| Total energy, T -> 0          :", &
                 (total_energy + 2 * &
                 (dble(n_methfessel_paxton + 1) / &
                 dble(n_methfessel_paxton + 2)) * entropy_correction), &
                 " Ha", &
                 (total_energy + 2 * &
                 (dble(n_methfessel_paxton + 1) / &
                 dble(n_methfessel_paxton + 2)) * entropy_correction) &
                 * hartree, " eV  <-- do not rely on this value for ", &
                 "anything but (periodic) metals"
	    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
!         end if

      else

!         if (myid.eq.0) then
            write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,2A)') &
                 "| Total energy, T -> 0          :", &
                 (total_energy + entropy_correction), " Ha", &
                 (total_energy + entropy_correction) * hartree, &
                 " eV  <-- do not rely on this value for anything but ", &
                 "(periodic) metals"
	    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
!         end if

      end if

         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Electronic free energy        :", &
              (total_energy + 2 * entropy_correction), " Ha", &
              (total_energy + 2 * entropy_correction) * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )


      if (use_embedding_potential.or.use_embedding_pp) then

        ! also add contributions including embedding
!        if (myid.eq.0) then

          write(info_str, '(A)') ''
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
          write(info_str,'(2X,A)') &
          "Total energy including external embedding contributions:"
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
          out_energy = (total_energy+en_density_embed+en_ion_embed+en_nonlocal+en_xc_nlcc)

          write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
          "| Total energy + embedding      :", &
          out_energy, " Ha", &
          out_energy* hartree, " eV"
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
          if (occupation_type .eq. 2) then
!           methfessel-paxton
            out_energy = (total_energy + 2 * &
            (dble(n_methfessel_paxton + 1) / &
            dble(n_methfessel_paxton + 2)) * entropy_correction) &
            +en_density_embed+en_ion_embed+en_nonlocal+en_xc_nlcc
          else
            out_energy = total_energy + entropy_correction &
            +en_density_embed+en_ion_embed+en_nonlocal+en_xc_nlcc
          end if
          write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
          "| Total energy + embedding, T->0:", &
          out_energy, " Ha", &
          out_energy* hartree, " eV  <-- do not rely on this value for anything but (periodic) metals"

          out_energy = (total_energy + 2 * entropy_correction) &
            +en_density_embed+en_ion_embed + en_nonlocal+ en_xc_nlcc
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

          write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
          "| Elec. free energy + embedding :", &
          out_energy, " Ha", &
          out_energy* hartree, " eV"
          call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
        end if
!      end if

         write(info_str, '(A)') ''
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write(info_str,'(2X,A)') "Derived energy quantities:"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         if (force_new_functional) then
           ! The formula below is correct only for the "new" (default!) functional
	   if (solvent_method.eq.SOLVENT_MPB) then
                   if (use_hartree_fock) then
                     write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
                      "| Kinetic energy                :", &
                       ev_sum - en_pot_xc -en_alpha_rho-en_eps_rho-en_pot_Gnonmf  -2* en_fock &
                      -(en_elec_free + en_elec_delta) - hartree_multipole_correction &
                       ,  " Ha", &
                       (ev_sum - en_pot_xc -2* en_fock &
                      -(en_elec_free + en_elec_delta) - hartree_multipole_correction) * hartree &
                     , " eV"
                     call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
                   else
                     write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
                      "| Kinetic energy                :", &
                       ev_sum - en_pot_xc-en_eps_rho-en_alpha_rho-en_pot_Gnonmf &
                      -(en_elec_free + en_elec_delta) - hartree_multipole_correction &
                       ,  " Ha", &
                       (ev_sum - en_pot_xc &
                      -(en_elec_free + en_elec_delta) - hartree_multipole_correction) * hartree &
                     , " eV"
                     call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
                   endif ! use_hartree_fock
          else
                   if (use_hartree_fock) then
                     write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
                      "| Kinetic energy                :", &
                       ev_sum - en_pot_xc  -2* en_fock &
                      -(en_elec_free + en_elec_delta) - hartree_multipole_correction &
                       ,  " Ha", &
                       (ev_sum - en_pot_xc -2* en_fock &
                      -(en_elec_free + en_elec_delta) - hartree_multipole_correction) * hartree &
                     , " eV"
                     call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
                   else
                     write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
                      "| Kinetic energy                :", &
                       ev_sum - en_pot_xc &
                      -(en_elec_free + en_elec_delta) - hartree_multipole_correction &
                       ,  " Ha", &
                       (ev_sum - en_pot_xc &
                      -(en_elec_free + en_elec_delta) - hartree_multipole_correction) * hartree &
                     , " eV"
                     call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
                   endif ! use_hartree_fock
	  end if

!  total electrostatic energy including electron-electron, electron-nuclear, and 
!  nuclear-nuclear interaction

           en_elecstat_tot = en_elec_free + en_elec_delta + hartree_multipole_correction &
                            - 0.5d0*(hartree_energy_free + hartree_delta_energy)
             write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Electrostatic energy          :", &
               en_elecstat_tot , " Ha", &
               en_elecstat_tot*hartree , " eV"
             call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         end if ! force_new_functional

         write (info_str,'(2X,A)') &
              "| Energy correction for multipole"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| error in Hartree potential    :", &
              0.5*hartree_multipole_correction, &
              " Ha", 0.5*hartree_multipole_correction * hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
         write(info_str,'(2X,A,F20.8,A)') &
              "| Sum of eigenvalues per atom                           : ", &
              ev_sum*hartree/dble(n_real_atoms), " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

!       if (occupation_type .eq. 2) then
! !     methfessel-paxton
! 
!          if (myid.eq.0) then
!             write(use_unit,'(2X,A,F20.8,A)') &
!                  "| Total energy (T->0) per atom                      : ", &
!                  (total_energy + 2 * &
!                  (dble(n_methfessel_paxton + 1) / &
!                  dble(n_methfessel_paxton + 2)) * entropy_correction)* &
!                  hartree/dble(n_atoms), " eV"
!          end if
! 
!       else
! 
!          if (myid.eq.0) then
!             write(use_unit,'(2X,A,F20.8,A)') &
!                  "| Total energy (T->0) per atom                      : ", &
!                (total_energy+entropy_correction)*hartree/dble(n_atoms), &
!                  " eV"
!          end if
! 
!       end if
! 
!       if (myid.eq.0) then
!          write(use_unit,'(2X,A,F20.8,A)') &
!               "| Free energy per atom                              : ", &
!               (total_energy + 2 * entropy_correction) * &
!               hartree/dble(n_atoms), " eV"
!       end if
      if (occupation_type .eq. 2) then
!     methfessel-paxton

!         if (myid.eq.0) then
            write(info_str,'(2X,A,F20.8,A)') &
                 "| Total energy (T->0) per atom                          : ", &
                 (total_energy + 2 * &
                 (dble(n_methfessel_paxton + 1) / &
                 dble(n_methfessel_paxton + 2)) * entropy_correction)* &
                 hartree/dble(n_real_atoms), " eV  <-- do not rely on this value for anything but (periodic) metals"
	     call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
!         end if

      else

!         if (myid.eq.0) then
            write(info_str,'(2X,A,F20.8,A)') &
                 "| Total energy (T->0) per atom                          : ", &
               (total_energy+entropy_correction)*hartree/dble(n_real_atoms), &
                 " eV  <-- do not rely on this value for anything but (periodic) metals"
	   call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
!         end if

      end if

!      if (myid.eq.0) then
         write(info_str,'(2X,A,F20.8,A)') &
              "| Electronic free energy per atom                       : ", &
              (total_energy + 2 * entropy_correction) * &
              hartree/dble(n_real_atoms), " eV"
	call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
!      end if


      if (solvent_method.eq.SOLVENT_MPE) then
         call localorb_info(" ", use_unit, '(A)', OL_norm)
         write (info_str,'(2X,A)') &
            "MPE Solvation Free Energy contributions:"
         call localorb_info(info_str, use_unit, '(A)', OL_norm)
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
            "| MPE electrostatic energy      :", &
            en_mpe_electrostatic , " Ha", &
            en_mpe_electrostatic*hartree , " eV"
         call localorb_info(info_str, use_unit, '(A)', OL_norm)
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
            "| MPE non-electrostatic energy  :", &
            en_mpe_nonelectrostatic , " Ha", &
            en_mpe_nonelectrostatic*hartree , " eV"
         call localorb_info(info_str, use_unit, '(A)', OL_norm)
         write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
            "| MPE total energy              :", &
            en_mpe_total , " Ha", &
            en_mpe_total*hartree , " eV"
         call localorb_info(info_str, use_unit, '(A)', OL_norm)
      end if

	if (solvent_method.eq.SOLVENT_MPB.and.surface_and_volume_calc) then
	  write (info_str,'(2X,A)') &
		" "
	  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
	  write (info_str,'(2X,A)') &
		"MPBE Solvation Additional Free Energies:"
	  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
	  if (.not.(solve_lpbe_only.and..not.use_mpbe_free_energy).and..not.mpbe_no_ion) then
             write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Free Energy in Electrolyte          :", &
               total_energy + en_Gnonmf + freeen_MPBE_nonelstat_solv_energy , " Ha", &
               (total_energy + en_Gnonmf + freeen_MPBE_nonelstat_solv_energy)*hartree , " eV"
             call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
	  else
             write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
              "| Free Energy in Electrolyte          :", &
               total_energy + en_Gnonmf + freeen_LPBE_nonelstat_solv_energy , " Ha", &
               (total_energy + en_Gnonmf + freeen_LPBE_nonelstat_solv_energy)*hartree , " eV"
             call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
	  end if
	  write (info_str,'(2X,A,1X,F20.8,A)') &
		"| Surface Area of Cavity            :", surface_mpb, " a.u."
	  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
	  write (info_str,'(2X,A,1X,F20.8,A)') &
		"| Volume of Cavity             :", volume_mpb, " a.u."
	  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
	  write (info_str,'(2X,A,1X,F20.8,A)') &
		"| Nonelectrostatic Free Energy             :", &
		en_Gnonmf , " Ha"
	  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
	  if (.not.(solve_lpbe_only.and..not.use_mpbe_free_energy).and..not.mpbe_no_ion) then
	    write (info_str,'(2X,A,F20.8,A)') &
		  "| Additional Nonelstatic MPBE Solvation Energy:",&
		  freeen_MPBE_nonelstat_solv_energy," Ha"
	    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
	  else if (solve_lpbe_only.and..not.use_mpbe_free_energy) then
	    write (info_str,'(2X,A,F20.8,A)') &
		  "| Additional Nonelstatic LPBE Solvation Energy:",&
		  freeen_LPBE_nonelstat_solv_energy," Ha"
	    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
	  end if
	end if

       ! Paula: The total_energy need to be the total energy, 
       ! so that relaxation work correctly.
       total_energy = total_energy  + en_density_embed + en_ion_embed + en_nonlocal + en_xc_nlcc

      end subroutine
!******	
