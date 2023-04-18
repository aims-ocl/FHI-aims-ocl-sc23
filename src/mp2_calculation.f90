!****s* FHI-aims/mp2_calculation
!  NAME
!   mp2_calculation
!  SYNOPSIS 

      subroutine mp2_calculation &
       ( )

!  PURPOSE
!    Calculation of the total energy perturbation terms up to the second order using Moller-Plesset theory (MP2) to the HF/DFT Energy.
!    A variant, called SCS-MP2, scales the two spin channels separately.
!    Frozen core approximation is available too. 
! USES

      use physics
      use runtime_choices
      use hartree_fock
      use mpi_tasks
      use prodbas
      use localorb_io, only: use_unit
      use dimensions, only: use_dftpt2, use_os_mp2

      implicit none

!  ARGUMENTS

!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society.
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


      if (myid.eq.0) then
          write(use_unit,'(A)') &
          &      "---------------------------------------------"
          write(use_unit,'(A)') &
          &      "---------------------------------------------"
          if (lrc_pt2_started) then
              write(use_unit,'(A)') " | lrc-PT2 calculation starts ..."
          else if (use_dftpt2) then
              write(use_unit,'(A)') " | PT2 calculation starts ..."
          else
              write(use_unit,'(A)') " | MP2 calculation starts ..."
          end if
          write(use_unit,'(A)') &
          &      "---------------------------------------------"
          write(use_unit,'(A)') &
          &      "---------------------------------------------"
          write(use_unit,'(A)') ""
      endif

      if (flag_en_shift .and. use_2d_corr) then
          write(use_unit,'(1X,A)')&
               "*--------------------------------------------"
        write(use_unit,'(1X,A)')&
               "* Error: Currently only lapack implementation is"
        write(use_unit,'(1X,A)')&
               "* avaliable for (s)BGE2. Aborting sBGE2 part"
        write(use_unit,'(1X,A)')&
               "*--------------------------------------------"
        stop

      endif
      ! en_shift_type = 1: for hole-hole EN2
      !                 2: for constant energy shift
      !                 3: for (screening) coupling PT2
      !                 4: for (screening) coupling PT2 + the first-order term
      if(flag_en_shift .and. (en_shift_type .eq. 1)) then
          ! the simplest level shift
          call evaluate_en_mp2_correlation_energy &
             (n_electrons,dftpt2_total_energy_dft,en_xc, &
              occ_numbers,KS_eigenvalue,KS_eigenvector,&
              post_scf_total_energy) 
      elseif(flag_en_shift .and. (en_shift_type .eq. 3)) then
          ! sBGE2, also called IEPA in quantum chemistry
          call evaluate_iepa_mp2_correlation_energy &
             (n_electrons,dftpt2_total_energy_dft,en_xc, &
              occ_numbers,KS_eigenvalue,KS_eigenvector,&
              post_scf_total_energy) 
      elseif(flag_en_shift .and. (en_shift_type .eq. 4)) then
          ! extended version of sBGE2 (still experimental)
          call evaluate_iepa_mp2_correlation_energy_v02 &
             (n_electrons,dftpt2_total_energy_dft,en_xc, &
              occ_numbers,KS_eigenvalue,KS_eigenvector,&
              post_scf_total_energy) 
      else
          if (use_dftpt2) then
            if(use_2d_corr) then
              call evaluate_dftpt2_correlation_2 &
                 (n_electrons,dftpt2_total_energy_dft,en_xc, &
                  occ_numbers,KS_eigenvalue,KS_eigenvector,&
                  post_scf_total_energy) 
            else
              call evaluate_dftpt2_correlation &
                 (n_electrons,dftpt2_total_energy_dft,en_xc, &
                  occ_numbers,KS_eigenvalue,KS_eigenvector,&
                  post_scf_total_energy) 
            endif
          else
            if(use_2d_corr) then
                call evaluate_mp2_correlation_energy_2 &
                   (n_electrons,total_energy,en_xc, &
                    occ_numbers,KS_eigenvalue,KS_eigenvector,&
                    post_scf_total_energy) 
            else
              if (use_os_mp2) then
                call evaluate_osmp2_correlation_energy &
                   (n_electrons,total_energy,en_xc, &
                    occ_numbers,KS_eigenvalue,KS_eigenvector,&
                    post_scf_total_energy) 
              else
                call evaluate_mp2_correlation_energy &
                   (n_electrons,total_energy,en_xc, &
                    occ_numbers,KS_eigenvalue,KS_eigenvector,&
                    post_scf_total_energy) 
              endif
            endif
          endif
      endif

      return

      end subroutine mp2_calculation
!****** 

