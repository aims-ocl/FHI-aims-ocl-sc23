!****s* FHI-aims/output_real_eigenfunctions
!  NAME
!   output_real_eigenfunctions
!  SYNOPSIS

      subroutine output_real_eigenfunctions &
      ( KS_eigenvalue, KS_eigenvector, occ_numbers &
      )

!  PURPOSE
!  Subroutines output_*_eigenfunctions() writes eigenvalues and eigenvectors
!  of the current Kohn-Sham eigenstates.
!
!  This is real version!
!
!  USES

      use dimensions
      use runtime_choices
      use basis
      use localorb_io
      use constants
      use pbc_lists, only : k_point_list

! >>> AB: feb 2012
      use timing,    only : number_of_loops
! <<< end with insert: AB, feb 2012
      use get_dielectric_function
      implicit none

!  ARGUMENTS

      real*8 KS_eigenvalue (n_states, n_spin,n_k_points)
      real*8 KS_eigenvector (n_basis, n_states, n_spin)
      real*8, dimension(n_states, n_spin,n_k_points) :: occ_numbers

!  INPUTS
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o KS_eigenvector -- Kohn-Sham eigenvectors
!   o occ_numbers -- occupation weights of eigenstates
!  OUTPUT
!    none
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

      character l_char

      character*22, dimension(n_spin) :: file_name

      character*150 :: info_str
      real*8, dimension(n_spin) :: n_elec_per_channel

      real*8 :: homo_level, lumo_level

!  counters
      integer i_basis, i_state, i_fn
      integer i_state_2
      integer i_spin, i_spin2
      integer i_k_point

! >>> AB: feb 2012
!     counters
      integer j, mincol, maxcol
!     output file
      integer mosfile
! <<< done with insert: AB, feb 2012

!  functions

      character l_to_str

      character(*), parameter :: func = 'output_real_eigenfunctions'

!  begin work

      write (info_str,'(2X,A,A)') &
           "Writing Kohn-Sham eigenvalues."
      call localorb_info(info_str,use_unit,'(A)', OL_low)

      n_elec_per_channel = 0.d0

      do i_k_point = 1, min(n_k_points,out_k_points_eigenvalues)

         if (spin_treatment.eq.1) then
            write(info_str,*)
            call localorb_info(info_str,use_unit,'(A)', OL_low)
            write(info_str,'(2X,A)') "Spin-up eigenvalues: "
            call localorb_info(info_str,use_unit,'(A)',OL_low)
         end if

         do i_spin = 1, n_spin, 1

            if (i_spin.eq.2) then
               write(info_str,'(2X,A)') "Spin-down eigenvalues: "
               call localorb_info(info_str,use_unit,'(A)',OL_low)

            end if

            if(n_k_points>1)then
               write(info_str,'(2X,A,I8,A,3(2X,F10.6), A)') &
                  'K-point:',i_k_point, ' at',k_point_list(i_k_point,1), &
                  k_point_list(i_k_point,2),k_point_list(i_k_point,3), &
                  ' (in units of recip. lattice)'
               call localorb_info(info_str,use_unit,'(A)',OL_low)
            end if

            write(info_str,*)
            call localorb_info(info_str,use_unit,'(A)',OL_low)

            write(info_str,'(2X,A,4X,A,4X,A,4X,A)') &
                 "State", "Occupation", "Eigenvalue [Ha]", "Eigenvalue [eV]"
            call localorb_info(info_str,use_unit,'(A)',OL_low)
            do i_state = 1, n_states, 1
               n_elec_per_channel(i_spin) =  n_elec_per_channel(i_spin) + &
                    occ_numbers(i_state, i_spin, i_k_point)

               write(info_str,'(2X,I5,6X,F8.5,5X,F14.6,4X,F15.5)') &
                    i_state, occ_numbers(i_state, i_spin, i_k_point), &
                    (KS_eigenvalue(i_state, i_spin, i_k_point)), &
                    (KS_eigenvalue(i_state, i_spin, i_k_point)*hartree)
               call localorb_info(info_str,use_unit,'(A)',OL_low)
            enddo
            write(info_str,*)
            call localorb_info(info_str,use_unit,'(A)',OL_low)

         enddo
      end do

      if (n_spin .eq. 2) then
         write(info_str,'(2X,A)') 'Current spin moment of the entire structure :'
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
        if (n_k_points.le.1) then

           write(info_str,'(2X,A,F8.2)') '| N = N_up - N_down :', &
                abs(n_elec_per_channel(2) - n_elec_per_channel(1))
           call localorb_info(info_str,use_unit,'(A)',OL_norm)
           write(info_str,'(2X,A,F8.2)') '| S                 :', &
                abs(n_elec_per_channel(2) - n_elec_per_channel(1)) * 0.5d0
           call localorb_info(info_str,use_unit,'(A)',OL_norm)
           write(info_str,'(2X,A,F8.2)') '| J                 :', &
                abs(n_elec_per_channel(2) - n_elec_per_channel(1)) + 1.d0
           call localorb_info(info_str,use_unit,'(A)',OL_norm)
           write(info_str,*)
           call localorb_info(info_str,use_unit,'(A)',OL_norm)

        else
        ! compute overall occupation of spin channels as sum over ALL k points

           n_elec_per_channel = 0.d0

           do i_k_point = 1, n_k_points
              do i_spin = 1, n_spin, 1
                 do i_state = 1, n_states, 1
                    n_elec_per_channel(i_spin) = n_elec_per_channel(i_spin) + &
                         occ_numbers(i_state, i_spin, i_k_point)
                 end do
              end do
           end do

           write(info_str,'(2X,A,F16.5)') '| N = N_up - N_down (sum over all k points):', &
                abs(n_elec_per_channel(2) - n_elec_per_channel(1)) / (n_k_points * 1.d0)
           call localorb_info(info_str,use_unit,'(A)',OL_norm)
           write(info_str,'(2X,A,F16.5)') '| S (sum over all k points)                :', &
                (abs(n_elec_per_channel(2) - n_elec_per_channel(1)) * 0.5d0) / (n_k_points * 1.d0)
           call localorb_info(info_str,use_unit,'(A)',OL_norm)
           write(info_str,*)
           call localorb_info(info_str,use_unit,'(A)',OL_norm)

        end if

      end if

      ! Find HOMO and LUMO levels, output results to screen
      call find_and_output_homo_lumo_gap &
           ( KS_eigenvalue, occ_numbers, homo_level, lumo_level )

      !>> KR setting of homo/lumo variables that have to be known in module "get_dielectric_function"
      gk_homo_level = homo_level
      gk_lumo_level = lumo_level
      !<< KR


      ! Output here works only for non-periodic systems in this way.
      ! (and also not for scalapack when eigenvectors are not collected.)
      if ( (out_eigenvec) .and. (n_periodic.eq.0) .and. (.not.(use_scalapack.and.(.not.collect_eigenvectors)) )) then

         write(info_str,'(2X,A,A)') &
              "Writing Kohn-Sham orbitals (eigenvectors) to file ", &
              "KS_eigenvectors.out ."
         call localorb_info(info_str,use_unit,'(A)',OL_low)

         if (n_spin.eq.1) then
            write(file_name(1),'(A)') "KS_eigenvectors.out"
         else if (n_spin.eq.2) then
            write(file_name(1),'(A)') "KS_eigenvectors_up.out"
            write(file_name(2),'(A)') "KS_eigenvectors_dn.out"
         end if

         do i_spin = 1, n_spin,1

            open (50, FILE=file_name(i_spin))

            i_state = 1

            do while (i_state.le.(n_states-4))

               write(info_str,'(30X,5(4X,I4,4X))') &
                    (i_state_2, i_state_2 = i_state, i_state+4, 1)
               call localorb_info(info_str,50,'(A)',OL_low)

               do i_basis = 1, n_basis, 1

                  i_fn = basis_fn(i_basis)
                  l_char = l_to_str(basis_l(i_basis))

                  write(info_str, &
                       '(I5,1X,I5,1X,A8,I3,1X,A1,1X,I3,3X,5(F15.9,2X))') &
                       i_basis, basis_atom(i_basis), basisfn_type(i_fn), &
                       basisfn_n(i_fn), l_char, basis_m(i_basis), &
                       (KS_eigenvector(i_basis, i_state_2,i_spin), &
                       i_state_2 = i_state, i_state+4, 1)
                  call localorb_info(info_str,50,'(A)',OL_low)

               enddo

               i_state = i_state+5
               write(info_str,*)
               call localorb_info(info_str,50,'(A)',OL_low)

            enddo

            write(info_str,'(30X,5(4X,I4,4X))') &
                 (i_state_2, i_state_2 = i_state, n_states, 1)
            call localorb_info(info_str,50,'(A)',OL_low)

            do i_basis = 1, n_basis, 1

               i_fn = basis_fn(i_basis)
               l_char = l_to_str(basis_l(i_basis))

               write(info_str, &
                    '(I5,1X,I5,1X,A8,I3,1X,A1,1X,I3,3X,5(F15.9,2X))') &
                    i_basis, basis_atom(i_basis), basisfn_type(i_fn), &
                    basisfn_n(i_fn), l_char, basis_m(i_basis), &
                    (KS_eigenvector(i_basis, i_state_2,i_spin), &
                    i_state_2 = i_state, n_states, 1)
               call localorb_info(info_str,50,'(A)',OL_low)

            enddo

            close(50)

            ! spin
         enddo

      end if

! >>> AB: feb 2012
!     output KS orbitals and energies to external file(s) in a TURBOMOLE
!     format compatible with "aitranss" module (Karlsruhe transport code)
!
!     Output here works only for non-periodic systems in this way.
!     (and also not for scalapack when eigenvectors are not collected.)
!
!     >>> output KS eigenfunctions only if number_of_loops >= 1 ! <<< AB: march 2012
      if ( (out_aitranss.and.(number_of_loops>=1)) .and. (n_periodic.eq.0) .and. &
           (.not.(use_scalapack.and.(.not.collect_eigenvectors)) ) ) then

          if (n_spin.eq.1) then
           write(info_str,'(2x,a)') &
             'Writing Kohn-Sham orbitals (eigenvectors) to a file "mos.aims" '
          else
           write(info_str,'(2x,a)') &
             'Writing spin-resolved Kohn-Sham orbitals to files "alpha.aims" and "beta.aims" '
          end if
          call localorb_info(info_str,use_unit,'(a)',OL_low)
          write(info_str,'(2x,a)') ''
          call localorb_info(info_str,use_unit,'(a)',OL_low)

          if (n_spin.eq.1) then
            write(file_name(1),'(A)') 'mos.aims'
          else if (n_spin.eq.2) then
            write(file_name(1),'(A)') 'alpha.aims'
            write(file_name(2),'(A)') 'beta.aims'
          end if

          mosfile = 50
          do i_spin = 1, n_spin

           open(mosfile,file=file_name(i_spin))
!          put header >>>
           if (n_spin.eq.1) then
             write(info_str,'(a,3x,a,i4,3x,a)') &
                   '$scfmo.aims', 'iteration: ', number_of_loops, 'format(4d20.14)'
             call localorb_info(info_str,mosfile,'(a)',OL_low)
           else if (i_spin.eq.1) then
             write(info_str,'(a,3x,a,i4,3x,a)') &
                   '$uhfmo_alpha.aims', 'iteration: ', number_of_loops, 'format(4d20.14)'
             call localorb_info(info_str,mosfile,'(a)',OL_low)
           else
             write(info_str,'(a,3x,a,i4,3x,a)') &
                   '$uhfmo_beta.aims', 'iteration: ', number_of_loops, 'format(4d20.14)'
             call localorb_info(info_str,mosfile,'(a)',OL_low)
           end if

           do i_state = 1, n_basis
!            output a KS orbital energy >>>
             write(info_str,'(i6,2x,a,6x,a,d20.14,3x,a,i6)') &
                   i_state, 'a', 'eigenvalue=', &
                   KS_eigenvalue(i_state, i_spin, n_k_points), 'nsaos=', n_basis
             call localorb_info(info_str,mosfile,'(a)',OL_low)
!            output an eigenvector >>>
             maxcol = 0
             do while (maxcol<n_basis)
               mincol = maxcol+1
               maxcol = min(maxcol+4,n_basis)
               write(info_str,'(4(d20.14))') &
                    (KS_eigenvector(j,i_state,i_spin),j=mincol,maxcol)
               call localorb_info(info_str,mosfile,'(a)',OL_low)
             end do
           end do ! i_state

!          all done: close file >>>
           write(info_str,'(a)') '$end'
           call localorb_info(info_str,mosfile,'(a)',OL_low)
           close(mosfile)

          end do ! i_spin

      end if
! <<< end with insert: AB, jan 2012

      return
    end subroutine output_real_eigenfunctions
!----------------------------------------------------------------------
!******

!****s* FHI-aims/find_and_output_homo_lumo_gap
!  NAME
!   find_and_output_homo_lumo_gap
!  SYNOPSIS
subroutine find_and_output_homo_lumo_gap &
      ( KS_eigenvalue, occ_numbers, homo_level, lumo_level )
!  PURPOSE
!    Find the HOMO and LUMO levels and return them, also write results to screen
!  USES
  use dimensions, only : n_states, n_spin, spin_degeneracy, n_periodic, use_gw,&
      n_k_points
  use pbc_lists, only : k_point_list
  use physics, only : estimate_low_gap
  use constants, only : hartree
  use localorb_io, only : localorb_info, use_unit, OL_norm
  implicit none

  real*8, dimension(n_states,n_spin,n_k_points), intent(in) :: KS_eigenvalue
  real*8, dimension(n_states,n_spin,n_k_points), intent(in) :: occ_numbers
  real*8, intent(out) :: homo_level
  real*8, intent(out) :: lumo_level

  real*8  :: homo_occ, lumo_occ
  integer :: i_kpt_homo, i_kpt_lumo
  integer :: i_spin_homo, i_spin_lumo
  logical :: found_min_direct_gap
  real*8  :: min_direct_gap
  integer :: min_direct_gap_k_point !k-point at which the direct gap is smallest
  integer :: min_direct_homo_spin, min_direct_lumo_spin

  character*150 :: info_str

  real*8 :: occupation_def
  logical :: fractionally_occupied

  ! Here's the "find" part of the subroutine
  call find_homo_lumo_gap &
      ( KS_eigenvalue, occ_numbers, homo_level, lumo_level, homo_occ, lumo_occ,&
        i_kpt_homo, i_kpt_lumo, i_spin_homo, i_spin_lumo, found_min_direct_gap,&
        min_direct_gap, min_direct_gap_k_point, min_direct_homo_spin, &
        min_direct_lumo_spin )

  ! And now for the "output" part of the subroutine

  if (n_k_points.gt.1) then

     write(info_str,'(2X,A)') &
        'What follows are estimated values for band gap, HOMO, LUMO, etc.'
     call localorb_info(info_str, use_unit, "(A)", OL_norm)
     write(info_str,'(2X,A)') &
        '| They are estimated on a discrete k-point grid and not necessarily exact.'
     call localorb_info(info_str, use_unit, "(A)", OL_norm)
     write(info_str,'(2X,A)') &
        '| For converged numbers, create a DOS and/or band structure plot on a denser k-grid.'
      call localorb_info(info_str, use_unit, "(A)", OL_norm)
     write(info_str,'(2X,A)') &
        ''
      call localorb_info(info_str, use_unit, "(A)", OL_norm)

  end if

  ! Write HOMO (VBM) information
  if (i_kpt_homo.le.0) then
    ! We never hit a single level that counts as occupied.
    ! Unlikely, but could happen for a positively charged hydrogen-like ion.
    write(info_str,'(2X,A)') "*** Could not determine highest occupied state (HOMO or VBM)."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A)') "*** This could happen, but only with less than one electron in the system."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
  else
    ! HOMO was found, output as planned.
    if (n_periodic.eq.0) then
      write(info_str,'(2X,A,F15.8,A)') "Highest occupied state (VBM) at ",homo_level*hartree, &
      " eV "
    else
      if(use_gw) then
        write(info_str,'(2X,A,F15.8,A)') "Highest occupied state (VBM) after the GW correction at ",homo_level*hartree, &
        " eV (relative to internal zero)"
      else
        write(info_str,'(2X,A,F15.8,A)') "Highest occupied state (VBM) at ",homo_level*hartree, &
        " eV (relative to internal zero)"
      endif
    end if
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,F15.8,A)') "| Occupation number: ",homo_occ
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    if (n_periodic.gt.0) then
      write(info_str,'(2X,A,I8,A,3(2X,F10.6), A)') &
          '| K-point:',i_kpt_homo, ' at',k_point_list(i_kpt_homo,1), &
          k_point_list(i_kpt_homo,2),k_point_list(i_kpt_homo,3), &
          ' (in units of recip. lattice)'
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if
    if (n_spin.gt.1) then
      write(info_str,'(2X,A,I8)') &
          '| Spin channel: ',i_spin_homo
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if
    if ((homo_level .gt. 0d0).and.(n_periodic .eq. 0)) then
      write(info_str,'(1X,A)') "* WARNING: The VBM is higher than the internal energy zero!"
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if
    write(info_str,*)
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
  end if

  ! Write LUMO (CBM) information
  if (i_kpt_lumo.le.0) then
    ! that should be a good enough indicator that we never hit a single level which was "unoccupied"
    ! This could happen for a minimal basis.
    write(info_str,'(2X,A)') "*** Could not determine lowest unoccupied state (LUMO or CBM)."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A)') "*** This may happen in a calculation with a minimal basis set."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
  else
    ! LUMO was found, output as planned.
    if (n_periodic.eq.0) then
      write(info_str,'(2X,A,F15.8,A)') "Lowest unoccupied state (CBM) at",lumo_level*hartree, &
      " eV "
    else
      if(use_gw) then
         write(info_str,'(2X,A,F15.8,A)') "Lowest unoccupied state (CBM) after the GW correction at",lumo_level*hartree, &
         " eV (relative to internal zero)"
      else
         write(info_str,'(2X,A,F15.8,A)') "Lowest unoccupied state (CBM) at",lumo_level*hartree, &
        " eV (relative to internal zero)"
      endif
    end if
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,F15.8,A)') "| Occupation number: ",lumo_occ
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    if (n_periodic.gt.0) then
      write(info_str,'(2X,A,I8,A,3(2X,F10.6), A)') &
          '| K-point:',i_kpt_lumo, ' at',k_point_list(i_kpt_lumo,1), &
          k_point_list(i_kpt_lumo,2),k_point_list(i_kpt_lumo,3), &
          ' (in units of recip. lattice)'
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if
    if (n_spin.gt.1) then
      write(info_str,'(2X,A,I8)') &
          '| Spin channel: ',i_spin_lumo
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if
    write(info_str,*)
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
  end if

  if ( (i_kpt_homo.gt.0) .and. (i_kpt_lumo.gt.0) ) then

       if (n_periodic.gt.0) then
         write(info_str,'(2X,A,F15.8,A,I0,A,I0)') 'ESTIMATED overall HOMO-LUMO gap: ', (lumo_level-homo_level)*hartree, &
                             & ' eV between HOMO at k-point ', i_kpt_homo, &
                             &  ' and LUMO at k-point ',  i_kpt_lumo
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
       else
         write(info_str,'(2X,A,F15.8,A)') 'Overall HOMO-LUMO gap: ', (lumo_level-homo_level)*hartree, &
                             & ' eV.'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
       end if

       if (i_kpt_homo .ne. i_kpt_lumo ) then
         write(info_str,'(2X,A)') '| This appears to be an indirect band gap.'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)

         if (found_min_direct_gap) then

           write(info_str,'(2X,A,F15.8,A,I0,A,3(2X,F10.6), A)') '| Smallest direct gap : ', min_direct_gap*hartree, &
                       & ' eV for k_point ', min_direct_gap_k_point, ' at',k_point_list(i_kpt_lumo,1), &
                       & k_point_list(i_kpt_lumo,2),k_point_list(i_kpt_lumo,3), &
          ' (in units of recip. lattice)'
           call localorb_info(info_str, use_unit, "(A)", OL_norm)

           if (n_spin.gt.1) then
             write(info_str,'(2X,A,I5)') '| Direct gap HOMO spin channel : ', min_direct_homo_spin
             call localorb_info(info_str, use_unit, "(A)", OL_norm)
             write(info_str,'(2X,A,I5)') '| Direct gap LUMO spin channel : ', min_direct_lumo_spin
             call localorb_info(info_str, use_unit, "(A)", OL_norm)
           end if

         else

           write(info_str,'(2X,A)') '| A valid direct gap (defined HOMO and LUMO at same k-point) could not be found.'
           call localorb_info(info_str, use_unit, "(A)", OL_norm)

         end if

       else

         if (n_k_points.gt.1) then
            write(info_str,'(2X,A)') '| This appears to be a direct band gap.'
            call localorb_info(info_str, use_unit, "(A)", OL_norm)
         end if

       end if
  end if

  ! Determine if there are unambiguously fractionally occupied orbitals in the system
  ! (not just a small effect)
  ! Define a tolerance that determines whether the states in a system are definitely
  ! fractionally occupied. Currently, this is hardwired here since it is not an obvious
  ! parameter with which a user needs to interact. If this ever changes, it should
  ! become configurable in control.in.
  occupation_def = 0.05d0
  if ( (lumo_occ.ge.occupation_def) .or. (homo_occ.le.(spin_degeneracy-occupation_def)) ) then
     fractionally_occupied = .true.
     ! no further action here for now
  else
     fractionally_occupied = .false.
  end if

    if (n_periodic.eq.0) then
       if (fractionally_occupied) then
         write(info_str,'(2X,A)') &
           'This system has fractional occupation numbers. Thus, it is a degenerate system'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         write(info_str,'(2X,A)') &
           'within the approximate finite broadening function (occupation_type)'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         write(info_str,'(2X,A)') &
           'applied to determine the occupation numbers.'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         write(info_str,'(2X,A)') &
           'This may not be harmful. However, for some systems (free atoms), a different, non-degenerate'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         write(info_str,'(2X,A)') &
           'ground state might exist by using a smaller broadening (see keyword occupation_type).'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         estimate_low_gap = .true.
       else if (((lumo_level-homo_level)*hartree) .lt.0.2) then ! 0.2 eV
         estimate_low_gap = .true.
       else
         estimate_low_gap = .false.
       end if
    end if

    if (n_periodic.gt.0) then
       ! In this case, we try to characterize the band structure a little bit.
       ! For molecules, the wording does not apply.
       if (fractionally_occupied) then
         write(info_str,'(2X,A)') &
           'However, this system has fractional occupation numbers. Since we use a finite k-point grid,'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         write(info_str,'(2X,A)') &
           'this material is metallic within the approximate finite broadening function (occupation_type)'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         write(info_str,'(2X,A)') &
           'applied to determine the occupation numbers.'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         estimate_low_gap = .true.
      else if (((lumo_level-homo_level)*hartree) .lt.0.1) then ! 0.1 eV
         write(info_str,'(2X,A)') &
           'However, the gap value is rather small. Since we use a finite k-point grid,'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         write(info_str,'(2X,A)') &
           'this material is most likely metallic (in the sense that there are states'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         write(info_str,'(2X,A)') &
           'at or near the Fermi level). A DOS and/or band structure plot is required'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         write(info_str,'(2X,A)') &
           'for an unambiguous distinction.'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         estimate_low_gap = .true.
       else if (((lumo_level-homo_level)*hartree) .lt.0.2) then ! 0.2 eV
         write(info_str,'(2X,A)') &
           'The gap value is rather small. This may be a metal or a small-gap semiconductor.'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         write(info_str,'(2X,A)') &
            'Since we use a finite k-point grid, a DOS and/or band structure plot is required'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         write(info_str,'(2X,A)') &
           'for an unambiguous assignment.'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         estimate_low_gap = .true.
       else ! gap value above 0.2 is probably semiconducting or insulating.
         write(info_str,'(2X,A)') &
           'The gap value is above 0.2 eV. Unless you are using a very sparse k-point grid,'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         write(info_str,'(2X,A)') &
           'this system is most likely an insulator or a semiconductor.'
         call localorb_info(info_str, use_unit, "(A)", OL_norm)
         estimate_low_gap = .false.
       end if
    end if
end subroutine find_and_output_homo_lumo_gap


!****s* FHI-aims/output_energy_levels
!  NAME
!   output_energy_levels
!  SYNOPSIS
      subroutine output_energy_levels &
      ( KS_eigenvalue, occ_numbers, chemical_potential, chemical_potential_spin, conv )
!  PURPOSE
!  Subroutine output_energy_levels writes energy levels (with respect to the vacuum level)
!  of the current Kohn-Sham eigenstates.
!  Attention:
!  please only use it if((n_periodic > 0) .and. (use_dipole_correction .or. calculate_work_function))
!
!
!  USES
      use dimensions
      use runtime_choices
      use basis
      use localorb_io
      use constants
      use pbc_lists, only : k_point_list
      use physics, only: vacuum_level_upper, vacuum_level_lower
      implicit none
!  ARGUMENTS
      real*8 KS_eigenvalue (n_states, n_spin,n_k_points)
      real*8, dimension(n_states, n_spin,n_k_points) :: occ_numbers
      real*8 :: chemical_potential
      real*8, dimension(n_spin) :: chemical_potential_spin
      ! MR: add the option to know whether calculation reached convergence
      logical,intent(in) :: conv
!  INPUTS
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o occ_numbers -- occupation weights of eigenstates
!   o chemical_potential, chemical_potential_spin
!  OUTPUT
!    none
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
      character l_char
      character*150 :: info_str
      real*8, dimension(n_spin) :: homo_level, lumo_level
      real*8, dimension(n_spin) :: homo_occ, lumo_occ, i_kpt_homo, i_kpt_lumo
      real*8 :: midpoint
!  counters
      integer i_basis, i_state, i_fn
      integer i_state_2
      integer i_spin, i_spin2
      integer i_k_point
!  others
      character l_to_str
      character*300 tmp_str
      logical :: passit ! MR: this is just checking whether it is converged or not
!  begin work

!      passit=.false.
      passit=conv
! TODO MR: Set optional argument. for some reason present() was always returning true
! even when it was not present (conv was an optional argument). Probably just my inability to code.
!      write(use_unit,*) "present?", present(conv), conv
!      if (present(conv)) then
!          passit=conv
!      endif

      write(info_str,*)
      call localorb_info(info_str,use_unit,'(A)', OL_low)
      write(info_str,'(2X,A)') "Writing energy levels: "
      call localorb_info(info_str,use_unit,'(A)',OL_low)


      ! Determine HOMO and LUMO values (VBM and CBM in case of periodic systems) for spins
      ! safe initial values. If we break these, we have a problem somwehere else.
      homo_level = -10000000.
      lumo_level = 10000000.
      ! Define the correct "half occupation" (with or without spin)
      midpoint = spin_degeneracy/2.0d0

      i_kpt_homo = 0
      i_kpt_lumo = 0
      do i_k_point = 1, n_k_points, 1
        do i_spin = 1, n_spin, 1
          do i_state = 1, n_states, 1
            if (occ_numbers(i_state,i_spin,i_k_point).ge.midpoint) then
              ! check if homo
              ! "HOMO" also includes Fermi level ("ge" above)
              if (KS_eigenvalue(i_state,i_spin,i_k_point).gt.homo_level(i_spin)) then
                homo_level(i_spin) = KS_eigenvalue(i_state,i_spin,i_k_point)
                homo_occ(i_spin) = occ_numbers(i_state,i_spin,i_k_point)
                i_kpt_homo(i_spin) = i_k_point
              end if
            end if
            if (occ_numbers(i_state,i_spin,i_k_point).le.midpoint) then
              ! check if lumo
              ! LUMO must also include Fermi level ("le" above), else we may get nonsensical gaps (i.e., a gap in a molecule with half-occupied orbitals)
              if (KS_eigenvalue(i_state,i_spin,i_k_point).lt.lumo_level(i_spin)) then
                lumo_level(i_spin) = KS_eigenvalue(i_state,i_spin,i_k_point)
                lumo_occ(i_spin) = occ_numbers(i_state,i_spin,i_k_point)
                i_kpt_lumo(i_spin) = i_k_point
              end if
            end if
          enddo
        enddo
      enddo

      !write vacuum level upper/lower
      write(info_str,'(2X,A,F15.8,A)') &
                               '| Potential vacuum level, "upper" slab surface:', &
                               vacuum_level_upper*hartree, ' eV'
      call localorb_info(info_str,use_unit,'(A)',OL_low)

      write(info_str,'(2X,A,F15.8,A)') &
                               '| Potential vacuum level, "lower" slab surface:', &
                               vacuum_level_lower*hartree, ' eV'
      call localorb_info(info_str,use_unit,'(A)',OL_low)

      !write work function upper/lower for both spins
      if (n_spin==1) then

         write(info_str,'(2X,A,F15.8,A)') &
                               '| Work function ("upper" slab surface)        :', &
                               - (chemical_potential- vacuum_level_upper)*hartree, ' eV'
         call localorb_info(info_str,use_unit,'(A)',OL_low)

         write(info_str,'(2X,A,F15.8,A)') &
                               '| Work function ("lower" slab surface)        :', &
                               - (chemical_potential- vacuum_level_lower)*hartree, ' eV'
         call localorb_info(info_str,use_unit,'(A)',OL_low)

         !warning for systems with finite band gap (> 0.1 eV):
         if (i_kpt_homo(1) > 0 .and. i_kpt_lumo(1) > 0) then
            if ((lumo_level(1)-homo_level(1))*hartree> 0.1) then
               write(info_str,'(2X,A)') &
               '*** Warning: band gap may exist, check the suitability to define work function using chemical potential.'
               call localorb_info(info_str,use_unit,'(A)',OL_low)
            end if
         end if

          if (passit) then
            if (use_pimd_wrapper .and. ipi_work) then
              write(tmp_str,'(2X,A,F15.8,A)') &
                               ' | Vacuum level up:', &
                                   vacuum_level_upper*hartree, ' eV'
             comm_string = trim(comm_string) // trim(tmp_str)

             write(tmp_str,'(2X,A,F15.8,A)') &
                               ' | Vacuum level down:', &
                                 vacuum_level_lower*hartree, ' eV'
             comm_string = trim(comm_string) // trim(tmp_str)

             write(tmp_str,'(2X,A,F15.8,A)') &
                               ' | Work function up:', &
                      - (chemical_potential- vacuum_level_upper)*hartree, ' eV'
             comm_string = trim(comm_string) // trim(tmp_str)

             write(tmp_str,'(2X,A,F15.8,A)') &
                               ' | Work function down:', &
                      - (chemical_potential- vacuum_level_lower)*hartree, ' eV'
             comm_string = trim(comm_string) // trim(tmp_str)
             write(tmp_str,'(2X,A,F15.8,A)') &
                               ' | Chempot:', &
                      - (chemical_potential)*hartree, ' eV'
             comm_string = trim(comm_string) // trim(tmp_str)
            endif ! ipi
          endif ! conv true

          !write HOMO(VBM)/LUMO(CBM) with repect to upper vacuum level
        if (i_kpt_homo(1) > 0) then
            write(info_str,'(2X,A,F15.8,A)') &
                               '| VBM (reference: upper vacuum level)         :', &
                               - (homo_level- vacuum_level_upper)*hartree, ' eV'
            call localorb_info(info_str,use_unit,'(A)',OL_low)
        end if

        if (i_kpt_lumo(1) > 0) then
            write(info_str,'(2X,A,F15.8,A)') &
                               '| CBM (reference: upper vacuum level)         :', &
                               - (lumo_level- vacuum_level_upper)*hartree, ' eV'
            call localorb_info(info_str,use_unit,'(A)',OL_low)
        end if

          if (passit) then
           if (use_pimd_wrapper .and. ipi_work) then
             write(tmp_str,'(2X,A,F15.8,A)') &
                               ' | VBM:', &
                               - (homo_level- vacuum_level_upper)*hartree, ' eV'
             comm_string = trim(comm_string) // trim(tmp_str)

             write(tmp_str,'(2X,A,F15.8,A)') &
                               ' | CBM:', &
                               - (lumo_level- vacuum_level_upper)*hartree, ' eV'
             comm_string = trim(comm_string) // trim(tmp_str)
           endif ! ipi
          endif


      else if (n_spin==2) then
         write(info_str,*)
         call localorb_info(info_str,use_unit,'(A)', OL_low)
         write(info_str,'(2X,A)') "Spin-up energy levels: "
         call localorb_info(info_str,use_unit,'(A)',OL_low)

         do i_spin=1,n_spin,1

            if (i_spin==2) then
               write(info_str,*)
               call localorb_info(info_str,use_unit,'(A)', OL_low)
               write(info_str,'(2X,A)') "Spin-down energy levels: "
               call localorb_info(info_str,use_unit,'(A)',OL_low)
            end if

            write(info_str,'(2X,A,F15.8,A)') &
                               '| Work function-spin ("upper" slab surface)   :', &
                               - (chemical_potential_spin(i_spin)- vacuum_level_upper)*hartree, ' eV'
            call localorb_info(info_str,use_unit,'(A)',OL_low)

            write(info_str,'(2X,A,F15.8,A)') &
                               '| Work function-spin ("lower" slab surface)   :', &
                               - (chemical_potential_spin(i_spin)- vacuum_level_lower)*hartree, ' eV'
            call localorb_info(info_str,use_unit,'(A)',OL_low)

            if (passit) then
                if (use_pimd_wrapper .and. ipi_work) then
                  write(tmp_str,'(2X,A,F15.8,A)') &
                               ' | WF sp:', &
                                - (chemical_potential_spin(i_spin)- vacuum_level_upper)*hartree, ' eV'
                  comm_string = trim(comm_string) // trim(tmp_str)

                  write(tmp_str,'(2X,A,F15.8,A)') &
                               ' | WF sp:', &
                               - (chemical_potential_spin(i_spin)- vacuum_level_lower)*hartree, ' eV'
                  comm_string = trim(comm_string) // trim(tmp_str)
                endif ! ipi
            endif

            if (i_kpt_homo(i_spin) > 0 .and. i_kpt_lumo(i_spin) > 0) then
              !warning for systems with finite band gap (> 0.1 eV):
              if ((lumo_level(i_spin)-homo_level(i_spin))*hartree> 0.1) then
                 write(info_str,'(2X,A)') &
                   '*** Warning: band gap may exist, check the suitability to define work function using chemical potential.'
                   call localorb_info(info_str,use_unit,'(A)',OL_low)
              end if
            end if

            !write HOMO(VBM)/LUMO(CBM) with repect to upper vacuum level
            if (i_kpt_homo(i_spin) > 0) then
               write(info_str,'(2X,A,F15.8,A)') &
                               '| VBM-spin (reference: upper vacuum level)    :', &
                               - (homo_level(i_spin)- vacuum_level_upper)*hartree, ' eV'
               call localorb_info(info_str,use_unit,'(A)',OL_low)

               if (passit) then
                  if (use_pimd_wrapper .and. ipi_work) then
                     write(tmp_str,'(2X,A,F15.8,A)') &
                               ' | VBM sp:', &
                                 - (homo_level(i_spin)- vacuum_level_upper)*hartree, ' eV'
                     comm_string = trim(comm_string) // trim(tmp_str)
                  end if !ipi
               endif

            end if

            if (i_kpt_lumo(i_spin) > 0) then
               write(info_str,'(2X,A,F15.8,A)') &
                               '| CBM-spin (reference: upper vacuum level)    :', &
                               - (lumo_level(i_spin)- vacuum_level_upper)*hartree, ' eV'
               call localorb_info(info_str,use_unit,'(A)',OL_low)

               if (passit) then
                  if (use_pimd_wrapper .and. ipi_work) then
                     write(tmp_str,'(2X,A,F15.8,A)') &
                               ' | CBM sp:', &
                                 - (homo_level(i_spin)- vacuum_level_upper)*hartree, ' eV'
                     comm_string = trim(comm_string) // trim(tmp_str)
                  end if !ipi
               endif

            end if

         end do !i_spin

      end if

      return
    end subroutine output_energy_levels
!----------------------------------------------------------------------
!******
