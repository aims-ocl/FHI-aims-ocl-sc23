!****s* FHI-aims/aufbau_principle.f90
!  NAME
!    aufbau_principle
!  SYNOPSIS

    subroutine aufbau_principle

!  PURPOSE
!  Sets the occupation numbers of a single atom so that it complies with 
!  the aufbau principle
!
!  USES
!TODO: Check and cut out what is not needed

      use localorb_io
      use dimensions
      use timing
      use runtime_choices
      use physics
      use species_data
  use synchronize_mpi


      implicit none

!  ARGUMENTS


!  INPUTS
!    o none
!  OUTPUTS
!    o occupation numbers are updated accordingly
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
!    Release version, FHI-aims (2008).
!  SOURCE
!

      integer :: i_state
      integer :: i_electron
      integer :: i_n, i_l, i_m, l_max, i_spin
      integer :: n_el

      character*100 :: info_str
      character*16 :: func = 'aufbau_principle'

      ! begin work

       call localorb_info('Single atom found. Occupation will be taken from hardwired values')
!Check if this call is valid
       if (n_atoms.ne.1) call aims_stop('Aufbau principle initilisation only allowed for a single atom')
       if (n_periodic.ne.0) call aims_stop('Aufbau principle initialisation only allowed for non-periodic calculations')
       if (flag_moment.eqv..true.) call aims_stop('Aufbau principle initialisation only allowed for Hund rule')
       if (n_spin.ne.2.and.modulo(int(n_electrons),2).ne.0)  & 
          call aims_stop('Single atom with odd number of electrons reqiures open shell calculation.') 
       n_el = nint(n_electrons)
       if (abs(n_el - n_electrons) > 1d-13) call aims_stop('No aufbau principle with fractional total charges', func)


       !clear occupation numbers
       occ_numbers = 0.d0
!Do the aufbau principle
        if (n_spin.eq.2) then
          if (n_el.le.18 .or. &  !straight forward until Ar
           & (n_el.ge.30.and.n_el.le.36)   ) then  !between Zn and Kr
                 i_state =0
                 i_electron=n_el !i_electron - number of electrons not yet assigned
                 do i_n = 0,3,1 !main quantum number
                    do i_l = 0,i_n,1 !angular quantum number
                     do i_spin = 1,2,1
                       do i_m = -i_l,i_l,1 !magnetic quantum number
                           i_state = i_state+1
                           if (i_electron.ge.1) then
                                 occ_numbers(i_state,i_spin,1) = 1.0
                                  i_electron = i_electron-1
                           endif
                        enddo
                        if (i_spin==1) i_state = i_state-(2*i_l+1) !reset state before filling the shell with second spin
                      enddo
                    enddo
                  enddo
          else
                if (n_el.gt.18.and.n_el.le.36) occ_numbers(1:9,1:2,1)= 1.0 !Ar core
                if (n_el.gt.36.and.n_el.le.54) occ_numbers(1:18,1:2,1)=1.0 !Kr core
                select case (int(n_el))
                       case (19) !K
                           occ_numbers(10,1,1) = 1.0  !4s below 3d
                       case (20) !Ca
                           occ_numbers(10,1:2,1) = 1.0 !4s below 3d
                       case (21) !Sc
                           occ_numbers(10,1:2,1)= 1.0 !4s below 3d
                           occ_numbers(11,1,1) = 1.0
                       case (22) !Ti
                            !strange..4s below 3d for spin=2 only
                             occ_numbers(15,1,1) = 1.0 !4s1
                             occ_numbers(10,2,1) = 1.0 !4s2
                             occ_numbers(10:11,1,1) = 1.0 !4d2
                       case (23) !V
                            !same as Ti
                             occ_numbers(15,1,1) = 1.0 !4s1
                             occ_numbers(10,2,1) = 1.0 !4s2
                             occ_numbers(10:12,1,1) = 1.0 !4d3
                       case (24) !Cr
                            !same as Ti and V 
                             occ_numbers(15,1,1) = 1.0 !4s1
                             occ_numbers(10:14,1,1) = 1.0 !4d3
                      case (25) !Mn
                             occ_numbers(15,1:2,1) = 1.0 !4s2
                             occ_numbers(10:14,1,1) = 1.0 !4d5
                      case (26) !Fe
                             occ_numbers(15,1:2,1) = 1.0 !4s2
                             occ_numbers(10:14,1,1) = 1.0 
                             occ_numbers(10,2,1) = 1.0
                      case (27) !Co
                             occ_numbers(15,1:2,1) = 1.0 !4s2
                             occ_numbers(10:14,1,1) = 1.0 
                             occ_numbers(10:11,2,1) = 1.0
                      case (28) !Ni
                             occ_numbers(15,1:2,1) = 1.0 !4s2
                             occ_numbers(10:14,1,1) = 1.0
                             occ_numbers(10:12,2,1) = 1.0
                      case (29) !Cu
                             occ_numbers(15,1,1) = 1.0 !4s1
                             occ_numbers(10:14,1:2,1) = 1.0
                      case (37) !Rb
                            occ_numbers(19,1,1) = 1.0 !5s below 4d
                      case (38) !Sr
                            occ_numbers(19,1:2,1) = 1.0 !5s2 
                      case (39) !Y
                            occ_numbers(19,1:2,1) = 1.0 !5s2
                            occ_numbers(20,1,1)= 1.0
                      case (40) !Zr
                            occ_numbers(19,1:2,1) = 1.0 !5s below 4d
                            occ_numbers(20:21,1,1) = 1.0
                      case (41) !Nb
                            occ_numbers(19,1,1) = 1.0 !5s below, but only s1 occ
                            occ_numbers(20:23,1,1) = 1.0
                      case (42) !Mo
                            occ_numbers(19:24,1,1) = 1.0 !4d54s1
                      case (43) !Tc !TODO: CHECK
                             occ_numbers(19:24,1,1) = 1.0
                             occ_numbers(20,2,1) = 1.0
                      case (44) !Ru
                            occ_numbers(24,1,1) = 1.0 !5s1
                            occ_numbers(19:23,1,1) = 1.0
                            occ_numbers(20:21,2,1) = 1.0
                       case (45) !Rh
                            occ_numbers(24,1,1) = 1.0
                            occ_numbers(19:23,1,1) = 1.0
                            occ_numbers(20:22,2,1) = 1.0
                       case (46) !Pd
                            occ_numbers(19:23,1:2,1) = 1.0
                       case (47) !Ag
                            occ_numbers(19:23,1:2,1) = 1.0
                            occ_numbers(24,1,1) = 1.0
                       case (48) !Cd
                            occ_numbers(19:23,1:2,1) = 1.0
                            occ_numbers(24,1:2,1) = 1.0
                       case (49) !In
                            occ_numbers(19:23,1:2,1) = 1.0
                            occ_numbers(24,1:2,1) = 1.0
                            occ_numbers(25,1,1) = 1.0 
                       case (50) !Sn
                            occ_numbers(19:23,1:2,1) = 1.0
                            occ_numbers(24,1:2,1) = 1.0
                            occ_numbers(25:26,1,1) = 1.0 
                       case (51) !Sb
                            occ_numbers(19:23,1:2,1) = 1.0
                            occ_numbers(24,1:2,1) = 1.0
                            occ_numbers(25:27,1,1) = 1.0 
                       case (52) !Te
                            occ_numbers(19:23,1:2,1) = 1.0
                            occ_numbers(24,1:2,1) = 1.0
                            occ_numbers(25:27,1,1) = 1.0
                            occ_numbers(25,2,1) = 1.0
                       case (53) !I
                            occ_numbers(19:23,1:2,1) = 1.0
                            occ_numbers(24,1:2,1) = 1.0
                            occ_numbers(25:27,1,1) = 1.0
                            occ_numbers(25:26,2,1) = 1.0
                       case (54) !Xe
                            occ_numbers(19:23,1:2,1) = 1.0
                            occ_numbers(24,1:2,1) = 1.0
                            occ_numbers(25:27,1:2,1) = 1.0
                      case default
                         call aims_stop('Atom not yet specified for aufbau principle initialisation. Sorry for the inconvenience.')
                end select

          endif
       elseif(n_spin.eq.2)  then
          !Check if applicable:
         if(n_el.eq.2 .or. n_el.eq.8 .or. n_el.eq.18 .or. n_el.eq.36 .or. n_el.eq.54 .or. n_el.eq.86 .or. & !Noble gas
         &  n_el.eq.4 .or. n_el.eq.12. .or. n_el.eq.20 .or. n_el.eq.38 .or. n_el.eq.56 .or. n_el.eq.88 .or. & !Earth alkali
         &  n_el.eq.30 .or. n_el.eq.48 .or. n_el.eq.80 .or. n_el.eq.112 .or. & !group 12 elements
         &  n_el.eq.70 .or. n_el.eq.102) then !Yb and Nb
            occ_numbers(1:(n_el/2),1,1)=2.0
          else 
           write(info_str,*) 'Closed shell calculation for ', species_name(1), ' atom not possible using aufbau principle.'
           call aims_stop(info_str)
         endif       
       endif
         
! Nonsense-check:
          if (abs(n_el - sum(occ_numbers)) > 1d-13) call aims_stop('Electrons not completely distributed', func)

                   
end subroutine aufbau_principle
