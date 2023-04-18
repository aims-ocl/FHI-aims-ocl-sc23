!****s* FHI-aims/read_species_data
!  NAME
!   read_species_data
!  SYNOPSIS

subroutine read_species_data &
   ( i_species, flag_radial, flag_angular, &
      flag_angular_min, flag_angular_acc, flag_cut_free_atom, &
      flag_prodbas_acc, flag_max_n_prodbas, flag_max_l_prodbas &
   )

!  PURPOSE
!  Subroutine read_species_data belongs to read_control.f;
!  reads electronic information and basis states for a species
!
!  USES

   use dimensions
   use grids
   use applicable_citations
   use species_data
   use mpi_tasks
   use localorb_io
   use constants
   use read_fixed_grid
   use runtime_choices
   use vdw_correction
   use pseudodata

   implicit none

!  ARGUMENTS

   integer i_species
   logical flag_radial
   logical flag_angular
   logical flag_angular_min
   logical flag_angular_acc
   logical, intent(OUT) :: flag_prodbas_acc
   logical, intent(OUT) :: flag_max_n_prodbas
   logical, intent(OUT) :: flag_max_l_prodbas

!  INPUTS
!  o i_species -- species index
!
!  OUTPUT
!  o flag_radial -- Was the radial grid read in?
!  o flag_angular-- Was the angular grid read in?
!  o flag_angular_min-- Was the angular grid minimum grid read in?
!  o flag_angular_acc -- Was the angular grid accuracy read in?
!  o flag_prodbas_acc -- Was the prodbas accuracy read in?
!  o flag_max_n_prodbas -- Was the maximum number of one-particle radial parts
!                          for prodbas read in?
!  o flag_max_l_prodbas -- Was the maximum angular momentum for prodbas read in?
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



!  local variables

!  desc_str : line identifier in control.in
!  i_code   : return status flag from read operation.
!  flag_eos      : end-of-species flag for read operation
!  flag_eof      : end-of-file flag
!  n_shell  : main quantum number of a shell
!  l_shell  : angular momentum quantum number of a shell
!  l_shell_str : s, p, d, f ... well, angular momentum identifiers.
!  occ_shell: occupation number of a shell
!  test_occ : checksum for total valence occupation
!  i_types : Counts number of requested basis types
!            (atomic, ionic, cutoff, ...); triggers warning if max_basis is
!            underdimensioned.

   character*132 inputline
   character*20 desc_str
   integer i_code,linecount
   logical flag_eos
   logical flag_eof
   integer n_shell, l_shell, kappa
   character l_shell_str
   real*8 occ_shell
   real*8 test_occ
   integer i_types
   real*8 radius
   real*8 z_eff
   integer cartesian_l
   integer min_gaussian_l
   integer n_contracted
   real*8 alpha(n_max_contracted)
   real*8 coeff(n_max_contracted)
   real*8 aux_alpha(n_max_aux_contracted)
   real*8 aux_coeff(n_max_aux_contracted)

   integer n_basis_sp
   integer n_aux_basis_sp
   integer new_int
   logical info_log

   integer radial_multiplier
   

!  flags to see which items were set explicitly

   logical flag_z
   logical flag_valence
   logical flag_ionic(0:l_wave_max, n_wave_max)
   logical flag_acc
   logical flag_inmax
   logical flag_core_pot
   logical flag_cut_pot
   logical flag_cutoff_type
   logical flag_cut_free_atom
   logical flag_cut_core
   logical flag_log
   logical flag_l_hartree
   logical flag_include_min_basis
   logical flag_radial_multiplier
   logical flag_specify_grids
   logical flag_pure_gauss
   logical :: flag_element

   logical auto_ionic(n_max_ind_fns)
   logical auto_conf(n_max_ind_fns)

   character*150 :: info_str

!  counters

   integer i_l, i_shell
   integer i_contracted
   integer :: i_basis
   integer :: i_shell_plus_u

   integer :: i_pp_species


!  functions

   integer str_to_l
   character l_to_str
   character(*), parameter :: func = 'read_species_data'


   
!  begin work



   if (myid.eq.0) then
      write(use_unit,*)
      write (use_unit,'(2X,A,A,A)') &
            "Reading configuration options for species ", &
            species_name(i_species), "."
   end if

   flag_eos = .false.
   flag_eof = .false.

!  initialize

   flag_z = .false.
   flag_element = .false.
   flag_acc = .false.
   flag_prodbas_acc = .false.
   flag_max_n_prodbas = .false.
   flag_max_l_prodbas = .false.
   flag_inmax = .false.
   flag_core_pot = .false.
   flag_cut_pot = .false.
   flag_cutoff_type = .false.
   flag_cut_free_atom = .false.
   flag_cut_core = .false.
   flag_radial = .false.
   flag_angular = .false.
   flag_angular_min = .false.
   flag_angular_acc = .false.
   flag_log = .false.
   flag_l_hartree = .false.
   flag_include_min_basis = .false.
   flag_radial_multiplier = .false.
   flag_specify_grids = .false.
   flag_pure_gauss = .true.

   l_shell_max(i_species) = 0

   n_atomic(i_species) = 0
   n_ionic(i_species) = 0
   n_conf (i_species) = 0
   n_hydro(i_species) = 0
!   n_hydro_aux(i_species) = 0
   n_gaussian(i_species) = 0
   n_aux_gaussian(i_species) = 0

   n_basis_sp = 0

   ! set the default just in case someone misses it later.
   cutoff_type(i_species) = 3

!     set default minimum number of angular grid points per radial shell directly,
!     rather than checking by flag angular_min.
   angular_min(i_species) = 0

   core_n_max(:,i_species) = -1
   core_fn(i_species,:) = .false.

   flag_valence = .false.

   do i_l = 0, l_wave_max, 1
      valence_n_max(i_l,i_species) = 0
      if (use_ionic) then
         ion_n_max(i_l,i_species) = 0
      end if
      do i_shell = i_l+1, n_wave_max,1
         valence_occ(i_shell, i_l, i_species) = 0.
         flag_ionic (i_l, i_shell) = .false.
      enddo
   enddo

   ! basis-dependent cutoff potential defaults:
   basis_dep_cutoff_thresh(i_species)    = basis_dep_cutoff_default
   cut_atomic_basis_functions(i_species) = .false.

   i_shell_plus_u = 0

!  read species information

   linecount = 0

   lineloop: do

      read(7,'(A)',iostat=i_code) inputline
      if(i_code<0) then
         backspace(7)
        exit lineloop        ! end of file reached
      end if
      if(i_code>0) then
         call aims_stop ("Error reading file 'control.in'", "read_species_data")
      end if

      linecount = linecount + 1

      read(inputline,*,iostat=i_code) desc_str
      if(i_code/=0) cycle lineloop             ! skip empty line

      if (desc_str(1:1).eq.'#') cycle lineloop ! skip comment

!DB: adding pseudpcore infrastructure here
      if (desc_str.eq."pseudo") then
         read(inputline,*,end=88,err=99) desc_str, pp_path_string(i_species)
         species_pseudoized(i_species) = .true. 

         n_pp_species_counter_read_in = n_pp_species_counter_read_in +1
         i_pp_species = n_pp_species_counter_read_in


     else  if (desc_str.eq."nonlinear_core") then
         read(inputline,*,end=88,err=99) desc_str, species_nonlinear_core(i_species)
      
     else  if (desc_str.eq."nucleus") then
         read(inputline,*,end=88,err=99) desc_str, species_z(i_species)

         if (species_z(i_species).le.0) then
           if (.not.override_warning_negative_nucleus) then
            if (myid.eq.0) then
               write(use_unit,*) "* Illegal nuclear charge set for species ", &
                  species_name(i_species), "."
            end if
           call aims_stop()
           endif
         end if

         flag_z = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,F8.4)') &
               "| Found nuclear charge : ", species_z(i_species)
         end if

      else if (desc_str.eq.'mass') then
         read(inputline,*,end=88,err=99) desc_str, species_m(i_species)
         if (myid.eq.0) then
            write(use_unit,*) ' | Found atomic mass : ',species_m(i_species), &
                  'amu'
         end if
      
      else if (desc_str == 'element') then
         read (inputline, *, end=88, err=99) desc_str, species_element(i_species)
         if (myid == 0) then
             write (use_unit,*) ' | Found element : ', species_element(i_species)
         end if
         flag_element = .true.

      else if (desc_str.eq.'hirshfeld_param') then
         read(inputline,*,end=88,err=99) desc_str, vdw_hirshfeld_C6(i_species), &
               vdw_hirshfeld_alpha(i_species), vdw_hirshfeld_R0(i_species)
         vdw_hirshfeld_data_external(i_species) = .true.
         write(info_str,'(2X,A)') "| Found explicit parameters for TS-vdw correction: "
         call localorb_info(info_str)
         write(info_str,'(2X,A,F10.6)') "|    C6    = ",vdw_hirshfeld_C6(i_species)
         call localorb_info(info_str)
         write(info_str,'(2X,A,F10.6)') "|    alpha = ",vdw_hirshfeld_alpha(i_species)
         call localorb_info(info_str)
         write(info_str,'(2X,A,F10.6)') "|    R0    = ",vdw_hirshfeld_R0(i_species)
         call localorb_info(info_str)

      else if (desc_str.eq."core") then
         read(inputline,*,end=88,err=99) desc_str, n_shell, l_shell_str

         l_shell = str_to_l(l_shell_str)

         core_n_max (l_shell,i_species) = n_shell

      else if (desc_str.eq."valence") then
         read(inputline,*,end=88,err=99) desc_str, n_shell, l_shell_str, occ_shell

         l_shell = str_to_l(l_shell_str)

         if (l_shell.gt.l_shell_max(i_species)) then
            l_shell_max(i_species) = l_shell
         end if

         valence_n_max (l_shell,i_species) = n_shell
         do i_shell = l_shell+1,n_shell-1,1
            valence_occ(i_shell, l_shell, i_species) = &
            2 * ( 2*l_shell + 1 )
         enddo
         valence_occ   (n_shell, l_shell, i_species) = occ_shell

         flag_valence = .true.

         if (allocated(n_fc_shell))then
            if ( n_fc_shell(i_species).lt.n_shell) then
               n_fc_shell(i_species)=n_shell
            endif
         endif

         if (myid.eq.0) then
            write(use_unit,'(2X,A,I2,1X,A,1X,F7.3)') &
               "| Found free-atom valence shell : ", &
               n_shell, l_shell_str, occ_shell
         end if

      else if (desc_str.eq."ion_occ") then

         if (use_ionic) then

            read(inputline,*,end=88,err=99) desc_str, n_shell, l_shell_str, occ_shell

            l_shell = str_to_l(l_shell_str)

            if (l_shell.gt.l_shell_max(i_species)) then
               l_shell_max(i_species) = l_shell
            end if

!           verify validity of main quantum number
            if (n_shell.le.l_shell) then

               if (myid.eq.0) then
                  write(use_unit,*) "* Ionic shell ",n_shell, l_shell_str, ":", &
                        " n is too low for this l."
               end if

              call aims_stop()
            end if

            ion_n_max(l_shell,i_species) = &
            max( n_shell,ion_n_max(l_shell,i_species) )

            ion_occ(n_shell,l_shell,i_species) = occ_shell

            flag_ionic(l_shell, n_shell) = .true.

            if (myid.eq.0) then
               write(use_unit,'(2X,A,I2,1X,A,1X,F7.3)') &
                     "| Found free-ion valence shell : ", &
                     n_shell, l_shell_str, occ_shell
            end if

         else

            if (myid.eq.0) then
               write(use_unit,'(2X,A)') &
                  "| No ionic wave fns used. Skipping ion_occ."
            end if

         end if

      else if (desc_str.eq."ionic") then

         read(inputline,*,end=88,err=99) desc_str, n_shell, l_shell_str, desc_str


         l_shell = str_to_l(l_shell_str)
         if (l_shell.gt.l_shell_max(i_species)) then
         l_shell_max(i_species) = l_shell
         end if

!         verify validity of main quantum number
         if (n_shell.le.l_shell) then
            if (myid.eq.0) then
               write(use_unit,*) "* Ionic shell ",n_shell, l_shell_str, ":", &
                  " n is too low for this l."
            end if
           call aims_stop()
         end if

         if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks ) then
            if (l_shell.gt.0)then
                n_ionic(i_species) = n_ionic(i_species)+1
                ionic_in_large_basis(i_species,n_ionic(i_species)) = .false.
                ionic_n(i_species,n_ionic(i_species)) = n_shell
                ionic_l(i_species,n_ionic(i_species)) = l_shell
                ionic_kappa(i_species,n_ionic(i_species)) = l_shell

                n_ionic(i_species) = n_ionic(i_species)+1
                ionic_in_large_basis(i_species,n_ionic(i_species)) = .false.
                ionic_n(i_species,n_ionic(i_species)) = n_shell
                ionic_l(i_species,n_ionic(i_species)) = l_shell
                ionic_kappa(i_species,n_ionic(i_species)) = -l_shell-1
            else
                n_ionic(i_species) = n_ionic(i_species)+1
                ionic_in_large_basis(i_species,n_ionic(i_species)) = .false.
                ionic_n(i_species,n_ionic(i_species)) = n_shell
                ionic_l(i_species,n_ionic(i_species)) = l_shell
                ionic_kappa(i_species,n_ionic(i_species)) = -l_shell-1
            endif
            n_basis_sp = n_basis_sp + (2*l_shell+1) ! What's this? Check further.
         else
            n_ionic(i_species) = n_ionic(i_species)+1
            ionic_in_large_basis(i_species,n_ionic(i_species)) = .false.
            ionic_n (i_species, n_ionic(i_species)) = n_shell
            ionic_l (i_species, n_ionic(i_species)) = l_shell
            n_basis_sp = n_basis_sp + (2*l_shell+1)
         endif

         if (desc_str.eq."auto") then
            auto_ionic(n_ionic(i_species)) = .true.
            if((flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks).and.l_shell.gt.0)then
              auto_ionic(n_ionic(i_species)-1) = .true.
            endif

            if (myid.eq.0) then
               write(use_unit,'(2X,A,I2,1X,A,1X,A)') &
                  "| Found ionic basis function : ", n_shell, l_shell_str, &
                  ", default cutoff radius."
            end if

         else
            auto_ionic(n_ionic(i_species)) = .false.

            read(inputline,*,end=88,err=99) desc_str, n_shell, l_shell_str, radius

            ionic_rad (i_species, n_ionic(i_species)) = radius

            if((flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks).and.l_shell.gt.0)then
              auto_ionic(n_ionic(i_species)-1) = .false.
              ionic_rad(i_species,n_ionic(i_species)-1) = radius
            endif

            if (myid.eq.0) then
               write(use_unit,'(2X,A,I2,1X,A,1X,F7.3)') &
                  "| Found ionic basis function : ", &
                  n_shell, l_shell_str, radius
            end if

         end if

      else if (desc_str.eq."confined") then

         read(inputline,*,end=88,err=99) desc_str, n_shell, l_shell_str, desc_str

         l_shell = str_to_l(l_shell_str)
         if (l_shell.gt.l_shell_max(i_species)) then
            l_shell_max(i_species) = l_shell
         end if

!         verify validity of main quantum number
         if (n_shell.le.l_shell) then
            if (myid.eq.0) then
               write(use_unit,*) "* Confined shell ",n_shell, &
                  l_shell_str, ":", &
                  " n is too low for this l."
            end if
           call aims_stop()
         end if

         if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks ) then
            if (l_shell.gt.0)then
                n_conf(i_species) = n_conf(i_species)+1
                conf_n(i_species,n_conf(i_species)) = n_shell
                conf_l(i_species,n_conf(i_species)) = l_shell
                conf_kappa(i_species,n_conf(i_species)) = l_shell

                n_conf(i_species) = n_conf(i_species)+1
                conf_n(i_species,n_conf(i_species)) = n_shell
                conf_l(i_species,n_conf(i_species)) = l_shell
                conf_kappa(i_species,n_conf(i_species)) = -l_shell-1
            else
                n_conf(i_species) = n_conf(i_species)+1
                conf_n(i_species,n_conf(i_species)) = n_shell
                conf_l(i_species,n_conf(i_species)) = l_shell
                conf_kappa(i_species,n_conf(i_species)) = -l_shell-1
            endif
            n_basis_sp = n_basis_sp + (2*l_shell+1) ! What's this? Check further.
         else
            n_conf(i_species) = n_conf(i_species)+1
            conf_n(i_species,n_conf(i_species)) = n_shell
            conf_l(i_species,n_conf(i_species)) = l_shell
            n_basis_sp = n_basis_sp + (2*l_shell+1)
         endif

         if (desc_str.eq."auto") then
            auto_conf(n_conf(i_species)) = .true.

            if (myid.eq.0) then
               write(use_unit,'(2X,A,I2,1X,A,1X,A)') &
                  "| Found confined basis function : ", n_shell, l_shell_str, &
                  ", default cutoff radius."
            end if

         else
            auto_conf(n_conf(i_species)) = .false.

            read(inputline,*,end=88,err=99) desc_str, n_shell, l_shell_str, radius

            conf_rad(i_species,n_conf(i_species)) = radius

            if (myid.eq.0) then
               write(use_unit,'(2X,A,I2,1X,A,1X,F7.3)') &
                  "| Found confined basis function : ", &
                  n_shell, l_shell_str, radius
            end if
         end if

      else if (desc_str.eq."hydro") then

         read(inputline,*,end=88,err=99) desc_str, n_shell, l_shell_str, z_eff

         l_shell = str_to_l(l_shell_str)

         if (l_shell.gt.l_shell_max(i_species)) then
            l_shell_max(i_species) = l_shell
         end if

!         verify validity of main quantum number
         if (n_shell.le.l_shell) then

            if (myid.eq.0) then
               write(use_unit,*) "* Hydrogenic shell ", &
                     n_shell,l_shell_str,":", &
                     " n is too low for this l."
            end if

           call aims_stop()
         end if

         if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks ) then
            if (l_shell.gt.0)then
                n_hydro(i_species) = n_hydro(i_species)+1
                hydro_n(i_species,n_hydro(i_species)) = n_shell
                hydro_l(i_species,n_hydro(i_species)) = l_shell
                hydro_kappa(i_species,n_hydro(i_species)) = l_shell
                hydro_scale(i_species,n_hydro(i_species)) = z_eff
                hydro_in_large_basis(i_species,n_hydro(i_species)) = .false.

                n_hydro(i_species) = n_hydro(i_species)+1
                hydro_n(i_species,n_hydro(i_species)) = n_shell
                hydro_l(i_species,n_hydro(i_species)) = l_shell
                hydro_kappa(i_species,n_hydro(i_species)) = -l_shell-1
                hydro_scale(i_species,n_hydro(i_species)) = z_eff
                hydro_in_large_basis(i_species,n_hydro(i_species)) = .false.
            else
                n_hydro(i_species) = n_hydro(i_species)+1
                hydro_n(i_species,n_hydro(i_species)) = n_shell
                hydro_l(i_species,n_hydro(i_species)) = l_shell
                hydro_kappa(i_species,n_hydro(i_species)) = -l_shell-1
                hydro_scale(i_species,n_hydro(i_species)) = z_eff
                hydro_in_large_basis(i_species,n_hydro(i_species)) = .false.
            endif
            n_basis_sp = n_basis_sp + (2*l_shell+1) ! What's this? Check further.
         else
            n_hydro(i_species) = n_hydro(i_species)+1
            hydro_n(i_species,n_hydro(i_species)) = n_shell
            hydro_l(i_species,n_hydro(i_species)) = l_shell
            hydro_scale(i_species,n_hydro(i_species)) = z_eff
            hydro_in_large_basis(i_species,n_hydro(i_species)) = .false.

            n_basis_sp = n_basis_sp + (2*l_shell+1)
         endif

         if (myid.eq.0) then
            write(use_unit,'(2X,A,I2,1X,A,1X,F7.3)') &
                  "| Found hydrogenic basis function : ", &
                  n_shell, l_shell_str, z_eff
         end if

      else if (desc_str.eq."gaussian") then

!         read entire description of possibly contracted Gaussian radial function
         read(inputline,*,end=88,err=99) desc_str, cartesian_l, n_contracted

         if (cartesian_l.gt.l_shell_max(i_species)) then
            l_shell_max(i_species) = cartesian_l
         end if

         if (n_contracted.lt.1) then

            if (myid.eq.0) then
               write(use_unit,'(1X,A,A,I10,A)') &
               "* Error: Illegal number of elementary Gaussians in ", &
               "Gaussian basis function: n = ", n_contracted, "."
               write(use_unit,'(1X,A)') &
               "* Please check and correct your control.in file."
            end if

           call aims_stop()
         end if

         if (n_contracted.eq.1) then

            coeff(1) = 1.d0

            read(inputline,*,end=88,err=99) desc_str, cartesian_l, n_contracted, &
                     alpha(1)

            if (alpha(1).le.0.d0) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A,1X,F10.5,1X,A)') &
                     "* Error: Found illegal Gaussian fn. ", &
                     "exponent alpha =", &
                     alpha(1), "."
                  write(use_unit,'(1X,A,A)') &
                     "* Please check and correct your ", &
                     "control.in file."
               end if

              call aims_stop()
            end if

            if (myid.eq.0) then
               write(use_unit,'(2X,A,A,I2,1X,E12.6)') &
                  "| Found primitive cartesian ", &
                  "Gaussian basis function : ", &
                  cartesian_l, alpha(1)
            end if

         else

            do i_contracted = 1, n_contracted, 1
               read(7,'(A)') inputline
               if(i_code<0)then
                  if (myid.eq.0) then
                     write(use_unit,*) "Error reading specied data from file 'control.in':"
                     write(use_unit,*) "EOF reached within 'gaussian' block"
                  endif
                  call aims_stop_coll('', func)
               endif
               if(i_code>0)then
                  call aims_stop_coll("Unknown error reading file species data from 'control.in'. (within gaussian block)", func)
               endif

               read(inputline,*,end=88,err=99) alpha(i_contracted), coeff(i_contracted)

               if (alpha(i_contracted).le.0.d0) then

                  if (myid.eq.0) then
                     write(use_unit,'(1X,A,A,1X,F10.5,1X,A)') &
                           "* Error: Found illegal Gaussian ", &
                           "fn. exponent alpha =", &
                           alpha(i_contracted), "."
                     write(use_unit,'(1X,A,A)') &
                           "* Please check and correct your ", &
                           "control.in file."
                  end if

                 call aims_stop()
               end if

            enddo

            if (myid.eq.0) then
               write(use_unit,'(2X,A,A,I2,1X,A,I3,A)') &
                  "| Found contracted cartesian Gaussian basis function : ", &
                  " L =", cartesian_l, ", ", n_contracted, &
                  " elementary Gaussians:"
               do i_contracted = 1, n_contracted, 1
                  write(use_unit,'(2X,A,E12.6,A,E12.6)') &
                     "|   alpha = ", &
                     alpha( i_contracted ), &
                     " weight = ", &
                     coeff( i_contracted )
               enddo
            end if

         end if

!         treat every single angular momentum instance of the same cartesian
!         Gaussian as a different function

!         This is mildly inefficient in terms of memory, but makes the code more
!         consistent with the other basis function types ...

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A)') &
                  "|   In terms of angular momentum, ", &
                  "this radial function adds: "
         end if

         if (.not.use_sph_gaussian) then
            do l_shell = cartesian_l, 0, -2

               n_gaussian(i_species) = n_gaussian(i_species) + 1

               gaussian_n(i_species,n_gaussian(i_species)) = cartesian_l
               gaussian_l(i_species,n_gaussian(i_species)) = l_shell
               gaussian_n_contr (i_species,n_gaussian(i_species)) = &
               n_contracted

               do i_contracted = 1, n_contracted, 1
                  gaussian_alpha( i_species, n_gaussian(i_species), &
                                 i_contracted ) = alpha(i_contracted)
                  gaussian_coeff( i_species, n_gaussian(i_species), &
                                 i_contracted ) = coeff(i_contracted)
               enddo

               n_basis_sp = n_basis_sp + (2*l_shell+1)

               l_shell_str = l_to_str(l_shell)
               if (l_shell.eq.0) then

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,I3,1X,A,A)') &
                        "|   ", (2*l_shell+1),l_shell_str, &
                        "-type basis function"
                  end if

               else

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,I3,1X,A,A)') &
                        "|   ", (2*l_shell+1),l_shell_str, &
                        "-type basis functions"
                  end if

               end if

            enddo

         else

            l_shell = cartesian_l
!         do l_shell = cartesian_l, 0 ,-3
!            if (use_sph_gaussian.and.l_shell.ne.cartesian_l) cycle
            n_gaussian(i_species) = n_gaussian(i_species) +1

            gaussian_n(i_species,n_gaussian(i_species)) = cartesian_l
!            if (.not.flag_pure_gauss ) then
            gaussian_l(i_species,n_gaussian(i_species)) = l_shell
!             else
!            gaussian_l(i_species,n_gaussian(i_species)) = 0.d0
!            endif
            gaussian_n_contr (i_species,n_gaussian(i_species)) = &
               n_contracted

            do i_contracted = 1, n_contracted, 1
               gaussian_alpha( i_species, n_gaussian(i_species), &
                              i_contracted ) = alpha(i_contracted)
               gaussian_coeff( i_species, n_gaussian(i_species), &
                              i_contracted ) = coeff(i_contracted)
            enddo

            if  (l_shell.gt.2) then
               n_basis_sp = n_basis_sp + (2*l_shell+1) + l_shell-1
            else
               n_basis_sp = n_basis_sp + (2*l_shell+1)
            endif

            l_shell_str = l_to_str(l_shell)
            if (l_shell.eq.0) then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,I3,1X,A,A)') &
                     "|   ", (2*l_shell+1),l_shell_str, &
                     "-type basis function"
               end if

            else

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,I3,1X,A,A)') &
                     "|   ", (2*l_shell+1),l_shell_str, &
                     "-type basis functions"
               end if

            end if
!           enddo
         endif

      else if (desc_str == 'sto') then
         read(inputline,*,end=88,err=99) desc_str, n_shell, l_shell, z_eff
         l_shell_max(i_species) = max(l_shell_max(i_species), l_shell)

         if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks ) then
            if (l_shell.gt.0)then
                n_sto(i_species) = n_sto(i_species)+1
                sto_n(i_species,n_sto(i_species)) = n_shell
                sto_l(i_species,n_sto(i_species)) = l_shell
                sto_k(i_species,n_sto(i_species)) = l_shell
                sto_zeta(i_species,n_sto(i_species)) = z_eff

                n_sto(i_species) = n_sto(i_species)+1
                sto_n(i_species,n_sto(i_species)) = n_shell
                sto_l(i_species,n_sto(i_species)) = l_shell
                sto_k(i_species,n_sto(i_species)) = -l_shell-1
                sto_zeta(i_species,n_sto(i_species)) = z_eff
            else
                n_sto(i_species) = n_sto(i_species)+1
                sto_n(i_species,n_sto(i_species)) = n_shell
                sto_l(i_species,n_sto(i_species)) = l_shell
                sto_k(i_species,n_sto(i_species)) = -l_shell-1
                sto_zeta(i_species,n_sto(i_species)) = z_eff
            endif
            n_basis_sp = n_basis_sp + (2*l_shell+1)
         else
            n_sto(i_species) = n_sto(i_species)+1
            sto_n(i_species,n_sto(i_species)) = n_shell
            sto_l(i_species,n_sto(i_species)) = l_shell
            sto_k(i_species,n_sto(i_species)) = -l_shell-1
            sto_zeta(i_species,n_sto(i_species)) = z_eff

            n_basis_sp = n_basis_sp + (2*l_shell+1)
         endif

      else if (desc_str.eq."pure_gauss") then

         read(inputline,*,end=88,err=99) desc_str, pure_gaussian(i_species)

         flag_pure_gauss = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,G5.2)') &
                  "| Found request to include pure gaussian fns. : ", &
                  pure_gaussian(i_species)
         end if

!            gaussian_l(:,:) = 0.d0


      else if (desc_str.eq."aux_gaussian" ) then

!         read entire description of possibly contracted Gaussian radial function
         read(inputline,*,end=88,err=99) desc_str, cartesian_l, n_contracted

         if (n_contracted.lt.1) then

            if (myid.eq.0) then
               write(use_unit,'(1X,A,A,I10,A)') &
               "* Error: Illegal number of elementary Gaussians in ", &
               "auxiliary Gaussian basis function: n = ", &
               n_contracted, "."
               write(use_unit,'(1X,A)') &
                     "* Please check and correct your control.in file."
            end if

           call aims_stop()
         end if

         if (n_contracted.eq.1) then

            aux_coeff(1) = 1.d0

            read(inputline,*,end=88,err=99) desc_str, cartesian_l, n_contracted, &
                     aux_alpha(1)

            if (aux_alpha(1).le.0.d0) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A,1X,F10.5,1X,A)') &
                     "* Error: Found illegal auxiliary Gaussian fn. ", &
                     "exponent alpha =", &
                     aux_alpha(1), "."
                  write(use_unit,'(1X,A,A)') &
                     "* Please check and correct your ", &
                     "control.in file."
               end if

              call aims_stop()
            end if

            if (myid.eq.0) then
               write(use_unit,'(2X,A,A,I2,1X,E12.6)') &
                  "| Found primitive auxiliary cartesian ", &
                  "Gaussian basis function : ", &
                  cartesian_l, aux_alpha(1)
            end if

         else
            write(use_unit,*) "Contracted auxiliary gaussian basis is not ", &
                        "allowed, please use elementary ones."
           call aims_stop()
         endif


!          if (myid.eq.0) then
!             write(use_unit,'(2X,A,A)')
!     +            "|   In terms of angular momentum, ",
!     +            "this radial function adds: "
!          end if

!          min_gaussian_l = 0
         min_gaussian_l = cartesian_l

         do l_shell = cartesian_l, min_gaussian_l, -2

            n_aux_gaussian(i_species) = n_aux_gaussian(i_species) + 1

            aux_gaussian_n(i_species,n_aux_gaussian(i_species)) &
                           = cartesian_l
            aux_gaussian_l(i_species,n_aux_gaussian(i_species)) &
                           = l_shell
            aux_gaussian_n_contr (i_species,n_aux_gaussian(i_species)) &
                           = n_contracted

            do i_contracted = 1, n_contracted, 1
               aux_gaussian_alpha( i_species, n_aux_gaussian(i_species), &
                              i_contracted ) = aux_alpha(i_contracted)
               aux_gaussian_coeff( i_species, n_aux_gaussian(i_species), &
                              i_contracted ) = aux_coeff(i_contracted)
            enddo

         enddo

      else if (desc_str.eq."for_aux") then

!   read whole line of aux_basis input 
          read(inputline,*,end=88,err=99) desc_str,desc_str
!          write(use_unit,*) "aux_basis: ", desc_str
!   check type of the aux_basis function         
          if (desc_str .eq. "hydro") then
            read(inputline,*,end=88,err=99) desc_str,desc_str, n_shell, l_shell_str, z_eff

            l_shell = str_to_l(l_shell_str)

            if (l_shell.gt.l_shell_max(i_species)) then
                l_shell_max(i_species) = l_shell
            end if

!         verify validity of main quantum number
            if (n_shell.le.l_shell) then

                if (myid.eq.0) then
                    write(use_unit,*) "* Hydrogenic shell ", &
                        n_shell,l_shell_str,":", &
                        " n is too low for this l."
                    end if

               call aims_stop()    
            end if

            n_hydro(i_species) = n_hydro(i_species)+1
            hydro_n(i_species,n_hydro(i_species)) = n_shell
            hydro_l(i_species,n_hydro(i_species)) = l_shell
            hydro_scale(i_species,n_hydro(i_species)) = z_eff
            hydro_in_large_basis(i_species,n_hydro(i_species)) = .true.
 
!   Auxiliary basis functions should not be used to count possible states         
!            n_basis_sp = n_basis_sp + (2*l_shell+1)


            if (myid.eq.0) then
                write(use_unit,'(2X,A,I2,1X,A,1X,F7.3)') &
                    "| Found hydrogenic extended basis function : ", &
                    n_shell, l_shell_str, z_eff
            end if
            
          else if (desc_str.eq."ionic") then

             read(inputline,*,end=88,err=99) desc_str, desc_str, n_shell, l_shell_str, desc_str

              n_ionic(i_species) = n_ionic(i_species)+1

              l_shell = str_to_l(l_shell_str)
              if (l_shell.gt.l_shell_max(i_species)) then
                  l_shell_max(i_species) = l_shell
              end if

!         verify validity of main quantum number
              if (n_shell.le.l_shell) then
                  if (myid.eq.0) then
                      write(use_unit,*) "* Ionic shell ",n_shell, l_shell_str, ":", &
                        " n is too low for this l."
                  end if
                 call aims_stop()
              end if

              ionic_in_large_basis(i_species,n_ionic(i_species)) = .true.
              ionic_n (i_species, n_ionic(i_species)) = n_shell
              ionic_l (i_species, n_ionic(i_species)) = l_shell
              n_basis_sp = n_basis_sp + (2*l_shell+1)

              if (desc_str.eq."auto") then
                  auto_ionic(n_ionic(i_species)) = .true.
                  if (myid.eq.0) then
                      write(use_unit,'(2X,A,I2,1X,A,1X,A)') &
                        "| Found ionic extended basis function : ", n_shell, l_shell_str, &
                        ", default cutoff radius."
                  end if

              else
                  auto_ionic(n_ionic(i_species)) = .false.
                  read(inputline,*,end=88,err=99) desc_str, n_shell, l_shell_str, radius

                  ionic_rad (i_species, n_ionic(i_species)) = radius
                  if (myid.eq.0) then
                      write(use_unit,'(2X,A,I2,1X,A,1X,F7.3)') &
                          "| Found ionic extended basis function : ", &
                          n_shell, l_shell_str, radius
                  end if
              end if   
          else
            write(use_unit,*) &
                "Error: Unsupported extended basis function!", &
                " Will exit now!"
           call aims_stop()
        end if




      else if (desc_str.eq."basis_acc") then

         read(inputline,*,end=88,err=99) desc_str, basis_acc(i_species)

         flag_acc = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,E10.4)') &
                  "| Found basis singularity cutoff : ", &
                  basis_acc(i_species)
         end if

      else if (desc_str.eq."prodbas_acc") then

         read(inputline,*,end=88,err=99) desc_str, prodbas_acc(i_species)

         flag_prodbas_acc = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,E10.4)') &
                  "| Found auxiliary basis singularity cutoff : ", &
                  prodbas_acc(i_species)
         end if


      else if (desc_str.eq."innermost_max") then

         read(inputline,*,end=88,err=99) desc_str, innermost_max(i_species)
         flag_inmax = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A,I5)') &
                  "| Found threshold radial shell", &
                  " for innermost basis fn maxima : ", &
                  innermost_max(i_species)
         end if

      else if (desc_str.eq."cut_pot") then

         read(inputline,*,end=88,err=99) desc_str, r_cutoff(i_species), w_cutoff(i_species), &
            scale_cutoff(i_species)

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A,F10.5,1X,F10.5,1X,F10.5)') &
                  "| Found cutoff potl. onset [A], width [A], ", &
                  "scale factor : ", &
                  r_cutoff(i_species), w_cutoff(i_species), &
                  scale_cutoff(i_species)
         end if

         flag_cut_pot = .true.
         r_cutoff(i_species) = r_cutoff(i_species)/bohr
         w_cutoff(i_species) = w_cutoff(i_species)/bohr

      else if (desc_str.eq."cutoff_type") then

         read(inputline,*,end=88,err=99) desc_str, desc_str

         if (desc_str.eq."x2_(1-x2)") then
            cutoff_type(i_species) = 1
         else if (desc_str.eq."junquera") then
            cutoff_type(i_species) = 2
         else if (desc_str.eq."exp(1_x)_(1-x)2") then
            cutoff_type(i_species) = 3
         else
            if (myid.eq.0) then
            write(use_unit,'(1X,A,A,A)') &
                  "* Unknown cutoff potential type ", &
                  desc_str," - please correct!"
            end if
           call aims_stop()
         end if

         flag_cutoff_type = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A)') &
                  "| Found cutoff potential type : ", &
                  desc_str
         end if

      else if (desc_str.eq."cut_free_atom") then
         ! We can (but do not have to) set a cutoff potential which bounds the
         ! free atom calculation. By not doing this, we get a superposition of free
         ! atom densities at the outset which "exact", but which can also reach quite
         ! far away from the nucleus. Since out integration and density partitioning
         ! functions are both tied to the free atom density, a cutoff radius for the
         ! free atom will also affect the long-range tails of those operations.

         read(inputline,*,end=88,err=99) desc_str, desc_str

         if (desc_str.eq.'infinite') then

            write(info_str,'(2X,A,A)') &
            "| Found cutoff specification for free atom: ", &
            "No cutoff will be applied."
            call localorb_info(info_str)

            cut_free_atom(i_species) = .false.

         else if (desc_str.eq.'finite') then

            read(inputline,*,end=88,err=99) desc_str, desc_str, free_r_cut(i_species)

            write(info_str,'(2X,A,A,F10.5,A)') &
            "| Found cutoff specification for free atom: ", &
            "Cutoff potential begins at ",free_r_cut(i_species)," A."
            call localorb_info(info_str)

            cut_free_atom(i_species) = .true.
            free_r_cut(i_species) = free_r_cut(i_species)/bohr

         else

            write(info_str,'(1X,A,A,A,A)') &
            "* Cutoff specification for free atom: ", &
            "Invalid specification ",desc_str, " - please correct."
            call localorb_info(info_str)
           call aims_stop()

         end if

         flag_cut_free_atom = .true.

      else if (desc_str.eq."cut_core") then
         ! We can (but do not have to) set a cutoff potential which bounds the
         ! atomic core basis functions.

         read(inputline,*,end=88,err=99) desc_str, desc_str

         if (desc_str.eq.'infinite') then

            write(info_str,'(2X,A,A)') &
            "| Found cutoff specification for core functions: ", &
            "No cutoff will be applied."
            call localorb_info(info_str)

            cut_core(i_species) = .false.

         else if (desc_str.eq.'finite') then

            read(inputline,*,end=88,err=99) desc_str, desc_str, core_r_cut(i_species)

            write(info_str,'(2X,A,A,F10.5,A)') &
            "| Found cutoff specification for core functions: ", &
            "Cutoff potential begins at ",core_r_cut(i_species)," A."
            call localorb_info(info_str)

            cut_core(i_species) = .true.
            core_r_cut(i_species) = core_r_cut(i_species)/bohr

         else

            write(info_str,'(1X,A,A,A,A)') &
            "* Cutoff specification for core functions: ", &
            "Invalid specification ",desc_str, " - please correct."
            call localorb_info(info_str)
           call aims_stop()

         end if

         flag_cut_core = .true.
      else if (desc_str.eq.'basis_dep_cutoff') then
         read(inputline,*,end=88,err=99) desc_str, desc_str
         if (desc_str.eq.'.false.') then
            write(info_str,'(2X,A)') '| Explicitly not using basis function dependent cutoff potential.'
            call localorb_info(info_str)
            basis_dep_cutoff_thresh(i_species) = 0d0
         else
            read(inputline,*,end=88,err=99) desc_str, basis_dep_cutoff_thresh(i_species)
            write(info_str,'(2X,A,E14.6)') '| Threshold for basis-dependent cutoff potential is ', & 
                  basis_dep_cutoff_thresh(i_species) 
            call localorb_info(info_str)
            if ( (basis_dep_cutoff_thresh(i_species).lt.0.d0) .or. & 
                 (basis_dep_cutoff_thresh(i_species).ge.1.d0) ) then
               write(info_str,'(1X,A)') & 
                 '*  Error: The basis_dep_cutoff keyword is intended to remove a small positive'
               call localorb_info(info_str)
               write(info_str,'(1X,A)') & 
                 '*  fraction of the tail of a specified radial function, in units of the norm of that function.'
               call localorb_info(info_str)
               write(info_str,'(1X,A)') & 
                 '*  You are either requesting a cutoff value below zero or equal to or above one.'
               call localorb_info(info_str)
               write(info_str,'(1X,A)') & 
                 '*  A requested value of 1, for instance, would cut away all of your radial functions.'
               call localorb_info(info_str)
               write(info_str,'(1X,A)') & 
                 '*  This is not possible. Likely, this is a typo in your input file control.in.'
               call localorb_info(info_str)
               write(info_str,'(1X,A)') & 
                 '*  Please check your input file.'
               call localorb_info(info_str)
               call aims_stop('Illegal value requested for basis_dep_cutoff in control.in.', func)
            end if
         end if
      else if (desc_str.eq.'cut_atomic_basis') then 
         read(inputline,*,end=88,err=99) desc_str, cut_atomic_basis_functions(i_species)
         if (cut_atomic_basis_functions(i_species)) then
            write(info_str,'(2X,A)') "| Using individual cutoff potential for atomic basis functions."              
         else
            write(info_str,'(2X,A)') "| Explicitly NOT cutting atomic basis individually. "
         end if
         call localorb_info(info_str)           
      else if (desc_str.eq."radial_base") then

!         We read the number of the grid points, and the outermost integration
!         radius, of each species. scale_radial then follows from that radius ...
         read(inputline,*,end=88,err=99) desc_str, n_radial (i_species), &
            radius

            
!test - to make sure all control files reflect the new style - can be taken out
!       by anyone later
         if ((radius.lt.5.d0) ) then

            if (myid.eq.0) then
               write(use_unit,'(1X,A,A,A)') &
                     "* Just a safety check: ", &
                     "Outer radius seems small - ", &
                     "does the control file reflect the new input yet?"
            end if

!           call aims_stop()
         end if

!test end

         radius = radius/bohr

         scale_radial(i_species) = - radius / &
         log( 1.d0 - &
               (n_radial(i_species)/(1.d0+n_radial(i_species)))**2.d0 )

         flag_radial = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,I5,A,F8.3,A)') &
                  "| Found data for basic radial integration grid : ", &
                  n_radial(i_species), " points, outermost radius = ", &
                  radius*bohr, " A"
         end if

         
!test
!          write(use_unit,*) scale_radial(i_species)
!test end

      else if (desc_str.eq."radial_multiplier") then

!         integer scaling factor for the actual number of radial grid points
!         compared to the basic radial grid of N grid points:
!         1 : n_radial = N
!         2 : n_radial = 2*N+1
!         4 : n_radial = 4*N+3
!         ...
!         This is in preparation of an adaptive radial grid: Each larger grid contains
!         all radial shells of each smaller grid, so that we could easily scale up step
!         by step until a required total grid accuracy is reached
         read(inputline,*,end=88,err=99) desc_str, radial_multiplier

!         Sanity check
         if (radial_multiplier.lt.1) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') &
                     "* radial_multiplier less than 1 requested."
               write(use_unit,'(1X,A)') &
                     "* Not possible - please fix your control.in input file."
            end if
           call aims_stop()
         end if

         flag_radial_multiplier = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,I3)') &
                  "| Found multiplier for basic radial grid : ", &
                  radial_multiplier
         end if

      else if (desc_str.eq."angular") then

!         FIXME: All this should be renamed to "angular_max" for consistency.

         read(inputline,*,end=88,err=99) desc_str, angular_limit ( i_species )

         flag_angular = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A,I5)') &
                  "| Found max. number of angular integration ", &
                  "points per radial shell : ", &
                  angular_limit(i_species)
         end if

      else if (desc_str.eq."angular_min") then

         read(inputline,*,end=88,err=99) desc_str, angular_min ( i_species )

         flag_angular_min = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A,I5)') &
                  "| Found min. number of angular integration ", &
                  "points per radial shell : ", &
                  angular_min(i_species)
         end if


      else if (desc_str.eq."angular_acc") then

         read(inputline,*,end=88,err=99) desc_str, angular_acc ( i_species )

         flag_angular_acc = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A,E10.4)') &
               "| Found accuracy criterion for ", &
               "angular integrations : ", &
               angular_acc(i_species)
         end if

         if (angular_acc(i_species).gt.0.d0) then
            if (n_periodic.eq.0) then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                     "| Will adapt angular grid densities ", &
                     "automatically."
               end if

            else
            ! Periodic boundary conditions.
            ! Automatically converged angular grids are not supported for periodic
            ! boundary conditions because all integration points are first mapped into the  
            ! central unit cell, breaking the spherical symmetry of each shell.
            ! This leads to overly dense grids (lack of error cancellation)
            ! that are ultimately unusable.

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                     "| Adapt angular grids are not supported ", &
                     "for periodic boundary conditions. "
                  write(use_unit,'(2X,A,A)') &
                     "| Please use specified angular grids ", &
                     "instead. "
               end if

              call aims_stop()

            end if
         end if

      else if (desc_str.eq."angular_grids") then

         read(inputline,*,end=88,err=99) desc_str, desc_str

         if (desc_str.eq."specified") then

            specified_grid(i_species) = .true.

            if (myid.eq.0) then
               write(use_unit,'(2X,A,A)') &
                     "| Found angular grid specification: ", &
                     "user-specified."
            end if

            call read_specified_grid &
            ( i_species )

         else if (desc_str.eq."auto") then

            specified_grid(i_species) = .false.

            if (myid.eq.0) then
               write(use_unit,'(2X,A,A)') &
                     "| Found angular grid specification: ", &
                     "automatic."
            end if

         else

            if (myid.eq.0) then
               write(use_unit,'(1X,A,A,A)') &
                  "* angular_grids: Unknown tag ", &
                  desc_str, ". Please correct."
            end if
           call aims_stop()

         end if

         flag_specify_grids = .true.

      else if (desc_str.eq."logarithmic") then

         read(inputline,*,end=88,err=99) desc_str, r_grid_min(i_species), &
            r_grid_max(i_species), r_grid_inc(i_species)

         ! Only for forces
         log_r_grid_inc(i_species) = log(r_grid_inc(i_species))

         flag_log = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,E10.4,1X,E10.4,1X,E10.4)') &
                  "| Found logarithmic grid data [bohr] : ", &
                  r_grid_min(i_species), r_grid_max(i_species), &
                  r_grid_inc(i_species)
         end if

      else if (desc_str.eq."l_hartree") then

         read(inputline,*,end=88,err=99) desc_str, l_hartree(i_species)

         flag_l_hartree = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,I3)') &
                  "| Found l_max for Hartree potential  : ", &
                  l_hartree(i_species)
         end if

         ! Safety checks to avoid unwanted exceptions
         if (l_hartree(i_species).lt.0) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') &
                  "* Warning: l_hartree cannot be less than zero."
               write(use_unit,'(1X,A)') &
                  "* Likely, this points to a bug in your control.in file. Please check."
               write(use_unit,'(1X,A)') &
                  "*call aims_stop()ping program execution for now."
            end if
           call aims_stop()
         end if

         if (l_hartree(i_species).gt.16) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') &
                  "* Warning: You requested l_hartree = 16 or larger."
               write(use_unit,'(1X,A)') &
                  "* While this is possible mathematically, we have not implemented such high"
               write(use_unit,'(1X,A)') &
                  "* multipole expansions so far."
               write(use_unit,'(1X,A)') &
                  "* If needed, this could be implemented, but would require careful testing to"
               write(use_unit,'(1X,A)') &
                  "* make sure that no numerical problems arise with terms like mu_lm/r^(16+1) ..."
               write(use_unit,'(1X,A)') &
                  "* With the present code version, mustcall aims_stop() the program execution."
            end if
           call aims_stop()
         end if


      else if (desc_str.eq."core_states") then

         read(inputline,*,end=88,err=99) desc_str, n_core_states_species(i_species)

         if (myid.eq.0) then
            write(use_unit,'(2X,A,I3)') &
                  "| Number of core states : ", &
                  n_core_states_species(i_species)
         end if

      else if (desc_str.eq."KH_core_states") then

         read(inputline,*,end=88,err=99) desc_str, n_KH_core_states_species(i_species)

         if (myid.eq.0) then
            write(use_unit,'(2X,A,I3)') &
                  "| Number of KH core states : ", &
                  n_KH_core_states_species(i_species)
         end if

         if(n_KH_core_states_species(i_species)>0)then
            flag_KH_core_states = .true.
         end if



      else if (desc_str.eq."include_min_basis") then

         read(inputline,*,end=88,err=99) desc_str, include_min_basis(i_species)

         flag_include_min_basis = .true.

         if (myid.eq.0) then
            write(use_unit,'(2X,A,G5.2)') &
                  "| Found request to include minimal basis fns. : ", &
                  include_min_basis(i_species)
         end if

      else if (desc_str.eq."max_n_prodbas") then
         if(use_prodbas) then
            read(inputline,*,end=88,err=99) desc_str, max_n_prodbas(i_species)
            call localorb_info('** Now ignoring max_n_prodbas.', use_unit, '(2X,A)')
            flag_max_n_prodbas = .true.
         endif

      else if (desc_str.eq."max_l_prodbas") then
         if(use_prodbas) then
            read(inputline,*,end=88,err=99) desc_str, max_l_prodbas(i_species)
            if (myid.eq.0) then
               write(use_unit,'(2X,A,A,I4)') &
                     "| Found request of the maximal angular quantum ", &
                     "number of the product basis. : ", &
                     max_l_prodbas(i_species)
            end if
            flag_max_l_prodbas = .true.
         endif

      else if ((desc_str .eq. "plus_u") .or. &
            (desc_str .eq. "plus_U")) then

         plus_u_flag(i_species) = .true.
         i_shell_plus_u = i_shell_plus_u + 1
         read(inputline,*,end=88,err=99) desc_str, plus_u_n(i_species, i_shell_plus_u), &
               plus_u_l_str(i_species, i_shell_plus_u), &
               plus_u_value(i_species, i_shell_plus_u)
         plus_u_l(i_species, i_shell_plus_u) =  str_to_l(plus_u_l_str(i_species, i_shell_plus_u))
         if (myid .eq. 0) then
            write(use_unit,'(2X,A,I2,A,A,A,F12.6,A)') &
            "| Found request for '+U' treatment for ", &
                  plus_u_n(i_species, i_shell_plus_u), " ", &
                  plus_u_l_str(i_species, i_shell_plus_u), " shell, U = ", &
                  plus_u_value(i_species, i_shell_plus_u), " eV."
         end if

      else if ((desc_str .eq. "plus_u_ramping") .or. &
              (desc_str .eq. "plus_U_ramping")) then
         read(inputline,*,end=88,err=99) desc_str, plus_u_ramping_increment(i_species)

      else if (desc_str .eq. "hubbard_coefficient") then
         read(inputline,*,end=88,err=99) desc_str, plus_u_hubc(i_species, 1), &
               plus_u_hubc(i_species, 2), &
               plus_u_hubc(i_species, 3), &
               plus_u_hubc(i_species, 4)
         if (myid .eq. 0) then
            write(use_unit,'(2X,A,I2,A,A,A,F12.6,A)') &
            "| Found coefficents for linear combination of hubbard projectors "
         end if

      else if ( (desc_str.eq."division") .or. &
               (desc_str.eq."outer_grid") ) then

         if (myid.eq.0) then
            write(use_unit,'(1X,A,A)') &
            "* Misplaced specified grid descriptor: ", &
            desc_str
            write(use_unit,'(1X,A,A)') &
            "* Can't be sure that your grid is ok - please correct."
            end if
        call aims_stop()

      else if (desc_str.eq."pp_charge") then
         read(inputline,*,end=88,err=99) desc_str, &
            pp_charge(i_pp_species)

         if (myid.eq.0) then
            write(use_unit,*) ' | Found pseudocore charge  : ',pp_charge(i_pp_species) 
         end if
      else if (desc_str.eq."pp_local_component") then
         read(inputline,*,end=88,err=99) desc_str, &
         pp_local_component(i_pp_species)

         if (myid.eq.0) then
            write(use_unit,*) ' | For this pseudo species ',pp_local_component(i_pp_species), &
                   ' is used for the local component!' 
         end if


      else if (desc_str.eq."cite_reference") then
         read(inputline,*,end=88,err=99) desc_str, desc_str

         call cite_reference(trim(desc_str))

      else
!         must have reached the end of the current species description.

         backspace(7)
         exit lineloop

      end if

   enddo lineloop

   if (linecount == 0) then
      if (myid.eq.0) then
         write(use_unit,*) "No data for species. Assuming defaults ", &
               "(occupied shells only)."
      end if
   end if


!  Anything which was not set should come from species defaults

!       were any valence electrons specified?


   if (.not.flag_valence) then
      if (myid.eq.0) then
         write(use_unit,'(1X,A,A,A,A)') &
         "* Species ", species_name(i_species), ": Missing ", &
         "valence electrons for neutral atomic configuration."
      end if

     call aims_stop()
   end if

!       were all other input values provided?

   if (.not.flag_z) then

      if (myid.eq.0) then
         write(use_unit,*) &
               "Species ", species_name(i_species), ": ", &
               "Missing nuclear charge."
      end if

     call aims_stop()
   end if

   if (.not. flag_element) then
       ! In the case of fractional charges, we need an unambiguous association
       ! of the species with a chemical element. This is done for
       ! the sake of van der Waals corrections, which rely on atomic data
       ! and must thus necessarily know what is the intended chemical element.
       if ( (abs(species_z(i_species)-nint(species_z(i_species))) < 0.34d0) .and. &
            (nint(species_z(i_species)).ge.1) .and. (nint(species_z(i_species)).le.102) ) then
           ! This case covers (n).66 < Z < (n+1).34 , for Z = 1 through 102.
           ! In this case, we use the nearest integer nuclear charge
           ! to determine which element we mean. 
           species_element(i_species) = periodic_table(nint(species_z(i_species)))
       else if (get_nucleus(species_name(i_species)) /= -1) then
           species_element(i_species) = species_name(i_species)(1:2)
       else if (nint(species_z(i_species)) == 0) then
           species_element(i_species) = 'X'
       else
           ! All other cases are treated as unknown and cause aims to abort.
           if (myid == 0) then
               write (use_unit, '(1X,A,F10.3)') &
                   "*** Error: This species has a 'nucleus' charge of ", &
                   species_z(i_species)
               write (use_unit, '(1X,2A)') &
                   '*** This nucleus is located sufficiently far from an existing chemical element'
               write (use_unit, '(1X,2A)') &
                   '*** that it becomes somewhat ambiguous to decide which chemical element is intended.'
               write (use_unit, '(1X,2A)') &
                   '*** For some internal decisions (related to van der Waals), FHI-aims needs to know '
               write (use_unit, '(1X,2A)') &
                   '*** the chemical element identity intended for this calculations.'
               write (use_unit, '(1X,2A)') &
                   "*** Please use the 'element' keyword in control.in, in addition to the 'species' keyword,"
               write (use_unit, '(1X,2A)') &
                   "*** to unambiguously designate the intended identity of this nucleus."
           end if
           call aims_stop_coll()
       end if
   end if

   if (.not.flag_cutoff_type) then

      ! Default to exp(1_x)_(1-x)2 type cutoff potential
      cutoff_type(i_species) = 3
      if (myid.eq.0) then
         write(use_unit,'(2X,A,A,A,A)') &
         "Species ", species_name(i_species), ": ", &
         "Missing cutoff potential type."
         write(use_unit,'(2X,A)') &
         "Defaulting to exp(1/x)/(1-x)^2 type cutoff potential."
      end if

   end if

   if (.not.flag_cut_pot) then
      ! use proven defaults

      if (myid.eq.0) then
         write(use_unit,'(2X,A,A,A,A)') &
         "Species ", species_name(i_species), ": ", &
         "Missing cutoff potential data."
      end if

      if (cutoff_type(i_species).eq.1) then

         r_cutoff(i_species) = 5.0   ! 5 A
         w_cutoff(i_species) = 2.5   ! 2.5 A
         scale_cutoff(i_species) = 1.d0

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A,F10.5,1X,F10.5,1X,F10.5)') &
               "| Defaulting to cutoff potl. onset [A], width [A], ", &
               "scale factor : ", &
               r_cutoff(i_species), w_cutoff(i_species), &
               scale_cutoff(i_species)
         end if

      else if (cutoff_type(i_species).eq.2) then

         r_cutoff(i_species) = 5.0   ! 5 A
         w_cutoff(i_species) = 2.0   ! 2 A
         scale_cutoff(i_species) = 1.d0

         if (myid.eq.0) then
         write(use_unit,'(2X,A,A,F10.5,1X,F10.5,1X,F10.5)') &
            "| Defaulting to cutoff potl. onset [A], width [A], ", &
            "scale factor : ", &
            r_cutoff(i_species), w_cutoff(i_species), &
            scale_cutoff(i_species)
         end if

      else if (cutoff_type(i_species).eq.3) then

         r_cutoff(i_species) = 5.0   ! 5 A
         w_cutoff(i_species) = 2.0   ! 2.0 A
         scale_cutoff(i_species) = 1.d0

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A,F10.5,1X,F10.5,1X,F10.5)') &
               "| Defaulting to cutoff potl. onset [A], width [A], ", &
               "scale factor : ", &
               r_cutoff(i_species), w_cutoff(i_species), &
               scale_cutoff(i_species)
         end if

      else
        call aims_stop()
      end if

      ! convert to bohrs
      r_cutoff(i_species) = r_cutoff(i_species)/bohr
      w_cutoff(i_species) = w_cutoff(i_species)/bohr

   end if

      ! Quick adjustment of grid parameters if not all were specified
   if (specified_grid(i_species)) then


!        Check that angular grid is big enough for l_hartree
      do i_shell = 1, n_ang_shells(i_species), 1

         call  verify_angular_grid &
               ( n_ang_points(i_shell,i_species), &
               l_hartree(i_species), new_int, info_log)

         if(info_log)then

            if (myid.eq.0) then
               write(use_unit,'(1X,A,A2,A,A,A,I5,A)') &
                  "* Species ", species_name(i_species), &
                  ": Specified number of angular ", &
                  "integration points", &
                  " is ",    n_ang_points(i_shell,i_species), ","
               write(use_unit,'(1X,A,I3,A,I5,A,A,I5,A)') &
                  "* but required value for l_hartree = ", &
                  l_hartree(i_species), &
                  " is ", new_int,". "

            end if
            n_ang_points(i_shell,i_species) =  new_int

         end if
      end do


      ! If other angular parameters were not set, adjust for consistency
      if (.not.flag_angular) then

         angular_limit(i_species) = 0
         do i_shell = 1, n_ang_shells(i_species), 1
            angular_limit(i_species) = &
               max( angular_limit(i_species), &
                     n_ang_points(i_shell,i_species) )
         enddo
         flag_angular = .true.

      end if

      if (.not.flag_angular_min) then

         angular_min(i_species) = angular_limit(i_species)
         do i_shell = 1, n_ang_shells(i_species), 1
            angular_min(i_species) = &
               min( angular_min(i_species), &
                  n_ang_points(i_shell,i_species) )
         enddo
         flag_angular_min = .true.

      end if

      if (.not.flag_angular_acc) then
         ! If grid was specified, do not require adaptation by default.

         angular_acc(i_species) = 0.d0
         flag_angular_acc = .true.

      end if


   end if

!       next, check logarithmic grid
   if (.not.flag_log) then

      r_grid_min(i_species) = 0.0001
      r_grid_max(i_species) = 100.
      r_grid_inc(i_species) = 1.0123

      ! Only for forces
      log_r_grid_inc(i_species) = log(r_grid_inc(i_species))         

      if (myid.eq.0) then
         write(use_unit,'(2X,A,A2,A,A)') &
               "Species ", species_name(i_species), ": ", &
               "No 'logarithmic' tag. Using default grid for free atom:"
         write(use_unit,'(2X,A,E10.4,1X,E10.4,1X,E10.4)') &
               "| Default logarithmic grid data [bohr] : ", &
               r_grid_min(i_species), r_grid_max(i_species), &
               r_grid_inc(i_species)
      end if

   end if

   if (.not.flag_l_hartree) then

      if (myid.eq.0) then
         write(use_unit,*) &
               "Species ", species_name(i_species), ": ", &
               "Missing Hartree l_max specification."
      end if

     call aims_stop()
   end if

!  now verify species data

!     complete and verify valence occupation

!     valence orbitals are a type of basis fns which is always present
   i_types = 1
   n_atomic(i_species) = 0

   test_occ = 0.

   do i_l = 0, l_shell_max(i_species), 1
      do i_shell = i_l+1, valence_n_max(i_l,i_species),1
         test_occ = test_occ + &
         valence_occ(i_shell, i_l, i_species)
         if (valence_occ(i_shell,i_l,i_species).gt.0.d0) then
         
! This should be change to use something else than .ge.                    
            if(use_relativistic_basis) then
               if(i_l.gt.0) then
                  n_atomic(i_species)=n_atomic(i_species)+1
                  atomic_n(i_species,n_atomic(i_species))=i_shell
                  atomic_l(i_species,n_atomic(i_species))=i_l
                  atomic_k(i_species,n_atomic(i_species))=i_l             
                  n_atomic(i_species)=n_atomic(i_species)+1
                  atomic_n(i_species,n_atomic(i_species))=i_shell
                  atomic_l(i_species,n_atomic(i_species))=i_l
                  atomic_k(i_species,n_atomic(i_species))=-(i_l+1)
               else
                  n_atomic(i_species)=n_atomic(i_species)+1
                  atomic_n(i_species,n_atomic(i_species))=i_shell
                  atomic_l(i_species,n_atomic(i_species))=i_l
                  atomic_k(i_species,n_atomic(i_species))=-(i_l+1)
               end if
            else
               n_atomic(i_species)=n_atomic(i_species)+1
               atomic_n(i_species,n_atomic(i_species))=i_shell
               atomic_l(i_species,n_atomic(i_species))=i_l
            end if
            if (i_shell.le.core_n_max(i_l,i_species)) then
               core_fn(i_species,n_atomic(i_species)) = .true.
            end if
         end if
      enddo
   enddo

!       original      
!      do i_l = 0, l_shell_max(i_species), 1
!        do i_shell = i_l+1, valence_n_max(i_l,i_species),1
!          test_occ = test_occ + &
!          valence_occ(i_shell, i_l, i_species)
!          if (valence_occ(i_shell,i_l,i_species).gt.0.d0) then
!            n_atomic(i_species)=n_atomic(i_species)+1
!            atomic_n(i_species,n_atomic(i_species))=i_shell
!            atomic_l(i_species,n_atomic(i_species))=i_l
!            if (i_shell.le.core_n_max(i_l,i_species)) then
!              core_fn(i_species,n_atomic(i_species)) = .true.
!            end if
!          end if
!        enddo
!      enddo

   if (abs(test_occ-species_z(i_species)).gt.1d-6 ) then
      if (myid.eq.0) then
         write(use_unit,*) "* Incorrect electron count for species."
         write(use_unit,*) "* Should be ", &
            species_z(i_species), ", is ", test_occ, " ."
      end if

     call aims_stop()

   end if

   if ((flag_include_min_basis).and.((n_pp_species.eq.0).and.(n_empty_atoms.eq.0)).and..not.use_sto)then
      if (.not.include_min_basis(i_species)) then
       if (.not.override_warning_nobasis) then
!         Check if there are enough other orbitals to treat
!         at least the occupied states ...
         if (dble(2*n_basis_sp).lt.species_z(i_species)) then

            if (myid.eq.0) then
               write(use_unit,'(1X,A)') &
                     "* Too few basis functions provided for species."
               write(use_unit,'(1X,A,A)') &
                     "* Please request additional basis functions", &
                     " in input file control.in ."
            end if

           call aims_stop()
         end if
       endif
      end if
   end if



   if (.not.flag_include_min_basis) then
!       quietly set default for use of minimal basis to true -
!       this is a test flag only and should not bother the user.
      include_min_basis(i_species) = .true.
   end if

   ! check whether i_species has basis functions at all

   if((n_basis_sp.eq.0).and.(.not.include_min_basis(i_species)) ) then
        no_basis(i_species) = .true.
   endif

!     Complete and verify ionic configuration, if ionic basis fns were requested.

   if ( n_ionic(i_species).gt.0 ) then

      test_occ = 0.

!       Fill up ionic occupation by atomic occupation, unless otherwise specified
      do i_l = 0, l_shell_max(i_species), 1
         do i_shell = i_l+1, ion_n_max(i_l,i_species), 1
            if (.not.flag_ionic(i_l, i_shell)) then
            ion_occ(i_shell,i_l,i_species) = &
            valence_occ(i_shell,i_l,i_species)
            end if
            test_occ = test_occ+ion_occ(i_shell,i_l,i_species)
         enddo
      enddo

!       Sanity check for ionic configuration
      test_occ = test_occ - species_z(i_species)
      if (test_occ.gt.0. ) then

         if (myid.eq.0) then
            write(use_unit,'(1X,A,A)') "* Prescribed ionic ", &
                  "occupation would be negative."
            write(use_unit,'(1X,A,A)') "* Please set if you wish to ", &
                  "use ionic functions."
         end if

        call aims_stop()
      else if (test_occ.eq.0. ) then

         if (myid.eq.0) then
            write(use_unit,*) "* Prescribed ionic occupation is neutral."
            write(use_unit,'(1X,A,A)') "* Please set if you wish to ", &
                  "use ionic functions."
         end if

        call aims_stop()
      else if (test_occ.le.(-species_z(i_species))) then

         if (myid.eq.0) then
            write(use_unit,'(1X,A,A)') "* Prescribed ionic ", &
                  "occupation has no electrons."
            write(use_unit,'(1X,A,A)') "* Please set if you wish to use ", &
                  "ionic functions."
         end if

        call aims_stop()
      else

         if (myid.eq.0) then
            write(use_unit,'(2X,A,F4.1,A,A,A)') &
                  "| Will include ionic basis functions of ", &
                  (-test_occ), "-fold positive ", &
                  species_name(i_species), " ion."
         end if

!         one more basis fn type
         i_types = i_types+1
      end if

      do i_basis = 1, n_ionic(i_species), 1
         if (auto_ionic(i_basis)) then
            ! requested automatic cutoff onset assignment for this function
            ionic_rad (i_species, i_basis) = r_cutoff(i_species)
         end if
      enddo

!     Ionic occupation complete.
   end if

!     Verify confined states.

   if (n_conf(i_species).gt.0) then

      i_types = i_types+1

      do i_basis = 1, n_conf(i_species), 1
         if (auto_conf(i_basis)) then
            ! requested automatic cutoff onset assignment for this function
            conf_rad (i_species, i_basis) = r_cutoff(i_species)
         end if
      enddo

   end if

!     Verify hydrogenic states.

   if (n_hydro(i_species).gt.0) then

      i_types = i_types+1

   end if

!     Verify cartesian Gaussians.

   if (n_gaussian(i_species).gt.0) then

      i_types = i_types+1

   end if

   if (n_sto(i_species).gt.0) then

      i_types = i_types+1

   end if

!     If there were extra basis functions, these must be orthogonalized. Do we
!     know the accuracy?

!      if (i_types.gt.1) then
   if (.not.flag_acc) then
      basis_acc(i_species) = 1.0d-4

      if (myid.eq.0) then
         write(use_unit,'(2X,A,A2,A,A,A)') &
            "Species ", species_name(i_species), ": ", &
            "On-site basis accuracy parameter (for Gram-Schmidt ", &
            "orthonormalisation) not specified."
         write(use_unit,'(2X,A,E14.7,A)') &
            "Using default value basis_acc = ", &
            basis_acc(i_species), "."
      end if

   end if

   if (.not.flag_inmax) then
!       use default value
      innermost_max(i_species) = 2

      if (myid.eq.0) then
         write(use_unit,'(2X,A,A,A,A,A,I3,A)') &
               "Species ", species_name(i_species), ": ", &
               "Using default innermost maximum threshold ", &
               "i_radial=", innermost_max(i_species), &
               " for radial functions."
      end if

      flag_inmax = .true.
   end if

!  Verify cutoff potential settings for the free atom

   if (.not.flag_cut_core) then
   ! Use default:
   ! do not yet use a free-atom cutoff potential yet
   ! This should be changed asap at least for
   ! * the periodic case
   ! * charged molecules

      cut_core(i_species) = .false.
      core_r_cut(i_species) = r_cutoff(i_species)

   end if

   if (.not.flag_cut_free_atom) then
   ! Use default:

      if ( (use_prodbas) .and. .not. ( use_lvl .and. use_logsbt )) then
         ! The product basis for the Coulomb potential is used and integrated
         ! for Hartree-Fock, hybrid functionals, MP2, GW etc.
         ! The required Coulomb integral is different from the integrals required in
         ! normal DFT-LDA/GGA in that it extends over a much longer range.
         !
         ! Among other things, cut_free_atom controls the extent of the partition function
         ! used in our integrations, and thus also the range that can maximally be
         ! integrated using our standard technology.
         !
         ! In the long run, the numerical integrals for the Coulomb matrix must 
         ! be separated from all other "standard DFT" integrals, but there is no quick 
         ! fix to do this (yet). Therefore, we extend the partition function as much as possible
         ! for now. This is a hack that happens to work - not a real fix.

         cut_free_atom(i_species) = .false.
         free_r_cut(i_species) = r_grid_max(i_species)

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A,A,A)') &
               "Species ", species_name(i_species), ": ", &
               "Default cutoff onset for free atom density etc. is infinite"
            write(use_unit,'(2X,A)') &
               "since the product basis is used (hybrid functionals, Hartree-Fock, GW etc.)."
         end if

      else
         ! no product basis - standard DFT-LDA/GGA only

         ! or, in fact, periodic Hartree-Fock based on LVL RI

         ! Set free_r_cut for rigorously bounded basis functions (only
         ! standard DFT-LDA / GGA). Choose cut_free_atom as tight as
         ! we can in this case.  In this case, default to cutoff onset
         ! of basis functions - no aspect of the final charge density
         ! can reach beyond this radius.

         ! Note: free_r_cut is also set for Gaussian orbitals here,
         ! but it may be modified later because Gaussian are not
         ! bounded in a transparent way and required special
         ! treatment.

         ! Default would be the species-wide onset radius specified by cut_pot
         cut_free_atom(i_species) = .true.
         free_r_cut(i_species) = r_cutoff(i_species)

         ! check confined and ionic basis functions separately
         ! must increase cut_free_atom if specified cutoff radius exceeds cut_pot
         do i_basis = 1, n_conf(i_species), 1
            if ( free_r_cut(i_species).lt. conf_rad(i_species,i_basis)) &
                 then
               free_r_cut(i_species) = conf_rad(i_species,i_basis)
            end if
         enddo
         do i_basis = 1, n_ionic(i_species), 1
            if ( free_r_cut(i_species).lt. ionic_rad(i_species,i_basis)) &
                 then
               free_r_cut(i_species) = ionic_rad(i_species,i_basis)
            end if
         enddo

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A,A,A,E15.8,A)') &
                 "Species ", species_name(i_species), ": ", &
                 "Default cutoff onset for free atom density etc. :", &
                 free_r_cut(i_species)*bohr, &
                 " AA."
         end if

      end if

   end if

   if (cut_free_atom(i_species)) then
      ! Next, check that the free atom cutoff potential at least covers all basis functions
      ! in the code
      !
      ! This is needed because the free-atom cutoff determines the extent of
      ! * the integral partition function
      ! * the density partition function for the Hartree potential

      ! check general confinement of minimal and hydrogenic basis functions
      if (free_r_cut(i_species).lt.r_cutoff(i_species)) then
         ! We could set the default reasonably here, but I prefer to simplycall aims_stop()
         ! for now - change later.
         write(info_str,'(2X,A,A,A,A,A)') &
               "* Species ", species_name(i_species), ": ", &
               "Warning - cutoff potential for free atom ", &
               "below that for basis functions."
         call localorb_info(info_str)
         write(info_str,'(2X,A)') &
               "* Please correct."
        call aims_stop()
      end if

      ! check confined and ionic basis functions separately
      do i_basis = 1, n_conf(i_species), 1
         if ( free_r_cut(i_species).lt. conf_rad(i_species,i_basis)) then
            write(info_str,'(2X,A,A,A,A,A,I5,A)') &
            "* Species ", species_name(i_species), ": ", &
            "Warning - cutoff radius for free atom ", &
            "below that for confined fn. number ", i_basis,"."
            call localorb_info(info_str)
            write(info_str,'(2X,A)') &
            "* Please correct."
           call aims_stop()
         end if
      enddo

      do i_basis = 1, n_ionic(i_species), 1
         if ( free_r_cut(i_species).lt. ionic_rad(i_species,i_basis)) then
            write(info_str,'(2X,A,A,A,A,A,I5,A)') &
            "* Species ", species_name(i_species), ": ", &
            "Warning - cutoff radius for free atom ", &
            "below that for ionic fn. number ", i_basis,"."
            call localorb_info(info_str)
            write(info_str,'(2X,A)') &
            "* Please correct."
           call aims_stop()
         end if
      enddo

      if (n_gaussian(i_species).gt.0) then
         ! give a warning message for now.
         ! We can only ensure proper treatment of free-atom cutoff radius
         ! with Gaussians later in the code, when we actually know the proper radial function.
         ! This is not implemented yet.

            write(info_str,'(2X,A,A,A,A)') &
            "* Species ", species_name(i_species), ": Warning - ", &
            "free-atom cutoff radius not yet verified for Gaussians. "
            call localorb_info(info_str)
            write(info_str,'(2X,A,A)') &
            "* This should be corrected by adding a check ", &
            "for each gaussian radial function later."
            write(info_str,'(2X,A,A)') &
            "* You may wish to check the convergence of ", &
            "your cut_free_atom setting by hand!"

      end if

   end if


!  rescale logarithmic grid

   r_grid_min(i_species) = &
      r_grid_min(i_species)/species_z(i_species)

!  increment basic radial grid acccording to radial_multiplier
   if (flag_radial_multiplier) then
      n_radial(i_species) = radial_multiplier * &
      (n_radial(i_species)+1) - 1

      if (myid.eq.0) then
         write(use_unit,'(2X,A,A,A,A,A,I3,A,I5,A)') &
               "Species ", species_name(i_species), ": Basic radial", &
               " grid will be enhanced according ", &
               "to radial_multiplier = ", &
               radial_multiplier, ", to contain ", n_radial(i_species), &
               " grid points."
      end if

   end if

   ! check whether i_species has basis functions at all

!   if((n_basis_sp.eq.0).and.(.not.flag_include_min_basis) ) then
!        no_basis(i_species) = .true.
!   endif


!   if(use_prodbas) then
!     if(myid.eq.0) then
!       write(use_unit, '(2X,A,A)') "Hartree-Fock or related methods is used here ", &
!        "and set the basis depth cutoff parameter (basis_dep_cutoff) to zero. "
!     endif
!     basis_dep_cutoff_thresh(i_species) = 0d0
!   endif
   return

88 continue
   write(use_unit,'(A)') & 
   "Syntax error reading data for species '"//trim(species_name(i_species))//"' from 'control.in' (missing arguments)"
   write(use_unit,'(A)') "in line: '"//trim(inputline)//"'"
  call aims_stop()

99 continue
   write(use_unit,'(A)') "Syntax error reading data for species '"//trim(species_name(i_species))//"' from 'control.in'"
   write(use_unit,'(A)') "in line: '"//trim(inputline)//"'"
  call aims_stop()

end subroutine read_species_data
