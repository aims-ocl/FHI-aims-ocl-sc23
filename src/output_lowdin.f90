!****s* FHI-aims/output_lowdin
!  NAME
!    output_lowdin
!  SYNOPSIS

subroutine output_lowdin &
     ( KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, occ_numbers, overlap_matrix, &
     filename,chemical_potential, n_electrons )

!  PURPOSE
!  The subroutine makes and prints out Loewdin analysis and l-projected density of states.
!  Lowdin weigth for orbital i  w_i = \sum_{jk} sqrt(S)_ij D_jk sqrt(S)_ki
!  where S is the overlap matrix for orbital basis, and D is the density matrix.
!
!  USES

  use dimensions
  use constants
  use mpi_tasks
  use localorb_io
  use basis
  use geometry
  use species_data
  use runtime_choices
  use pbc_lists
  use lapack_wrapper
  use synchronize_mpi
!  use density_matrix_evaluation
  use scalapack_wrapper, only: my_k_point, construct_lowdin_decomp_scalapack
  use generate_aims_uuid, only: write_aims_uuid
  use arch_specific, only: arch_erf
  implicit none

!  ARGUMENTS

  real*8, dimension(n_states, n_spin, n_k_points)               :: KS_eigenvalue 
  real*8, dimension(n_basis, n_states, n_spin, n_k_points_task) :: KS_eigenvector
  complex*16, dimension(n_basis, n_states, n_spin, n_k_points_task) :: KS_eigenvector_complex
  real*8, dimension(n_states, n_spin,n_k_points)                :: occ_numbers
  real*8, dimension( n_hamiltonian_matrix_size )                :: overlap_matrix
  real*8:: chemical_potential
  character*40 :: filename
  real*8 :: n_electrons
  character(LEN=80) :: uuid_str


!  INPUTS
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o KS_eigenvector -- Kohn-Sham eigenvectors if real format
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors if complex format
!   o occ_numbers -- occupations of the states
!   o overlap_matrix -- overlap matrix
!   o chemical_potential -- chemical potential
!   o filename -- file where results are printed
!   o n_electrons -- number of electrons in systems
!
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
!

  ! local variables

  ! We will do a Loewdin charge decomposition by the following criteria:
  ! Number of electrons in each KS state, per atom, spin channel, and 
  ! angular momentum component
  ! lowdin_decomp will contain this decomposition. All derived quantities
  ! are then accessible by way of appropriate sums.
  real*8, dimension( 0:l_wave_max, n_atoms, n_states, n_spin,n_k_points ) :: lowdin_decomp
  character*40 :: proj_dos_filename
  real*8, dimension( 0:l_wave_max, n_atoms, n_spin) :: l_projected_charge
  real*8, dimension( n_atoms, n_spin ) :: at_projected_charge
  real*8, dimension( n_atoms ) :: spin_per_atom
  real*8, dimension( n_spin ) :: total_charge
  real*8 :: charge_difference

  integer, dimension(n_spin,n_k_points) :: max_occ_number
  integer :: n_nonsingular 

  character*100 :: info_str, outputformat
  character*100 :: outputformat_header

  real*8  :: diff_electrons
  real*8  :: de, occu, E1, E2, buf
  real*8  :: en
  real*8, dimension(:,:,:,:), allocatable :: KS_dos
  real*8, dimension(:,:,:),   allocatable :: KS_sum_dos
  real*8, dimension(:),       allocatable :: proj_dos_erfs
  real*8, dimension(:),       allocatable :: ovlp_eigenvalues
  real*8, dimension(:,:),     allocatable :: ovlp_transform
  real*8, dimension(:,:),     allocatable :: sqrt_overlap
  complex*16:: mul_temp

  ! darn, need an extra aux array of characters
  character*3,dimension(0:l_wave_max) :: l_channel

  ! counters

  integer :: i_spin, i_state
  integer :: i_basis_1, i_basis_2
  integer :: i_index
  integer :: i_info
  integer :: i_atom, i_l, i_k, i_k_point, i_size, i_cell, i_e,i_species

  ! begin work

  write(info_str,'(2X,A)') &
       ' '
  call localorb_info(info_str)
  write(info_str,'(2X,A)') &
       'Performing Loewdin charge analysis on all atoms.'
  call localorb_info(info_str)


  lowdin_decomp = 0.d0

  ! Notice we do not sum up the density matrix because want a per-state 
  ! projection analysis first - Mulliken charges are summed up only thereafter.

if(use_scalapack) then

  call construct_lowdin_decomp_scalapack(lowdin_decomp(0,1,1,1,my_k_point))

else

  ! First obtain the square root of the overlap matrix -- a crucial step for
  ! Loewdin analysis.

  if(n_periodic == 0 .and. packed_matrix_format==PM_none)then


  ! Note the format of overlap_matrix is packed - hence we have to set up the 
     if (myid.eq.0) then

        allocate( ovlp_eigenvalues(n_basis), stat=i_info)
        call check_allocation(i_info, 'ovlp_eigenvalues            ')
        allocate( ovlp_transform(n_basis,n_basis), stat=i_info)
        call check_allocation(i_info, 'ovlp_transform              ')
        allocate( sqrt_overlap(n_basis,n_basis), stat=i_info)
        call check_allocation(i_info, 'sqrt_overlap                ')

        call diagonalize_overlap &
          ( n_basis, overlap_matrix, safe_minimum, basis_threshold, &
            n_nonsingular, ovlp_eigenvalues, ovlp_transform, 1 &
          )

        sqrt_overlap = 0.d0
        do i_state = 1, n_nonsingular, 1
          do i_basis_1 = 1, n_basis, 1
            do i_basis_2 = 1, n_basis, 1
                sqrt_overlap(i_basis_2,i_basis_1) = &
                  sqrt_overlap(i_basis_2,i_basis_1) + &
                  ovlp_transform(i_basis_2,i_state) *  &
                  ovlp_transform(i_basis_1,i_state) *  &
                  sqrt(ovlp_eigenvalues(i_state)) 
            enddo
          enddo
        enddo

        do i_spin = 1, n_spin, 1

           call dgemm('N','N',n_basis,n_states,n_basis,1.d0, &
                      sqrt_overlap,n_basis,KS_eigenvector(1,1,i_spin,1), &
                      n_basis, 0.d0,  &
                      ovlp_transform,n_basis)

           do i_state = 1, n_states, 1

              i_index = 0
              do i_basis_1 = 1, n_basis, 1
                    lowdin_decomp( basis_l(i_basis_1), basis_atom(i_basis_1), i_state, i_spin,1 ) = &
                         lowdin_decomp( basis_l(i_basis_1), basis_atom(i_basis_1), i_state, i_spin,1 ) + &
                         ovlp_transform(i_basis_1, i_state) * ovlp_transform(i_basis_1, i_state)
              enddo

           enddo
        enddo
     end if
  else

     ! RJ: When I got it right, the complete code below is for Mulliken, not for Loewdin,
     ! i.e. output_lowdin just doesn't work in this case:

     call aims_stop('output_lowdin needs either scalapack or packed_matrix_format==none!')

     ! RJ: Remove the aims_stop above if the code below makes sense, otherwise the code
     ! should be cleaned up (removed)!

     if(packed_matrix_format /= PM_index)then
        write(use_unit,*) 'Error: periodic Mulliken supports only packed matrix format index'
        return
     end if


     if(real_eigenvectors)then

        do i_spin = 1, n_spin, 1
           i_k = 0           
           do i_k_point = 1, n_k_points
              if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
                 i_k = i_k + 1

                 do i_state = 1, n_states, 1
                    do i_cell = 1,n_cells_in_hamiltonian-1

                       do i_basis_2 = 1, n_basis

                          if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then

                             i_index = index_hamiltonian(1,i_cell, i_basis_2)-1

                             do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)

                                i_index = i_index + 1
                                i_basis_1 =  column_index_hamiltonian(i_index)


    ! 1st pass over all matrix elements
    lowdin_decomp( basis_l(Cbasis_to_basis(i_basis_1)), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k_point ) = &
         lowdin_decomp( basis_l(Cbasis_to_basis(i_basis_1)), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k_point ) + &
         KS_eigenvector(Cbasis_to_basis(i_basis_1), i_state, i_spin,i_k) * &
         overlap_matrix(i_index) * &
         KS_eigenvector(Cbasis_to_basis(i_basis_2), i_state, i_spin,i_k) * dble(k_phase(i_cell,i_k_point))

    ! 2nd pass: must average all off-diagonal matrix elements (but not the diagonal)
    if (i_basis_1.ne.i_basis_2) then
       lowdin_decomp( basis_l(Cbasis_to_basis(i_basis_2)), Cbasis_to_atom(i_basis_2), i_state, i_spin,i_k_point ) = &
            lowdin_decomp( basis_l(Cbasis_to_basis(i_basis_2)), Cbasis_to_atom(i_basis_2), i_state, i_spin,i_k_point ) + &
            KS_eigenvector(Cbasis_to_basis(i_basis_2), i_state, i_spin,i_k) * &
            overlap_matrix(i_index) * &
            KS_eigenvector(Cbasis_to_basis(i_basis_1), i_state, i_spin,i_k) * dble(k_phase(i_cell,i_k_point))

                                end if
                             end do
                          end if
                       end do
                    end do
                 end do
              end if
           end do
        end do


     else


        do i_spin = 1, n_spin, 1
           i_k = 0           

           do i_k_point = 1, n_k_points
             ! write(use_unit,*) i_k_point

              if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
                 i_k = i_k + 1

                 do i_state = 1, n_states, 1
                    do i_cell = 1,n_cells_in_hamiltonian-1

                       do i_basis_2 = 1, n_basis

                          if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then

                             i_index = index_hamiltonian(1,i_cell, i_basis_2)-1

                             do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)

                                i_index = i_index + 1
                                i_basis_1 =  column_index_hamiltonian(i_index)


                                mul_temp =  KS_eigenvector_complex(i_basis_1, i_state, i_spin,i_k) * &
                                     dconjg(KS_eigenvector_complex(i_basis_2, i_state, i_spin,i_k)) &
                                     * dconjg(k_phase(i_cell,i_k_point)) &
                                     * overlap_matrix(i_index)

                                lowdin_decomp( basis_l(i_basis_1), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k_point ) = &
                                  lowdin_decomp( basis_l(i_basis_1), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k_point ) + &
                                  dble(mul_temp)

                                ! 2nd pass: must average all off-diagonal matrix elements (but not the diagonal)
                                if (i_basis_1.ne.i_basis_2) then
                                  lowdin_decomp( basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin,i_k_point ) = &
                                    lowdin_decomp( basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin,i_k_point ) + &
                                    dble(mul_temp)

                                end if
                             end do
                          end if
                       end do
                    end do
                 end do



              end if
           end do

        end do

     end if
  end if

end if ! use_scalapack

  
  call sync_vector( lowdin_decomp(0,1,1,1,1),(1+l_wave_max)*n_atoms*n_states*n_spin*n_k_points)

  if( out_lowdin )then

  ! We have the basic decomposition. Now sum up individual contributions ahead of time.

  call check_occs('output_lowdin', occ_numbers, .false.)

  l_projected_charge = 0.d0
  at_projected_charge = 0.d0
  spin_per_atom = 0.d0

  total_charge = 0.d0

  do i_spin = 1, n_spin, 1
     i_k = 0           
     do i_k_point = 1, n_k_points

        if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
           i_k = i_k + 1


           i_state = n_states
           do while ( (occ_numbers(i_state,i_spin,i_k_point).eq.0.d0) .and. (i_state.ge.1) )
              i_state = i_state - 1
           enddo
           max_occ_number(i_spin,i_k_point) = i_state

           do i_atom = 1, n_atoms, 1

              ! sum up angular-momentum resolved charge projected on atoms in each spin channel
              ! sum up total charge projected on atoms in each spin channel
              do i_l = 0, l_shell_max(species(i_atom))

                 do i_state = 1, max_occ_number(i_spin,i_k_point), 1
                    l_projected_charge( i_l, i_atom, i_spin) = &
                         l_projected_charge( i_l, i_atom, i_spin) + &
                         occ_numbers(i_state,i_spin,i_k_point) * & 
                         lowdin_decomp(i_l, i_atom, i_state, i_spin,i_k_point)*k_weights(i_k_point)

                 enddo

              enddo

           enddo

        end if
     end do
  enddo

  call sync_vector( l_projected_charge(0,1,1),(1+l_wave_max)*n_atoms*n_spin)

  write(info_str,'(2X,A)') &
       "Full analysis will be written to separate file '"//trim(filename)//"'."
  call localorb_info(info_str)
  write(info_str,'(2X,A)') &
       "Summary of the per-atom charge analysis:"
  call localorb_info(info_str)

  do i_spin = 1, n_spin, 1
     do i_atom = 1, n_atoms, 1
        at_projected_charge(i_atom, i_spin) =  sum(l_projected_charge( :, i_atom, i_spin ))
     end do
     total_charge(i_spin) = sum(at_projected_charge(:, i_spin))
  end do

  if (myid.eq.0) then

     charge_difference = 0.d0
     do i_spin = 1, n_spin, 1
        charge_difference = charge_difference + total_charge(i_spin)
     enddo
     do i_atom = 1, n_atoms, 1
        if(species_pseudoized(species(i_atom))) cycle
        if (spin_treatment .eq. 1) then
           spin_per_atom(i_atom) = at_projected_charge(i_atom, 1) - at_projected_charge(i_atom, 2)
        end if
        charge_difference = charge_difference - species_z(species(i_atom))
     enddo

     ! At this point, write detailed Mulliken analysis output

     ! prepare angular momentum output strings
     do i_l = 0, l_wave_max, 1
        write(l_channel(i_l),'(A,I1)') "l=",i_l
     enddo

     ! First, write summary output into standard output file
     write(use_unit,'(2X,A)') "|"
     write(use_unit,'(2X,A,A,3X,A,6X,A,10(9X,A))') &
          "| "," atom","electrons","charge",( l_channel(i_l), i_l=0,l_wave_max,1 )

     do i_atom = 1, n_atoms, 1
        if(species_pseudoized(species(i_atom))) cycle

        write(use_unit,'(2X,A,I5,2X,F10.6,2X,F10.6,10(2X,F10.6))') &
             "| ", i_atom, sum(at_projected_charge(i_atom,1:n_spin)), &
             -( sum(at_projected_charge(i_atom,1:n_spin))-species_z(species(i_atom)) ), &
             (sum(l_projected_charge(i_l,i_atom,1:n_spin)), i_l=0,l_shell_max(species(i_atom)),1 )

     enddo

     write(use_unit,'(2X,A)') "|"

     write(use_unit,'(2X,A,2X,F10.6,2X,F10.6)') "| Total", sum(total_charge(1:n_spin)), -charge_difference

     write(use_unit,*)

     if (spin_treatment .eq. 1) then

        write(info_str,'(2X,A)') &
             "Summary of the per-atom spin analysis:"
        call localorb_info(info_str)

        write(use_unit,'(2X,A)') "|"
        write(use_unit,'(2X,A,A,3X,A,6X,A,10(9X,A))') &
             "| "," atom","spin",( l_channel(i_l), i_l=0,l_wave_max,1 )

        do i_atom = 1, n_atoms, 1
        if(species_pseudoized(species(i_atom))) cycle

           write(use_unit,'(2X,A,I5,2X,F10.6,2X,F10.6,10(2X,F10.6))') &
                "| ", i_atom, spin_per_atom(i_atom), &
                (l_projected_charge(i_l,i_atom,1) - l_projected_charge(i_l,i_atom,2), i_l=0,l_shell_max(species(i_atom)),1 )

        enddo

        write(use_unit,'(2X,A)') "|"

        write(use_unit,'(2X,A,2X,F10.6)') "| Total", sum(spin_per_atom(1:n_atoms))

        write(use_unit,*)
     end if

     ! Now write spin-/state-resolved projection to file
     open(50, file=filename)
     call write_aims_uuid(uuid_str)
     write(50,'(A,2X,A)') '#', uuid_str

     do i_atom = 1, n_atoms, 1
        if(species_pseudoized(species(i_atom))) cycle

        write(50,*)
        write(50,'(A,I5,A)') "Atom number ",i_atom, ":"
        write(50,*)

        do i_spin = 1, n_spin, 1

           if (n_spin.gt.1) then
              write(50,*)
              if (i_spin.eq.1) then
                 write(50,'(2X,A,A)') "Spin channel: ", "up"
              else 
                 write(50,'(2X,A,A)') "Spin channel: ", "down"
              end if
              write(50,*)
           end if

           write(50,'(4X,A,7X,A,2X,A,7X,A,10(9X,A3))') & 
                "State", "eigenvalue", "occ.number", "total", (l_channel(i_l),i_l=0,l_shell_max(species(i_atom)))


           do i_k_point = 1, n_k_points

              do i_state=1,n_states,1
                 write(50,'(2X,I7,2X,F15.5,2X,F10.7,2X,F10.5,10(2X,F10.5))') &
                      i_state, KS_eigenvalue(i_state,i_spin,i_k_point)* hartree, occ_numbers(i_state,i_spin,i_k_point), &
                      sum( lowdin_decomp(0:l_shell_max(species(i_atom)),i_atom,i_state,i_spin,i_k_point) ), &
                      ( lowdin_decomp(i_l,i_atom,i_state,i_spin,i_k_point), i_l=0,l_shell_max(species(i_atom)) )
              enddo
           end do
        enddo

     enddo

     close(50)


  end if      ! end exclusion of all threads other than #0

end if


if(out_l_proj_dos)then

!   VB: In my view, this call to check_norm is incorrect - at the very least,
!   it will now no longer work when two chemical potentials (spin constraint) are employed.
!   call check_norm_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_spin)

   if(.not. allocated (KS_dos)) then
      allocate(KS_dos(0:l_wave_max,n_spin,l_proj_dos_n_en_points, n_species),stat=i_state) 
      call check_allocation(i_state, 'KS_dos                        ')
   endif
   if(.not. allocated (KS_sum_dos)) then
      allocate(KS_sum_dos(n_spin,l_proj_dos_n_en_points, n_species),stat=i_state) 
      call check_allocation(i_state, 'KS_sum_dos                    ')
   endif
   if (.not.allocated(proj_dos_erfs)) then
      allocate(proj_dos_erfs(l_proj_dos_n_en_points),stat=i_state)
      call check_allocation(i_state, 'proj_dos_erfs                 ')
   end if
   
   write(info_str,'(2X,A)') 'Calculating angular momentum projected density of states ...'
   call localorb_info(info_str)
   write(info_str,'(2X,A,F10.6,A)') '| Chemical potential is ', chemical_potential* hartree,' eV.'
   call localorb_info(info_str)
   
   de= (l_proj_dos_high_energy - l_proj_dos_low_energy)/dble(l_proj_dos_n_en_points-1)
   KS_dos =0.d0 
   KS_sum_dos =0.d0 
   
   ! compute DOS for angular momentum projection. Make sure that the mulliken decomposition is properly normalized!!!
   do i_k_point = 1, n_k_points
      if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
         do i_spin = 1, n_spin
            do i_state = 1, n_states
               ! calculate state- and k-dependent prefactors
               do i_e = 1, l_proj_dos_n_en_points
                  en= l_proj_dos_low_energy + dble(i_e-1)*de   
                  E1 = en - de/2d0
                  E2 = en + de/2d0
                  proj_dos_erfs(i_e) = &
                       arch_erf((E2-(KS_eigenvalue(i_state,i_spin,i_k_point))*hartree)*one_over_sqrt2/l_proj_dos_alpha)  &
                       -arch_erf((E1-(KS_eigenvalue(i_state,i_spin,i_k_point))*hartree)*one_over_sqrt2/l_proj_dos_alpha)
               end do
               proj_dos_erfs(:) = proj_dos_erfs(:)*k_weights(i_k_point)/(2d0*dE)
               ! calculate dos projection 
               do i_e = 1, l_proj_dos_n_en_points
                  do i_atom = 1, n_atoms, 1
                     do i_l = 0, l_shell_max(species(i_atom))
                        KS_dos(i_l, i_spin,i_e,species(i_atom)) =  KS_dos (i_l,i_spin,i_e,species(i_atom)) + &
                             lowdin_decomp(i_l,i_atom,i_state,i_spin,i_k_point) * proj_dos_erfs(i_e)
                     end do
                  end do
               end do
            end do
         end do
      end if
   end do
   
   call sync_vector( KS_dos(0,1,1,1),(1+l_wave_max)*n_spin*l_proj_dos_n_en_points*n_species)
   call sync_vector( KS_sum_dos(1,1,1),n_spin*l_proj_dos_n_en_points*n_species)

   do i_spin = 1, n_spin, 1
      do i_e = 1, l_proj_dos_n_en_points, 1
         do i_species = 1, n_species
            do i_l = 0, l_shell_max(i_species)
               KS_sum_dos(i_spin,i_e,i_species) =  KS_sum_dos(i_spin,i_e,i_species) + KS_dos(i_l, i_spin,i_e,i_species)
            end do
         end do
      end do
   end do

   ! output on thread 0 only 
   if (myid.eq.0) then
      if (n_spin.eq.1) then
         do i_species = 1, n_species
            write(proj_dos_filename,'(2A)') trim(species_name(i_species)), '_l_proj_dos.dat'
            write(outputformat,'("(",I2,"(2X,F16.8))")') l_shell_max(i_species)+3
            write(use_unit,'(2X,5A)') '| writing projected DOS (shifted by the chemical potential) for species ',& 
                 trim(species_name(i_species)), &
                 ' to file ',trim(proj_dos_filename),'.'
            open(88, file=proj_dos_filename)
            call write_aims_uuid(uuid_str)
            write(88,'(A,2X,A)') '#', uuid_str

            write(88,'(3A)') '# Angular momentum resolved density for species ',trim(species_name(i_species)), &
                 'as calculated by FHI-aims'
            write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', & 
                 chemical_potential*hartree, ' eV'
            
            write(proj_dos_filename,'(2A)') trim(species_name(i_species)), '_l_proj_dos_raw.dat'
            write(use_unit,'(2X,5A)') '| writing projected DOS (raw data) for species ',trim(species_name(i_species)),& 
                 ' to file ',trim(proj_dos_filename),'.'
            open(89, file=proj_dos_filename)
            call write_aims_uuid(uuid_str)
            write(89,'(A,2X,A)') '#', uuid_str

            write(89,'(3A)') '# Angular momentum resolved density for species ',& 
                 trim(species_name(i_species)), 'as calculated by FHI-aims'
            write(89,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level '
            
            write(outputformat_header,'("(A18,2X,A16,",I2,"(A7,I2,9X))")') l_shell_max(i_species)+1
            write(88,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))
            write(89,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))   
            
            do i_e = 1, l_proj_dos_n_en_points ,1
               en = l_proj_dos_low_energy + dble(i_e -1) *de 
               write(88,outputformat) en-chemical_potential*hartree, KS_sum_dos(1,i_e,i_species), &
                    (KS_dos(i_l, 1,i_e,i_species), i_l=0,l_shell_max(i_species))
               write(89,outputformat) & 
                    en, KS_sum_dos(1,i_e,i_species), (KS_dos(i_l, 1,i_e,i_species), i_l=0,l_shell_max(i_species))
            enddo
            close(88)
            close(89)
         end do
      else 
         do i_species = 1, n_species
            write(proj_dos_filename,'(2A)') trim(species_name(i_species)), '_l_proj_dos_spin_up.dat'
            write(use_unit,'(2X,5A)') '| writing spin-up projected DOS (shifted by the chemical potential) for species ',&
                 trim(species_name(i_species)),' to file ',trim(proj_dos_filename),'.'
            open(88, file=proj_dos_filename)
            call write_aims_uuid(uuid_str)
            write(88,'(A,2X,A)') '#', uuid_str

            write(proj_dos_filename,'(2A)') trim(species_name(i_species)), '_l_proj_dos_spin_down.dat'
            write(use_unit,'(2X,5A)') '| writing spin-down projected DOS (shifted by the chemical potential) for species ',&
                 trim(species_name(i_species)),' to file ',trim(proj_dos_filename),'.'
            open(89, file=proj_dos_filename)           
            call write_aims_uuid(uuid_str)
            write(89,'(A,2X,A)') '#', uuid_str

            write(88,'(3A)') '# spin-up angular momentum projected density of states for species ',& 
                 species_name(i_species), ' as calculated by FHI-aims '
            write(89,'(3A)') '# spin-down angular momentum projected density of states for species ',& 
                 species_name(i_species), ' as calculated by FHI-aims '
            
            write(proj_dos_filename,'(2A)') trim(species_name(i_species)), '_l_proj_dos_spin_up_raw.dat'
            write(use_unit,'(2X,5A)') '| writing spin-up projected DOS (raw data) for species ',trim(species_name(i_species)),&
                 ' to file ',trim(proj_dos_filename),'.'
            open(90, file=proj_dos_filename)
            call write_aims_uuid(uuid_str)
            write(90,'(A,2X,A)') '#', uuid_str

            write(proj_dos_filename,'(2A)') trim(species_name(i_species)), '_l_proj_dos_spin_down_raw.dat'
            write(use_unit,'(2X,5A)') '| writing spin-down projected DOS (raw data) for species ',trim(species_name(i_species)),&
                 ' to file ',trim(proj_dos_filename),'.'
            open(91, file=proj_dos_filename)                      
            call write_aims_uuid(uuid_str)
            write(91,'(A,2X,A)') '#', uuid_str

            write(outputformat,'("(",I2,"(2X,F16.8))")') l_shell_max(i_species)+3
            write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ',& 
                 chemical_potential*hartree,' eV'
            write(89,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ',& 
                 chemical_potential*hartree,' eV'
            write(90,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
            write(91,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
            
            write(outputformat_header,'("(A18,2X,A16,",I2,"(A7,I2,9X))")') l_shell_max(i_species)+1
            write(88,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))
            write(89,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))   
            write(90,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))
            write(91,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))   
            
            do i_e = 1, l_proj_dos_n_en_points ,1
               en = l_proj_dos_low_energy + dble(i_e -1) *de 
               write(88,outputformat) en-chemical_potential*hartree, KS_sum_dos(1,i_e,i_species), &
                    (KS_dos(i_l, 1,i_e,i_species), i_l=0,l_shell_max(i_species))
               write(89,outputformat) en-chemical_potential*hartree, KS_sum_dos(2,i_e,i_species), &
                    (KS_dos(i_l, 2,i_e,i_species), i_l=0,l_shell_max(i_species))
               write(90,outputformat) & 
                    en, KS_sum_dos(1,i_e,i_species), (KS_dos(i_l, 1,i_e,i_species), i_l=0,l_shell_max(i_species))
               write(91,outputformat) & 
                    en, KS_sum_dos(2,i_e,i_species), (KS_dos(i_l, 2,i_e,i_species), i_l=0,l_shell_max(i_species))
            enddo
            close(89)
            close(88)  
            close(90)
            close(91)
         end do
      end if
      
   end if
   if (allocated(KS_dos)       ) deallocate(KS_dos)        
   if (allocated(KS_sum_dos)   ) deallocate(KS_sum_dos)    
   if (allocated(proj_dos_erfs)) deallocate(proj_dos_erfs)
end if

! use same multipole decomposition to calculate atom-projected density of states if desired ... 
if(out_atom_dos)then

!   VB: In my view, this call to check_norm is incorrect - at the very least,
!   it will now no longer work when two chemical potentials (spin constraint) are employed.
!   call check_norm_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_spin)
   
   call localorb_info(" ")
   write(info_str,'(2X,A)') 'Calculating atom-projected density of states ...'
   call localorb_info(info_str)
   write(info_str,'(2X,A,F10.6,A)') '| Chemical potential is ', chemical_potential* hartree,' eV.'
   call localorb_info(info_str)
   
   if(.not. allocated (KS_dos)) then
      allocate(KS_dos(0:l_wave_max,n_spin,atom_dos_n_en_points, n_atoms),stat=i_state) 
      call check_allocation(i_state, 'KS_dos                        ')
   endif
   if(.not. allocated (KS_sum_dos)) then
      allocate(KS_sum_dos(n_spin,atom_dos_n_en_points, n_atoms),stat=i_state) 
      call check_allocation(i_state, 'KS_sum_dos                    ')
   endif
   if (.not. allocated (proj_dos_erfs)) then
      allocate(proj_dos_erfs(atom_dos_n_en_points),stat=i_state)
      call check_allocation(i_state, 'proj_dos_erfs                 ')
   end if
   
   de= (atom_dos_high_energy - atom_dos_low_energy)/dble(atom_dos_n_en_points-1)
   KS_dos =0.d0 
   do i_k_point = 1, n_k_points
      if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
         do i_spin = 1, n_spin
            do i_state = 1, n_states
               do i_e = 1, atom_dos_n_en_points
                  en= atom_dos_low_energy + dble(i_e-1)*de   
                  E1 = en - de/2d0
                  E2 = en + de/2d0
                  proj_dos_erfs(i_e) = &
                        arch_erf((E2-(KS_eigenvalue(i_state,i_spin,i_k_point))*hartree)*one_over_sqrt2/atom_dos_alpha)  &
                       -arch_erf((E1-(KS_eigenvalue(i_state,i_spin,i_k_point))*hartree)*one_over_sqrt2/atom_dos_alpha)  
               end do
               proj_dos_erfs(:) = proj_dos_erfs(:)*k_weights(i_k_point)/(2d0*dE)
               do i_e = 1, atom_dos_n_en_points
                  do i_atom = 1, n_atoms
        if(species_pseudoized(species(i_atom))) cycle
                     do i_l = 0, l_shell_max(species(i_atom))
                        KS_dos(i_l, i_spin,i_e,i_atom) =  KS_dos (i_l,i_spin,i_e,i_atom) + &
                             lowdin_decomp(i_l,i_atom,i_state,i_spin,i_k_point) * proj_dos_erfs(i_e)
                     end do
                  end do
               end do
            end do
         end do
      end if
   end do
  

   call sync_vector( KS_dos(0,1,1,1),(1+l_wave_max)*n_spin*atom_dos_n_en_points*n_atoms)
   call sync_vector( KS_sum_dos(1,1,1),n_spin*atom_dos_n_en_points*n_atoms)

   KS_sum_dos =0.d0 
   do i_spin = 1, n_spin
      do i_e = 1, atom_dos_n_en_points
         do i_atom = 1, n_atoms
            if(species_pseudoized(species(i_atom))) cycle
            do i_l = 0, l_shell_max(species(i_atom))
               KS_sum_dos(i_spin,i_e,i_atom) =  KS_sum_dos(i_spin,i_e,i_atom) + KS_dos(i_l, i_spin,i_e,i_atom)
            end do
         end do
      end do
   end do
   
   if (myid.eq.0) then
      if (n_spin.eq.1) then
         do i_atom = 1, n_atoms
        if(species_pseudoized(species(i_atom))) cycle
            
            if (i_atom.lt.10) then
               write(proj_dos_filename,'(3A,I1,A4)') 'atom_projected_dos_',trim(species_name(species(i_atom))),'000',i_atom,'.dat'
            else if (i_atom.lt.100) then
               write(proj_dos_filename,'(3A,I2,A4)') 'atom_projected_dos_',trim(species_name(species(i_atom))),'00',i_atom,'.dat'
            else if (i_atom.lt.1000) then
               write(proj_dos_filename,'(3A,I3,A4)') 'atom_projected_dos_',trim(species_name(species(i_atom))),'0',i_atom,'.dat'
            else
               write(proj_dos_filename,'(2A,I4,A4)') 'atom_projected_dos_',trim(species_name(species(i_atom))),i_atom,'.dat'
            end if
            
            write(outputformat,'("(",I2,"(2X,F16.8))")') l_shell_max(species(i_atom))+3
            write(use_unit,'(2X,5A)') '| writing projected DOS (shifted by the chemical potential) for species ',& 
                 trim(species_name(species(i_atom))),' to file ',trim(proj_dos_filename),'.'
            open(88, file=proj_dos_filename)
            call write_aims_uuid(uuid_str)
            write(88,'(A,2X,A)') '#', uuid_str

            write(88,'(3A)') '# Angular momentum resolved density for species ',& 
                 trim(species_name(species(i_atom))), 'as calculated by FHI-aims'
            write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', & 
                 chemical_potential*hartree, ' eV'
            
            
            if (i_atom.lt.10) then
               write(proj_dos_filename,'(3A,I1,A8)') 'atom_proj_dos_',trim(species_name(species(i_atom))),'000',i_atom,'_raw.dat'
            else if (i_atom.lt.100) then
               write(proj_dos_filename,'(3A,I2,A8)') 'atom_proj_dos_',trim(species_name(species(i_atom))),'00',i_atom,'_raw.dat'
            else if (i_atom.lt.1000) then
               write(proj_dos_filename,'(3A,I3,A8)') 'atom_proj_dos_',trim(species_name(species(i_atom))),'0',i_atom,'_raw.dat'
            else
               write(proj_dos_filename,'(2A,I4,A8)') 'atom_proj_dos_',trim(species_name(species(i_atom))),i_atom,'_raw.dat'
            end if
            
            write(use_unit,'(2X,5A)') '| writing projected DOS (raw data) for species ',& 
                 trim(species_name(species(i_atom))),' to file ',trim(proj_dos_filename),'.'
            open(89, file=proj_dos_filename)
            call write_aims_uuid(uuid_str)
            write(89,'(A,2X,A)') '#', uuid_str

            write(89,'(3A)') '# Angular momentum resolved density for species ',& 
                 trim(species_name(species(i_atom))), 'as calculated by FHI-aims'
            write(89,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level '
            
            write(outputformat_header,'("(A18,2X,A16,",I2,"(A7,I2,9X))")') l_shell_max(species(i_atom))+1
            write(88,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))
            write(89,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))   
            
            do i_e = 1, atom_dos_n_en_points ,1
               en = atom_dos_low_energy + dble(i_e -1) *de 
               write(88,outputformat) en-chemical_potential*hartree, KS_sum_dos(1,i_e,i_atom), &
                    (KS_dos(i_l, 1,i_e,i_atom), i_l=0,l_shell_max(species(i_atom)))
               write(89,outputformat) en, KS_sum_dos(1,i_e,i_atom), (KS_dos(i_l, 1,i_e,i_atom), & 
                    i_l=0,l_shell_max(species(i_atom)))
            enddo
            close(88)
            close(89)
         end do
      else 
         do i_atom = 1, n_atoms
            if(species_pseudoized(species(i_atom))) cycle
            
            if (i_atom.lt.10) then
               write(proj_dos_filename,'(3A,I1,A4)') 'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),'000',i_atom,'.dat'
            else if (i_atom.lt.100) then
               write(proj_dos_filename,'(3A,I2,A4)') 'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),'00',i_atom,'.dat'
            else if (i_atom.lt.1000) then
               write(proj_dos_filename,'(3A,I3,A4)') 'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),'0',i_atom,'.dat'
            else
               write(proj_dos_filename,'(2A,I4,A4)') 'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),i_atom,'.dat'
            end if
            
            write(use_unit,'(2X,5A)') '| writing spin-up projected DOS (shifted by the chemical potential) for species ',&
                 trim(species_name(species(i_atom))),' to file ',trim(proj_dos_filename),'.'
            open(88, file=proj_dos_filename)
            call write_aims_uuid(uuid_str)
            write(88,'(A,2X,A)') '#', uuid_str

            
            if (i_atom.lt.10) then
               write(proj_dos_filename,'(3A,I1,A4)') 'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),'000',i_atom,'.dat'
            else if (i_atom.lt.100) then
               write(proj_dos_filename,'(3A,I2,A4)') 'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),'00',i_atom,'.dat'
            else if (i_atom.lt.1000) then
               write(proj_dos_filename,'(3A,I3,A4)') 'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),'0',i_atom,'.dat'
            else
               write(proj_dos_filename,'(2A,I4,A4)') 'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),i_atom,'.dat'
            end if
            
            write(use_unit,'(2X,5A)') '| writing spin-down projected DOS (shifted by the chemical potential) for species ',&
                 trim(species_name(species(i_atom))),' to file ',trim(proj_dos_filename),'.'
            open(89, file=proj_dos_filename)           
            call write_aims_uuid(uuid_str)
            write(89,'(A,2X,A)') '#', uuid_str

            write(88,'(3A)') '# spin-up angular momentum projected density of states for species ',& 
                 species_name(species(i_atom)), ' as calculated by FHI-aims '
            write(89,'(3A)') '# spin-down angular momentum projected density of states for species ',& 
                 species_name(species(i_atom)), ' as calculated by FHI-aims '
            
            if (i_atom.lt.10) then
               write(proj_dos_filename,'(3A,I1,A8)') & 
                    'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),'000',i_atom,'_raw.dat'
            else if (i_atom.lt.100) then
               write(proj_dos_filename,'(3A,I2,A8)') & 
                    'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),'00',i_atom,'_raw.dat'
            else if (i_atom.lt.1000) then
               write(proj_dos_filename,'(3A,I3,A8)') & 
                    'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),'0',i_atom,'_raw.dat'
            else
               write(proj_dos_filename,'(2A,I4,A8)') & 
                    'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),i_atom,'_raw.dat'
            end if
            
            write(use_unit,'(2X,5A)') &
               '| writing spin-up projected DOS (raw data) for species ', &
               trim(species_name(species(i_atom))), &
               ' to file ',trim(proj_dos_filename),'.'
            open(90, file=proj_dos_filename)
            call write_aims_uuid(uuid_str)
            write(90,'(A,2X,A)') '#', uuid_str

            if (i_atom.lt.10) then
               write(proj_dos_filename,'(3A,I1,A8)') & 
                    'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),'000',i_atom,'_raw.dat'
            else if (i_atom.lt.100) then
               write(proj_dos_filename,'(3A,I2,A8)') & 
                    'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),'00',i_atom,'_raw.dat'
            else if (i_atom.lt.1000) then
               write(proj_dos_filename,'(3A,I3,A8)') & 
                    'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),'0',i_atom,'_raw.dat'
            else
               write(proj_dos_filename,'(2A,I4,A8)') &
                    'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),i_atom,'_raw.dat'
            end if
            
            write(use_unit,'(2X,5A)') &
               '| writing spin-down projected DOS (raw data) for species ', &
               trim(species_name(species(i_atom))), &
               ' to file ',trim(proj_dos_filename),'.'
            open(91, file=proj_dos_filename)                      
            call write_aims_uuid(uuid_str)
            write(91,'(A,2X,A)') '#', uuid_str

            write(outputformat,'("(",I2,"(2X,F16.8))")') l_shell_max(species(i_atom))+3
            write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', &
                 chemical_potential*hartree, ' eV'
            write(89,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', &
                 chemical_potential*hartree, ' eV'
            write(90,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
            write(91,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
            
            write(outputformat_header,'("(A18,2X,A16,",I2,"(A7,I2,9X))")') l_shell_max(species(i_atom))+1
            write(88,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))
            write(89,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))   
            write(90,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))
            write(91,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))   
            
            do i_e = 1, atom_dos_n_en_points ,1
               en = atom_dos_low_energy + dble(i_e -1) *de 
               write(88,outputformat) en-chemical_potential*hartree, KS_sum_dos(1,i_e,i_atom), (KS_dos(i_l, 1,i_e,i_atom), &
                    i_l=0,l_shell_max(species(i_atom)))
               write(89,outputformat) en-chemical_potential*hartree, KS_sum_dos(2,i_e,i_atom), (KS_dos(i_l, 2,i_e,i_atom), &
                    i_l=0,l_shell_max(species(i_atom)))
               write(90,outputformat) en, KS_sum_dos(1,i_e,i_atom), (KS_dos(i_l, 1,i_e,i_atom), i_l=0,l_shell_max(species(i_atom)))
               write(91,outputformat) en, KS_sum_dos(2,i_e,i_atom), (KS_dos(i_l, 2,i_e,i_atom), i_l=0,l_shell_max(species(i_atom)))
            enddo
            close(89)
            close(88) 
            close(90)
            close(91)
         end do
      end if
      
   end if  ! output on master thread

  if (allocated (KS_dos)        ) deallocate(KS_dos)
  if (allocated (KS_sum_dos)    ) deallocate(KS_sum_dos)
  if (allocated (proj_dos_erfs) ) deallocate(proj_dos_erfs)

end if   ! out_atom_dos

end subroutine output_lowdin
!******
