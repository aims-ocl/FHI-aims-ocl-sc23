!****s* FHI-aims/RRS-PBC/parse_rrs_pbc()
!  NAME
!    parse_rrs_pbc() 
!  SYNOPSIS

    subroutine parse_rrs_pbc()

!  PURPOSE
!  High-level wrapper around the post RRS-PBC projection depending on the
!  cluster calculation.
!  This routine can only be called after a succeful scf procedure. Necessary
!  conditions for running are:
!  * an converged overlap matrix
!  * an converged hamiltonian matrix,
!  and all of that with the correct array dimensions
!
!  USES

      use localorb_io
      use species_data
      use dimensions
      use physics
      use geometry
      use basis
      use bravais
      use runtime_choices
      use numerical_utilities
      use mpi_tasks
      implicit none

!  ARGUMENTS

!  INPUTS
!    none
!  OUTPUT
!    none
!  COPYRIGHT
!   
!   
!   
!  HISTORY
!   
!  SOURCE
!

      ! imported variables

      ! local variables
      character*132 :: rrs_pbc_info
      character*132 :: func = 'parse_rrs_pbc()'
      integer :: i,i_center
      integer :: info, rank
      real*8, dimension(3)   :: tmp_vec
      real*8  :: tmp_n_electron

      ! parse the lattice vector
      call get_cell_volume(rrs_pbc_lattice_vector, rrs_pbc_cell_volume)
      call get_reciprocal_vectors(3, rrs_pbc_lattice_vector, &
                                  rrs_pbc_recip_lattice_vector)
      call get_length_of_lv(3, rrs_pbc_lattice_vector, &
                            rrs_pbc_length_of_lattice_vector)
      call pseudo_inverse(func, 3, 3, rrs_pbc_lattice_vector, &
                          rrs_pbc_inv_lattice_vector, 1d-10, rank)
      if(rank /= 3) then
         call aims_stop('ERROR: rrs_pbc_lattice_vector is singular!')
      end if
      call output_rrs_pbc_lattice_vector()

      ! parse electronic structure info. in the unit cell
      write(rrs_pbc_info,'(4X,A,I12)') &
          '| Atom number in unit cell             :',rrs_pbc_n_center_atom
      call localorb_info(rrs_pbc_info)

      ! now count rrs_pbc_n_electron
      tmp_n_electron = 0.d0
      do i_center = 1, rrs_pbc_n_center_atom, 1
          tmp_n_electron = tmp_n_electron + &
              species_z(species(rrs_pbc_center_atom(1,i_center)))
      enddo

      rrs_pbc_n_electron(1)     = (tmp_n_electron + n_spin - 1) / 2
      rrs_pbc_n_electron(2)     = (tmp_n_electron - n_spin + 1) / 2
      rrs_pbc_n_electron_int(1) = int((tmp_n_electron + n_spin - 1) / 2)
      rrs_pbc_n_electron_int(2) = int((tmp_n_electron - n_spin + 1) / 2)

      write(rrs_pbc_info,'(4X,A,2(F12.3))') &
          '| Associated electron number           :',&
          (rrs_pbc_n_electron(i), i=1,2,1)
      call localorb_info(rrs_pbc_info)


      ! parse the basis table in rrs_pbc_center_atom and rrs_pbc_equal_atom
      call get_rrs_pbc_atom2basis()

      call output_rrs_pbc_unit_cell()

      do i=1,rrs_pbc_n_equal,1
          call get_rrs_pbc_cell_index(1,i,tmp_vec)
      enddo

      ! very rough algorithm, might not be suitable for complecated spin state,
      ! make note here for futhur improvement.
      rrs_pbc_occ_num(:,:) = 0
      if (n_spin .eq. 1) then
          do i = 1, int(tmp_n_electron)/2, 1
              rrs_pbc_occ_num(i,1) = 2.0d0
          enddo
      endif

      ! parse the k_point_list
      call get_rrs_pbc_k_point_list()

      ! parse the rationality of the running job for RRS-PBC scheme
      call parse_rrs_pbc_rationality()
      
      end subroutine

