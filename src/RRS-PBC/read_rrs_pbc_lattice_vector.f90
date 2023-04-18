!****s* FHI-aims/RRS-PBC/read_rrs_pbc_lattice_vector
!  NAME
!   read_rrs_pbc_lattice_vector
!  SYNOPSIS

    subroutine read_rrs_pbc_lattice_vector(fn,iop)
!  USES
    use dimensions
    use runtime_choices
    use localorb_io
    use geometry
    use bravais
    use numerical_utilities
    use mpi_tasks, only: aims_stop, aims_stop_coll

    implicit none

!  PURPOSE
!   The subroutine reads the lattice vector from control.in for the 
!   RRS-PBC method.
!   This subroutine is called from read_control and dimension
!  INPUTS
!    fn   :: the file flow id of control.in
!    iop  :: 0 for actual cell_rrs_pbc loading
!            other value for check only
!  OUTPUT
!    none
!  AUTHOR
!                                                                  
!  SEE ALSO
!    
!  COPYRIGHT
!   
!  HISTORY
!    
!  SOURCE

    integer :: fn, iop
    integer :: i, j, rank
    character*130 :: inputline
    character*30 desc_str
    integer :: i_code
    character(*), parameter :: func = 'read_rrs_pbc_unit_cell'

     do i=1, 3, 1
         read(fn,'(A)',iostat=i_code) inputline
         if(i_code .ne. 0) then
            call aims_stop_coll("Unknown error reading file 'control.in'...", func)
         endif
         read(inputline,*,iostat=i_code) desc_str
         if(i_code/=0 .or. desc_str(1:1).eq.'#') then
            call aims_stop_coll( &
                "Error: empty and comment lines are forbidden following the keyword of rrs_pbc_lattice_vector" &
                , func)
         endif

         if (iop.eq.0) then
             read(inputline,*) (rrs_pbc_lattice_vector(j,i), j=1,3,1)
             do j = 1,3,1
                 rrs_pbc_lattice_vector(j,i) = &
                 rrs_pbc_lattice_vector(j,i) / bohr
             enddo
         endif
     enddo

     ! calculate rrs_pbc_recip_lattice_vector, rrs_pbc_inv_lattice_vector,
     ! rrs_pbc_cell_volume rrs_pbc_length_of_lattice_vector and so on.
     if (iop.eq.0) then
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
     endif
     end subroutine
