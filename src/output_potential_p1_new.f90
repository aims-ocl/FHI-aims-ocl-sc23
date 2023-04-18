!****s* FHI-aims/output_potential_p1
!  NAME
!    output_potential_p1
!  SYNOPSIS

subroutine output_potential_p1_new( local_potential_parts, mode )

!  PURPOSE
!  Subroutine output_potential
!
!  Writes the effective potential on the integration grid. In practise it prints out the Hartree potential part.
!
!  Note: This is not the effective potential for  a real fix, we must evaluate the 
!  xc potential right here so that we can actually write out the relevant bits.
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use localorb_io, only: use_unit
  use mpi_tasks
  use mpi_utilities
  use species_data

  use directories
  implicit none

!  ARGUMENTS

  real*8, dimension(n_spin, n_full_points) :: local_potential_parts

!  INPUTS
!   o local_potential_parts -- potential which is printed out.
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

  real*8  :: grid_coord(3)

!  counters

  integer              :: i_coord
  integer              :: i_spin
  integer              :: i_point, i_my_batch, i_index, i_atom
  character*(*)        :: mode
  character( len = 8 ) :: myid_string, number_string
  integer, save        :: consistent_counter, iteration_counter, &
                          initialize_counter
  integer              :: write_once(n_atoms)
  real*8               :: x_atom, y_atom, z_atom
  integer              :: n_radial_of_atom, n_angular_of_atom
  integer              :: actuall_index_radial, actuall_index_angular
  logical              :: file_exists
  integer              :: n_fields =1 

  logical, parameter   :: debug=.false.
  
  logical              :: dirSuccess, doWrite
  integer,save         :: record, i
  integer              :: iolengthInt, iolengthChar, iolengthReal8, &
                          reclen, recbytes, ndoubles
  character            :: c
  real*8               :: r

  call get_my_task()

  doWrite = .false.

  ! create directory structure (if it does not allready exist)
  dirSuccess = create_dir("output/")
  print *,"creating output/"
  
  if (.not.(dirSuccess)) then
     print *,"could not create directory output"
     stop
  endif
  
  dirSuccess = create_dir("output/"//mode)
  
  if (.not.(dirSuccess)) then
     print *,"could not create directory output/"//mode
     stop
  endif
  

  print *,"creating output/"//mode


  !  begin work

  write_once(:) = 0
  
  ! convert mpi task id to a string
  write(unit=myid_string, fmt='(i8.8)') myid 
  
  ! set counters; convert them to string
  if (mode .eq. "consistent") then
     consistent_counter = consistent_counter + 1
     write(unit=number_string, fmt='(i8.8)') consistent_counter
  endif
  
  if (mode .eq. "iteration") then
     iteration_counter = iteration_counter + 1
     write(unit=number_string, fmt='(i8.8)') iteration_counter
  endif
  
  if (mode .eq. "initialize") then
     initialize_counter = initialize_counter + 1
     write(unit=number_string, fmt='(i8.8)') initialize_counter
  endif
  

  ! check whether the respective counter satisfies the writting condition

  if (mode .eq. "consistent") then
     
     if (out_nconsistent .eq. -1) then
        ! never write consistent steps
        doWrite = .false.
     else
        if (mod(consistent_counter,out_nconsistent) .eq. 0) then
           ! write consistent step
           doWrite = .true.
        endif
     endif
  endif

  if (mode .eq. "iteration") then
     
     if (out_niteration .eq. -1) then
        ! never write iteration steps
        doWrite = .false.
     else
        if (mod(iteration_counter,out_niteration) .eq. 0) then
           ! write consistent step
           doWrite = .true.
        endif
     endif
  endif

  
  if (mode .eq. "initialize") then
     
     if (out_ninitialize .eq. -1) then
        ! never write initialization steps
        doWrite = .false.
     else
        if (mod(initialize_counter,out_ninitialize) .eq. 0) then
           ! write consistent step
           doWrite = .true.
        endif
     endif
  endif

  
  ! message for debug
  if (debug) then
     write(use_unit,'(2X,A,A)') &
          "Writing the effective potential at each grid point", &
          " to file"," ./output/"//mode//"/v_eff_"//number_string//"_rank"//myid_string//".dat"
  endif



  ! only do sth. if flag is set
  if (doWrite) then

     inquire(file="./output/"//mode//"/v_eff_"//number_string//"_rank"//myid_string//".dat",exist=file_exists)

     ! write ascii output
     if(out_ascii) then
        
        if (.not.file_exists) then
           
           open (50, file="./output/"//mode//"/v_eff_"//number_string//"_rank"//myid_string//".dat", &
                form='formatted', status='new')
           write(50,'(a,i0)') "ATOMS ",n_atoms
           write(50,'(a,i0)') "FIELDS ",n_fields
        else
           open (50, file="./output/"//mode//"/v_eff_"//number_string//"_rank"//myid_string//".dat", &
                form='formatted', position="append")
        endif
        
     ! write binary output
     else

        ! number of doubles per line/record (at least 4)
        ndoubles = 6
        
        ! ### determination of the record length based on actual data type sizes ###
        inquire(iolength=iolengthInt) i
        inquire(iolength=iolengthChar) c
        inquire(iolength=iolengthReal8) r
        !
        ! units are bytes for the GNU compiler, and others for the Intel Compiler!!!
        ! write(use_unit,*) iolengthInt, iolengthChar, iolengthReal8
        !
        ! ==> the following should be portable
        reclen = iolengthReal8*ndoubles
        recbytes = 8*ndoubles
        
        if (.not.file_exists) then
           
           open(50, file="./output/"//mode//"/v_eff_"//number_string &
              //"_rank"//myid_string//".dat", status='new', &
              form='unformatted', access='direct', recl=reclen)
           record=1    ! first record contains the record length in bytes
           write(50, rec=record) "RECL ", recbytes
           record=record+1  ! all data which follows is written in the same way as it is done using ASCII files
           write(50, rec=record) "ATOMS ", n_atoms
           record=record+1 
           write(50, rec=record) "FIELDS ", n_fields
        else
           open(50, file="./output/"//mode//"/v_eff_"//number_string &
                 //"_rank"//myid_string//".dat", status='old', &
                 form='unformatted', access='direct', recl=reclen)
        endif
     endif

     
        
     !     calculate potential, atom by atom, and point by point 
     do i_atom = 1, n_atoms
        
        i_point = 0
        
        ! loop over all batches
        do i_my_batch = 1, n_my_batches, 1
           
           ! loop over all points in a batch
           do i_index = 1, batches(i_my_batch)%size, 1
              
              i_point = i_point + 1
              
              ! calculate grid point coordinate
              grid_coord(:) = batches(i_my_batch) % points(i_index) % coords(:)*bohr
              x_atom = coords(1,i_atom)*bohr
              y_atom = coords(2,i_atom)*bohr
              z_atom = coords(3,i_atom)*bohr
              
              ! we want to know how many radial and angular points
              ! are in this grid for this atom
              
              actuall_index_radial = batches(i_my_batch)%points(i_index)%index_radial
              actuall_index_angular= batches(i_my_batch)%points(i_index)%index_angular
              
              n_radial_of_atom = n_radial(species(i_atom))
              n_angular_of_atom = n_angular(actuall_index_radial,species(i_atom))
              
              ! write potential
              
              
              ! write ascii output
              if(out_ascii) then
                 
                 ! first write the header of the actual species
                 if (i_atom .eq. batches(i_my_batch)%points(i_index)%index_atom .and. &
                      write_once(i_atom) .eq. 0) then
                    
                    write(50,'(a,i0)') "ATOM ",i_atom
                    write(50,'(a,a)') "SPECIES ",trim(species_name(species(i_atom)))
                    write(50,'(a,3(F20.10,1x))') "CENTER ",x_atom, y_atom, z_atom
                    write_once(i_atom) = 1

                    write(50,'(a)') "FIELDS potential"
                 endif
                 
                 if (i_atom .eq. batches(i_my_batch)%points(i_index)%index_atom) then
                    
                    write(50,'(3(F20.10,1x),i7,1x,i7,1x,F20.10)') &
                         !!'(7(F20.10,1X))') &
                         (grid_coord(i_coord), i_coord=1,3,1), &
                         batches(i_my_batch)%points(i_index)%index_radial, &
                         batches(i_my_batch)%points(i_index)%index_angular, &
                         ( local_potential_parts(i_spin,i_point), &
                         i_spin = 1, n_spin, 1 )
!, ( local_potential_parts(i_spin,i_point), &
!                         i_spin = 1, n_spin, 1 )
                 endif

 
              else
                 ! first write the header of the actual species
                 if (i_atom .eq. batches(i_my_batch)%points(i_index)%index_atom .and. &
                      write_once(i_atom) .eq. 0) then
                    record = record+1
                    write(50, rec=record) "ATOM ", i_atom
                    record = record+1
                    write(50, rec=record) "SPECIES ", trim(species_name(species(i_atom)))
                    record = record+1
                    write(50, rec=record) "CENTER ",x_atom, y_atom, z_atom
                    record = record+1
                    write(50, rec=record) "FIELDS potential"
                    write_once(i_atom) = 1
                 endif
                 
                 
                 if (i_atom .eq. batches(i_my_batch)%points(i_index)%index_atom) then
                    record = record+1
!                    write(50, rec=record) grid_coord(1),grid_coord(2),grid_coord(3),batches(i_my_batch)%points(i_index)%index_radial,batches(i_my_batch)%points(i_index)%index_angular,local_potential_parts(1,i_point),local_potential_parts(1,i_point)
                    write(50, rec=record) (grid_coord(i_coord), i_coord=1,3,1),  &
                         batches(i_my_batch)%points(i_index)%index_radial, &
                         batches(i_my_batch)%points(i_index)%index_angular,&
                         ( local_potential_parts(i_spin,i_point), &
                         i_spin = 1, n_spin, 1 )
                    
                 endif
                 
                 
              endif
              
              
              ! end loop over batch points
           enddo
           ! end loop over batches
        enddo
        ! end loop over atoms   
     enddo
     
     close(50)
     
     
  endif ! doWrite
  

  return
end subroutine output_potential_p1_new
!----------------------------------------------------------------------
!******	
