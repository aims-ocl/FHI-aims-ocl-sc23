!****s* FHI-aims/output_bands_during_scf
!  NAME
!   read_plot_band
!  SYNOPSIS

  subroutine output_bands_during_scf
!  USES

    use dimensions
    use runtime_choices
    use localorb_io
    use mpi_tasks
    use pbc_lists, only : k_point_list
    use physics

!  PURPOSE
!
!  The purpose of this routine is to 
!  - take the normal k-point grid that we use for Brillouin zone integrals in the s.c.f. cycle
!  - identify any and all k-points that are located on a specific line in k-space 
!    (for example, a high-symmetry line for band structure output)
!  - and output the Kohn-Sham eigenvalues at these k-points (the band structure) in the same format
!    which we would use for the band structure output.
!
!  By default, we only do this after the s.c.f. cycle is converged and after any post-processing
!  corrections have been added (scaled ZORA!), but the same functionality can also be used
!  during the s.c.f. cycle, for example to monitor the convergence of a given calculation.
!
!  INPUTS
!    none
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
!    Release version, FHI-aims (2012).
!  SOURCE

  implicit none

  integer :: i_band, i_coord, i_k_point
  real*8, dimension(3) :: band_dir
  real*8 :: band_length
  real*8, dimension(3) :: k_direction
  real*8 :: k_distance, scalar_product, difference

  integer, dimension(n_plot_band_scf) :: n_points_on_band
  integer, dimension(n_k_points,n_plot_band_scf) :: k_point_index
  real*8,  dimension(n_k_points,n_plot_band_scf) :: k_distance_on_band

  integer, dimension(n_k_points) :: old_k_index
  real*8,  dimension(n_k_points) :: new_k_distance_on_band
  integer, dimension(:),allocatable :: temp_k_index
  integer :: i_point_on_band

  character*50 :: file_name, file_name2
  character*80 :: info_str

  integer :: i_state

   ! This routine is not yet functional. What needs to be done:

   ! Loop over the requested bands - their number is given by
   ! n_plot_band_scf (determined in dimensions.f90)

   ! For each band i_band = 1, n_plot_band_scf
   ! (beginning and end points are specified by band_scf_begin and band_scf_end,
   ! declared in runtime_choices and read in read_plot_band_during_scf)

     ! Go through all k-points (i_k_point = 1, n_k_points). 

       ! Check if this k point lies on the line between band_scf_begin(i_band) and band_scf_end(i_band)

       ! If so, calculate its (one-dimensional!) distance from band_scf_begin by projection
       ! write this k point to a list
       ! store the eigenvalues at this k point in a list also

     ! after we are done with all k points for this band:

     ! From all CPUs, collect all k-points in this direction in a single array on
     ! CPU number zero (Note: k points are distributed on different CPUs. Once we
     ! are done finding them on every single CPUs, we only need to output them on
     ! a single CPU, CPU number zero. This CPU is identified by (myid == 0) .

     ! sort all k points in this band according to their distance from band_scf_begin

     ! write the list of identified k-points into a file in the same format as that
     ! used by the normal band structure output. The file name should be different, though.

   ! next band

   do i_band = 1, n_plot_band_scf, 1

      band_dir(:) = band_scf_end(i_band,:) - band_scf_begin(i_band,:)

      band_length = 0.d0
      do i_coord = 1,3,1
         band_length = band_length + band_dir(i_coord)*band_dir(i_coord)
      enddo
      band_length = sqrt(band_length)

      n_points_on_band (i_band) = 0

      do i_k_point = 1, n_k_points, 1

         k_direction(:) = k_point_list(i_k_point,:) - band_scf_begin(i_band,:) 

         k_distance = 0.d0
         do i_coord = 1,3,1
            k_distance = k_distance + k_direction(i_coord)*k_direction(i_coord)
         enddo
         k_distance = sqrt(k_distance)

         scalar_product = 0.d0
         do i_coord = 1,3,1
            scalar_product = scalar_product + band_dir(i_coord)*k_direction(i_coord)
         enddo

         ! could compute division (cosine of angle between vectors) instead but subtraction is faster
         difference = scalar_product - k_distance*band_length

         if (dabs(difference).le.1d-9) then
            ! k-point on band

            if ( (scalar_product.ge.0.d0).and.(k_distance.le.band_length)) then
               ! k-point also inside band_begin, band_end
 
               if (myid.eq.0) then ! Output only needed on first processor as k-points are known to all.
                  write(use_unit,'(2X,A,I8,A,I4,A,A,F15.8,A,F15.8,A,F15.8,A)') & 
                    "myid = ", myid, ", band number ", i_band, ": ", "Found matching k-point ( ", & 
                    k_point_list(i_k_point,1), ",", k_point_list(i_k_point,2), ",", k_point_list(i_k_point,3), " )"
               end if

               ! memorize information for this k-point

               n_points_on_band(i_band) = n_points_on_band(i_band) + 1
               k_point_index(n_points_on_band(i_band),i_band) = i_k_point
               k_distance_on_band(n_points_on_band(i_band),i_band) = k_distance

            end if

         end if
        
      enddo

      ! Found all k points on this line. Now sort according to k_distance_on_band.

      if (n_points_on_band(i_band).gt.0) then

        ! use the existing insertionsort subroutine here

        call insertionsort ( & 
          k_distance_on_band(1,i_band),n_points_on_band(i_band), &
          new_k_distance_on_band, old_k_index &
          )

        ! And reinsert permutation into previous arrays to avoid huge confusion
        k_distance_on_band(1:n_points_on_band(i_band),i_band) = new_k_distance_on_band(1:n_points_on_band(i_band)) 

        allocate(temp_k_index(n_points_on_band(i_band)))
        temp_k_index(:)=k_point_index(1:n_points_on_band(i_band),i_band)
        do i_point_on_band = 1, n_points_on_band(i_band), 1
           k_point_index(i_point_on_band,i_band) = temp_k_index(old_k_index(i_point_on_band))
        enddo
        deallocate(temp_k_index)

        if (myid.eq.0) then
           write(use_unit,*)
           write(use_unit,'(2X,A,I5,A,I5,A)') "Band number: ", i_band, ": Found ", n_points_on_band(i_band), " k-points."
           do i_point_on_band = 1, n_points_on_band(i_band), 1
              write(use_unit,'(2X,I8,A,F15.8,A,F15.8,A,F15.8,A)') & 
                i_point_on_band, " : ( ", k_point_list(k_point_index(i_point_on_band,i_band),1), ",", & 
                k_point_list(k_point_index(i_point_on_band,i_band),2), ",", & 
                k_point_list(k_point_index(i_point_on_band,i_band),3), " )"
           enddo
        end if

        ! Now write this band
        if(myid==0)then

          ! File name first

          if(i_band < 10)then
             write(file_name, '(A10,I1,A4)') 'scfband100', i_band,'.out'
          else if(i_band < 100)then
             write(file_name, '(A9,I2,A4)') 'scfband10', i_band,'.out'
          else if(i_band < 1000)then
             write(file_name, '(A8,I3,A4)') 'scfband1', i_band,'.out'
          else
             write(use_unit,*) 'Band output error: automatic file name does not work with more than 999 bands!'
             stop
          end if
          open(88, file=file_name)

          if(n_spin >1)then

             if(i_band < 10)then
                write(file_name2, '(A10,I1,A4)') 'scfband200', i_band,'.out'
             else if(i_band < 100)then
                write(file_name2, '(A9,I2,A4)') 'scfband20', i_band,'.out'
             else if(i_band < 1000)then
                write(file_name2, '(A8,I3,A4)') 'scfband2', i_band,'.out'
             else
                write(use_unit,*) 'Band output error: automatic file name does not work with more than 999 bands!'
                stop
             end if
             open(89, file=file_name2)

          end if

          do  i_point_on_band = 1,  n_points_on_band(i_band)

             write(88,'(I4,2X,3F15.7)',ADVANCE='NO') i_point_on_band,  & 
                k_point_list(k_point_index(i_point_on_band,i_band),1), & 
                k_point_list(k_point_index(i_point_on_band,i_band),2), & 
                k_point_list(k_point_index(i_point_on_band,i_band),3)

             do  i_state = 1,  n_states

                write(88,'(F12.5,F15.5)',ADVANCE='NO') & 
                     occ_numbers(i_state,1,k_point_index(i_point_on_band,i_band)), &
                     (KS_eigenvalue(i_state,1,k_point_index(i_point_on_band,i_band))-chemical_potential)* hartree
             end do
             write(88,'()') 

             if(n_spin ==2)then

                write(89,'(I4,2X,3F15.7)',ADVANCE='NO') i_point_on_band,  & 
                k_point_list(k_point_index(i_point_on_band,i_band),1), & 
                k_point_list(k_point_index(i_point_on_band,i_band),2), & 
                k_point_list(k_point_index(i_point_on_band,i_band),3)

                do  i_state = 1,  n_states

                   write(89,'(F12.5,F15.5)',ADVANCE='NO') & 
                        occ_numbers(i_state,2,k_point_index(i_point_on_band,i_band)), &
                        (KS_eigenvalue(i_state,2,k_point_index(i_point_on_band,i_band))-chemical_potential)* hartree

                end do
                write(89,'()') 
             end if

          end do

          close(88)
          if(n_spin==2) close(89)

        end if ! myid == 0

      else

        if (myid.eq.0) then
           write(use_unit,*)
           write(use_unit,'(2X,A,I5,A,I5,A)') "Band number: ", i_band, ": Found ", n_points_on_band(i_band), " k-points."
           write(use_unit,'(2X,A)') "Not writing any band structure file for this band."
        end if

      end if

   enddo

  end subroutine output_bands_during_scf
!******
