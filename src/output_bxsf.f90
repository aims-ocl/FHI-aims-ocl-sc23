!****s* FHI-aims/output_bxsf
!  NAME
!  output_bxsf
!  SYNOPSIS
! Purpose of this routine is to create a band-xsf file (bxsf), which
! can be visualized with XCrysden, e.g., to visualize Fermi surfaces
subroutine output_bxsf ()
use dimensions
use runtime_choices
use geometry
use physics
use constants
use localorb_io
use mpi_tasks, only: myid, aims_stop

implicit none

!real*8 :: KS_eigenvalue(n_states,n_spin,n_k_points)
character*100 :: filename = 'Test.bxsf' !FIXME: Name shoud be custimizable
integer :: ikx, iky, ikz, i_kpoint, i_state
 character(LEN=100) :: info_str


!Safeguards
if (n_spin.ne.1) then
   write(info_str,*) 'Right now, this routine only works for spin-unrestricted, though that would be easy to change'
   call aims_stop(info_str)
   stop
endif
if (.not.collect_eigenvectors) call aims_stop('output of fermisurface works only with collect_eigenvectors .true.')

!Do only  on single processor
if (myid.eq.0) then
   write(use_unit,*) 'Starting output of eigenvalues on grid'
   open(unit=20,file=trim(filename))
   
   !Header
   write(20,*) 'BEGIN_INFO'
   write(20,fmt='(2X,A,10X,E15.5)') 'Fermi Energy: ', chemical_potential*hartree !TODO: Add actual chemical potential
   write(20,*) 'END_INFO'
   write(20,*) ''
   write(20,*) "BEGIN_BLOCK_BANDGRID_3D"
   write(20,*) "Eigenvalues_for_bands "
   write(20,*) "BEGIN_BANDGRID_3D_simple_example"
   write(20,fmt='(4X,I0)') n_states  !FIXME: make custimizable. Currently: all
   write(20,fmt='(2X,I0,2X,I0,2X,I0)') n_k_points_xyz(1), n_k_points_xyz(2), n_k_points_xyz(3)
   write(20,fmt='(2X,I0,2X,I0,2X,I0)') 0,0,0 !FIXME: Check for actual origin and stop / warn for off-center 
   write(20,*) recip_lattice_vector(1:3,1) !FIXME: Formatting
   write(20,*) recip_lattice_vector(1:3,2) !FIXME: Formatting
   write(20,*) recip_lattice_vector(1:3,3) !FIXME: Formatting
   !End header
   
   do i_state=1,n_states,1 !go ever each state
      write(20,fmt='(4X,A,4X,I0)') 'BAND:', i_state
      i_kpoint=0
      do ikx=1,n_k_points_xyz(1),1
         do iky=1,n_k_points_xyz(2),1
            do ikz=1,n_k_points_xyz(3),1
                i_kpoint=i_kpoint+1
                write(20,fmt='(E15.5,2X,$)') KS_eigenvalue(i_state,1,i_kpoint)*hartree !TODO: Check ordering of k-points
            enddo
            write(20,*) '' !Linebreak
         enddo
         write(20,*) '' !Linebreak
      enddo
   
   enddo
   
   
   write(20,*) '   END_BANDGRID_3D'
   write(20,*) 'END_BLOCK_BANDGRID_3D'
   
   close(20)
endif

write(info_str,*) ' ouptut of eigenvalues on grid finished'
call localorb_info(info_str)




end subroutine output_bxsf

