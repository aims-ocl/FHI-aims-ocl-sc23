
!****s* FHI-aims/write_preamble
!  NAME
!    write_preamble
!  SYNOPSIS

    subroutine write_preamble()

      use generate_aims_uuid, only: write_aims_uuid
      use localorb_io
      use timing
      use mpi_utilities, only: initial_mpi_report

      implicit none

      character*150 info_str
      character(LEN=80) ::  uuid_str


      write(info_str,'(A)') "------------------------------------------------------------"
      call localorb_info( info_str )
      write(info_str,'(10X,A)') "Invoking FHI-aims ..."
      call localorb_info( info_str )

      call write_version_stamp()

      call localorb_info( " " )
      write(info_str,'(10X,A)') "When using FHI-aims, please cite the following reference:"
      call localorb_info( info_str )
      call localorb_info( " " )
      write(info_str,'(10X,A)') "Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, "
      call localorb_info( info_str )
      write(info_str,'(10X,A)') "Ville Havu, Xinguo Ren, Karsten Reuter, and Matthias Scheffler,"
      call localorb_info( info_str )
      write(info_str,'(10X,A)') "'Ab Initio Molecular Simulations with Numeric Atom-Centered Orbitals',"
      call localorb_info( info_str )
      write(info_str,'(10X,A)') "Computer Physics Communications 180, 2175-2196 (2009)"
      call localorb_info( info_str )
      call localorb_info( " " )
      write(info_str,'(10X,A)') "For any questions about FHI-aims, please visit the aimsclub website"
      call localorb_info( info_str )
      write(info_str,'(10X,A)') "with its forums and wiki. Contributions to both the forums and the"
      call localorb_info( info_str )
      write(info_str,'(10X,A)') "wiki are warmly encouraged - they are for you, and everyone is welcome there."
      call localorb_info( info_str )
      call localorb_info( " " )
      write(info_str,'(A)') "------------------------------------------------------------"
      call localorb_info( info_str )

      call localorb_info( " " )


      call localorb_info( " " )

      ! report initial information on time & MPI parallelism (if relevant)
      call initial_timings ( )

      write(info_str,'(2X,A)')"FHI-aims created a unique identifier for this run for later identification"
      call localorb_info( info_str )
      call write_aims_uuid(uuid_str)
      write(info_str,'(2X,A)') uuid_str
      call localorb_info( info_str )
      call localorb_info( " " )

      call cmake_info()

      call initial_mpi_report ( )

    end subroutine write_preamble
!******
