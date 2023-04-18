module directories


  implicit none

  public create_dir

  private

  contains

    function create_dir(dirName) result(success)

      implicit none

      character*(*), intent(in) :: dirName
      logical                   :: success

      integer                   :: istat
      integer, parameter        :: dirUnit = 34
      logical                   :: fileExists

      integer, dimension(8)     :: dateAndTime 
      success = .false.
      
      ! check whether the directory exits
      ! use a dirty hack: try to open the file ./dirName/log

!      open(dirUnit,file=trim(dirName)//"/log",status="old",iostat=istat)

      inquire(file=trim(dirName)//"/log",exist=fileExists)


      if (fileExists) then
         ! could open file
         success = .true.
         
         ! close file again

         close(dirUnit)

      else
         ! could not open file

         ! silently create the directory and the log file

         call system('mkdir ' // trim(dirName) )
         
         
         success = .true.

  
         open(dirUnit,file=trim(dirName)//"/log",status="new",iostat=istat)
         
         call date_and_time(values=dateAndTime(:))


         write(dirUnit,'(6i4)') dateAndTime(1),dateAndTime(2),dateAndTime(3),dateAndTime(5),dateAndTime(6),dateAndTime(7)

         close(dirUnit)
      
      endif


      return
    end function create_dir


end module directories
