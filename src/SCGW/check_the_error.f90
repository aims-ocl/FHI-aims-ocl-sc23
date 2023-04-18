      subroutine check_the_error (polar_freq, polar_green_freq, &
          n_basbas, n_loc_prodbas, omega, nomega, womega, name_of_quantity,&
                 average_error, max_error, printout )


      use localorb_io, only: use_unit
      use mpi_tasks


      implicit none
!ARGUMENTS
      integer n_basbas, n_loc_prodbas, nomega
      real*8 omega(-nomega:nomega)
      real*8 womega(nomega)
      real*8 polar_freq       (n_basbas, n_loc_prodbas,1:nomega)
      real*8 polar_green_freq (n_basbas, n_loc_prodbas,1:nomega)
      logical printout

! statistical stuff      
      real*8 error (n_basbas, n_loc_prodbas)
      real*8 average_polar_matr (n_basbas, n_loc_prodbas)
      real*8 error_row(n_basbas)
      real*8 average_polar
      real*8 average_error
!      real*8 average_relative_error
      real*8 max_error
      real*8 hystogram_interval
      integer n_hystogram_columns
      real*8, dimension(:,:), allocatable :: hystogram

      real*8 reference_value (n_basbas, n_loc_prodbas)
 
!counters
      integer i_freq
      integer i_basbas
      integer j_basbas
      integer i_interval
      integer counter
      character*5 name_of_quantity
      character*40 file_name
      character*4 extension     
      logical end_of_hystogram_loop          
      real*8 threshold 

      extension = '.dat'
      error(:,:) = 0.d0
      max_error = 0.d0
      reference_value (:,:) = 0.d0
      average_polar = 0.d0
      average_polar_matr (:,:)= 0.d0
      do i_freq=1,nomega, 1
!        do i_basbas = 1, n_basbas,1 
!          do j_basbas = 1, n_basbas,1

!            write(use_unit,*) i_freq
!            error(i_basbas,j_basbas) = (polar_freq(i_basbas,j_basbas,i_freq)-&
!            polar_green_freq (i_basbas,j_basbas,i_freq))**2 *womega(i_freq)

           error(:,:) = error(:,:)+ (polar_freq(:,:,i_freq)-&
           polar_green_freq (:,:,i_freq))**2 *womega(i_freq)
           
           average_polar_matr (:,:) =  average_polar_matr (:,:)+&
           polar_freq (:,:, i_freq)*womega(i_freq) 
!           reference_value (:,:)= reference_value (:,:)+ &
!           polar_freq(:,:,i_freq)**2 *womega(i_freq)

!          enddo
!        enddo
      enddo

!  some statistics on the error

! average error----------------------------------------------------------
      average_error = 0.d0
      !naverage_relative_error = 0.d0
      threshold = 1.d-6
      do i_basbas = 1, n_basbas,1 
        do j_basbas = 1, n_loc_prodbas,1

          if (printout)then
           if (error(i_basbas,j_basbas).gt.threshold) then
            if(myid.eq.0)then

             write(use_unit,*) " *** matrix element ", i_basbas,j_basbas

             if(i_basbas.lt.10)then
                if(j_basbas .lt.10)then 
                   write(extension,'(A,I1,A,I1)') "0",i_basbas,"0",j_basbas
                else 
                   write(extension,'(A,I1,I2)') "0",i_basbas,j_basbas
                endif
             else
                if(j_basbas .lt.10)then
                   write(extension,'(I2,A,I1)') i_basbas,"0",j_basbas
                else 
                   write(extension,'(I2,I2)') i_basbas,j_basbas
                endif
             endif

             if(.false.)then
                file_name = 'green'//extension//'.dat'
                open(44,file=file_name)
                  do i_freq=1,nomega, 1
                    write(44,*) omega(i_freq),&
                polar_green_freq (i_basbas,j_basbas,i_freq),&
                polar_freq (i_basbas,j_basbas,i_freq) 
                  enddo 
                close(44)
             endif
 
            endif
           endif
          endif

          average_error = average_error + &
           error(i_basbas,j_basbas)

          average_polar = average_polar +&
             average_polar_matr(i_basbas,j_basbas)

         

          if (error(i_basbas,j_basbas).gt.max_error)then
             max_error = error(i_basbas,j_basbas)
          endif

        enddo
      enddo    
 
      max_error = sqrt(max_error)
      average_polar = average_polar / (n_basbas*n_loc_prodbas)
      average_error = sqrt(average_error /(n_basbas*n_loc_prodbas))
!      average_relative_error = average_error / average_polar

      if(myid.eq.0)then
        write(use_unit,'(A,A,E14.3)') "          | Average Deviation on  ",&
                   name_of_quantity, average_error
!        write(use_unit,'(A,A,E14.3)') "          | Relative Deviation on  ",&
!                   name_of_quantity , average_relative_error
        !write(use_unit,'(A,A,E14.3)') "          | Maximum Deviation on  ",&
        !        name_of_quantity , max_error
      endif
!------------------------------------------------------------------------


! do the average of the errors of each row --------------------------
!      error_row(:) = 0.d0
!      file_name = 'error_row_'//name_of_quantity//extension
!
!      if (myid.eq.0)then      
!      open(33,file=file_name)
!      do i_basbas = 1, n_basbas,1
!        do j_basbas = 1, n_loc_prodbas,1
!          error_row(i_basbas) =  error_row(i_basbas) + &
!            error(i_basbas, j_basbas)
!        enddo
!        write(33,*) i_basbas , error_row(i_basbas)/n_basbas
!      enddo
!      close(33)
!      endif
!--------------------------------------------------------------------

! do an hystograms of the errors

! number of suddivisions of the interval
!      n_hystogram_columns = 200
!      hystogram_interval = max_error / n_hystogram_columns

! the hystograms is defined in the interval [ 0 : n_hystogram_columns*hystogram_interval]
!      hystogram_interval = 2*average_error / n_hystogram_columns/5.
!      allocate(hystogram(n_hystogram_columns,2))
!!      hystogram(:,:) = 0
!
!      end_of_hystogram_loop = .false.
!      counter = 0 

!      file_name = 'error_hystogram_'//name_of_quantity//extension
!      if(myid.eq.0)then
!      open(67,file=file_name)

      !refine the definition of the hystogram_interval until the first interval of
      ! the hystogram contains a number of points minor than the desired rate of total point
      !e.g. it stops reiterating only if hystogram(1,2).lt.(n_basbas**2)/4.
      
!      do while (.not. end_of_hystogram_loop) 
!        hystogram(:,:) = 0.d0
!        counter = counter + 1
!        do i_interval = 1, n_hystogram_columns, 1
!   
!          hystogram(i_interval,1 ) = i_interval * hystogram_interval
!   
!          do i_basbas = 1, n_basbas,1
!            do j_basbas = 1, n_loc_prodbas,1
!
!!here it counts how many points stays in each interval         
!
!         if(error(i_basbas,j_basbas).gt.((i_interval-1)*hystogram_interval) &
!          .and. error(i_basbas,j_basbas).le.((i_interval)*hystogram_interval))then
!           hystogram(i_interval,2) = hystogram(i_interval,2) + 1 
!         endif
!         
!            enddo
!          enddo
!writes down the result of the first iteration
!         if(counter .eq. 1)then
!           write(67,*) hystogram(i_interval,1 ), hystogram(i_interval,2 )
!         endif

!        enddo
!        if(hystogram(1,2).lt.(n_basbas**2)/5.)then
!          end_of_hystogram_loop = .true.
!        else 
!          hystogram_interval = hystogram_interval/3
!        endif
!        if (counter .gt. 100) end_of_hystogram_loop = .true.
!      enddo
!     close(67)

!      file_name = 'refined_error_hystogram_'//name_of_quantity//extension
!      open(35,file=file_name)
!      do i_interval = 1, n_hystogram_columns,1  
!          if(hystogram(i_interval,2 ).gt.0)then
!            write(35,*) hystogram(i_interval,1 ), hystogram(i_interval,2 )
!          endif
!      enddo
!      close(35)
!      endif
 
      end subroutine check_the_error
