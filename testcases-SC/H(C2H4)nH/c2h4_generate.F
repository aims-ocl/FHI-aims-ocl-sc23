!-------------------------------------------------------------------------
! a small code to generate H-[c2h4]n-H  mol and [c2h4]n extended system. 
! use exp data as input from P.W. Teare, Acta Cryst. 12,294 (1959). 
! you can control to make mol or extended systems by using add_H keyword. 
!                                              shanghui@Berlin, 2015.10.24
!------------------------------------------------------------------------
      program c2h4_generate
      implicit none
      
      real*8, parameter :: pi = 3.14159265358979323846264338327950288419716939937510d0
      real*8, parameter :: deg = pi / 180.0d0

      real*8, parameter :: r_CC=1.53d0      ! r_CC: C-C single bond
      real*8, parameter :: r_CH=1.07d0       ! r_CH: C-H bond
      real*8, parameter :: alpha_CCC=112.0d0     ! C-C-C 
      real*8, parameter :: alpha_HCH=107.5d0     ! H-C-H
      real*8, parameter :: b=40.0d0       ! 40 ang in b and c. 
  
!-------------add H-----------------------------------------------
      real*8, parameter :: alpha_CCH = 110.0d0 
      logical, parameter :: add_H = .true.  ! mol=add_H
      ! integer, parameter :: n = 640 !5000 !10000 !36667 ! 19600  !6144 ! 1667
      integer :: n
      integer, parameter :: max_num_of_atoms = 2000000

  ! 1667 -> 10000
  ! 3333 -> 20000
  ! 6667 -> 40000

  ! 1000 -> 6002

!-----------for mol------------------------
      integer i_atom,i_cell
      real*8 coords(3,max_num_of_atoms), phi_2
      character*10 species(max_num_of_atoms) 
      character(len=64) :: arg

!----------for extended system-------------
      integer i_coord,i_periodic
      real*8 c1(3),c2(3),h1_c1(3),h2_c1(3),h2_c2(3),h1_c2(3), a,phi
      real*8 lattice_vector(3,3)     
 
      if (command_argument_count() /= 1) then
        write(6,'(A)')   "Error: This program requires exactly one "
        write(6,'(A)')   "       integer argument as the n of "
        write(6,'(A)')   "       H (C_2 H_4)_n H ."
        stop
      endif
    
      call get_command_argument(1, arg)
      read(arg, *) n
      ! print *, "The integer value is:", n

        a=dsqrt(r_CC**2+r_CC**2-2.0d0*r_CC*r_CC*dcos(alpha_CCC*deg))
        phi=(180.0d0-alpha_CCC)/2.0d0*deg
        !dacos( (r_CC**2+a**2-r_CC**2)/(2.0d0*r_CC*a) )

      if(add_H) then ! mol system 


       do i_cell = 1,n 

        !---C1----
         coords(1,1+(i_cell-1)*6)=-0.5d0*r_CC*dcos(phi) + a*(i_cell-1) 
         coords(2,1+(i_cell-1)*6)=-0.5d0*r_CC*dsin(phi)
         coords(3,1+(i_cell-1)*6)= 0.0d0
         species(1+(i_cell-1)*6) = ' C '

        !---C2---
         coords(1,2+(i_cell-1)*6)= 0.5d0*r_CC*dcos(phi) + a*(i_cell-1)
         coords(2,2+(i_cell-1)*6)= 0.5d0*r_CC*dsin(phi)
         coords(3,2+(i_cell-1)*6)= 0.0d0
         species(2+(i_cell-1)*6) = ' C '
       
        !---H1_C1---
         coords(1,3+(i_cell-1)*6)= coords(1,1+(i_cell-1)*6)
         coords(2,3+(i_cell-1)*6)= coords(2,1+(i_cell-1)*6)
     .                           - r_CH*dcos(alpha_HCH/2*deg)
         coords(3,3+(i_cell-1)*6)= coords(3,1+(i_cell-1)*6) 
     .                           +r_CH*dsin(alpha_HCH/2*deg)
         species(3+(i_cell-1)*6) = ' H '

        !---H2_C1---
         coords(1,4+(i_cell-1)*6)= coords(1,1+(i_cell-1)*6)
         coords(2,4+(i_cell-1)*6)= coords(2,1+(i_cell-1)*6)
     .                           - r_CH*dcos(alpha_HCH/2*deg)
         coords(3,4+(i_cell-1)*6)= coords(3,1+(i_cell-1)*6)
     .                            -r_CH*dsin(alpha_HCH/2*deg)
         species(4+(i_cell-1)*6) = ' H '

        !---H1_C2---
         coords(1,5+(i_cell-1)*6)= coords(1,2+(i_cell-1)*6)
         coords(2,5+(i_cell-1)*6)= coords(2,2+(i_cell-1)*6)
     .                           + r_CH*dcos(alpha_HCH/2*deg)
         coords(3,5+(i_cell-1)*6)= coords(3,2+(i_cell-1)*6) 
     .                           +r_CH*dsin(alpha_HCH/2*deg)
        
         species(5+(i_cell-1)*6) = ' H '

        !---H2_C2---
         coords(1,6+(i_cell-1)*6)= coords(1,2+(i_cell-1)*6)
         coords(2,6+(i_cell-1)*6)= coords(2,2+(i_cell-1)*6)
     .                           + r_CH*dcos(alpha_HCH/2*deg)
         coords(3,6+(i_cell-1)*6)= coords(3,2+(i_cell-1)*6) 
     .                           - r_CH*dsin(alpha_HCH/2*deg)
          
         species(6+(i_cell-1)*6) = ' H '
 
       enddo     

         phi_2= (180.0d0 - alpha_CCH - (180.0d0-alpha_CCC)/2.0d0)*deg

         coords(1,1 + n*6)= coords(1,1) - r_CH*dcos(phi_2) 
         coords(2,1 + n*6)= coords(2,1) + r_CH*dsin(phi_2)
         coords(3,1 + n*6)= coords(3,1) 
         species(1  + n*6) = ' H '
  
        
         coords(1,2 + n*6)= coords(1,n*6-4) + r_CH*dcos(phi_2)  
         coords(2,2 + n*6)= coords(2,n*6-4) - r_CH*dsin(phi_2)
         coords(3,2 + n*6)= coords(3,n*6-4)
         species(2  + n*6) = ' H '

         do i_atom = 1, 6*n+2 
            write(6,'(12X,A,3(2X,F16.8),2X,A)') ,  "atom ",
     .      (coords(i_coord,i_atom), i_coord=1,3,1),
     .        species(i_atom)
         enddo




      else          ! extended system 

      c1(1)=-0.5d0*r_CC*dcos(phi)
      c1(2)=-0.5d0*r_CC*dsin(phi)+b/2.0d0
      c1(3)= b/2.0d0

      c2(1)= 0.5d0*r_CC*dcos(phi)
      c2(2)= 0.5d0*r_CC*dsin(phi)+b/2.0d0
      c2(3)= b/2.0d0
       

      h1_c1(1)=c1(1)
      h1_c1(2)=c1(2)-r_CH*dcos(alpha_HCH/2*deg)
      h1_c1(3)=c1(3)+r_CH*dsin(alpha_HCH/2*deg)

      h2_c1(1)=c1(1)
      h2_c1(2)=c1(2)-r_CH*dcos(alpha_HCH/2*deg)
      h2_c1(3)=c1(3)-r_CH*dsin(alpha_HCH/2*deg)


      h1_c2(1)=c2(1)
      h1_c2(2)=c2(2)+r_CH*dcos(alpha_HCH/2*deg)
      h1_c2(3)=c2(3)+r_CH*dsin(alpha_HCH/2*deg)

      h2_c2(1)=c2(1)
      h2_c2(2)=c2(2)+r_CH*dcos(alpha_HCH/2*deg)
      h2_c2(3)=c2(3)-r_CH*dsin(alpha_HCH/2*deg)


      lattice_vector = 0.0d0
      lattice_vector(1,1) = a 
      lattice_vector(2,2) = b
      lattice_vector(3,3) = b

      do i_periodic = 1, 3, 1
          write(6,'(2X,A,3(2X,F16.8),2X)') 
     .         "lattice_vector ", 
     .          (lattice_vector(i_coord,i_periodic), i_coord=1,3,1)
      enddo


       write(6,'(12X,A,3(2X,F16.8),2X,A)')   "atom ", 
     .       (c1(i_coord), i_coord=1,3,1), 'C'
       write(6,'(12X,A,3(2X,F16.8),2X,A)')   "atom ", 
     .       (c2(i_coord), i_coord=1,3,1), 'C'
       write(6,'(12X,A,3(2X,F16.8),2X,A)')   "atom ", 
     .       (h1_c1(i_coord), i_coord=1,3,1), 'H'
       write(6,'(12X,A,3(2X,F16.8),2X,A)')   "atom ", 
     .       (h2_c1(i_coord), i_coord=1,3,1), 'H'
       write(6,'(12X,A,3(2X,F16.8),2X,A)')   "atom ", 
     .       (h1_c2(i_coord), i_coord=1,3,1), 'H'
       write(6,'(12X,A,3(2X,F16.8),2X,A)')   "atom ", 
     .       (h2_c2(i_coord), i_coord=1,3,1), 'H'
             
      endif


      end
