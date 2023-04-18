!****s* FHI-aims/nsolve_iepa_mp2
!  NAME
!   nsolve_iepa_mp2
!  SYNOPSIS

subroutine nsolver_iepa_mp2(numerator,denominator,&
     i_spin, i_spin2, i_state, i_state2, n_lumo,&
     initial_energy,initial_energy_ss,initial_energy_os,&
     n_lumo_min,enhanced_factor,screening_factor,shift_factor, &
     threshold,final_energy,final_energy_ss,final_energy_os)

!  PURPOSE
!  The subroutine solve equation x = \sum_{ij} A_{ij}/(B_{ij}-x).
!  USES

  use dimensions, only: n_states, n_spin
  use mpi_tasks
  use mpi_utilities
  use synchronize_mpi
  use arch_specific
  use localorb_io, only: use_unit
  implicit none

!  ARGUMENTS
   integer :: i_spin, i_spin2, i_state, i_state2, i_index
   integer :: n_lumo_min
   integer :: n_lumo(n_spin)
   real*8  :: initial_energy
   real*8  :: initial_energy_ss
   real*8  :: initial_energy_os
   real*8  :: final_energy
   real*8  :: final_energy_ss
   real*8  :: final_energy_os
   real*8  :: threshold
   real*8  :: enhanced_factor, screening_factor, shift_factor
   !real*8  :: numerator_a(n_lumo_min:n_states,n_lumo_min:n_states)
   !real*8  :: numerator_b(n_lumo_min:n_states,n_lumo_min:n_states)
   real*8  :: numerator(n_lumo_min:n_states,n_lumo_min:n_states)
   real*8  :: denominator(n_lumo_min:n_states,n_lumo_min:n_states)


! INPUTS
! o numerator -- A_{ij}
! o denominator -- B_{ij}
! o initial_energy -- Ec[MP2] term
! o n_unoccupied -- the number of virtual states
! o threshold -- convergence threshold
!
! OUTPUT
! o final_energy -- Ec[IEPA] term


! Error variables
  integer mpierr
  integer :: info
  character*150 :: info_str
  character(*), parameter :: func = 'nsolver_iepa_mp2'
 
! Local variables
  integer :: i,j,a,b
  integer :: max_scf
  integer :: n_scf
  real*8  :: tmp_energy
  real*8  :: tmp_energy_ss
  real*8  :: tmp_energy_os
  real*8  :: tmp_energy_term
  real*8  :: tmp_energy_term_ss
  real*8  :: tmp_energy_term_os
  real*8  :: tmp_energy_mpi
  real*8  :: tmp_energy_mpi_ss
  real*8  :: tmp_energy_mpi_os
  real*8  :: tmp_denominator

  real*8  :: delta_energy_0
  real*8  :: delta_energy_1
  real*8  :: previous_energy
  real*8  :: previous_energy_ss
  real*8  :: previous_energy_os

  max_scf           = 100
  n_scf             = 0

  tmp_energy        = 0.0d0
  tmp_energy_ss     = 0.0d0
  tmp_energy_os     = 0.0d0

  tmp_energy_mpi    = 0.0d0
  tmp_energy_mpi_ss = 0.0d0
  tmp_energy_mpi_os = 0.0d0


  ! ==============================================================
  !  In spin-unpolarized cases:
  !    e_{ij,ss} = \sum_{a,b}^{vir}    (ia|jb)^2-(ia|jb)(ib|ja)
  !                                 x ---------------------------
  !                                    e_i+e_j-e_a-e_b+e_{ij,ss}
  !    e_{ij,os} = \sum_{a,b}^{vir}    (ia|jb)^2
  !                                  x --------------------------
  !                                     e_i+e_j-e_a-e_b+e_{ij,os}
  ! ==============================================================
  if (screening_factor .eq. 0.0d0) then
      if (n_spin.eq.1) then
          i_index = 0
          do i=n_lumo(i_spin), n_states, 1
            do j=n_lumo(i_spin2), n_states, 1
                i_index = i_index + 1
                if(myid.eq.MOD(i_index,n_tasks)) then
                   tmp_energy_term_ss = &
                       (numerator(i,j)-numerator(j,i)) * & 
                       numerator(i,j) / &
                       (denominator(i,j) + enhanced_factor*initial_energy_ss)
                   tmp_energy_term_os = &
                       numerator(i,j)**(2.d0)/ &
                       (denominator(i,j) + enhanced_factor*initial_energy_os)
                   tmp_energy_term = tmp_energy_term_os + tmp_energy_term_ss
                   tmp_energy    = tmp_energy    + tmp_energy_term   
                   tmp_energy_ss = tmp_energy_ss + tmp_energy_term_ss
                   tmp_energy_os = tmp_energy_os + tmp_energy_term_os
              endif ! if(myid.eq.MOD(i_index,n_tasks))
            enddo ! j
          enddo ! i
      else ! if (n_spin.ne.1)
          if (i_spin .eq. i_spin2) then
              i_index = 0
              do i=n_lumo(i_spin), n_states, 1
                do j=i+1, n_states, 1
                  i_index = i_index + 1
                  if(myid.eq.MOD(i_index,n_tasks)) then
                      tmp_energy_term_ss = &
                          (numerator(i,j)-numerator(j,i))**(2.0d0) / &
                          (denominator(i,j)+ enhanced_factor*initial_energy_ss)
                      tmp_energy      = tmp_energy + tmp_energy_term_ss
                      tmp_energy_ss   = tmp_energy_ss + tmp_energy_term_ss
                  endif
                enddo ! do j=n_lumo(i_spin2), n_states, 1
              enddo
          else ! if (i_spin .ne. i_spin2)
              i_index = 0
              do i=n_lumo(i_spin), n_states, 1
                do j=n_lumo(i_spin2), n_states, 1
                  i_index = i_index + 1
                  if(myid.eq.MOD(i_index,n_tasks)) then
                      tmp_energy_term_os = &
                          numerator(i,j)**(2.0d0) / &
                          (denominator(i,j)+ enhanced_factor*initial_energy_os)
                      tmp_energy      = tmp_energy + tmp_energy_term_os
                      tmp_energy_os   = tmp_energy_os + tmp_energy_term_os
                  endif
                enddo ! do j=n_lumo(i_spin2), n_states, 1
              enddo
          endif ! if(i_spin.ne.i_spin2)
      endif !if (n_spin.ne.1)
  else
      if (myid.eq.0) then
          write(use_unit,'(2X,A,F12.6,F12.6,F12.6)') &
              'Screened coupling is invoked with', &
              enhanced_factor, screening_factor, shift_factor
      endif
      if (n_spin.eq.1) then
          i_index = 0
          do i=n_lumo(i_spin), n_states, 1
            do j=n_lumo(i_spin2), n_states, 1
              if(myid.eq.MOD(i_index,n_tasks)) then
                   tmp_denominator  = enhanced_factor * &
                       arch_erfc(screening_factor*(denominator(i,j)))
                   tmp_energy_term_ss = &
                       (numerator(i,j)-numerator(j,i)) * & 
                       numerator(i,j) / (denominator(i,j) + &
                       tmp_denominator * (initial_energy_ss + shift_factor))
                   tmp_energy_term_os = &
                       numerator(i,j)*numerator(i,j)/ (denominator(i,j) + &
                       tmp_denominator * (initial_energy_os + shift_factor))
                   tmp_energy_term = tmp_energy_term_os + tmp_energy_term_ss
                   tmp_energy    = tmp_energy    + tmp_energy_term   
                   tmp_energy_ss = tmp_energy_ss + tmp_energy_term_ss
                   tmp_energy_os = tmp_energy_os + tmp_energy_term_os
              endif ! if(myid.eq.MOD(i_index,n_tasks))
            enddo ! j
          enddo ! i
      else ! if (n_spin.ne.1)
          if (i_spin .eq. i_spin2) then
              i_index = 0
              do i=n_lumo(i_spin), n_states, 1
                do j=i+1, n_states, 1
                  i_index = i_index + 1
                  if(myid.eq.MOD(i_index,n_tasks)) then
                       tmp_denominator  = enhanced_factor * &
                           arch_erfc(screening_factor*(denominator(i,j)))
                       tmp_energy_term_ss = &
                           (numerator(i,j)-numerator(j,i))**(2.0d0) / &
                           (denominator(i,j) + &
                           tmp_denominator*(initial_energy_ss + shift_factor))
                       tmp_energy      = tmp_energy + tmp_energy_term_ss
                       tmp_energy_ss   = tmp_energy_ss + tmp_energy_term_ss
                  endif
                enddo ! do j=n_lumo(i_spin2), n_states, 1
              enddo
          else ! (i_spin .ne. i_spin2)
              i_index = 0
              do i=n_lumo(i_spin), n_states, 1
                do j=n_lumo(i_spin2), n_states, 1
                  i_index = i_index + 1
                  if(myid.eq.MOD(i_index,n_tasks)) then
                       tmp_denominator  = enhanced_factor * &
                           arch_erfc(screening_factor*(denominator(i,j)))
                        tmp_energy_term_os = &
                            numerator(i,j)**(2.0d0) / (denominator(i,j) +&
                            tmp_denominator*(initial_energy_os+shift_factor))
                        tmp_energy      = tmp_energy + tmp_energy_term_os
                        tmp_energy_os   = tmp_energy_os + tmp_energy_term_os
                  endif
                enddo ! do j=n_lumo(i_spin2), n_states, 1
              enddo
          endif ! (i_spin.eq.i_spin2) or (i_spin.ne.i_spin2)
      endif ! (n_spin.eq.1) or (n_spin.ne.1)
  endif

  if(use_mpi) then
      call MPI_ALLREDUCE(tmp_energy, tmp_energy_mpi, 1, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      tmp_energy=tmp_energy_mpi 
      call MPI_ALLREDUCE(tmp_energy_ss, tmp_energy_mpi_ss, 1, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      tmp_energy_ss=tmp_energy_mpi_ss 
      call MPI_ALLREDUCE(tmp_energy_os, tmp_energy_mpi_os, 1, &
           MPI_DOUBLE_PRECISION, &
           MPI_SUM, mpi_comm_global, mpierr)
      tmp_energy_os=tmp_energy_mpi_os 
  endif

  if (myid.eq.0) then
      write(use_unit,'(2X,A,E10.3,E10.3,E10.3,A,E10.3,E10.3,E10.3)') &
          'Initial energy ',initial_energy,initial_energy_ss,&
          initial_energy_os, &
          ', generated energy ', tmp_energy, tmp_energy_ss, &
          tmp_energy_os
  endif

  !if (initial_energy .lt. -3.0d0) initial_energy = -1.0d0
  delta_energy_0     = tmp_energy - initial_energy
  delta_energy_1     = tmp_energy - initial_energy
  previous_energy    = tmp_energy
  previous_energy_ss = tmp_energy_ss
  previous_energy_os = tmp_energy_os
  SCF_LOOP: do while (abs(delta_energy_0) .gt. threshold &
                      .and. n_scf .lt. max_scf)
            n_scf = n_scf + 1
            tmp_energy      = 0.0d0
            tmp_energy_ss   = 0.0d0
            tmp_energy_os   = 0.0d0
            if (screening_factor .eq. 0.0d0) then
                if (n_spin.eq.1) then
                    i_index = 0
                    do i=n_lumo(i_spin), n_states, 1
                      do j=n_lumo(i_spin2), n_states, 1
                        i_index = i_index + 1
                        if(myid.eq.MOD(i_index,n_tasks)) then
                             tmp_energy_term_ss = &
                                 (numerator(i,j)-numerator(j,i)) * & 
                                 numerator(i,j) / &
                                 (denominator(i,j)+enhanced_factor*previous_energy_ss)
                             tmp_energy_term_os = &
                                 numerator(i,j)**(2.d0)/ &
                                 (denominator(i,j)+enhanced_factor*previous_energy_os)
                             tmp_energy_term = tmp_energy_term_os + tmp_energy_term_ss
                             tmp_energy    = tmp_energy    + tmp_energy_term   
                             tmp_energy_ss = tmp_energy_ss + tmp_energy_term_ss
                             tmp_energy_os = tmp_energy_os + tmp_energy_term_os
                        endif
                      enddo
                    enddo
                else
                    if (i_spin.eq.i_spin2) then
                        i_index = 0
                        do i=n_lumo(i_spin), n_states, 1
                          do j=i+1, n_states, 1
                            i_index = i_index + 1
                            if(myid.eq.MOD(i_index,n_tasks)) then
                                tmp_energy_term_ss = &
                                    (numerator(i,j)-numerator(j,i))**(2.0d0) / &
                                    (denominator(i,j)+ enhanced_factor*previous_energy_ss)
                                tmp_energy      = tmp_energy + tmp_energy_term_ss
                                tmp_energy_ss   = tmp_energy_ss + tmp_energy_term_ss
                            endif
                          enddo ! do j=n_lumo(i_spin2), n_states, 1
                        enddo
                    else ! if (i_spin.ne.i_spin2)
                        i_index = 0
                        do i=n_lumo(i_spin), n_states, 1
                          do j=n_lumo(i_spin2), n_states, 1
                            i_index = i_index + 1
                            if(myid.eq.MOD(i_index,n_tasks)) then
                                tmp_energy_term_os = &
                                    numerator(i,j)**(2.0d0) / &
                                    (denominator(i,j)+ enhanced_factor*previous_energy_os)
                                tmp_energy      = tmp_energy + tmp_energy_term_os
                                tmp_energy_os   = tmp_energy_os + tmp_energy_term_os
                            endif
                          enddo ! do j=n_lumo(i_spin2), n_states, 1
                        enddo
                    endif
                endif
            else
                if (n_spin.eq.1) then
                    i_index = 0
                    do i=n_lumo(i_spin), n_states, 1
                      do j=n_lumo(i_spin2), n_states, 1
                        i_index = i_index + 1
                        if(myid.eq.MOD(i_index,n_tasks)) then
                             tmp_denominator  = enhanced_factor * &
                                 arch_erfc(screening_factor*(denominator(i,j)))
                             tmp_energy_term_ss = &
                                 (numerator(i,j)-numerator(j,i)) * & 
                                 numerator(i,j) / (denominator(i,j) + &
                                 tmp_denominator * (previous_energy_ss + shift_factor))
                             tmp_energy_term_os = &
                                 numerator(i,j)**(2.d0)/ (denominator(i,j) + &
                                 tmp_denominator * (previous_energy_os + shift_factor))
                             tmp_energy_term = tmp_energy_term_os + tmp_energy_term_ss
                             tmp_energy    = tmp_energy    + tmp_energy_term   
                             tmp_energy_ss = tmp_energy_ss + tmp_energy_term_ss
                             tmp_energy_os = tmp_energy_os + tmp_energy_term_os
                        endif
                      enddo ! j
                    enddo ! i
                else !if (n_spin.ne.1)
                    if (i_spin .eq. i_spin2) then
                        i_index = 0
                        do i=n_lumo(i_spin), n_states, 1
                          do j=i+1, n_states, 1
                            i_index = i_index + 1
                            if(myid.eq.MOD(i_index,n_tasks)) then
                                 tmp_denominator  = enhanced_factor * &
                                     arch_erfc(screening_factor*(denominator(i,j)))
                                 tmp_energy_term_ss = &
                                     (numerator(i,j)-numerator(j,i))**(2.0d0) / &
                                     (denominator(i,j) + &
                                     tmp_denominator*(previous_energy_ss + shift_factor))
                                 tmp_energy      = tmp_energy + tmp_energy_term_ss
                                 tmp_energy_ss   = tmp_energy_ss + tmp_energy_term_ss
                            endif
                          enddo ! j
                        enddo ! i
                    else ! if (i_spin .ne. i_spin2)
                        i_index = 0
                        do i=n_lumo(i_spin), n_states, 1
                          do j=n_lumo(i_spin2), n_states, 1
                            i_index = i_index + 1
                            if(myid.eq.MOD(i_index,n_tasks)) then
                                 tmp_denominator  = enhanced_factor * &
                                     arch_erfc(screening_factor*(denominator(i,j)))
                                 tmp_energy_term_os = &
                                     numerator(i,j)**(2.0d0) / (denominator(i,j) +&
                                     tmp_denominator*(previous_energy_os+shift_factor))
                                 tmp_energy      = tmp_energy + tmp_energy_term_os
                                 tmp_energy_os   = tmp_energy_os + tmp_energy_term_os
                            endif
                          enddo ! j
                        enddo ! i
                    endif
                endif ! if (n_spin.eq.1)
            endif
            if(use_mpi) then
                call MPI_ALLREDUCE(tmp_energy, tmp_energy_mpi, 1, &
                     MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_comm_global, mpierr)
                tmp_energy=tmp_energy_mpi 
                call MPI_ALLREDUCE(tmp_energy_ss, tmp_energy_mpi_ss, 1, &
                     MPI_DOUBLE_PRECISION, &
                     MPI_SUM, mpi_comm_global, mpierr)
                tmp_energy_ss=tmp_energy_mpi_ss 
                call MPI_ALLREDUCE(tmp_energy_os, tmp_energy_mpi_os, 1, &
                     MPI_DOUBLE_PRECISION, &
                     MPI_SUM, mpi_comm_global, mpierr)
                tmp_energy_os=tmp_energy_mpi_os 
            endif

            delta_energy_1 = delta_energy_0
            delta_energy_0 = tmp_energy - previous_energy
            if (delta_energy_0*delta_energy_1.lt.0.0d0) then
                tmp_energy     = (tmp_energy + previous_energy) * 0.5d0
                tmp_energy_ss  = (tmp_energy_ss + previous_energy_ss) * 0.5d0
                tmp_energy_os  = (tmp_energy_os + previous_energy_os) * 0.5d0
                delta_energy_0 = tmp_energy - previous_energy
            endif
            previous_energy    = tmp_energy
            previous_energy_ss = tmp_energy_ss
            previous_energy_os = tmp_energy_os
            if (myid.eq.0) then
                write(use_unit,'(2X,A,I3,A,E12.3,E12.3,E12.3,A,F16.8)') &
                    'The ',n_scf,'th iteration yields ',&
                    tmp_energy,tmp_energy_ss,tmp_energy_os,&
                    ' with the energy change of ',delta_energy_0
            endif
  end do SCF_LOOP

  final_energy = tmp_energy
  final_energy_ss = tmp_energy_ss
  final_energy_os = tmp_energy_os

end subroutine nsolver_iepa_mp2
!****** 
