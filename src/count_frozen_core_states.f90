!****s* FHI-aims/count_frozen_core_states()
!  NAME
!    count_frozen_core_states() 
!  SYNOPSIS

    subroutine count_frozen_core_states(n_low_state)

!  PURPOSE
!  This routine counts the number of frozen core states which do not take into
!  account in the post-KS calculation. 
!  "n_low_state" is returned as the lowest state taken effect on the post-KS 
!  calculation
!
!  USES

      use dimensions
      use geometry
      use species_data
      use runtime_choices
      use mpi_tasks
      use localorb_io,only : use_unit
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
      integer :: n_low_state 
      ! local variables
      character*20 :: species_temp
      integer :: i_atom
      logical :: flag_first_row, flag_second_row, flag_third_row, flag_four_row
      logical :: flag_five_row, flag_six_row
      logical :: flag_tm
      integer :: n_frozen_shell_1, n_frozen_shell_2, n_frozen_shell_curr
      integer :: n_tm
      real*8, dimension(3)   :: tmp_vec, tmp_vec_2
      flag_first_row  = .false.
      flag_second_row = .false.
      flag_third_row  = .false.
      flag_four_row   = .false.
      flag_five_row   = .false.
      flag_six_row    = .false.
      flag_tm         = .false.
      n_tm            = 0
      n_low_state     = 1  ! default one

      ! =========================================================================
      ! decouple the n_frozen_shell for main-group and transition-metal elements
      ! n_frozen_shell_1 for main-group ones
      ! n_frozen_shell_2 for transition-metals
      ! =========================================================================
      if (n_frozen_shell .gt. 10) then
          n_frozen_shell_1 = mod(n_frozen_shell, 10)
          n_frozen_shell_2 = n_frozen_shell/10
      else
          n_frozen_shell_1 = n_frozen_shell
          n_frozen_shell_2 = n_frozen_shell
      endif

      do i_atom=1, n_occ_atoms, 1
          species_temp = species_name(species(i_atom))
          flag_first_row  = (species_temp.eq.'H' ) .or. (species_temp.eq.'He')
          flag_second_row = (species_temp.eq.'Li') .or. (species_temp.eq.'Be') &
                       .or. (species_temp.eq.'B' ) .or. (species_temp.eq.'C' ) &
                       .or. (species_temp.eq.'N' ) .or. (species_temp.eq.'O' ) &
                       .or. (species_temp.eq.'F' ) .or. (species_temp.eq.'Ne')
          flag_third_row  = (species_temp.eq.'Na') .or. (species_temp.eq.'Mg') &
                       .or. (species_temp.eq.'Al') .or. (species_temp.eq.'Si') &
                       .or. (species_temp.eq.'P' ) .or. (species_temp.eq.'S' ) &
                       .or. (species_temp.eq.'Cl') .or. (species_temp.eq.'Ar')
          flag_four_row   = (species_temp.eq.'K')  .or. (species_temp.eq.'Ca') &
                       .or. (species_temp.eq.'Ga') .or. (species_temp.eq.'Ge') &
                       .or. (species_temp.eq.'As') .or. (species_temp.eq.'Se') &
                       .or. (species_temp.eq.'Br') .or. (species_temp.eq.'Kr') &
                       .or. (species_temp.eq.'Se') .or. (species_temp.eq.'Ti') &
                       .or. (species_temp.eq.'V')  .or. (species_temp.eq.'Cr') &
                       .or. (species_temp.eq.'Mn') .or. (species_temp.eq.'Fe') &
                       .or. (species_temp.eq.'Co') .or. (species_temp.eq.'Ni') &
                       .or. (species_temp.eq.'Cu') .or. (species_temp.eq.'Zn') 
          flag_five_row   = (species_temp.eq.'Rb') .or. (species_temp.eq.'Sr') &
                       .or. (species_temp.eq.'In') .or. (species_temp.eq.'Sn') &
                       .or. (species_temp.eq.'Sb') .or. (species_temp.eq.'Te') &
                       .or. (species_temp.eq.'I')  .or. (species_temp.eq.'Xe') &
                       .or. (species_temp.eq.'Y')  .or. (species_temp.eq.'Zr') &
                       .or. (species_temp.eq.'Nb') .or. (species_temp.eq.'Mo') &
                       .or. (species_temp.eq.'Tc') .or. (species_temp.eq.'Ru') &
                       .or. (species_temp.eq.'Rh') .or. (species_temp.eq.'Pd') &
                       .or. (species_temp.eq.'Ag') .or. (species_temp.eq.'Cd') 
          flag_six_row   =  (species_temp.eq.'Cs') .or. (species_temp.eq.'Ba') &
                       .or. (species_temp.eq.'Hf') .or. (species_temp.eq.'Ta') &
                       .or. (species_temp.eq.'W')  .or. (species_temp.eq.'Re') &
                       .or. (species_temp.eq.'Os') .or. (species_temp.eq.'Ir') &
                       .or. (species_temp.eq.'Pt') .or. (species_temp.eq.'Au') &
                       .or. (species_temp.eq.'Hg') .or. (species_temp.eq.'Tl') &
                       .or. (species_temp.eq.'Pb') .or. (species_temp.eq.'Bi') &
                       .or. (species_temp.eq.'Po') .or. (species_temp.eq.'At') &
                       .or. (species_temp.eq.'Rn')
          flag_tm        =  (species_temp.eq.'Se') .or. (species_temp.eq.'Ti') &
                       .or. (species_temp.eq.'V')  .or. (species_temp.eq.'Cr') &
                       .or. (species_temp.eq.'Mn') .or. (species_temp.eq.'Fe') &
                       .or. (species_temp.eq.'Co') .or. (species_temp.eq.'Ni') &
                       .or. (species_temp.eq.'Cu') .or. (species_temp.eq.'Zn') &
                       .or. (species_temp.eq.'Y')  .or. (species_temp.eq.'Zr') &
                       .or. (species_temp.eq.'Nb') .or. (species_temp.eq.'Mo') &
                       .or. (species_temp.eq.'Tc') .or. (species_temp.eq.'Ru') &
                       .or. (species_temp.eq.'Rh') .or. (species_temp.eq.'Pd') &
                       .or. (species_temp.eq.'Ag') .or. (species_temp.eq.'Cd') &
                       .or. (species_temp.eq.'Hf') .or. (species_temp.eq.'Ta') &
                       .or. (species_temp.eq.'W')  .or. (species_temp.eq.'Re') &
                       .or. (species_temp.eq.'Os') .or. (species_temp.eq.'Ir') &
                       .or. (species_temp.eq.'Pt') .or. (species_temp.eq.'Au') &
                       .or. (species_temp.eq.'Hg')

          if (flag_tm) then
              n_tm = n_tm + 1
              n_frozen_shell_curr = n_frozen_shell_2
          else
              n_frozen_shell_curr = n_frozen_shell_1
          endif

          if (flag_first_row) then
              n_low_state = n_low_state + 0
          else if (flag_second_row) then
              if (n_frozen_shell_curr .eq. 0) then
                  n_low_state = n_low_state + 0
              else if (n_frozen_shell_curr .eq. 1 .or. n_frozen_shell_curr .le. -1) then
                  n_low_state = n_low_state + 1
              else
                  n_low_state = n_low_state + 0
              endif
          else if (flag_third_row) then
              if (n_frozen_shell_curr .eq. 0) then
                  n_low_state = n_low_state + 0
              else if (n_frozen_shell_curr .eq. 1 .or. n_frozen_shell_curr .le. -2) then
                  n_low_state = n_low_state + 5
              else if (n_frozen_shell_curr .eq. 2 .or. n_frozen_shell_curr .eq. -1) then
                  n_low_state = n_low_state + 1
              else
                  n_low_state = n_low_state + 0
              endif
          else if (flag_four_row) then
              if (n_frozen_shell_curr .eq. 0) then
                  n_low_state = n_low_state + 0
              else if (n_frozen_shell_curr .eq. 1 .or. n_frozen_shell_curr .le. -3) then
                  n_low_state = n_low_state + 9
              else if (n_frozen_shell_curr .eq. 2 .or. n_frozen_shell_curr .eq. -2) then
                  n_low_state = n_low_state + 5
              else if (n_frozen_shell_curr .ge. 3 .or. n_frozen_shell_curr .eq. -1) then
                  n_low_state = n_low_state + 1 
              else
                  n_low_state = n_low_state + 0
              endif
          else if (flag_five_row) then
              if (n_frozen_shell_curr .eq. 0) then
                  n_low_state = n_low_state + 0
              else if (n_frozen_shell_curr .eq. 1 .or. n_frozen_shell_curr .le. -4) then
                  n_low_state = n_low_state + 18
              else if (n_frozen_shell_curr .eq. 2 .or. n_frozen_shell_curr .eq. -3) then
                  ! NOTE: for 4d-block elements, 3d orbitals are frozen as well for n_frozen_shell=2.
                  !       It leads to a core shell with 28 electrons
                  if (flag_tm) then
                      n_low_state = n_low_state + 14
                  else
                      n_low_state = n_low_state + 9
                  endif
              else if (n_frozen_shell_curr .eq. 3 .or. n_frozen_shell_curr .eq. -2) then
                  n_low_state = n_low_state + 5
              else if (n_frozen_shell_curr .ge. 4 .or. n_frozen_shell_curr .eq. -1) then
                  n_low_state = n_low_state + 1 
              else
                  n_low_state = n_low_state + 0
              endif
          else if (flag_six_row) then
              if (n_frozen_shell_curr .eq. 0) then
                  n_low_state = n_low_state + 0
              else if (n_frozen_shell_curr .eq. 1 .or. n_frozen_shell_curr .le. -5) then
                  n_low_state = n_low_state + 34 
              else if (n_frozen_shell_curr .eq. 2 .or. n_frozen_shell_curr .le. -4) then
                  ! NOTE: for 5d-block elements, 4d and 4f orbitals are frozen as well for n_frozen_shell=2.
                  !       It leads to a core shell with 60 electrons
                  if (flag_tm) then
                      n_low_state = n_low_state + 30
                  else
                      n_low_state = n_low_state + 18
                  endif
              else if (n_frozen_shell_curr .eq. 3 .or. n_frozen_shell_curr .eq. -3) then
                  n_low_state = n_low_state + 9
              else if (n_frozen_shell_curr .eq. 4 .or. n_frozen_shell_curr .eq. -2) then
                  n_low_state = n_low_state + 5
              else if (n_frozen_shell_curr .ge. 5 .or. n_frozen_shell_curr .eq. -1) then
                  n_low_state = n_low_state + 1 
              else
                  n_low_state = n_low_state + 0
              endif
          endif
      enddo
      if (myid .eq. 0) then
          write(use_unit,'(2X,A,I12)') &
              'First valence state in the frozen-core algorithm :', &
              n_low_state
          write(use_unit,'(2X,A20,A20,A20)') &
              '|                   ','Main-group', 'Transition-metal'
          write(use_unit,'(2X,A20,I20,I20)') &
              '| Atom number     ::', &
              n_occ_atoms-n_tm, n_tm
          write(use_unit,'(2X,A20,I20,I20)') &
              '| n_frozen_shell  ::', &
              n_frozen_shell_1, n_frozen_shell_2
      endif
end subroutine
