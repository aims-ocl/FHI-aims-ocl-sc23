!  module conjugate gradient to encapsulate all
!  necessary subroutines for a local relaxation
!  of the structure, energies and forces provided
!
!  R.Gehrke (2005)
!

module conjugate_gradient

  use vector_array_class
  use map_to_one_dimension, allocate_maptoone => allocate, deallocate_maptoone => deallocate, &
       coords_maptoone => coords
  use cluster
  use lennard_jones

  implicit none

  private

  type (vector_array) :: old_forces
  type (vector_array) :: forces
  type (vector_array) :: forces_try
  type (vector_array) :: F0

  ! variables
  public :: &
       trial_step, &
       trial_step_initial, &
       max_force_allowed, &
       max_cosine, &
       max_line_ts, &
       max_line_steps, &
       zeroin_tolerance, &
       max_n_cg_loops, &
       min_progress, &
       max_progress, &
       sw_prerelax, &
       max_pr_force, &
       damping_factor, &
       step_increase, &
       step_decrease, &
       max_increase_ts, &
       pr_step_decrease, &
       max_pr_ts

  ! subroutines
  public :: &
       optimize, &
       allocate, &
       deallocate

  ! maximum force component
  real*8 :: max_force_allowed
  ! force threshold for prerelaxation
  real*8 :: max_pr_force
  ! max cosine between forces and search direction for line-minimization
  real*8 :: max_cosine

  ! trial step in line minimization
  ! trial_step is self-adapting during one conjugate-gradient scheme
  ! trial_step_initial is the initial guess which stays constant
  real*8 :: trial_step
  real*8 :: trial_step_initial
  ! maximum displacement in line minimization for first trial step
  real*8 :: max_line_ts
  ! maximum displacement in line minimization when stepwidth is increased
  real*8 :: max_increase_ts
  ! minimal progress of line force in line minimization
  real*8 :: min_progress
  ! maximum progress of line force in line minimizatoin
  real*8 :: max_progress
  ! step increase/decrease for line minimization
  real*8 :: step_increase
  real*8 :: step_decrease

  ! step width for prerelaxation
  real*8 :: sw_prerelax
  ! maximum displacement for prerelaxation
  real*8 :: max_pr_ts

  real*8 :: damping_factor
  real*8 :: pr_step_decrease

  ! maximum number of cg-loops
  integer :: max_n_cg_loops
  ! maximum number of steps in line minimization
  integer :: max_line_steps
  
  real*8 :: zeroin_tolerance
  
contains

  subroutine allocate()
    call allocate_array(old_forces, n_atoms)
    call allocate_array(forces, n_atoms)
    call allocate_array(forces_try, n_atoms)
    call allocate_array(F0, n_atoms)
    ! set zeroin_tolerance so high, that it's never fulfilled => max_line_force is the determining
    ! criterium also for zeroin-method
    zeroin_tolerance = 1e-20
    trial_step = trial_step_initial
    open (25, FILE="optimization.log")

  end subroutine allocate

  subroutine deallocate()
    call deallocate_array(old_forces)
    call deallocate_array(forces)
    call deallocate_array(forces_try)
    call deallocate_array(F0)
    
    close (25)

  end subroutine deallocate
  
  subroutine damped_newton(coords, forces, energy, status)
     
     implicit none
     
     ! imported variables
     
     ! input
     type (vector_array), intent(inout)  :: coords
     
     ! output
     type (vector_array), intent(inout) :: forces
     real*8, intent(out) :: energy
     integer, intent(out) :: status
    
     ! local variables
     real*8 :: max_force_try
     real*8 :: max_force
     type (vector_array) :: coords_try
     type (vector_array) :: coords_old
     type (vector_array) :: forces_try
     real*8 :: sw
     real*8 :: line_force
     real*8 :: line_force_old
     real*8 :: energy_try
     
     ! functions
     real*8 :: get_max_force
     real*8 :: get_line_force
     
     ! counter
     integer :: i_step
     integer :: i_atom
    
     call allocate_array(coords_try, n_atoms)
     call allocate_array(coords_old, n_atoms)
     call allocate_array(forces_try, n_atoms)
     
     forces_try = forces
     max_force = get_max_force(forces, n_atoms)
     line_force = get_line_force(forces, forces, n_atoms)
     line_force_old = 0.5d0 * line_force
     coords_old = coords
     status = 0
     ! do prerelaxation until forces are small AND curvature of PES has the right sign
     do while ((max_force .gt. max_pr_force) .or. (abs(line_force) .gt. abs(line_force_old)))
        max_force = get_max_force(forces, n_atoms)
        sw = max_pr_ts / max_force
        
        if (sw .gt. sw_prerelax) then
           sw = sw_prerelax
        end if
        do i_atom = 1, n_atoms, 1
           coords_try%coords(i_atom) = coords%coords(i_atom) + forces%coords(i_atom) * sw + &
                (coords%coords(i_atom) - coords_old%coords(i_atom)) * damping_factor
        end do
        call get_forces_and_energy(coords_try, energy_try, forces_try, status)
        if (status .gt. 0) then
           coords = coords_try
           forces = forces_try
           energy = energy_try
           return
        end if
        max_force_try = get_max_force(forces_try, n_atoms)
        write (25,'(A,1X,F18.6,1X,A,1X,E14.6,1X,A,1X,E14.6)') "energy=", energy_try, "[eV] max_force=", max_force_try, &
             "[eV/Ang] sw=", sw  
        do while (energy_try .gt. energy) 
           sw = sw * pr_step_decrease
           do i_atom = 1, n_atoms, 1
              coords_try%coords(i_atom) = coords%coords(i_atom) + forces%coords(i_atom) * sw + &
                   (coords%coords(i_atom) - coords_old%coords(i_atom)) * damping_factor
           end do
           call get_forces_and_energy(coords_try, energy_try, forces_try, status)
           if (status .gt. 0) then
              coords = coords_try
              return
           end if
           max_force_try = get_max_force(forces_try, n_atoms)
           write (25,'(A,1X,F18.6,1X,A,1X,E14.6,1X,A,1X,E14.6)') "energy=", energy_try, "[eV] max_force=", max_force_try, "[eV/Ang] sw=", sw  
        end do
        
        coords_old = coords
        coords = coords_try
        forces = forces_try
        energy = energy_try
        line_force_old = line_force
        line_force = get_line_force(forces, forces, n_atoms)
        max_force = max_force_try
        write (25,'(2(A,1X,E14.6))') "line_force= ", line_force, " line_force_old= ", line_force_old
     end do
     if (max_force .lt. max_force_allowed) then
        status = 3
     end if

     call deallocate_array(coords_try)
     call deallocate_array(coords_old)
     call deallocate_array(forces_try)
     
   end subroutine damped_newton

   subroutine line_minimization(coords, forces, energy, i_counter, status)
     
     implicit none
     
     ! imported variables
     
     ! input
     type (vector_array), intent(inout) :: coords
     type (vector_array), intent(inout) :: forces
     
     ! output
     real*8 :: energy
     real*8 :: lineforce
     integer, intent(inout) :: i_counter
     integer, intent(out) :: status
     
     ! local variables
     logical :: success
     logical :: bracketed
     real*8 :: a
     real*8 :: b
     real*8 :: Fa
     real*8 :: Fb
     real*8 :: l
     real*8 :: l0
     real*8 :: l1
     real*8 :: lmin
     real*8 :: Emin
     real*8 :: FLin0
     real*8 :: FLin1
     real*8 :: FLin_try
     real*8 :: energy_try
     real*8 :: cosine_try
     real*8 :: max_force_try
     real*8 :: max_force
     real*8 :: cosine
     real*8 :: cosine_old
     real*8 :: E0
     real*8 :: E_zero
     real*8 :: ts
     real*8 :: max_search_direction
     
     ! functions
     real*8 :: get_max_force
     
     ! counter
     integer :: i_atom
     integer :: i_loop
     
     zero_coords = coords
     success   = .true.
     bracketed = .false.
     status = 0
     
     ! just to make sure that first energy and forces are not lost if first
     ! get_forces_and_energy_call leads to status = 0
     energy_try = energy
     forces_try = forces
     
     l0 = 0
     a  = 0
     F0 = forces
     E0 = energy
     E_zero = E0
     i_loop = 0
     ! !!
     ! !! step ONE of line minimization (adaptive trial step)
     ! !!
     FLin0 = 0
     do i_atom = 1, n_atoms, 1
        FLin0 = FLin0 - forces%coords(i_atom)%x * search_direction%coords(i_atom)%x
        FLin0 = FLin0 - forces%coords(i_atom)%y * search_direction%coords(i_atom)%y
        FLin0 = FLin0 - forces%coords(i_atom)%z * search_direction%coords(i_atom)%z
     end do
     
     Fa = FLin0
     
     i_loop   =  i_loop    + 1
     i_counter = i_counter + 1
     
     max_search_direction = get_max_force(search_direction, n_atoms)
     l = max_line_ts / max_search_direction
     
     if (l .gt. trial_step) then
        l = trial_step
     end if
     ts = l
     call evaluate(l, forces_try, FLin_try, energy_try, cosine_try, max_force_try, status)
     if (status .eq. 0) then
        write (25,'(A,1X,F18.6,1X,A,1X,F10.4,1X,A,1X,E14.6,1X,A,1X,E14.6,1X,A)') "energy=", energy_try, "[eV] cosine=", cosine_try, &
             "max_force=", max_force_try, "[eV/Ang] l=", l, "(ts)" 
     end if
     ! trial step, step_width is adaptive
     do while ((status .eq. 0) .and. (dabs(cosine_try) .gt. max_cosine) .and. (max_force_try .gt. max_force_allowed) &
          .and. (energy_try .gt. energy))
        ! if forces increased too much, decrease stepwidth
        i_loop = i_loop + 1
        i_counter = i_counter + 1
        if (i_loop .gt. max_line_steps) then
           status = 2
           exit
        end if
        l = l * step_decrease
        trial_step = l
        ts = ts * step_decrease
        call evaluate(l, forces_try, FLin_try, energy_try, cosine_try, max_force_try, status)
        write (25,*) energy_try, energy
        if (status .gt. 0) then
           exit
        end if
        write (25,'(A,1X,F18.6,1X,A,1X,F10.4,1X,A,1X,E14.6,1X,A,1X,E14.6,1X,A)') "energy=", energy_try, "[eV] cosine=", cosine_try, &
             "max_force=", max_force_try, "[eV/Ang] l=", l, "(ts)"
     end do
     
     FLin1 = FLin_try
     Fb    = FLin1
     energy = energy_try
     forces = forces_try
     cosine = cosine_try
     max_force = max_force_try
     lmin = l
     l1   = l
     b    = l
     ! !!
     ! !! step TWO of line minimization (bracketing of minimum)
     ! !!
     if ((FLin0 * FLin1) .lt. 0) then
        bracketed = .true.
     else
        bracketed = .false.
     end if
     ! first, interval has to be found, so make trial steps or cubic interpolations until
     ! root is bracketed
     do while ((success .or. .not.(bracketed)) .and. (abs(cosine) .gt. max_cosine) .and. (max_force .gt. max_force_allowed) .and. (status .eq. 0)) 
        
        ! first find a trial lmin
        
        ! as long as cubic_interpolations make sense, do them
        if (success) then
          call cubic_interpolation(l0, l1, E0, energy, FLin0, FLin1, lmin, Emin, success)
          if (success) then
             i_loop = i_loop + 1
             i_counter = i_counter + 1
             if (i_loop .gt. max_line_steps) then
                status = 2
                exit
             end if
             call evaluate(lmin, forces_try, FLin_try, energy_try, cosine_try, max_force_try, status)
             if (status .gt. 0) then
                exit
             end if
             write (25,'(A,1X,F18.6,1X,A,1X,F10.4,1X,A,1X,E14.6,1X,A,1X,E14.6,1X,A)') "energy=", energy_try, "[eV] cosine=", cosine_try, "max_force=", &
                  max_force_try, "[eV/Ang] l=", lmin, "(ce)"
          end if
       end if
       ! if interval not yet found and cubic interpolation useless, make further trial steps
       if ((.not.bracketed) .and. (.not.success)) then
          ! since cubic interpolation didn't make sense => a further trial step
          lmin = l1 + ts
          i_loop = i_loop + 1
          i_counter = i_counter + 1
          if (i_loop .gt. max_line_steps) then
             status = 2
             exit
          end if
          call evaluate(lmin, forces_try, FLin_try, energy_try, cosine_try, max_force_try, status)
          if (status .gt. 0) then
             exit
          end if
          write (25,'(A,1X,F18.6,1X,A,1X,F10.4,1X,A,1X,E14.6,1X,A,1X,E14.6,1X,A)') "energy=", energy_try, "[eV] cosine=", cosine_try, "max_force=", &
               max_force_try, "[eV/Ang] l=", lmin, "(ts)"
       end if
       ! energy increases though force decreases => far away from harmonic region =>
       ! cubic extrapolation useless
       if ((dabs(FLin_try) .lt. dabs(FLin1)) .and. (energy_try .gt. energy)) then
          if (success) then
             write (25,*) "too far from harmonic region!"
             success = .false.
             cycle
          end if
       end if
       ! forces increase too much => cubic interpolation not successful or decrease trial step
       if (dabs(FLin_try / FLin1) .gt. max_progress) then
          write (25,*) "forces increase too much!"
          if (success) then
             success = .false.
          else
             ts = ts * step_decrease
          end if
          ! discard that step
          cycle 
       end if
       ! forces change too slowly => cubic interpolation not efficient enough or increase trial step
       ! if ((dabs(FLin1 / FLin_try) .lt. min_progress) .and. (dabs(FLin_try / FLin1) .lt. min_progress)) then
       ! in harmonic region ...                            ... away from harmonic region ...   
       if (dabs(FLin1 / FLin_try) .lt. min_progress) then
 !         write (25,*) "progress too slow!"
          if (success) then
             success = .false.
          else
             ts = ts * step_increase
             ! beware that you don't displace too much
             if (ts .gt. max_increase_ts / max_force_try) then
                ts = max_increase_ts / max_force_try
             end if
          end if
       end if
       l0 = l1
       l1 = lmin
       FLin0 = FLin1
       FLin1 = FLin_try
       E0 = energy
       energy = energy_try
       cosine_old = cosine
       cosine = cosine_try
       F0 = forces
       forces = forces_try
       max_force = max_force_try
       ! keep track of interval for the root
       call check_bound(a, b, Fa, Fb, l1, FLin1, bracketed, success)
    end do
    ! !!
    ! !! step THREE of line minimization (zeroin)
    ! !!
    ! cubic interpolations doesn't make sense any more, root is bracketed, if the lineforce is still
    ! big, let's start zeroin...
    if (status .eq. 0) then
       ! there are still forces in the force-file, so let's go on if necessary
       if ((abs(cosine) .gt. max_cosine) .and. (max_force .gt. max_force_allowed)) then
          call zeroin(a, b, Fa, Fb, lmin, FLin1, forces, energy, cosine, max_force, i_loop, status)
       else
          ! regularly optimized
          status = 3
       end if
    else
       ! no forces left, so new localorb-call required
       ! store previous forces and energy for opt_data.out
       energy = energy_try
       forces = forces_try 
    end if
    
    do i_atom = 1, n_atoms, 1
       coords%coords(i_atom)%x = coords%coords(i_atom)%x + lmin * search_direction%coords(i_atom)%x
       coords%coords(i_atom)%y = coords%coords(i_atom)%y + lmin * search_direction%coords(i_atom)%y
       coords%coords(i_atom)%z = coords%coords(i_atom)%z + lmin * search_direction%coords(i_atom)%z
    end do
    
  end subroutine line_minimization
  
  subroutine cubic_interpolation(l0, l1, E0, E1, F0, F1, lmin, Emin, success)

!  imported variables

!  input
    real*8, intent(in) :: l0 ! first function argument
    real*8, intent(in) :: l1 ! second function argument
    real*8, intent(in) :: E0 ! first function value
    real*8, intent(in) :: E1 ! second function value
    real*8, intent(in) :: F0 ! first derivative
    real*8, intent(in) :: F1 ! second derivative

!  output
    real*8, intent(out) :: lmin
    real*8, intent(out) :: Emin
    logical, intent(out) :: success
    
!  local variables
    real*8 :: Det
    real*8 :: Da
    real*8 :: Db
    real*8 :: a
    real*8 :: b
    real*8 :: c
    real*8 :: d
    real*8 :: s
    real*8 :: l

    l = l1 - l0

    success = .true.

    Det  = - l*l*l*l
    Da   = 2*l*(E1 - E0 - l*F0) - l*l*(F1 - F0)
    Db   = l*l*l*(F1 - F0) - 3*l*l*(E1 - E0 - l*F0)

    a = Da / Det
    b = Db / Det
    c = F0
    d = E0

    if (abs(a) .gt. 1E-10) then
       s = b*b / (9*a*a) - c / (3*a)

       if (s .gt. 0) then
          lmin = - b / (3*a) + sqrt(s)
          if ( (6*a*lmin + b) .lt. 0) then   ! f''(x) < 0 -> maximum instead of minimum was found
             lmin = - b / (3*a) - sqrt(s)
          end if
       else  ! no minimum could be found
          success = .false.
       end if
    else   ! quadratic fit was found
       if (2*b .gt. 0) then
          lmin = - c / (2*b) 
       else   ! no minimum could be found
          success = .false.
       end if
    end if
    if (success) then
       Emin = a*lmin**3 + b*lmin**2 + c*lmin + d 
       lmin = lmin + l0
    end if
    
  end subroutine cubic_interpolation

  subroutine zeroin(ax, bx, fax, fbx, lmin, line_force, forces, energy, cosine, max_force, i_counter, status)

    implicit none

    ! imported variables
    
    ! input
    real*8, intent(in) :: ax
    real*8, intent(in) :: bx
    real*8, intent(in) :: fax
    real*8, intent(in) :: fbx 

    ! output
    real*8, intent(out) :: lmin
    type (vector_array), intent(inout) :: forces
    real*8, intent(out) :: line_force
    integer, intent(inout) :: i_counter
    real*8, intent(inout) :: energy
    integer, intent(out) :: status
    real*8, intent(out) :: cosine
    real*8, intent(inout) :: max_force

    !
    ! a zero of the function  f(x)  is computed in the interval ax,bx .
    !
    ! input..
    !
    !  ax     left endpoint of initial interval
    !  bx     right endpoint of initial interval
    !  f      function subprogram which evaluates f(x) for any x in
    !         the interval  ax,bx
    !  zeroin_tolerance  desired length of the interval of uncertainty of the
    !         final result (.ge.0.)
    !
    !  output..
    !
    !  zeroin abscissa approximating a zero of  f  in the interval ax,bx
    !
    !      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
    !  this is checked, and an error message is printed if this is not
    !  satisfied.   zeroin  returns a zero  x  in the given interval
    !  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
    !  the  relative machine precision defined as the smallest representable
    !  number such that  1.+macheps .gt. 1.
    !      this function subprogram is a slightly  modified  translation  of
    !  the algol 60 procedure  zero  given in  richard brent, algorithms for
    !  minimization without derivatives, prentice-hall, inc. (1973).
    !
    real*8 ::  a, b, c, d, e, eps, fa, fb, fc, tol1, xm, p, q, r, s
    real*8 ::  dabs, d1mach
    status = 0
10  eps = d1mach(4)
    tol1 = eps+1.0d0
    a  = ax
    b  = bx
    fa = fax
    fb = fbx
    ! check that f(ax) and f(bx) have different signs
    if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
    if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20
    return
20  c  = a
    fc = fa
    d  = b - a
    e  = d
30  if (dabs(fc).ge.dabs(fb)) go to 40
    a=b
    b=c
    c=a
    fa=fb
    fb=fc
    fc=fa
40  tol1=2.0d0*eps*dabs(b)+0.5d0 * zeroin_tolerance
    xm = 0.5d0*(c-b)
!    if ((dabs(xm).le.tol1).or.(fb.eq.0.0d0)) go to 150
    if ((abs(cosine) .le. max_cosine) .or. (max_force .le. max_force_allowed)) then
       ! regularly optimized
       status = 3
       go to 150
    end if
    if (i_counter .gt. max_line_steps) then
       status = 2
       go to 150
    end if
!
! see if a bisection is forced
!
    if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50
    d=xm
    e=d
    go to 110
50  s=fb/fa
    if (a.ne.c) go to 60
!
! linear interpolation
!
    p=2.0d0*xm*s
    q=1.0d0-s
    go to 70
!
! inverse quadratic interpolation
!
60  q=fa/fc
    r=fb/fc
    p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
    q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
70  if (p.le.0.0d0) go to 80
    q=-q
    go to 90
80  p=-p
90  s=e
    e=d
    if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.(p.ge.dabs(0.5d0*s*q))) go to 100
     d=p/q
     go to 110
100  d=xm
     e=d
110  a=b
     fa=fb
     if (dabs(d).le.tol1) go to 120
     b=b+d
     go to 140
120  if (xm.le.0.0d0) go to 130
     b=b+tol1
     go to 140
130  b=b-tol1
140  i_counter = i_counter + 1
     call evaluate(b, forces, fb, energy, cosine, max_force, status)
     if (status .gt. 0) then
        ! store previous line-force
        line_force = fb
        lmin = b
        return
     end if
     write (25,'(A,1X,F18.6,1X,A,1X,F10.4,1X,A,1X,E14.6,1X,A,1X,E14.6,1X,A)') "energy=", energy, "[eV] cosine=", cosine, "max_force=", max_force, &
          "[eV/Ang] l=", b, "(z)"
     if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20
     go to 30
150  lmin = b
     line_force = fb
     return

   end subroutine zeroin

   subroutine check_bound(l0, l1, F0, F1, l, F, bracketed, success)

!  imported variables

!  input
     real*8, intent(inout) :: l0
     real*8, intent(inout) :: l1
     real*8, intent(inout) :: F0
     real*8, intent(inout) :: F1
     real*8, intent(in) :: l
     real*8, intent(in) :: F
     logical, intent(inout) :: bracketed
     logical, intent(inout) :: success
     
     if (bracketed) then
! new value inside hitherto found interval?
        if ((l0 - l) * (l - l1) .gt. 0) then
           if ((F0 * F) .lt. 0) then
              F1 = F
              l1 = l
           else 
              F0 = F
              l0 = l
           end if
        else 
           success = .false.
        end if
     else 
! not yet bracketed
        if ((F0 * F) .lt. 0) then
! zero bracketed!
           if (abs(F0 * F) .lt. abs(F1 * F)) then
! zero between F0 and F
              F1 = F
              l1 = l
           else 
! zero between F and F1
              F0 = F
              l0 = l
           end if
           bracketed = .true.
        else 
           F0 = F1
           l0 = l1
           F1 = F
           l1 = l
        end if
     end if

   end subroutine check_bound
   
   subroutine optimize(try_coords, coords, energy, forces, max_force, i_line, i_evaluations, status)

     ! imported variables

     ! input
     type (vector_array), intent(in) :: try_coords

     ! output
     type (vector_array), intent(out) :: coords
     real*8, intent(out) :: energy
     real*8, intent(out) :: max_force
     integer, intent(inout) :: i_evaluations
     integer, intent(inout) :: i_line
     integer, intent(inout) :: status
     type (vector_array), intent(out) :: forces

     ! local variables
     type (vector_array) :: old_forces
     
     real*8 :: c1
     real*8 :: c2
     real*8 :: gamma
     
     ! functions
     real*8 :: get_max_force

     ! counter
     integer :: i_atom
     
     i_line        = 0
     i_evaluations = 0
     status = 0
     coords = try_coords
     
     call allocate_maptoone(n_atoms)
     call allocate_array(old_forces, n_atoms)
     
     write (25,*) "start prerelaxation with damped newton ..."
     call get_forces_and_energy(coords, energy, forces, status)
     if (status .eq. 0) then
        max_force = get_max_force(forces, n_atoms)
        write (25,'(A,1X,F18.6,1X,A,1X,E14.6,1X,A)') "energy=", energy, "[eV] max_force=", max_force, "[eV/Ang]" 
        call damped_newton(coords, forces, energy, status)
        max_force = get_max_force(forces, n_atoms)
     end if
     ! search direction in map_to_one_dimension
     if (status .eq. 0) then
        search_direction = forces
        max_force = get_max_force(forces, n_atoms)
     end if
     write (25,*) "start conjugate-gradient ..."
     trial_step = trial_step_initial
     do while ((max_force .gt. max_force_allowed) .and. ((status .eq. 0) .or. (status .eq. 3)))
        write (25,'(A,1X,I4)') "line minimization #", i_line
        old_forces = forces
        i_line = i_line + 1
        if (i_line .gt. max_n_cg_loops) then
           status = 2
           exit
        end if
        call line_minimization(coords, forces, energy, i_evaluations, status)
        max_force = get_max_force(forces, n_atoms)
        if ((status.eq.0) .or. (status .eq. 3)) then
           write (25,'(A,1X,I4)') "after line minimization #", i_line - 1
           write (25,'(A,1X,F18.6,1X,A,1X,E14.6,1X,A)') "energy=", energy, "[eV] max_force=", max_force, "[eV/Ang]" 
        end if
        c1 = 0.d0
        c2 = 0.d0
        do i_atom = 1, n_atoms, 1
           ! polak-ribiere
           c1 = c1 + (forces%coords(i_atom)%x - old_forces%coords(i_atom)%x) * forces%coords(i_atom)%x
           c1 = c1 + (forces%coords(i_atom)%y - old_forces%coords(i_atom)%y) * forces%coords(i_atom)%y
           c1 = c1 + (forces%coords(i_atom)%z - old_forces%coords(i_atom)%z) * forces%coords(i_atom)%z

           ! fletcher-reeves
           !c1 = c1 + forces%coords(i_atom)%x * forces%coords(i_atom)%x
           !c1 = c1 + forces%coords(i_atom)%y * forces%coords(i_atom)%y
           !c1 = c1 + forces%coords(i_atom)%z * forces%coords(i_atom)%z
           
           c2 = c2 + old_forces%coords(i_atom)%x * old_forces%coords(i_atom)%x
           c2 = c2 + old_forces%coords(i_atom)%y * old_forces%coords(i_atom)%y
           c2 = c2 + old_forces%coords(i_atom)%z * old_forces%coords(i_atom)%z
        end do
        
        gamma = c1 / c2

        do i_atom = 1, n_atoms, 1
           search_direction%coords(i_atom)%x = forces%coords(i_atom)%x + gamma * search_direction%coords(i_atom)%x
           search_direction%coords(i_atom)%y = forces%coords(i_atom)%y + gamma * search_direction%coords(i_atom)%y
           search_direction%coords(i_atom)%z = forces%coords(i_atom)%z + gamma * search_direction%coords(i_atom)%z
        end do
     end do
     if ((status.eq.0) .or. (status .eq. 3)) then
        write (25,'(A)') "done."
     end if
     call deallocate_maptoone()
     call deallocate_array(old_forces)

   end subroutine optimize

end module conjugate_gradient
