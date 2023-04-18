!   VB 2005:
!   Routine taken from http://www.netlib.org/go/
!
!   Copyright: My impression is that these routines come with a liberal (if nonexistant)
!   copyright provision, but netlib states that if in doubt, ask the original
!   author. In any case, it is publicly available for anyone to use.
!
!   Adapted by Ralf Gehrke 2005 to work for our purposes
!
!   Description from netlib
!   lib      ../go
!   for      Golden Oldies:  widely used,  but not in standard libraries.
!   #      Nominations welcome!
!   rel      excellent
!   age      old
!   editor      Eric Grosse
!   master      netlib.bell-labs.com
!
!  To get d1mach, mail netlib
!       send d1mach from core
!  VB: In my opinion, we already have the machine accuracy in question, but we can always get it
!      from netlib by downloading the dependencies also.
!
      subroutine zeroin_v2(ax, bx, fax, fbx, tol, &
           occupation_acc, max_zeroin, &
           n_states, KS_eigenvalue, &
           n_electrons, chemical_potential, diff_electrons, &
           occ_numbers, i_counter)


      use localorb_io, only: use_unit
      implicit none

! imported variables

! input
      real*8, intent(in) :: ax
      real*8, intent(in) :: bx
      real*8, intent(in) :: fax
      real*8, intent(in) :: fbx
      real*8, intent(in) :: tol

      real*8, intent(in) :: occupation_acc
      integer, intent(in) :: max_zeroin

      integer, intent(in) :: n_states
!      integer, intent(in) :: n_spin
      real*8, dimension(n_states), intent(in) :: KS_eigenvalue
      real*8, intent(in) :: n_electrons

! output
      real*8, dimension(n_states), intent(out) :: occ_numbers
      real*8, intent(out) :: chemical_potential
      real*8, intent(out) :: diff_electrons
      integer, intent(inout) :: i_counter



!
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
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
      real*8  a, b, c, d, e, eps, fa, fb, fc, tol1, xm, p, &
           q, r, s
      real*8  dabs, d1mach

 10   eps = d1mach(4)
      tol1 = eps+1.0d0
!
      a  = ax
      b  = bx
      fa = fax
      fb = fbx

!test
!      write(use_unit,*) "After 10: "
!      write(use_unit,*) " a = ", a
!      write(use_unit,*) " b = ", b
!      write(use_unit,*) "fa = ", fa
!      write(use_unit,*) "fb = ", fb
!test end

!     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
      if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20
         write(use_unit,2500)
2500     format(1x,'f(ax) and f(bx) do not have different signs,', &
                   ' zeroin is aborting')
         return
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if (dabs(fc).ge.dabs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0d0*eps*dabs(b)+0.5d0*tol
      xm = 0.5d0*(c-b)
      if (abs(fb) .le. occupation_acc) go to 150
!
! see if a bisection is forced
!
      if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
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
   60 q=fa/fc
      r=fb/fc
      p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
      q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
   70 if (p.le.0.0d0) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.(p.ge. &
      dabs(0.5d0*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if (dabs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0d0) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
!  140 fb=f(b)
 140  continue
      call check_norm_v2(b, KS_eigenvalue, n_electrons, &
           occ_numbers, fb, i_counter)
!test
!      write(use_unit,*) "i_counter = ", i_counter
!      write(use_unit,*) "n_electrons: ", n_electrons
!      write(use_unit,*) "occ_numbers: ", occ_numbers
!test end

      if (i_counter .gt. max_zeroin) then
!       force end despite non-convergence
        go to 150
      end if
      if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20
      go to 30
  150 chemical_potential = b
      diff_electrons     = fb
      return
      end
