!***h* FHI-aims/statistics
!  NAME
!    statistics
!  SYNOPSIS
module statistics
!  PURPOSE
!    This module is used to generate a set of statistical averages for a
!    provided array, suitable for condensing large arrays into smaller,
!    easier-to-parse summaries.
!
!    Currently, only 1D integer arrays are supported, but the module has been
!    deliberately structured to liberally use module procedures, so extending
!    to other data types should be as simple as copy/pasting and modifying a
!    few lines.
!
!    I've tested this module out versus R and numpy for sample sizes of 1, 2,
!    200, 2000, and 20000 drawn from linear, Gaussian (via Box-Muller), and
!    Dirac probability distributions.  Still, this is custom-rolled code by a
!    non-statistician, so...
!  USES
  implicit none
!  AUTHOR
!    William Huhn (Duke University)
!  HISTORY
!    January 2018 - Created
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  SOURCE
  private

  public :: reset_stats
  public :: calc_stats_array
  public :: bubblesort_array
  public :: calc_std_dev_array
  public :: calc_mode_array
  public :: calc_median_array
  public :: calc_fractiles_array
  public :: calc_histogram_array

  public :: stats_int
  type stats_int
    logical  :: is_valid   ! Whether the current instance is initialized
    integer  :: num_values ! Number of values in data set
    integer  :: sum_values ! Sum of values in data set
    integer  :: min_value  ! Minimum value in data set
    integer  :: max_value  ! Maximum value in data set
    real*8   :: mean_value ! Arthimetic mean value in data set
    real*8   :: median_value  ! Median value in data set
    integer  :: mode_value ! Mode in data set
                           ! If there's more than one mode, the value reported
                           ! will be the lowest value.
    integer  :: mode_count ! Number of times mode appeared in data set
    real*8   :: std_dev    ! Standard deviations of data set
    real*8, allocatable :: fractiles(:) ! Fractiles for data set
    integer  :: n_fractiles ! # divisions used when calculating fractile
                            ! The fractiles array will have one additional
                            ! element, because it includes the 0th fractile,
                            ! i.e. the minimum value
    integer, allocatable :: histogram(:) ! Histogram for data set
    real*8, allocatable :: histogram_bounds(:) ! Bounds of histograms
    integer  :: n_bins ! # bins used when calculating histogram
  end type stats_int

  interface reset_stats
       module procedure reset_stats_int
  end interface

  interface calc_stats_array
       module procedure calc_stats_1D_int_array
  end interface

  interface bubblesort_array
       module procedure bubblesort_1D_int_array
  end interface

  interface calc_std_dev_array
       module procedure calc_std_dev_1D_int_array
  end interface

  interface calc_mode_array
       module procedure calc_mode_1D_int_array
  end interface

  interface calc_median_array
       module procedure calc_median_1D_int_array
  end interface

  interface calc_fractiles_array
       module procedure calc_fractiles_1D_int_array
  end interface

  interface calc_histogram_array
       module procedure calc_histogram_1D_int_array
  end interface

contains

  subroutine reset_stats_int(stats)
    implicit none

    type(stats_int), intent(out) :: stats

    stats%is_valid = .false.
    stats%num_values = -1
    stats%sum_values = 0
    stats%min_value = 0
    stats%max_value = 0
    stats%mean_value = 0.0d0
    stats%median_value = 0.0d0
    stats%mode_value = 0
    stats%mode_count = -1
    stats%std_dev = 0.0d0
    stats%n_fractiles = -1
    if (allocated(stats%fractiles)) deallocate(stats%fractiles)
    stats%n_bins = -1
    if (allocated(stats%histogram)) deallocate(stats%histogram)
    if (allocated(stats%histogram_bounds)) deallocate(stats%histogram_bounds)
  end subroutine reset_stats_int

!-------------------------------------------------------------------------------
!****s* statistics/calc_stats_1D_int_array
!  NAME
!    calc_stats_1D_int_array
!  SYNOPSIS
  subroutine calc_stats_1D_int_array(array, stats, n_fractiles, n_bins, &
       hist_start_in, hist_end_in, is_array_sorted)
!  PURPOSE
!    Takes in an array of integers and outputs the Stat 101 suite of values for
!    them.  Calculates histograms and fractiles if user supplies the
!    appropriate values.
!  USES
    use mpi_tasks, only: aims_stop
    implicit none
!  ARGUMENTS
    integer, intent(in)  :: array(:)
    type(stats_int), intent(out)   :: stats
    integer, intent(in),  optional :: n_fractiles
    integer, intent(in),  optional :: n_bins
    real*8,  intent(in),  optional :: hist_start_in
    real*8,  intent(in),  optional :: hist_end_in
    logical, intent(in),  optional :: is_array_sorted
!  INPUTS
!    o array - the array of values to find statistical averages on
!    o n_fractiles - the number of fractiles to compute (optional).  Required to
!      compute fractiles.
!    o n_bins - the number of bins to use in histograms (optional).  Required to
!      compute histograms.  If both hist_start_in and hist_end_in are not
!      specified, the starting and ending values used will be the minimum and
!      maximum values of the array.
!    o hist_start_in - the starting value for the histogram (optional)
!    o hist_end_in - the ending value for the histogram (optional)
!    o is_array_sorted - whether the array is already sorted, and we can skip
!      the creation of a temporary sorted array
!  OUTPUT
!    o stats - a derived type containing statistical averages for array
! SOURCE
    integer, allocatable :: sorted_array(:)
    integer :: mode_value_temp, mode_count_temp
    logical :: array_is_sorted
    real*8  :: hist_start, hist_end

    character(*), parameter :: func = "calc_stats_1D_int_array"

    call reset_stats(stats)

    ! Calculate statistics that don't require the array be sorted
    stats%num_values = SIZE(array,1)
    stats%sum_values = SUM(array)
    stats%min_value  = MINVAL(array)
    stats%max_value  = MAXVAL(array)
    stats%mean_value = dble(SUM(array)) / dble(SIZE(array,1))
    call calc_std_dev_array(array, stats%std_dev)

    ! Sort the array, if it isn't already sorted.
    array_is_sorted = .false.
    if (present(is_array_sorted)) array_is_sorted = is_array_sorted

    allocate( sorted_array(SIZE(array,1)) )
    if (.not.array_is_sorted) then
      call bubblesort_array(array, sorted_array)
    else
      sorted_array = array
    end if

    ! Now calculate statistics that require a sorted version of the array
    call calc_median_array(sorted_array, stats%median_value)
    call calc_mode_array(sorted_array, stats%mode_value, stats%mode_count)

    ! Calculate fractiles if n_fractiles is present
    if (present(n_fractiles)) then
      stats%n_fractiles = n_fractiles
      allocate( stats%fractiles(n_fractiles + 1) )
      call calc_fractiles_array(sorted_array, stats%n_fractiles, &
           stats%fractiles)
    end if

    ! Calculate the histogram if n_bins is present
    if (present(n_bins)) then
      stats%n_bins = n_bins
      ! If user specified both starting and ending values, use those
      if (present(hist_start_in) .and. present(hist_end_in)) then
        hist_start = hist_start_in
        hist_end = hist_end_in
      ! If they specified one or the other, exit
      else if ((present(hist_start_in) .and. .not. present(hist_end_in)) .or. &
               (.not. present(hist_start_in) .and. present(hist_end_in))) then
         call aims_stop("When specifying the range for a histogram, both the &
                        &minimum and maximum boundaries must be specified.", &
                        func)
      ! If they specified neither, use the edges of the distribution
      else
        hist_start = stats%min_value
        hist_end = stats%max_value
      end if

      ! Edge case; for zero interval size, we override the interval and use a
      ! small interval around the original values instead
      if (hist_start == hist_end) then
        hist_start = hist_start - 0.5
        hist_end   = hist_end + 0.5
      end if

      allocate( stats%histogram(n_bins) )
      allocate( stats%histogram_bounds(n_bins + 1) )
      call calc_histogram_array(sorted_array, stats%n_bins, hist_start, &
              hist_end, stats%histogram, stats%histogram_bounds)
    end if

    stats%is_valid = .true.

  end subroutine calc_stats_1D_int_array
!******

  ! Bubblesort subroutines

  subroutine bubblesort_1D_int_array(array, sorted_array)
    implicit none
    integer, intent(in)  :: array(:)
    integer, intent(out) :: sorted_array(:)

    integer :: i, j, array_size, temp

    character(*), parameter :: func = "bubblesort_1D_int_array"

    array_size = size(array,1)
    sorted_array = array
    do i = 1, array_size, 1
      do j = i + 1, array_size, 1
        if (sorted_array(i) > sorted_array(j)) then
          temp = sorted_array(i)
          sorted_array(i) = sorted_array(j)
          sorted_array(j) = temp
        end if
      end do
    end do
  end subroutine bubblesort_1D_int_array

  ! Standard deviation subroutines

  subroutine calc_std_dev_1D_int_array(array, std_dev)
    implicit none
    integer, intent(in)  :: array(:)
    real*8,  intent(out) :: std_dev

    real*8 :: mean_value, variance
    integer :: i, array_size

    character(*), parameter :: func = "calc_std_dev_1D_int_array"

    array_size = size(array,1)
    mean_value = dble(sum(array))/dble(array_size)
    variance = 0.0d0
    if (array_size > 1) then
      do i = 1, array_size, 1
        variance = variance + (dble(array(i)) - mean_value)**2
      end do
      variance = variance / dble(array_size-1)
      std_dev = sqrt(variance)
    else
      std_dev = 0.0d0
    end if
  end subroutine calc_std_dev_1D_int_array

  ! Mode subroutines

  subroutine calc_mode_1D_int_array(array, mode_value, mode_count)
    use mpi_tasks, only: aims_stop
    implicit none
    integer, intent(in)  :: array(:)
    integer, intent(out) :: mode_value
    integer, intent(out) :: mode_count

    integer :: i, array_size, this_value, this_count

    character(*), parameter :: func = "calc_mode_1D_int_array"

    array_size = size(array,1)
    mode_value = array(1)
    this_value = array(1)
    mode_count = 1
    this_count = 1
    do i = 2, array_size
      ! Check that the array is sorted
      if (array(i-1) .gt. array(i)) &
           call aims_stop("Array is not sorted.", func)

      if (this_value == array(i)) then
        ! We haven't reached a new value yet, increment
        this_count = this_count + 1
        if (this_count > mode_count) then
          ! The current value is the new mode.
          ! If there's more than one mode, the smallest will be returned.
          mode_value = array(i)
          mode_count = this_count
        end if
      else
        ! This is a new value, reset the count
        this_value = array(i)
        this_count = 1
      end if
    end do
  end subroutine calc_mode_1D_int_array

  ! Median subroutines

  subroutine calc_median_1D_int_array(array, median)
    use mpi_tasks, only: aims_stop
    implicit none
    integer, intent(in)  :: array(:)
    real*8,  intent(out) :: median

    integer :: i, array_size

    character(*), parameter :: func = "calc_median_1D_int_array"

    array_size = size(array,1)
    ! Check that the array is sorted
    do i = 2, array_size
      if (array(i-1) .gt. array(i)) &
           call aims_stop("Array is not sorted.", func)
    end do

    if (MOD(array_size,2) /= 0) then
      i = array_size/2 + 1
      median = dble(array(i))
    else
      i = array_size/2
      median = (dble(array(i)) + dble(array(i+1)))/2.0d0
    end if
  end subroutine calc_median_1D_int_array

  subroutine calc_fractiles_1D_int_array(array, n_fractiles, fractiles)
    use mpi_tasks, only: aims_stop
    implicit none
    integer, intent(in)  :: array(:)
    integer, intent(in)  :: n_fractiles
    real*8,  intent(out) :: fractiles(:)

    real*8 :: remain
    integer :: i, quot, i_fractile, array_size

    character(*), parameter :: func = "calc_fractiles_1D_int_array"

    if (n_fractiles < 1) &
         call aims_stop("Please specify a positive number for n_fractiles.", &
                        func)

    if (size(fractiles,1) < n_fractiles+1) &
         call aims_stop("fractiles array must have at least n_fractiles+1 many &
                        &elements.", func)

    array_size = size(array,1)
    ! Check that the array is sorted
    do i = 2, array_size
      if (array(i-1) .gt. array(i)) &
           call aims_stop("Array is not sorted.", func)
    end do

    fractiles = 0.0d0
    fractiles(1) = dble(array(1))

    do i_fractile = 1, n_fractiles - 1
      if (array_size == 1) then
        ! Essentially meaningless, but mathematically correct...
        fractiles(i_fractile+1) = dble(array(1))
      else
        ! We here use the C=1 variant of linear interpolation between adjacent
        ! ranks (c.f. the Wikipedia article on percentiles)
        quot = ((array_size-1)*i_fractile)/(n_fractiles)
        remain = dble((array_size-1)*i_fractile)/dble(n_fractiles) - dble(quot)
        fractiles(i_fractile+1) = dble(array(quot+1)) + &
                                  remain*dble( array(quot+2) - array(quot+1) )
      end if
    end do

    fractiles(n_fractiles+1) = dble(array(array_size))
  end subroutine calc_fractiles_1D_int_array

  subroutine calc_histogram_1D_int_array(array, n_bins, hist_start, &
                                         hist_end, histogram, &
                                         histogram_bounds)
    use mpi_tasks, only: aims_stop
    implicit none
    integer, intent(in)  :: array(:)
    integer, intent(in)  :: n_bins
    real*8,  intent(in)  :: hist_start
    real*8,  intent(in)  :: hist_end
    integer, intent(out) :: histogram(:)
    real*8,  intent(out) :: histogram_bounds(:)

    real*8 :: bin_size, bin_start, bin_end
    integer :: i_bin, i, array_size

    character(*), parameter :: func = "calc_histogram_1D_int_array"

    if (n_bins < 1) &
         call aims_stop("Please specify a positive number for n_bins.", &
                        func)

    if (size(histogram,1) < n_bins) &
         call aims_stop("histogram array must have at least n_bins many &
                        &elements.", func)

    if (size(histogram_bounds,1) < n_bins+1) &
         call aims_stop("histogram_bounds array must have at least n_bins + 1 &
                        &many elements.", func)

    if (hist_start >= hist_end) &
         call aims_stop("The starting value for the histogram must be smaller &
                        &than the final value.", func)

    array_size = size(array,1)
    ! Check that the array is sorted
    do i = 2, array_size
      if (array(i-1) .gt. array(i)) &
           call aims_stop("Array is not sorted.", func)
    end do

    histogram = 0
    bin_size = (hist_end - hist_start) / dble(n_bins)

    i = 1 ! Start with first element of array
    do i_bin = 1, n_bins
      ! Define the current bin
      ! Need to recalculate the values each time to avoid rounding errors
      bin_start = hist_start + dble(i_bin-1) * bin_size
      bin_end = hist_start + dble(i_bin) * bin_size

      histogram_bounds(i_bin) = bin_start

      ! Now loop over elements of the array for the current bin
      do while (.true.)
        if (i > array_size) exit ! Reach end of array, all bins are 0 afterwards

        ! We use the same criterion for setting up bins as histogram() in the
        ! numpy package, namely that the all bins are half-open, with the
        ! leading edge belonging to the bin and the trailing edge not belonging
        ! to the bin, except for the final bin where both edges reside in the
        ! bin
        if (bin_start <= dble(array(i)) .and. dble(array(i)) < bin_end) then
          ! Element belongs in current bin, add to bin and move to next element
          histogram(i_bin) = histogram(i_bin) + 1
          i = i + 1
        else if ((i_bin == n_bins) .and. (array(i) == hist_end)) then
          ! Element falls on the trailing edge of the final bin (which is the
          ! end of the histogram), add this point to the final bin, and move to
          ! next element
          histogram(i_bin) = histogram(i_bin) + 1
          i = i + 1
        else if (dble(array(i)) < bin_start) then
          ! We're below the bounds of the bin, move on to the next element
          ! This can happen when we start with a lower bound for the
          ! histogram larger than the lowest values in the sorted data set
          i = i + 1
        else
          ! Element belongs to a future bin, move on to next bin
          exit
        end if
      end do
    end do

    histogram_bounds(1) = hist_start
    histogram_bounds(n_bins+1) = hist_end
  end subroutine calc_histogram_1D_int_array

end module
!******
