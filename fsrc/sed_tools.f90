module sed_tools
  ! Fortran tools for working with Spectral Energy Distributions (SEDs).
  !
  ! This code is in parts inspired by the code distributed alongside GALAXEV
  ! (Bruzual & Charlot, Mon.Not.R.Astron.Soc.344,1000–1028, 2003).
  !
  ! Martin Glatzle, mglatzle@mpa-garching.mpg.de, 2019
  use, intrinsic :: iso_fortran_env, only: REAL64
  implicit none

contains

  subroutine bin_sed(n,m,x,y,bin_edges,z)
    ! Bin an SED to given frequency/wavelength values conserving luminosity.
    !
    ! x is nu or lambda at which SED is known, assumed to be sorted in
    ! ascending order.
    !
    ! y(x) is L_nu or L_lambda.
    !
    ! bin_edges are assumed to be sorted in ascending order with
    ! diff(bin_edges) > 0, bin_edges(1) >= x(1) and bin_edges(n) <= x(n).
    !
    ! The new L_nu or L_lambda values at each bin are stored in z. These values
    ! are averaged over the bin so as to conserve luminosity. Averages are
    ! computed using trapezoidal integration.
    integer     , intent(in) :: n, m
    real(REAL64), intent(in) :: x(n), y(n), bin_edges(m+1)

    real(REAL64), intent(inout) :: z(m)

    integer :: istart, j

    istart = 1

    ! Build binned SED.
    do j=1, m
       z(j) = trapzq_sorted( &
            n,x,y,bin_edges(j),bin_edges(j+1),istart &
            )/(bin_edges(j+1)-bin_edges(j))
    enddo

  contains

    function trapzq_sorted(n,x,y,x1,x2,istart)
      ! Finds area below y(x) from x=x1 to x=x2 using trapezoid rule. y is
      ! interpolated at x1 and x2 so these are not necessarily values of x.
      !
      ! x is assumed to be sorted in ascending order.
      !
      ! x1 and x2 must fulfill >=x(1) and <=x(n). Moreover, x1<=x2 must hold.
      !
      ! Starts search for x1 at istart if provided and if x1>=x(istart), else
      ! at one. Sets istart to the last index at which x<=x2.
      integer     , intent(in) :: n
      real(REAL64), intent(in) :: x(n), y(n), x1, x2

      integer, intent(inout), optional :: istart

      real(REAL64) :: trapzq_sorted, y1, y2
      integer      :: i, i1, i2

      if(present(istart)) then
         if(x1 < x(istart)) then
            i1 = 1
         else
            i1 = istart
         end if
      else
         i1 = 1
      endif

      ! Find last element <x1 in x.
      do while(x1 > x(i1+1))
         i1 = i1+1
      end do
      i2 = i1
      ! Find last element <x2 in x.
      do while(x2 > x(i2+1))
         i2 = i2+1
      end do

      ! Find y1=y(x1) by interpolation.
      y1 = ( (x(i1+1)-x1)*y(i1) + (x1-x(i1))*y(i1+1) ) &
           / (x(i1+1)-x(i1))

      ! Find y2=y(x2) by interpolation.
      y2 = ( (x(i2+1)-x2)*y(i2) + (x2-x(i2))*y(i2+1) ) &
           / (x(i2+1)-x(i2))

      if(i1 == i2) then
         ! Have only area from x1 to x2.
         trapzq_sorted = (x2-x1)*(y1+y2)/2
      else
         ! Add area from x1 to x(i1+1).
         trapzq_sorted = (x(i1+1)-x1)*(y1+y(i1+1))/2
         ! Add area from x(i1+1) to x(i2-1).
         do i=i1+1, i2-1
            trapzq_sorted = trapzq_sorted+(x(i+1)-x(i))*(y(i)+y(i+1))/2
         end do
         ! Add area from x(i2) to x2.
         trapzq_sorted = trapzq_sorted+(x2-x(i2))*(y2+y(i2))/2
      end if

      if(present(istart)) then
         istart = i2
      end if
    end function trapzq_sorted
  end subroutine bin_sed
end module sed_tools
