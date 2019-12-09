module sed_tools
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
    ! bin_edges are the edges of the bins. Assumed to be sorted in ascending
    ! order with bin_edges(1) >= x(1) and bin_edges(n) <= x(n).
    !
    ! The new L_nu or L_lambda values at each bin are stored in z. These values
    ! are averaged over the bin so as to conserve luminosity. Averages are
    ! computed using trapezoidal integration.
    integer     , intent(in) :: n, m
    real(REAL64), intent(in) :: x(n), y(n), bin_edges(m+1)

    real(REAL64), intent(inout) :: z(m)

    integer :: istart, j

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
      ! at one.
      integer     , intent(in) :: n
      real(REAL64), intent(in) :: x(n), y(n), x1, x2

      integer, intent(inout), optional :: istart

      real(REAL64) :: trapzq_sorted, y1, y2
      integer      :: i, i1, i2

      if(present(istart)) then
         if(x1.lt.x(i1)) then
            i1 = 1
         else
            i1 = istart
         end if
      else
         i1 = 1
      endif

      trapzq_sorted = 0.0_REAL64

      ! Find first element >=x1 in x.
      do i=i1, n
         if(x(i).ge.x1) then
            i1 = i
            if(i1==1) then
               y1 = y(1)
            else
               ! Find y1=y(x1) by interpolation.
               y1 = ( (x(i1)-x1)*y(i1-1) + (x1-x(i1-1))*y(i1) ) &
                    / (x(i1)-x(i1-1))
            endif
            ! Add area from x1 to x(i1).
            trapzq_sorted = trapzq_sorted+(x(i1)-x1)*(y1+y(i1))/2
            exit
         endif
      enddo

      ! Find last element <=x2 in x.
      do i=i1, n
         if(x(i).ge.x2) then
            if(i==n) then
               i2 = n
               y2 = y(n)
            else
               i2 = i-1
               ! Find y2=y(x2) by interpolation.
               y2 = ( (x(i2+1)-x2)*y(i2) + (x2-x(i2))*y(i2+1) ) &
                    / (x(i2+1)-x(i2))
            endif
            if(i1.eq.i2+1) then
               ! Just have area from x1 to x2.
               trapzq_sorted = (x2-x1)*(y2+y1)/2
               return
            else
               ! Add area from x(i2) to x2.
               trapzq_sorted = trapzq_sorted+(x2-x(i2))*(y2+y(i2))/2
            endif
            exit
         endif
      enddo

      ! Compute area from x(i1) to x(i2)
      do i=i1+1, i2
         trapzq_sorted = trapzq_sorted+(x(i)-x(i-1))*(y(i)+y(i-1))/2
      enddo
    end function trapzq_sorted
  end subroutine bin_sed
end module sed_tools
