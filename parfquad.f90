module parflib

implicit none

contains

function parallel_quadrature(order, coefficients, bounds, num_pts_in, num_cpu_in) result(quad)

	use omp_lib

	! Input variables
	integer, intent(in) :: order
	real, dimension(order+1), intent(in) :: coefficients ! e.g., c_0 + c_1 * x + c_2 * x^2 + ...
	real, dimension(2), intent(in) :: bounds ! assume evaluating at bounds(1) and bounds(2)
	integer, intent(in), optional :: num_pts_in, num_cpu_in

	! Internal variables
	integer :: num_pts, num_cpu, i, j
	integer(8) :: wall_clock1, wall_clock2, wc_rate
	real :: dx, xi, xf, quad, x, y

	call SYSTEM_CLOCK(count=wall_clock1, count_rate=wc_rate)

	! Sanity checks
	if(order < 0) then
		write(*, "(' Invalid order (< 0)')")
		stop
	endif

	if(bounds(2) <= bounds(1)) then
		write(*, "(' Invalid bounds, must be strictly ascending')")
		stop
	endif

	! Deal with optional arguments
	if(present(num_pts_in)) then
		num_pts = num_pts_in
	else 
		num_pts = 100
	end if

	if(num_pts < 2) then
		write(*, "(' Invalid number of evaluation points (< 2)')")
		stop
	endif

	! Provide user confirmation of polynomial and bounds
	write(*, "(' Integrating polynomial: ')", advance="no") 
	write(*, "(ES9.2)", advance="no") coefficients(0)
	if(order > 1) then
		do j = 2,order+1
			write(*, "(' + (',ES9.2,')x^',i0)", advance="no") coefficients(j), j-1
		end do
	end if
	write(*, "()") ! Finally advance
	write(*, "(' Between bounds: [',ES9.2,ES9.2,']')") bounds(1), bounds(2)

	if(present(num_cpu_in)) then
		num_cpu = num_cpu_in
	else 
		num_cpu = 2
	end if

	call OMP_set_num_threads(num_cpu)

	xf = bounds(2)
	xi = bounds(1)
	dx = (xf - xi) / (num_pts - 1) ! find delta_x

	! Equal spaced trapezoid rule: 
	! sum = dx / 2 * (y_0 + 2*y_1 + 2*y_2 + ... + y_n)
	quad = 0.0
	! Handle first and last points
	y = 0.0
	do j = 1,order+1
		y = y + coefficients(j) * xi**(j-1)
	end do
	quad = quad + 0.5 * y
	y = 0.0
	do j = 1,order+1
		y = y + coefficients(j) * xf**(j-1)
	end do
	quad = quad + 0.5 * y

	! Main loop, parallel
!$OMP PARALLEL
!$OMP MASTER
	write(*, "(' Using ',i0,' thread(s)...')") OMP_get_num_threads()
!$OMP END MASTER
! Use an OMP reduction to avoid using a CRITICAL section around the summation
!$OMP DO REDUCTION(+:quad)
	do i = 2, num_pts-1
		x = xi + dx*(i-1)
		y = 0.0
		do j = 1,order+1
			y = y + coefficients(j) * x**(j-1)
		end do
		quad = quad + y
	end do
!$OMP END DO
!$OMP END PARALLEL
	! scale
	quad = dx * quad

	call SYSTEM_CLOCK(count=wall_clock2)

	write(*, "(' Elapsed time (s): ',ES14.7)") real(wall_clock2 - wall_clock1) / real(wc_rate)

end function parallel_quadrature

function analytical(order, coefficients, bounds) result(ana)

	! Input variables
	integer, intent(in) :: order
	real, dimension(order+1), intent(in) :: coefficients ! e.g., c_0 + c_1 * x + c_2 * x^2 + ...
	real, dimension(2), intent(in) :: bounds ! assume evaluating at bounds(1) and bounds(2)

	! Internal variables
	real :: xi, xf, yi, yf, ana
	integer :: i

	! Analytical integration is sum( c_i / (i+1) * x^(i+1)) evaluated at end points
	xi = bounds(1)
	xf = bounds(2)
	yi = 0.0
	yf = 0.0
	do i = 1,order+1
		yi = yi + coefficients(i) * xi**i / real(i)
		yf = yf + coefficients(i) * xf**i / real(i)
	end do
	ana = yf - yi

end function analytical

end module parflib

program main

use parflib

integer :: order, num_pts_in, num_cpu_in, i
real :: quad, ana
real, dimension(:), allocatable :: coefficients
real, dimension(2) :: bounds
integer :: num_pts, num_cpu

write(*,"('  _____        _____  ______ ____  _    _         _____  ')")
write(*,"(' |  __ \ /\   |  __ \|  ____/ __ \| |  | |  /\   |  __ \ ')")
write(*,"(' | |__) /  \  | |__) | |__ | |  | | |  | | /  \  | |  | |')")
write(*,"(' |  ___/ /\ \ |  _  /|  __|| |  | | |  | |/ /\ \ | |  | |')")
write(*,"(' | |  / ____ \| | \ \| |   | |__| | |__| / ____ \| |__| |')")
write(*,"(' |_| /_/    \_\_|  \_\_|    \___\_\\____/_/    \_\_____/ ')")
write(*,"('')")
write(*,"(' Parallel Fortran Quadrature of Polynomials ')")

write(*,"(' Enter the polynomial order: ')")
read(*,*) order

allocate(coefficients(order+1))
write(*,"(' Enter the polynomial coefficients from zeroth to nth, press enter between:')")
do i = 1,order+1
	read(*,*) coefficients(i)
end do

write(*,"(' Enter the lower bound then upper bound, press enter between:')")
read(*,*) bounds(1)
read(*,*) bounds(2)

write(*,"(' Enter the number of equally-spaced integration points:')")
read(*,*) num_pts

write(*,"(' Enter the number of threads:')")
read(*,*) num_cpu

quad = parallel_quadrature(order, coefficients, bounds, num_pts, num_cpu)

write(*, "(' Quadrature result: ',ES14.7)") quad

ana = analytical(order, coefficients, bounds)

write(*, "(' Analytical result: ',ES14.7)") ana

end program main