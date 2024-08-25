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
	real :: dx, xi, xf, quad, x, y

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
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, x, y)
!$OMP MASTER
	write(*, "(' Using ',i0,' threads...')") OMP_get_num_threads()
!$OMP END MASTER
!$OMP DO
	do i = 2, num_pts-1
		x = xi + dx*(i-1)
		y = 0.0
		do j = 1,order+1
			y = y + coefficients(j) * x**(j-1)
		end do
!$OMP CRITICAL
		quad = quad + y
!$OMP END CRITICAL
	end do
!$OMP END DO
!$OMP END PARALLEL
	! scale
	quad = dx * quad

end function parallel_quadrature

end module parflib

program main

use parflib

integer :: order, num_pts_in, num_cpu_in, i
real :: quad
real, dimension(:), allocatable :: coefficients
real, dimension(2) :: bounds
integer :: num_pts, num_cpu

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

end program main