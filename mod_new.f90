
!>Module non utilisé
module mod_new
	use mod_var
	
	implicit none
		

	contains


	function f(x)
	implicit none
	real(rp), intent(in) :: x
	real(rp)             :: f

	f = (x-1)*(x+1)*(x-3)	

	end function f

	function fprim(x)
	implicit none
	real(rp),intent(in) :: x
	real(rp)            :: fprim

	fprim = x*(3*x-6) -1

	end function fprim

	subroutine m_newton(x_0)
	implicit none
	real(rp),intent(inout) :: x_0
	real(rp)               :: x 
	real(rp)               :: tol
	real(rp)               :: err
	integer                :: iter_max , i

! initialisation
	err = one
	i   = 0
	
	print*,"tolérance = "," iter_max = "
	read*,tol,iter_max	

	do while (err > tol .and. i<iter_max )
		x = x_0 - f(x_0)/fprim(x_0)
		err = abs(x-x_0)
		x_0 = x
		i = i+1
	end do
	
	end subroutine m_newton


	end module mod_new
