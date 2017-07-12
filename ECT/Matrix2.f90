program test_call_C
!gcc -c get_filled_ar.c
!gfortran get_filled_ar.o Matrix2.f90 -o test
use iso_c_binding
implicit none
interface c_interface
   subroutine get_filled_ar (ar) bind (C, name = "get_filled_ar")
   use iso_c_binding
   implicit none
   integer (c_int), intent (out), dimension (*) :: ar
   end subroutine get_filled_ar
end interface c_interface


integer (c_int), dimension (0:3) :: ar
call get_filled_ar (ar)
write (*, *) "Fortran: ar:", ar
stop
end program test_call_C
