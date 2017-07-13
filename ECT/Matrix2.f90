PROGRAM WAFFLE

use iso_c_binding
implicit none
interface c_interface
   subroutine get_filled_ar (nb,ar) bind (C, name = "get_filled_ar")
   use iso_c_binding
   implicit none
   real (c_float), intent (out), dimension (*) :: ar
   integer (c_int), value ::nb
   end subroutine get_filled_ar
   subroutine open_file(stringArray) bind (C, name = "open")
   use iso_c_binding
   implicit none
   CHARACTER(LEN=1), TARGET :: stringArray
   end subroutine open_file   
   subroutine close_file() bind (C, name = "close")
   use iso_c_binding
   implicit none
   end subroutine close_file
end interface c_interface


real:: g
real,allocatable:: matrix(:,:),eig(:),work(:),vmatrix(:,:,:,:),h_0(:)
integer:: i,j,k,l,k2,l2,m,p,info,p1,p2,p3,p4
integer:: nP, nS, nC,deg
integer,allocatable::states(:,:)
character:: unb
real (c_float), dimension (0:5) :: ar

open(unit=2,file="eigen.dat")

!stringPtrs(ns) = C_LOC(stringArray(ns))

call open_file("data.dat")
nS=12
nP=8
nC=INT(factorial(nS)/(factorial(nP)*factorial(nS-nP)))


allocate(states(nC,nP),vmatrix(nS,nS,nS,nS),h_0(nS))
k=1
do i=1,nP
	states(1,i)=k
	k=k+1
enddo


do i=2,nC
	states(i,:)=states(i-1,:)
	do j=0,nP-1
		states(i,nP-j)=states(i,nP-j)+1
		if (states(i,nP-j).gt.nS-j) then
			states(i,nP-j)=0
		else
        		exit
		endif
	enddo
	do j=1,nP
		if(states(i,j).eq.0) states(i,j)=states(i,j-1)+1
	enddo
enddo


allocate (matrix(1:nC,1:nC),eig(1:nC),work(1:3*nC))


vmatrix=0.d0
do p=1,343
	call get_filled_ar (1,ar)
	write (*, *) ar
	i=int(ar(0))
	j=int(ar(1))
	k=int(ar(2))
	l=int(ar(3))
	vmatrix(i,j,k,l)=ar(4)
	!write(*,*) i,j,k,l,vmatrix(i,j,k,l)
enddo
call close_file()

call open_file("data2.dat")
h_0(12)=0.d0
do p=1,12
	call get_filled_ar (2,ar)
	!write (*, *) ar
	h_0(p) = ar(1)
	write (*, *) h_0(p)	
enddo


return


!interaction
matrix=0.d0
do i=1,nC  
	do j=1,nC
		do k=1,nP
		do k2=k,nP
			do l=1,nP
			do l2=l,nP
				p1=states(i,k)
				p2=states(i,k2)
				p3=states(j,l)
				p4=states(j,l2)
				matrix(i,j)=matrix(i,j)+vmatrix(p1,p2,p3,p4)
			enddo
			enddo
		enddo
		enddo
	enddo
enddo

do i=1,nC  
	do k=1,nP
		p1=states(i,k)
		matrix(i,i)=matrix(i,i)+h_0(p1)
	enddo
enddo
work=0
info=0

call ssyev('v', 'u', nC , matrix, nC ,eig, work,3*nC, info )
write(2,*) g, eig(1:nC)
write(*,*) 'Yeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee', nP, nC, nS
write(*,*) g, eig(1:nC)

!enddo
close(2)
deallocate(states)
RETURN


CONTAINS

real FUNCTION factorial(n)

integer::n,i
factorial=1
do i=1,n
	factorial=factorial*i
enddo
END FUNCTION factorial

END PROGRAM WAFFLE
