PROGRAM WAFFLE
real:: g
real,allocatable:: matrix(:,:),eig(:),work(:)
integer:: i,j,k,l,m,p,info
integer:: nP, nS, nC,deg
integer,allocatable::states(:,:)
character:: unb

open(unit=2,file="eigen.dat")
write(*,*) 'Write the number of particles'
read(*,*) nP
write(*,*) 'Write the number of states'
read(*,*) nS
if(nS.lt.nP) then
	write(*,*) 'Learn to count'
	return
endif
write(*,*) 'Write the degeneracy of the orbitals:'
read(*,*) deg
if(nS.lt.deg) then
	write(*,*) 'Very funny'
	return
endif
write(*,*) 'Unbroken pairs? (y/n)'
read(*,*) unb


if (unb.eq.'y') then
	nC=INT(factorial(nS/2)/(factorial(nP/2)*factorial((nS-nP)/2)))
else if (unb.eq.'n') then
	nC=INT(factorial(nS)/(factorial(nP)*factorial(nS-nP)))
else
write(*,*) "Dumb idiot"
	return
endif


allocate(states(nC,nP))
k=1
do i=1,nP
	states(1,i)=k
	k=k+1
enddo


if(unb.eq.'n') then
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
else
do i=2,nC
	states(i,:)=states(i-1,:)
	!write(*,*) states(i,:)
	do j=1,nP-1,2
		states(i,nP-j)=states(i,nP-j)+2
		states(i,nP-j+1)=states(i,nP-j+1)+2
		if (states(i,nP-j).gt.nS-j) then
			states(i,nP-j)=0
			states(i,nP-j+1)=0
		else
        		exit
		endif
	enddo
	do j=1,nP
		if(states(i,j).eq.0) states(i,j)=states(i,j-1)+1
	enddo
	!write(*,*) states(i,:)
enddo
endif





allocate (matrix(1:nC,1:nC),eig(1:nC),work(1:3*nC))

do p=1,200
	g=(p-100)/200.d0
	matrix=0.d0

	do j=1,nC
	!interaction  
		do k=1,nP
	        	if(mod(states(j,k),deg).eq.1.and.states(j,k+1).eq.states(j,k)+1) then
	        	do i=1,nC
			do m=1,nP
	        	if(mod(states(i,k),deg).eq.1.and.states(i,m+1).eq.states(i,m)+1.and.states(i,m).eq.states(j,k)) matrix(i,j)=&
	        	                                                                                           matrix(i,j)-0.5*g
	        	enddo
	        	enddo
	        	endif
		enddo
	       ! write(*,*) states(j,:)
	!h_0
		do i=1,nC
		        if(i.eq.j) then
		        do k=1,nP
		        	matrix(i,j)=matrix(i,j)+ceiling(real(states(i,k))/real(deg))-1
			enddo
			endif
		enddo
	!write(*,*) matrix(1:nC,j)
	enddo

	work=0
	info=0

	call ssyev('v', 'u', nC , matrix, nC ,eig, work,3*nC, info )
	write(2,*) g, eig(1:nC)

enddo
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























