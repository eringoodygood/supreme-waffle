PROGRAM WAFFLE
real:: g
real,allocatable:: matrix(:,:),eig(:),work(:)
integer:: i,j,k,l,m,p,info
integer:: nP, nS, nC,deg
integer,allocatable::states(:,:),lastc(:),vvec(:)
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


allocate(states(nC,nS),lastc(nP),vvec(nS))
states=0

k=1
do i=1,nP
	states(1,i)=1
	lastc(i)=k
	k=k+1
enddo


do i=2,nC
	if(unb.eq.'n') then
		do j=0,nP-1
			lastc(nP-j)=lastc(nP-j)+1
			if (lastc(nP-j).gt.nS-j) then
				lastc(nP-j)=0
			else
				exit
			endif
		enddo
		do j=1,nP
			if(lastc(j).eq.0) lastc(j)=lastc(j-1)+1
		enddo
	else
		do j=1,nP-1,2
			lastc(nP-j)=lastc(nP-j)+2
			lastc(nP-j+1)=lastc(nP-j+1)+2
			if (lastc(nP-j).gt.nS-j) then
				lastc(nP-j)=0
				lastc(nP-j+1)=0
			else
				exit
			endif
		enddo
		do j=1,nP
			if(lastc(j).eq.0) lastc(j)=lastc(j-1)+1
		enddo
	endif
	do j=1,nP
		k=lastc(j)
		states(i,k)=1
	enddo
enddo





allocate (matrix(1:nC,1:nC),eig(1:nC),work(1:3*nC))

do l=1,200
	g=(l/100.d0)-1.d0
	matrix=0.d0

	do j=1,nC
	!interaction  
		do k=1,nS-1
	        	if(mod(k+1,deg).ne.1.and.states(j,k).eq.1.and.states(j,k+1).eq.1) then
	        	do p=1,nS-1
	        	vvec(:)=states(j,:)
	       		vvec(k)=0
	       		vvec(k+1)=0
       			if(mod(p+1,deg).ne.1.and.vvec(p).eq.0.and.vvec(p+1).eq.0) then
				vvec(p)=1
				vvec(p+1)=1
				do i=1,nC
					if(all(states(i,:).eq.vvec)) matrix(i,j)=matrix(i,j)-0.5*g
				enddo
       			endif
        		enddo
	        	endif
		enddo
	       ! write(*,*) states(j,:)
	!h_0
		do i=1,nC
		        if(i.eq.j) then
		        do k=1,nS
		        	if(states(i,k).eq.1) matrix(i,j)=matrix(i,j)+ceiling(real(k)/real(deg))-1
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






















