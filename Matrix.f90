PROGRAM WAFFLE

!Compile with gfortran Matrix.f90 -llapack

real:: g
real:: matrix(1:6,1:6),eig(1:6)
real:: work(100)
integer:: i,j,k,p,states(1:2,1:6),info

open(unit=2,file="eigen.dat")
k=1
do i=1,4
        do j=i+1,4
            states(1,k)=i
            states(2,k)=j
            k=k+1
        enddo
enddo


do p=1,200
g=(p/100.d0)-1.d0
matrix=0.d0

do i=1,6
        do j=1,6
                k=0
!interaction  
                if(states(1,i).eq.states(1,j)) k=k+1
                if(states(1,i).eq.states(2,j)) k=k+1
                if(states(2,i).eq.states(1,j)) k=k+1
                if(states(2,i).eq.states(2,j)) k=k+1
                matrix(i,j)=-k*g     
!h_0
                if(i.eq.j) matrix(i,j)=matrix(i,j)+2.d0*(states(1,i)-1)+2.d0*(states(2,i)-1)
        enddo
enddo

work=0
info=0

call ssyev('v', 'u', 6 , matrix, 6 ,eig, work,100, info )
write(2,*) g, eig(1:6)

enddo
close(2)

RETURN
END PROGRAM WAFFLE
