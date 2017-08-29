        program optsub

	Implicit none

        real*8 dt1(16,20),dt2(16,20),dt(16,20),tempvect(20)
        integer i,j

        do j=1,20
         tempvect(j) = 0.d0
        enddo

        open(unit=21, file='dtopt1.dat',status='unknown') 
        do i=1,16
         read(21,*) tempvect
         do j=1,20
         dt1(i,j) = tempvect(j)
         enddo
c         write(6,*) tempvect
        enddo
        close(21)

        do j=1,20
         tempvect(j) = 0.d0
        enddo

        open(unit=22, file='dtopt2.dat',status='unknown') 
        do i=1,16
         read(22,*) tempvect
         do j=1,20
          dt2(i,j) = tempvect(j)
         enddo
c         write(6,*) tempvect
        enddo
        close(22)

        do j=1,20
         tempvect(j) = 0.d0
        enddo

        open(unit=23, file='optcor.dat',status='unknown') 

        do i=1,16
         do j=1,20
          dt(i,j) = dt1(i,j)-dt2(i,j)
          if(abs(dt(i,j)).GT.2.) dt(i,j) = 0.
c          write(6,*) dt1(i,j),dt2(i,j),dt(i,j)
          tempvect(j) = dt(i,j)
         enddo
         write(23,8000) tempvect
        enddo

 8000 format(20f7.3) 

        end



