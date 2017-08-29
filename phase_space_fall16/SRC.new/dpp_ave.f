        Program dpp_ave

CCCCC   Reads in cross section files and averages data/mc rations for each dp/p bin  CCCCC
	implicit none


	real*4 rat(100,16),erat(100,16),temp(3),dpp(16),ratav(16)
        real*4 eratav(16),eave(16)
        integer i,j,num,runs(100),endof
        character*80 infile

        open(unit=35,file='dp_runlist.dat',status='old') 
        i = 1
        endof = 1
        dowhile(endof.NE.-1)
          read(35,*,END=1001) runs(i)
          write(6,*) i, runs(i)
          i=i+1
        enddo
        

 1001   endof = -1
	
        close(35)

        num = i-1

c        write(6,*) "num = ",num

        do i=1,16
          ratav(i) = 0.0
          eave(i) = 0.0
        enddo

        do i=1,num
          write(infile,'("rat",i5,".dat")')runs(i)
          open(unit=45,file=infile,status='old')
          do j=1,16
            read(45,*) temp
            dpp(j) = temp(1)
            rat(i,j) = temp(2)
            erat(i,j) = temp(3)
            dpp(i) = temp(1)
            write(6,*) i,j,dpp(i),rat(i,j),erat(i,j)
            ratav(j) = ratav(j) + rat(i,j) 
            eave(j) = eave(j) + 1./erat(i,j)/erat(i,j)
          enddo
          close(45)
        enddo

        open(25)

        do i=1,16
          ratav(i) = ratav(i)/num
          eave(i) = 1./sqrt(eave(i))
          dpp(i) = -8.5+float(i)
          write(25,*) dpp(i),ratav(i),eave(i)
        enddo
	
        close(25)

c 1001   endof = -1
	
	end




