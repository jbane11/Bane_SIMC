        SUBROUTINE poz_tar_5(ydar,ygol,target,factor) 
*       program poz
	IMPLICIT NONE
        
        INCLUDE 'logicals.cmn'
        INCLUDE 'cs.cmn'
	 
	real*4 ydar,ygol,yg1,yg2,factor
	real*8 dif,gipomal,gipo,cosygla,diffact,difygol
	real*4 e1,e2,e3,e23,em,ema,emin,emax,fact1,fact2,e
	real*4 e16,e17,e18,e19,e20,e21,e22
        character csdatfile 
	integer k,i,j
        real*4 str
 
        if(firstcs) then
         firstcs = .false.
 
         do k=1,51
           yd_min_t(k)=0.
           yd_max_t(k)=0.
           fact10_t(k)=0. 
           fact14_t(k)=0. 
           fact18_t(k)=0. 
           fact22_t(k)=0. 
         enddo
 
********* READING MIN VALUES FOR MOMENTA        

         if(target.EQ.15) then
          csdatfiel='theta_h_poz.dat'
         elseif (target.EQ.11) then 
          csdatfile='theta_h_poz.dat'
         endif 
                   
         open(10,file=csdatfile)
         read(10,1001)str,e
         yd_min_t(1)=e
 1001    format(a1,F7.4)
	 
*	write(*,1002)yd_min(1)
 1002    format(F7.4)
 
          do i=2,50
	   read(10,1003)str,e1
	   yd_min_t(i)=e1
 1003      format(a1,F7.4)
*           write(*,1002)yd_min(i)
	  enddo 
	  close(10)
	 
********* READING MAX VALUES FOR MOMENTA        

          open(10,file=csdatfile) 	          
	  read(10,1004)str,e2
	  yd_max_t(1)=e2
 1004    format(a9,F7.4)
	 
*	write(*,1002)yd_max(1)
 
          do i=2,50
	   read(10,1005)str,e3
	   yd_max_t(i)=e3
 1005       format(a9,F7.4)
*           write(*,1002)yd_max(i)
	  enddo 
	  close(10)
          
********* READING Factors [N(-) - N(+)]/[N(-)]        

**** 10.6 degree

	  open(10,file=csdatfile) 	          
	  read(10,1019)str,e16
	  fact10_t(1)=e16
 1019   format(a16,F8.5)
*        WRITE(*,1007) fact10(1)  

	  do i=2,50
	   read(10,1020)str,e17 
           fact10_t(i)=e17
 1020       format(a16,F8.5)
*           WRITE(*,1007) fact10(i)	  
	   enddo
          close(10) 

 1007     format(F7.5)
       

**** 14.6 degree

	  open(10,file=csdatfile) 	          
	  read(10,1021)str,e18
	  fact14_t(1)=e18
 1021    format(a24,F8.5)
*        WRITE(*,1007) fact14(1)  

	  do i=2,50
	   read(10,1022)str,e19 
           fact14_t(i)=e19
 1022       format(a24,F8.5)
*           WRITE(*,1007) fact14(i)	  
	  enddo
         close(10) 
 
**** 18.6 degree

	  open(10,file=csdatfile) 	          
	  read(10,1023)str,e20
	  fact18_t(1)=e20
 1023    format(a32,F8.5)
*        WRITE(*,1007) fact18(1)  

	  do i=2,50
	   read(10,1024)str,e21 
           fact18_t(i)=e21
 1024       format(a32,F8.5)
*           WRITE(*,1007) fact18(i)	  
	  enddo
         close(10) 
         
**** 22.6 degree

	 open(10,file=csdatfile) 	          
	 read(10,1025)str,e22
	 fact22_t(1)=e22
 1025    format(a40,F8.5)
*        WRITE(*,1007) fact22(1)  

	 do i=2,50
	  read(10,1026)str,e23 
          fact22_t(i)=e23
 1026       format(a40,F8.5)
*           WRITE(*,1007) fact22(i)	  
          enddo
         close(10) 

        endif          !!!!!    READ from file    !!!!!

       
*###############################################################*
*****************************************************************
*###############################################################*

*** OPREDELENIE BINA PO P

           em   = 0.
	   emin = 0. 
	   ema  = 0.
	   emax = 0.
	   difygol=0.
	   diffact=0.
	   dif=0.
	   gipo=0.
	   gipomal=0.
	   cosygla=0.
	   factor=0.
	    
           DO j=1,50
            if(ydar.ge.yd_min_h(j)) then
             emin=yd_min_h(j)
	     emax=yd_max_h(j)
             if(emin.gt.em) then
              em=emin
	      ema=emax
	      
	       if(ygol.lt.14.6) then
	        fact1 = fact10_h(j)
		fact2 = fact14_h(j)
		yg1=10.6
		yg2=14.6
	       endif

	       if(ygol.ge.14.6.and.ygol.le.18.6) then
	        fact1 = fact14_h(j)
		fact2 = fact18_h(j)
		yg1=14.6
		yg2=18.6
	       endif
	       
	       if(ygol.gt.18.6) then
	        fact1 = fact18_h(j)
		fact2 = fact22_h(j)
		yg1=18.6
		yg2=22.6
	       endif
              
	     endif
            endif
           enddo
 
*            write(*,1111) em,ydar,ema,fact1,fact2
 1111        format(F7.4,1x,F7.4,1x,F7.4,2x,F8.5,1x,F8.5)
	 
	
           difygol=yg2-ygol
           diffact=abs(fact2-fact1)         
	   dif=yg2-yg1	   
           gipo=sqrt(dif**2+diffact**2)
           cosygla=dif/gipo
           gipomal=difygol/cosygla	   
           factor=sqrt(gipomal**2-difygol**2)
 
          if(ygol.le.22.60) then
           if(fact1.gt.fact2) then
            factor=factor+fact2
            else
            factor=factor+fact1
           endif
	  endif
	  
	  if(ygol.gt.22.60) then
	   if(fact1.lt.fact2) then
	    factor=fact1-factor
	     else
	    factor=fact2-factor
	   endif
	  endif   
	   
    	
    	    if(ydar.ge.2.0) then
    	     factor=1.
    	    endif 
    	   	
*	write(8,1112) ydar,factor
 1112    format(F5.3,2x,F8.5)	 
	 
	 
	end	
		










