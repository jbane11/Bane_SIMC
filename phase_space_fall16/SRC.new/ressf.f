CCC  Version 050530  -  Author:  M.E. Christy                        CCC
CCC  Subroutine to get F1 and FL from fits to L/T cross sections.    CCC
CCC  The subroutine resmod.f is required.                            CCC


      SUBROUTINE ressf(W2,Q2,F1,FL)

      IMPLICIT NONE

      INCLUDE 'logicals.cmn'
      INCLUDE 'resparms.cmn'

      real*8 w2,q2,fn,fnerr,w2inc,temp(4)
      real*8 mp,mp2,pi,alpha,xb,f1,fl
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.036

      if(first_ressf) then
        first_ressf = .false.
        open(unit=32,file='f1parms.dat',status='old')
        do i=1,50
          read(32,*) temp(1)  ! par #
          read(32,*) xval1(i) ! starting value
          read(32,*) temp(2)   ! initial step (0 means fixed parm)
          read(32,*) temp(3) ! low limit
          read(32,*) temp(4) ! high limit 
        enddo
        close(32)

        open(unit=33,file='flparms.dat',status='old')
        do i=1,50
          read(33,*) temp(1)  ! par #
          read(33,*) xvall(i) ! starting value
          read(33,*) temp(2)   ! initial step (0 means fixed parm)
          read(33,*) temp(3) ! low limit
          read(33,*) temp(4) ! high limit 
        enddo
        close(33)
      endif


      xb = q2/(w2+q2-mp2)
      call resmod(1,w2,q2,xval1,f1)
      call resmod(2,w2,q2,xvall,fl)

      f1 = f1*(w2-0.9382727*0.9382727)/8./pi/pi/alpha/0.3894e3
      fl = fl*2*xb*(w2-0.9382727*0.9382727)/8./pi/pi/alpha/0.3894e3

      if(w2.LT.1.159) then
       f1 = 0.
       fl = 0.
      endif
    
      end



