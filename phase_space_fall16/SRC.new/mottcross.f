      SUBROUTINE MOTTCROSS(E,EP,THETA,SIG)

      implicit none
        

********************************************************************************
*   Calculates and returns Mott cross section.                                 *
********************************************************************************  


      real*4 E,EP,THETA,SIG
      real*8 TH,WW,QQ,x,nu,q4
      real*8 MP,MP2,mott,a1,a2
      real*8 sin2,cos2,tan2,alpha,alpha2
      integer tar
      
      character*1 tarf2

      if(tar.EQ.1) tarf2='h' 
      if(tar.EQ.2) tarf2='d'
      th = theta  
      mp = .9382727
      mp2 = mp*mp
      alpha = 1/137.03599
      alpha2 = alpha*alpha
      sin2 = sin(th/2.)*sin(th/2.)
      cos2 = cos(th/2.)*cos(th/2.)
      tan2 = tan(th/2.)*tan(th/2.)        
      qq = 4.*e*ep*sin2 
      q4 = qq*qq  
      nu = e-ep
      ww = mp2 + 2.*mp*nu -qq
      x = qq/2./mp/nu
      mott = 4.*alpha2*ep*ep*cos2/q4*0.389e6  !! In nb/GeV/sr   
      a1 = .3*nu/qq
      a2 = .3/nu    

      mott = mott*(a2 + 2.* a1*tan2)
        
      sig = mott

c      write(6,*) "mott: ", theta,mott    
    
      end
