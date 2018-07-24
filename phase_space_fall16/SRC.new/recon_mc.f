 	program recon_mc

	implicit none

        INCLUDE 'imodel.cmn'
        INCLUDE 'rad.cmn'
        INCLUDE 'logicals.cmn'
        INCLUDE 'flags.cmn'
       
	integer	nwpawc
	real	memor
        integer lrecl,nevt,ievt,ierr,istat
        real*4 ntuple_contents(17),xfoc,dydz,ypcor
        real*4 hse,hsp,hsec,hsev,eb,thetactemp,ebeam,radcon,thetac
        real*4 thetacrad,mp,mp2,nu,sin2,q2,w2,rcic,rcib,rcec
        real*4 rceb,trad,targetdata(6),eff_cal,eff_cer,emc
        real*4 pie,bcm1charge,bcm2charge,gbcm1charge,hpre,edt
        real*4 bcmavecharge,bmcur1,bmcur2,bmcur,eltime,cltime2
        real*4 positron_weight,posal,born,rci,rce,delup,deldown
        real*4 delini,ytarini,yptarini,xptarini,zini,thetaini
        real*4 delrec,ytarrec,yptarrec,xptarrec,zrec,thetarec
        real*4 w,normfac,charge,ydata,lumfract,sigmacent,sigmac
        real*4 lumdata,lummc,fract,sigave,sigtot,dxp,dyp,dep
        real*4 prescale,cltime,trackeff,trigeff,xpup,ypup,dtup,hstheta
        real*4 dt,phasespcor,denscor,delcor,phase_space
        real*4 hmsprlo,hmstof,hms34,rate,poscs,t1,t2,t3,t4
        real*8 fitemc,xb,AA

        real*4 ex,z,thick,emean,de

	character*80 infile,outfile,title,ishit,new
        character*80 directory,icycle
        integer nentries,target,tartemp,fail_id,runnum,maxev
	integer	i, m, NtupleSize, bank,ngen,ngentot,ntrecl,cs_flag,rc_flag
        logical firstr,ron,goodfit,newrc,cryo,docs,dorc,use_rcmod

        parameter (nwpawc=500000)
        parameter (nentries = 27)        
        parameter(bank = 1000)
        parameter(title = 'RECONTUPLE')
        character*80 NtupleTag(nentries)
        real*4 ntu(nentries)
        common /pawc/ memor(nwpawc)

        ntrecl = 4096

        first = .true.
        firstr = .true.
        firstcs = .true.
        firstcsdum = .true.
        firstres = .true.
        first_ressf = .true.
        ron = .true.
        newrc = .true.
        use_rcmod = .true.

c       use_rcmod = .false.

!!!!!       Options         !!!!!

        newrc = .true.
        usenmc = .true.
        cryo = .true.

        docs = .true.        !   if true then include charge-symmetric contributions  !
        dorc = .true.        !   if true then include radiative contributions         !


        ngen = 0        

        mp = .9382723
        mp2 = mp*mp

CCC c   read(5,*) runnum
        write(outfile,'("mcrun.rzdat")')
        open(unit=18, file='input.dat',status='old') 
        read(18,*) infile,ngentot,maxev,dxp,dyp,delup,deldown,cs_flag,
     &              rc_flag

        if(cs_flag.EQ.0) docs = .false.
        if(rc_flag.EQ.0) dorc = .false.
          

        radcon = 180./3.141592654
 
        sigave = 0.0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCC      Read run info from database for data           CCCCCCC

        open(unit=16,file='reconmc.in',status='old')

        read(16,*) target,ebeam,hsec,thetac,prescale,bmcur1,
     &            bcm1charge,cltime,eltime,trackeff,
     &            trigeff,rate

c        hmsprlo = 1.0   !!!!*****  TEST  *****!!!!

C        trigeff = hmsprlo*hmstof + (1. - hmsprlo)*hms34
          

        if(target.NE.11.OR.target.NE.15) cryo = .false.

        bcmavecharge = (bcm1charge+bcm1charge)/2.         !!! Average over BCMs !!!
        bmcur = (bmcur1+bmcur1)/2.
        close(16)
        

        write(6,*) 'Run#        = ', runnum
        write(6,*) 'Target#     = ', target
        write(6,*) 'Ebeam       = ', ebeam  
        write(6,*) 'Phms        = ', hsec 
        write(6,*) 'Theta       = ', thetac
        write(6,*) 'Prescale    = ', prescale 
        write(6,*) 'Comp. Ltime = ', cltime
        write(6,*) 'Elect. Ltime = ', eltime        
        write(6,*) 'Track Eff.  = ', trackeff
        write(6,*) 'Trig Eff.   = ', trigeff
        write(6,*) 'Charge BCM1 = ', bcm1charge
        write(6,*) 'Charge BCM2 = ', bcm2charge
        write(6,*) 'Charge      = ', bcmavecharge
        write(6,*) 'Current     = ', bmcur
        write(6,*) 'Rate (kHz)  = ', rate
        write(6,*)
        write(6,*) 'Including  Radiative Contributions? ', dorc
        write(6,*) 'Including Charge Symmetric Contributions?', docs
        write(6,*) 'Using Model from RC file?', use_rcmod
          

c        write(6,*) 'Charge      = ', charge
 
c        hsec = hsec*(1.0-.0022749)       !!!!!  Put in kinematic offsets  !!!!!        

        thetacrad = thetac/radcon
        dep = (delup-deldown)/100.*hsec

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCC        Read in target data from file                CCCCCCC

        open(unit=17, file='targetdata.dat',status='old') 
        do i=1,target
          read(17,*) targetdata
        enddo        

CCCCCC              Initialize radcor arrays               CCCCCCC
 
        if(newrc) then
         firstr = .true. 
         trad = thetac/radcon
CCCCC        write(6,*) 'debug::::  calling rc_mod ::'
         call rc_mod(firstr,trad,trad,hsec,1,rcic,t1)
CCCCC        write(6,*) 'debug:: value for rad table (1,1)::', exttab(1,1)
CCCCC        write(6,*) 'debug:: value for rad table (5,1)::', exttab(1,5)
        write(6,*) 'debug:: value for rad table (13,1)::', exttab(1,13)

         write(6,*) firstr,trad,hsec,rcic
         write(6,*) "test t1 value: " , t1


         firstr = .false.
        else
         firstr = .true.
         trad = thetac/radcon
         call rcint(firstr,trad,trad,hsec,target,rcic)
         trad = thetac/radcon
         firstr = .true. 
         call rcext(firstr,trad,trad,hsec,target,rcec)
         firstr = .false.
        endif
 
c        write(6,*) "rcic: ", rcic, "rcec: ",rcec

        rcic = 1.
        rcec = 1.
  
        lumdata = targetdata(6)*6.022137e-10/targetdata(4)
     &        *bcmavecharge/1.602177e-13
c     &       *targetdata(5)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c        ebeam = ebeam*(1.+.001)
        eb = ebeam

        call hlimit(nwpawc)
       
        lrecl = 0
        ishit = 'U'
        new = 'N'
        istat = 0

	call hropen (1, 'MCntuple', infile, ' ',ntrecl,istat) 
        call HCDIR(directory,'R')

        call hgnpar (1, 'readdat') 
        call hnoent (1, nevt)  
        call HCDIR(directory,'R') 
        write(6,*) '               Old ntuple is: ', infile
        call HLDIR(' ',' ')

                 
 	m = 0
        m = m+1
        NtupleTag(m) = 'xfoc'
        m = m+1
        NtupleTag(m) = 'yfoc'
        m = m+1 
        NtupleTag(m) = 'dx_fp'
        m = m+1
        NtupleTag(m) = 'dy_fp'
        m = m+1
        NtupleTag(m) = 'ytar_init'
        m = m+1
        NtupleTag(m) = 'dppi'
        m = m+1
        NtupleTag(m) = 'dthi'
        m = m+1
        NtupleTag(m) = 'dphi'
        m = m+1
        NtupleTag(m) = 'ytar_recon'
        m = m+1
        NtupleTag(m) = 'hsdelta'  
        m = m+1
        NtupleTag(m) = 'dydz_recon' 
        m = m+1
        NtupleTag(m) = 'dxdz_recon'        !dphr
        m =m+1
        NtupleTag(m) =  'fry'
        m = m+1
        NtupleTag(m) =  'frx'
        m = m+1
        NtupleTag(m) =  'good_evt'
        m = m+1
        NtupleTag(m) =  'stop_where'
        m = m+1
        NtupleTag(m) = 'x_stop'
        m = m+1 
        NtupleTag(m) = 'y_stop' 
        m = m+1          
        NtupleTag(m) = 'born'
        m = m+1
        NtupleTag(m) = 'rci'
        m = m+1
        NtupleTag(m) = 'rce' 
        m = m+1
        NtupleTag(m) = 'positron_weight'
        m = m+1
        NtupleTag(m) = 'hse'
        m = m+1
        NtupleTag(m) = 'hstheta'
        m = m+1
        NtupleTag(m) = 'sigmac'
        m = m+1
        NtupleTag(m) = 'q2'
        m = m+1
        NtupleTag(m) = 'w2'
               
        ntuplesize = m

        call hropen (2, 'reconmc', outfile,'N',ntrecl,istat) 
        call hrin (2, 99999,0)

        call HBOOKN(9040,'reconmc',NtupleSize,'reconmc',1000,NtupleTag)
        call hcdir('//reconmc',' ')
        write(6,*) '                    New ntuple is:  ', outfile
       
        write(6,*)
        write(6,*) 'Number of events in master file =',nevt  
        if(maxev.GT.nevt) maxev = nevt
        write(6,*)
        write(6,*) 'Number of events analyzing = ',maxev
        write(6,*)
 
       
        ebeam = eb
        do ievt = 1, maxev          
         if(ievt.EQ.10000) write(6,*) ' analyzed 10000 events'
         if(ievt.EQ.50000) write(6,*) ' analyzed 50000 events'
         if(mod(ievt,100000).EQ.0.)write(6,*) ' analyzed',ievt,' events'
         
         call hcdir('//MCntuple',' ')
         call hgnf(1, ievt, ntuple_contents, ierr)
         if (ierr .ne. 0) then
          write (6,*) 'hgnf err:', ierr
         endif 
 
         xfoc = ntuple_contents(1)  !??
         dydz = ntuple_contents(4)  !??
         delini = ntuple_contents(6)          
         yptarini = ntuple_contents(7)*1000.
         xptarini = ntuple_contents(8)*1000.  !convert to mrad
c         zrec = ntuple_contents(9)  !ask him later ??
         delrec = ntuple_contents(10)
         yptarrec = ntuple_contents(11)*1000. 
         xptarrec = ntuple_contents(12)*1000.
         fail_id = ntuple_contents(15)
         ytarini = ntuple_contents(5)
         ytarrec = ntuple_contents(9) 


         hse = hsec*(1.+delini/100.)
         if (targetdata(4).GT.1.1) then
           z = targetdata(4)/2.
         else
           z = 1.
         endif
         thick = targetdata(6)/2.
         ex = 21.8/1.0e6
         de = 0.0 
         emean = 0.0

c         call samp_eloss(1.0,1.0,0.511,thick,ex,ebeam*1000.,emean,de)
c         ebeam = eb + (emean - de)/1000.
c         call samp_eloss(1.0,1.0,0.511,thick,ex,hse*1000.,emean,de)
c         hsev hse + de/1000.         
c         hse = hse + emean/1000.  
          hsev = hse

c         write(6,*) hse,emean,de

         hsp = hse

c         yptarrec = yptarrec - 0.90    !!!!!!   Test
c         xptarrec = xptarrec - 0.69 
c         yptarini = yptarini - 0.90    !!!!!!   Test
c         xptarini = xptarini - 0.69

         call yp_optcor(xfoc,dydz*1000.,ypcor)

         ypcor = 0.0

c         write(6,*) xfoc,dydz,yptarrec,ypcor

         yptarrec = yptarrec-ypcor 

         thetaini = acos(cos(thetacrad+yptarini/1000.)   !initial theta
     &              *cos(xptarini/1000.))             
      
         hstheta = acos(cos(thetacrad+yptarrec/1000.)    !reconstruction theta
     &              *cos(xptarrec/1000.))
c          zrec  = ytarrec/cos(thetacrad+yptarrec/1000.)

CCCCCCC    Calculate the vertex kinematics for the event    CCCCCCCCC

 
         sin2 = sin(thetaini/2.)*sin(thetaini/2.)
         nu = ebeam - hsev
         q2 = 4.*hsev*ebeam*sin2
         w2 = mp2 + 2.*mp*nu-q2
      
ccccc     I need to set condition to get out the nonsense events
        if (w2.LT.0) then 
          born =0
          rci =1

          goto 777
        endif

         w = sqrt(w2) 
         xb = q2/2./mp/nu  
            


CCCCCCC           Get Model Cross section in nb/SR/GeV           CCCCCCCCC

         AA = targetdata(4)
	 

         if(use_rcmod) then
           call rc_mod(firstr,thetaini,thetacrad,hsev,1,rci,born)    !every event has theta and e' = hsev
           call rc_mod(firstr,thetacrad,thetacrad,hsev,1,t1,sigmac) !every event with theta center and e'
         else
           call model_new(ebeam,hsev,thetaini,target,born)      !!!   For actual angle     !!!
           call model_new(ebeam,hsev,thetacrad,target,sigmac)   !!!   For central angle    !!!

	endif
         dt = thetaini - trad
         phasespcor = 1./cos(dt)/cos(dt)/cos(dt)

!!!!!!!     Changed on 6/19/01  !!!!!!!!
c          phasespcor = 1.

c          write(6,*) dt,phasespcor
 
c         write(6,*) ebeam,hse,thetaini,target,born

c         write(6,*) xb,AA,emc

         if(abs(xptarini).LT.dxp.AND.abs(yptarini).LT.dyp.AND.
     &      delini.GT.deldown.AND.delini.LT.delup.AND.born.GE.0.) then
          sigave = sigave + born  
          ngen = ngen + 1
         endif

!!!!!!!   Apply delcor :  7/09/02
       
          delcor =  1.009+.4260E-02*delrec-.8603E-03*delrec*
     &       delrec-.10942E-03*delrec*delrec*delrec+.12697E-04*
     &       delrec*delrec*delrec*delrec+.1094E-07*delrec*
     &       delrec*delrec*delrec*delrec

c          delcor =  1.0077+.4707E-02*delrec-.84067E-03*delrec*
c     &       delrec-.16119E-03*delrec*delrec*delrec+.13636E-04*
c     &       delrec*delrec*delrec*delrec+.64174E-06*delrec*
c     &       delrec*delrec*delrec*delrec



          if(delrec.GT.-10.0.AND.delrec.LT.-9.0) delcor=delcor/1.04
c          if(delrec.GT.-10.0.AND.delrec.LT.-9.0) delcor=1.00
          if(delrec.GT.-9.0.AND.delrec.LT.-8.0) delcor=delcor/0.99
c          if(delrec.GT.-9.0.AND.delrec.LT.-8.0) delcor=1.00
c          if(delrec.GT.-6.0.AND.delrec.LT.-5.0) delcor=delcor/1.012
          if(delrec.GT.9.0.AND.delrec.LT.10.0) delcor=delcor/0.99

c          write(6,*) delrec,delcor

CCCCC      Turn delta correction off    CCCCC

c          delcor = 1. 


CCCCCCC   Now Get RC corrections for event    CCCCCCCC

         if(newrc) then
c          call rc_mod(firstr,thetaini,thetacrad,hsev,1,rci,t1)
          rce = 1.
          if(target.EQ.1.AND.w2.LT.1.18) then
c           rci = 1.
c           rce = 1.
          endif
         else
          call rcint(firstr,thetaini,thetacrad,hsev,target,rci)
          call rcext(firstr,thetaini,thetacrad,hsev,target,rce)
          if(w2.LT.1.18) then 
           rci = 1.
           rce = 1.
          endif
         endif    
     
         if(.not.dorc) then
           rci = 1.
           rce = 1.
         endif
        

c         write(6,*) hse,thetaini,w2,rci,rce
          
c         if(fail_id.EQ.0) then
 
CCCCCCC          Calculate the efficiencies for Calorimeter and Cerenkov          CCCCCCCCC
 
c         eff_cal = 1.*exp(-0.0192/hse)
c         eff_cal = .999*exp(-0.0017114*hse**(-2.73))
c         eff_cer = .99622     


          eff_cer = 0.998                                      !!!!  Cer Efficiency  !!!!
c         eff_cal = .999*exp(-0.0039258*hse**(-1.4930))        !!!!  Cal Efficiency  !!!!
          eff_cal = 0.96503+0.75590E-01*hse-0.65283E-01*hse**2+
     &       0.26938E-01*hse**3-0.53013E-02*hse**4+0.39896E-03*hse**5  !!!! For Jan05 !!!!
          eff_cal = eff_cal*(1.-0.000018*rate)



CCCCCCC    Calculate positron background weight.     CCCCCCCCC

 
         positron_weight = 1.                                !!!!  reset just in case   !!!!
         tartemp = target
         call csb_jan05(ebeam,hsev,thetaini,tartemp,poscs)
         poscs = targetdata(5)*poscs       !!!  put in ub/sr/GeV  !!!                                      
         positron_weight = 1.+poscs/(born/rci/rce)
         if(.not.docs) positron_weight = 1.

CCCCCCC            Fill new Ntuple                   CCCCCCCCC

 777	 do j = 1, 18
          ntu(j) = ntuple_contents(j)
         enddo
           
         j = 19
         ntu(j) = born*phasespcor*delcor
         j = j+1
         ntu(j) = rci
         j = j+1
         ntu(j) = rce
         j = j+1
         ntu(j) = positron_weight
         j = j+1
         ntu(j) = hse
         j = j+1
         ntu(j) = hstheta
         j = j+1
         ntu(j) = sigmac
         j = j+1
         ntu(j) = q2
         j = j+1
         ntu(j) = w2
        
         call hcdir('//reconmc',' ') 
         call HFN(9040,ntu)
        
c        endif 

c        write(6,*)  sigave,trackeff,prescale,cltime,
c     &              eff_cer,eff_cal     

       enddo

       if(use_rcmod) then
         call rc_mod(firstr,thetacrad,thetacrad,hsev,1,t1,sigmacent)

       write(6,*) "test: ",thetacrad,hsev,t1  
       else
         call model_new(ebeam,hsec,thetacrad,target,sigmacent)  
       endif

       write(6,*) "TEST, SIGMACENT = ",sigmacent

c       write(6,*) "NEW: ",dxp,dyp,dep,maxev,sigave

cc       sigave = sigave/maxev
       sigave = sigave/ngen
       phase_space = 4.0*dxp*dyp*dep/1000.0
       sigtot = sigave*phase_space 
       lummc = ngen/sigtot*1000.0
       lumfract = lumdata/lummc

       fract = lumdata*phase_space/ngentot/1000.00   !!! scale factor for MC (data must be      !!!
                                                  !!! corrected for eff. and prescale.       !!!   
       write(6, *) "fract befor corrected: ", fract
       denscor = 1.0 - bmcur/100.0*0.002
	if(target.GT.2) denscor = 1.0
      
c       denscor = denscor/1.015  !!!!****  TEST DG Corr  ****!!!! 

       write(6,*) "H2 density correction:  ", denscor

       fract = fract*trackeff*trigeff*cltime*eltime*eff_cer*eff_cal/
     &          prescale*denscor


       write(6,*)"** ", eff_cer,eff_cal


CCCC   The factor denscor is for target boiling which reduces the      CCCC
CCCC   # of MC events which should be observed for cryotargts.         CCCC

       write(6,*) nevt,sigave,lumdata,lummc,lumfract
       write(6,*) "MC scale factor is:  ",fract      
       write(6,*) "Average Sigma is :  ", sigave     
       write(6,*) "Sigma at central kinematics is:  ",sigmacent 
       write(6,*) "test the sigmac : ", sigmac  
       write(6,*) "test rci at center: ", rcic  
       call hcdir('//reconmc',' ')
       call HROUT(9040,ICYCLE,' ')
       call HREND('reconmc')         !CERNLIB close file   
       close (1)
       close (2)
       close(19)
 
      end





















