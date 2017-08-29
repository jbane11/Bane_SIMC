 	program opt

	Implicit none

	integer	nwpawc
	real	memor
        integer lrecl,nevt,ievt,ierr,istat
        real*4 ntuple_contents(17)
      
	character*80 infile,outfile,title,ishit,new
        character*80 directory,icycle
        real*4 radcon,delini,yptarini,xptarini,zrec,delrec
        real*4 yptarrec,xptarrec,ytarini,ytarrec,dt 
        real*4 thetacrad,dxp,dyp,thetaini,thetarec,dtheta
        real*4 yield(16,20),ytheta(16,20),thbin,thinc,delinc
        real*4 dtave(16,20),tempvect(20),dt1(16,20),dt2(20)
        integer nentries,target,fail_id,runnum,maxev
	integer	i,j,m,NtupleSize, bank,ngen,ntrecl
        integer nthbins,ndelbins,delbin
        logical first,ron,goodfit,newrc,cryo,fixopt


        parameter (nwpawc=5000000)
        parameter (nentries = 18)        
        parameter(bank = 1000)
        parameter(title = 'RECONTUPLE')
        character*80 NtupleTag(nentries)
        real*4 ntu(nentries)
        common /pawc/ memor(nwpawc)

        ntrecl = 2048

        first = .true.

        fixopt = .true.     !!! set to true to apply correction to optics  !!!

        do i=1,16
         dt2(j) = 0.0
         do j=1,20
          dt1(i,j) = 0.0
         enddo
        enddo
        if(fixopt) then 
c         open(unit=23, file='optcor.dat',status='old')
         open(unit=24, file='dtcor2.dat',status='old')
         read(24,*) dt2
         do j=1,20     
          write(6,*) j,dt2(j)         
         enddo
c         do i=1,16
c          read(23,*) tempvect
c          do j=1,20
c           dt1(i,j) = tempvect(j)
c          enddo
c          write(6,*) tempvect
c         enddo
        endif

        thetacrad = 12.974*3.14159/180.
        nthbins = 20
        ndelbins = 16
        delinc = 1.
        thinc = 70./nthbins

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        call hlimit(nwpawc)
       
        lrecl = 0
        ishit = 'U'
        new = 'N'
        istat = 0

        infile = 'opt.rzdat'
        outfile = 'optcor.rzdat'

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
        NtupleTag(m) = 'dydz'
        m = m+1
        NtupleTag(m) = 'dxdz'
        m = m+1
        NtupleTag(m) = 'zini'
        m = m+1
        NtupleTag(m) = 'dppi'
        m = m+1
        NtupleTag(m) = 'dthi'
        m = m+1
        NtupleTag(m) = 'dphi'
        m = m+1
        NtupleTag(m) = 'zrec'
        m = m+1
        NtupleTag(m) = 'dppr' 
        m = m+1
        NtupleTag(m) = 'dthr' 
        m = m+1
        NtupleTag(m) = 'dphr' 
        m = m+1
        NtupleTag(m) = 'fail_id'
        m = m+1  
        NtupleTag(m) = 'x_stop'
        m = m+1 
        NtupleTag(m) = 'y_stop' 
        m = m+1 
        NtupleTag(m) = 'yini'
        m = m+1
        NtupleTag(m) = 'yrec' 
        m = m+1
        NtupleTag(m) = 'dtcor' 


        ntuplesize = m

        call hropen (2, 'reconmc', outfile,'N',ntrecl,istat) 
        call hrin (2, 99999,0)

        call HBOOKN(9040,'reconmc',NtupleSize,'reconmc',1000,NtupleTag)
        call hcdir('//optcor',' ')
        write(6,*) '                    New ntuple is:  ', outfile



        write(6,*)
        write(6,*) 'Number of events in master file =',nevt  
        maxev = nevt
        write(6,*)
        write(6,*) 'Number of events analyzing = ',maxev
        write(6,*)
 
       
        do ievt = 1, maxev          
c        do ievt = 1, 500000        
         if(ievt.EQ.10000) write(6,*) ' analyzed 10000 events'
         if(ievt.EQ.50000) write(6,*) ' analyzed 50000 events'
         if(mod(ievt,100000).EQ.0.)write(6,*) ' analyzed',ievt,' events'
         
         call hcdir('//MCntuple',' ')
         call hgnf(1, ievt, ntuple_contents, ierr)
         if (ierr .ne. 0) then
          write (6,*) 'hgnf err:', ierr
         endif 

         delini = ntuple_contents(6)          
         yptarini = ntuple_contents(7)
         xptarini = ntuple_contents(8)
         zrec = ntuple_contents(9)
         delrec = ntuple_contents(10)
         yptarrec = ntuple_contents(11) 
         xptarrec = ntuple_contents(12)
         fail_id = ntuple_contents(13)
         ytarini = ntuple_contents(16)
         ytarrec = ntuple_contents(17) 

         thetaini = acos(cos(thetacrad-yptarini/1000.)
     &              *cos(xptarini/1000.))   
         thetarec = acos(cos(thetacrad-yptarrec/1000.)
     &              *cos(xptarrec/1000.))
         dt = 1000.*(thetarec - thetacrad)
         dtheta = thetarec-thetaini

         delbin = int((delrec +8.)/delinc) + 1
         thbin = int((dt + 35.)/thinc) + 1
          
c         if(fail_id.EQ.0.AND.abs(delrec).LT.8.0.AND.abs(thetarec).
c     &     LT.35.0) write(6,*) dt,thbin,delrec,delbin

         if(fail_id.EQ.0.AND.abs(delrec).LT.8.0.AND.abs(thetarec).
     &     LT.35.0) then
          yield(delbin,thbin) = yield(delbin,thbin) + 1.
          ytheta(delbin,thbin) = ytheta(delbin,thbin) + 1.*dtheta
         endif


CCCCCCC            Fill new Ntuple                   CCCCCCCCC

         do j = 1, 17
          ntu(j) = ntuple_contents(j)
         enddo
        
CCCCCCC   Apply theta optics correction    CCCCCCC

c         ntu(18) = dt1(delbin,thbin)
          ntu(18) = dt2(thbin) 

CCCCCCC

         call hcdir('//reconmc',' ') 
         call HFN(9040,ntu)       
        enddo   

        if(.not.fixopt) then
         open(unit=20, file='dtopt.dat',status='new')

         do i=1,16
          do j=1,20
           tempvect(j) = 0.
           if(yield(i,j).EQ.0) yield(i,j) = 1
           dtave(i,j) = ytheta(i,j)/yield(i,j)
           if(abs(1000.*dtave(i,j)).GT.100.) dtave(i,j) = 0.
           write(6,*) i,j,dtave(i,j)
           tempvect(j) = 1000.*dtave(i,j)
          enddo
          write(20,8000) tempvect
         enddo
        endif

       call hcdir('//reconmc',' ')
       call HROUT(9040,ICYCLE,' ')
       call HREND('reconmc')         !CERNLIB close file   


       close (1)
       close (2)
       close(19)
 
 8000 format(20f7.3)

      end





















