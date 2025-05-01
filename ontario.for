      Subroutine FirstLine(UseDatabase,SiDatabase,SIDataFileName,isupcx,
     +    iwsx,nrcx,iadjx,iasix,plevx,dlevx,limplnx,limdisx,sympx,msi)
c
c This subroutine is only called once to initialize Lake Ontario forecast
c from "OntarioInit.for"
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      logical UseDatabase,SiDatabase
      integer msi(1)
      integer*4 day,yr,day1,yr1,ksi(8),jsi(13),md(8),jd(8)
      integer*4 isea(12,31),iwt(12,31)
      integer*4 isup,iflo,qdev,prevstlou
C     Some compilers squawk about common blocks that don't total a multiple of
C     8 bytes, so use ipad to satisfy them when necessary.
      integer*4 ipd, ipe
      character*2 sympx,symdx
      character month(12)*3,monts(12)*3,dumtxt*100,SIDataFileName*(*)
      integer GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +        SupVerbosity, MLVerbosity, OntVerbosity
      COMMON/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
      common/c1/month,monts
      common/d1/mon,day,yr,md,jd
      common/d2/xdislev,AveOntLev
      common/p1/ipd,nm1,nd1,ny1,plevc,dlevc,itotdev,nwt1785,nwts,limplan
      common/q1/isea,iwt,isup,iflo,qdev
      common/r1/iprnt,kcrt,iforce,isetflow
      common/r2/prevstlou
      common/t1/ipe,nm,nd,ny,plev,dlev,iws,iws1,nrcp,nrcd,lasi
      common/t2/iadj,n13,ksi,jsi
      common/ss/iswitch
c
c Getting Friday's date (since the Starting date is always a Saturday
      ny=yr
      call date10(mon,day,yr,jd,md,2,'b')
      nm=md(1)
      nd=jd(1)
      if(mon.eq.1.and.nm.eq.12)ny=yr-1
      call quarter(nm,nd,mmm,mqt)
      mqtr=(mmm-1)*4+mqt
      iswitch=0
c
c If "UseDatabase" or "SiDatabase" variable is true, then it will read 
c the supply indicators (ie. SI) from the database "wk5mSI.txt"
c
      if(UseDatabase.or.SiDatabase)then
c
c UseDatabase = true or SiDatabase = true
c
c Getting the date 13 weeks previous from start date
c
          open(7,file=SIDataFileName)
          read(7,100)dumtxt
          im=nm
          id=nd
          iy=yr
          do i=1,13
              call date10(im,id,iy,jd,md,8,'b')
              if(im.eq.1.and.md(1).eq.12)then
                  iy=iy-1
              endif
              im=md(1)
              id=jd(1)
          end do
c
c Setting the data to read the 12 previous week's supply indicators
c
          jer=1
          do while (jer.eq.1)
              read(7,117)im1,id1,iy1
              if(iy1.eq.iy)then
C change by MMM 12Aug2002     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C                 if(id1.eq.id.and.im1.eq.im)then
C                     exit
C                 endif
                  IF(id1 .EQ. id .AND. im1 .EQ. im) GOTO 4444
C end change 12Aug2002        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              endif
          end do
C
C change by MMM 12Aug2002     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 4444 CONTINUE
C end change 12Aug2002        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c Reading the Supply Indicators & Chg-in-Sup Indicators,
c Reading the starting week's data:
c Reading from file "SIDataFileName"
c
          iqtdone=0
          do i=13,1,-1
              call quarter(im,id,jm,jqt)
              jqtr=(jm-1)*4+jqt
              if(i.eq.1)then
                  if(UseDatabase)then
c
c Note that this reading only takes place is "UseDatabase" is true.
c     Reading the Entire line
c
                      read(7,111)isupcx,iwsx,iwsx1,nntsx,jsi(1),ksi(1),
     +                    iadjx,iasix,krcx,isadx,nrcx,sympx,limplnx,
     +                    plevx,dlevx,limdisx,itotdev,prevstlou
                  else
c
c Reading up to the Adjustments ONLY
c
                      read(7,111)isupcx,iwsx,iwsx1,nntsx,jsi(1),ksi(1),
     +                    iadjx,iasix
                  endif
                  if(iqtdone .eq. 0)then
                      if(jqtr .ge. 47 .OR. jqtr .le. 13)then
                          ipreadjx=iadjx
                          iqtdone=1
                      end if
                  endif
              else
                  if(i.gt.8)then
                      read(7,120)jsi(i),idum,kadjx
                  else
                      read(7,120)jsi(i),ksi(i),kadjx
                  endif
                  if(iqtdone .eq. 0)then
                      if(jqtr .ge. 47 .OR. jqtr .le. 13)then
                          ipreadjx=kadjx
                          iqtdone=1
                      end if
                  endif
              endif
              call date10(im,id,iy,jd,md,8,'f')
              if(im.eq.12.and.md(8).eq.1)then
                  iy=iy+1
              endif
              im=md(8)
              id=jd(8)
c              
          end do
          close(7)
c         
          if(mqtr .le. 17 .or. mqtr .eq. 48)then
              if((mqtr.gt.13.or.mqtr.eq.48) .and. ipreadjx.ge.0)then
                  iswitch=1
              endif
          endif
      else
c
c UseDatabase = false or SiDatabase = false
c
          do i=2,13
              jsi(i)=msi(i)
              if(i.le.8)then
                  ksi(i)=msi(i)-msi(i+13)
              endif
          end do
          xx=iwsx/17.85
          call round1(xx,x1,1)
          iwsx1=x1
          nntsx=iwt(nm,nd)
          jsi(1)=iwsx1-nntsx
C   XXX Error-array MSI is only dimensioned to 1, so changed
C          ksi(1)=jsi(1)-msi(14)
C   to
          ksi(1)=jsi(1)-msi(1)
          krcx=0
          isadx=0
          if(mqtr .le. 17 .or. mqtr .eq.48)then
              if((mqtr.gt.13.or.mqtr.eq.48) .and. iasix .ge. 0)then
                  iswitch=1
              endif
          endif
      endif
c
c plan computations
c
      jlimpx=isupcx-limplnx
      as=jlimpx*100./3228.94
      call round1(as,aout,1)
      astagpx=aout*.01
c
c Discretion computations
c
      jlimdx=limdisx-limplnx
      as=itotdev*100./3228.94
      call round1(as,aout,1)
      astagdx=aout*.01
c
c Getting next week's date and the quarter it lies in
c
      call date10(nm,nd,ny,jd,md,8,'f')
      mon1=md(8)
      day1=jd(8)
      if(nm.eq.12.and.mon1.eq.1)yr1=ny+1
      call quarter(mon1,day1,mqb,iqtb)
      ipb=(mqb-1)*4+iqtb
      symdx=' '
      plev=plevx
      dlev=dlevx
      iws=iwsx
      iws1=iwsx1
      nrcp=limplnx
      nrcd=limdisx
      lasi=iasix
      iadj=iadjx
      n13=jsi(13)
      xdislev=dlevx
      if(ontverbosity.ge.5)then
          do i=13,2,-1
              if(i.gt.8)then
                  write(753,102)jsi(i)
              else
                  write(753,102)jsi(i),ksi(i)
              endif
          end do
          if(UseDatabase.or.SiDatabase)then
              write(753,103)nm,nd,ny,isupcx,iwsx,iwsx1,
     +            nntsx,jsi(1),ksi(1),iadjx,iasix,krcx,isadx,nrcx,
     +            limplnx,sympx,jlimpx,astagpx,plevx,limdisx,symdx,
     +            jlimdx,itotdev,astagdx,dlevx
          else 
              write(753,104)nm,nd,ny,isupcx,iwsx,iwsx1,
     +            nntsx,jsi(1),ksi(1),iadjx, lasi,           nrcx,
     +            limplnx,sympx,jlimpx,astagpx,plevx,limdisx,symdx,
     +            jlimdx,itotdev,astagdx,dlevx
          endif
      endif
      return
  100 format(a)
  102 format(30x,2i5)
  103 format(i3,'/',i2.2,'/',i4,i4,i6,i5,i4,4i5,i6,i5,2i5,1x,
     +    a,i4,f5.2,f8.2,i5,1x,a,i4,i5,f6.2,2f7.2)
  104 format(i3,'/',i2.2,'/',i4,i4,i6,i5,i4,4i5,  11x,2i5,1x,
     +    a,i4,f5.2,f8.2,i5,1x,a,i4,i5,f6.2,2f7.2)
C     Apparently unused line commented by mmm 13Mar09 to hush compiler
C  106 format(i3,'/',i2.2,'/',i4,f6.2)
  111 format(10x,11i6,6x,a,i6,2f6.0,i6,6x,i6,6x,i6)
  117 format(2(i2,1x),i4)
  120 format(34x,3i6)
      end
      Subroutine Plan58d(Onnts, OntDeviation, LouisOntFlowDiff, 
     +   StLouisRetard, DesPrFlow, MilleFlow, OntForceFlow,
     +   OntMeanLev, LakeOntarioFlow, Qdpmi, Qstfran, Qstmau,
     +   QRich, RJetty1, RVarennes, RSorel, RPierre, R3Rivers, 
     +   RBatiscan, Tidal)
      implicit integer*4 (i-n)
      implicit real*8 (a-h)
      implicit real*8 (o-z)
C     Apparently unused variables removed by mmm 13Mar09 to hush compiler
C     integer*4 fmon,fday,fyr,totday,elim
      integer*4 fmon,fday,fyr,       elim
C     Apparently unused variables removed by mmm 13Mar09 to hush compiler
C     integer*4 day,yr,daye,yre,md(8),jd(8)
      integer*4 day,yr,          md(8),jd(8)
      integer*4 isea(12,31),iwt(12,31),ksi(8),jsi(13)
      integer*4 q_stlaw_ont_diff,q_stlouis,q_desprairies
      integer*4 q_milleiles,ice_retard,prevstlou
C     Apparently unused variables removed by mmm 13Mar09 to hush compiler
C     integer*4 isup,iflo,qdev,erieout,Locsup
      integer*4 isup,iflo,qdev
C     Some compilers squawk about common blocks that don't total a multiple of
C     8 bytes, so use ipad to satisfy them when necessary.
      integer*4 ipd, ipe
c Inputs:
c Onnts- The sum of NiagaraConverted, SupplyON, and FlowFromWelland.
c     This includes all Lake Ontario water balance terms, except
c     outflow.  
c
      DOUBLE PRECISION Onnts, DesPrFlow, MilleFlow, OntDeviation
C     Apparently unused variables removed by mmm 13Mar09 to hush compiler
C     DOUBLE PRECISION StLouisRetard, StLouisFlow, OntForceFlow  
      DOUBLE PRECISION StLouisRetard,              OntForceFlow  
      DOUBLE PRECISION LouisOntFlowDiff
c
c Output:
c OntMeanLev      - Average Lake Ontario level for period
c LakeOntarioFlow - The Lake Ontario outflow for the routing period.
      DOUBLE PRECISION OntMeanLev, LakeOntarioFlow
c
      common/az/fmon,fday,fyr
      common/ck1/alevmon,alevtot,totv
      common/d1/mon,day,yr,md,jd
      common/d2/xdislev,AveOntLev
      common/p1/ipd,nm1,nd1,ny1,plevc,dlevc,itotdev,nwt1785,nwts,limplan
      common/p2/iasi,isi,icsi,elim
      common/q1/isea,iwt,isup,iflo,qdev
      common/r1/iprnt,kcrt,iforce,isetflow
      common/r2/prevstlou
      common/rr1/q_stlaw_ont_diff,q_stlouis,q_desprairies
      common/rr2/ice_retard,q_milleiles
      common/t1/ipe,nm,nd,ny,plev,dlev,iws,iws1,nrcp,nrcd,lasi
      common/t2/iadj,n13,ksi,jsi
c
c Setting the proper input data for model
c
      xx=Onnts*.1
      call round1(xx,x1,1)
      isup=x1
c      
      xx=OntForceFlow*.1
      call round1(xx,x1,1)
      iflo=x1
c      
      xx=OntDeviation*.1
      call round1(xx,x1,1)
      qdev=x1
c      
      xx=LouisOntFlowDiff*.1
      call round1(xx,x1,1)
      q_stlaw_ont_diff=x1
c
c Rounding QDpmi, QStFrancois, QstMaurice, QRichelieu (keeping them in cms)
c
      call round1(QDpmi,x1,1)
      IQDpmi=x1
c      
      call round1(QStfran,x1,1)
      IQstfran=x1
c      
      call round1(Qstmau,x1,1)
      IQstmau=x1
c      
      call round1(Qrich,x1,1)
      IQrich=x1
c
c Converting StLouisRetard, DesPrFlow & MilleFlow into 10cms (from cms).
c
      xx=StLouisRetard*.1
      call round1(xx,x1,1)
      ice_retard=x1
c      
      xx=DesPrFlow*.1
      call round1(xx,x1,1)
      q_desprairies=x1
c      
      xx=MilleFlow*.1
      call round1(xx,x1,1)
      q_milleiles=x1
c
c start of weekly computations
c
      call ComputeWeek(IQdpmi, IQstfran, IQstmau, IQRich, RJetty1, 
     +   RVarennes, RSorel, RPierre, R3Rivers, RBatiscan, Tidal)
c
c re-initialize the current week's values to next week's values
c
      nm=nm1
      nd=nd1
      ny=ny1
      plev=plevc
      dlev=dlevc
      iws=nwt1785
      iws1=nwts
      nrcp=limplan
      nrcd=limplan
      lasi=iasi
      do m=13,2,-1
          j=m-1
          jsi(m)=jsi(j)
          if(m.le.8)then
              ksi(m)=ksi(j)
          endif
      end do
      jsi(1)=isi
      ksi(1)=icsi
      n13=jsi(13)
      xdislev=dlevc
      OntMeanLev=AveOntLev
      LakeOntarioFlow=q_stlouis-q_stlaw_ont_diff
      prevstlou=q_stlouis
      return
      end
      Subroutine ComputeWeek(IQdpmi,IQstfran,IQstmau,IQRich,RJetty1, 
     +    RVarennes,RSorel,RPierre,R3Rivers,RBatiscan,Tidal)
      implicit integer*4 (i-n)
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      integer*4 md(8),jd(8)
      integer*4 isea(12,31),iwt(12,31),ksi(8),jsi(13)
      integer*4 elim,onoff
      integer*4 q_stlaw_ont_diff,q_stlouis,q_desprairies
      integer*4 isup,iflo,qdev
      character*2 symp,symd
      character month(12)*3,monts(12)*3
      integer*4 mtnn(12)
      integer*4 nlimMon(20),nlimDay(20),nlimYear(20),nlim
      INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +        SupVerbosity, MLVerbosity, OntVerbosity
C     Some compilers squawk about common blocks that don't total a multiple of
C     8 bytes, so use ipad to satisfy them when necessary.
      integer*4 ipd, ipe

      logical DCompute
      common/DComp/DCompute
      COMMON/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
      common/ae/nlimMon,nlimDay,nlimYear,nlim
      common/b1/mtnn
      common/c1/month,monts
      common/ck1/alevmon,alevtot,totv
      common/d2/xdislev,AveOntLev
      common/p1/ipd,nm1,nd1,ny1,plevc,dlevc,itotdev,nwt1785,nwts,limplan
      common/p2/iasi,isi,icsi,elim
      common/q1/isea,iwt,isup,iflo,qdev
      common/r1/iprnt,kcrt,iforce,isetflow
      common/rr1/q_stlaw_ont_diff,q_stlouis,q_desprairies
      common/t1/ipe,nm,nd,ny,plev,dlev,iws,iws1,nrcp,nrcd,lasi
      common/t2/iadj,n13,ksi,jsi
      common/ss/iswitch
      common/pd1/onoff
c
c Begin computations:
c -------------------
c
      call date10(nm,nd,ny,jd,md,8,'f')
      nm1=md(8)
      nd1=jd(8)
      ny1=ny
      if(nm.eq.12.and.md(8).eq.1)ny1=ny1+1
      call quarter(nm1,nd1,mq,iqt)
      ipn=(mq-1)*4+iqt
      nwt1785=iws-iws1+isup
      wtx=nwt1785/17.85
      call round1(wtx,wt,1)
      nwts=wt
      nnts=iwt(nm1,nd1)
      isi=nwts-nnts
      icsi=isi-n13
c
c check winter months to see if sea. adj. is to remain constant:
c    "iadj" is the counter to keep track of the value of the change-
c    in-supply indicator during mid-december.  If "iadj" is negative,
c    it will remain constant until 1st Qtr April.  If "iadj" is zero
c    or positive, it will remain constant until 1st Qtr May.
c If the current week is the 46th quarter and the following week is computed
c    to be the 47th quarter (ie. "ipn"), then the constant "iadj" is computed.
c    The 48th qtr is kept constant as is the 1st 17 qtrs of the following year.
c
      if(ontverbosity.ge.9)then
          write(812,201)nm1,nd1,ny1,ipn
          write(812,105)
      endif
      if(ipn.gt.13.and.ipn.ne.48)then
          if(ipn.le.17.and.iswitch.eq.1)go to 10
          if(ipn.eq.47)then
c
c Compute last week to see if it was also the 47th quarter.  If it is
c then skip since the constant was already computed.  (There is a chance 
c that you may have 2 consecutive weeks ending in the 3rd qtr December).
c
              call date10(nm1,nd1,ny1,jd,md,8,'b')
              call quarter(md(1),jd(1),ibm,ibq)
              ipc=(ibm-1)*4+ibq
              if(ipc.eq.47)then
                  if(ontverbosity.ge.9)then
                      write(812,202)
                  endif
                  goto 10
              endif
              jsum=0
              do j=1,8
                  jsum=jsum+ksi(j)
              end do
              xt=jsum*.111
              if(ontverbosity.ge.9)then
                  write(812,203)(ksi(l),l=1,8),xt
              endif
              call round1(xt,vt,1)
              iadj=vt
              if(iadj.lt.-20)iadj=-20
              if(iadj.gt.31) iadj= 31
              if(iadj.lt.0)then
                  iswitch=0
              else
                  iswitch=1
              endif
          else
              isum=0
              do j=1,3
                  isum=isum+ksi(j)
              end do
              xt=(isum+icsi)*.222
              if(ontverbosity.ge.9)then
                  write(812,204)icsi,(ksi(l),l=1,3),xt
              endif
              call round1(xt,vt,1)
              iadj=vt
              if(iadj.lt.-20)iadj=-20
              if(iadj.gt.31) iadj= 31
          endif
      else
          if(ontverbosity.ge.9)then
              write(812,206)
          endif
      endif
   10 if(ontverbosity.ge.9)then
          write(812,205)iadj
      endif
      iasi=isi+iadj
c
c Computing the rule curve, use last week's adj. supply indicator
c and last week's date for the month of computation (ie. "lasi")
c
      ir=nm
      if((nm.eq.1.or.nm.eq.8).and.nd.gt.27)ir=ir+1
      if(ontverbosity.ge.9)then
          write(812,106)
          write(812,210)plev,ir,lasi
      endif
      call rulecur(plev,ir,lasi,krc)
      isad=isea(nm,nd)
      nrc=krc+isad
c
c Compute the date for the 2nd week - this is to compute the "p"
c and "p*" values to be added to the current week's supply indicator
c and the quarter it lies in
c
      call date10(nm1,nd1,ny1,jd,md,8,'f')
      nm2=md(8)
      nd2=jd(8)
      ny2=ny1
      if(nm1.eq.12.and.md(8).eq.1)ny2=ny2+1
      call quarter(nm2,nd2,mq2,iqt2)
      ipn2=(mq2-1)*4+iqt2
c
c Compute limits for PLAN level and flow
c
      onoff=0
      call limits(mq,iqt,ipn,nrc,nrcp,jsi(1),plev,limplan,symp)
      ajl=isup-limplan
      jlimp=ajl
      as=ajl*100./3228.94
      call round1(as,aout,1)
      astagp=aout*.01
      plevc=plev+astagp
c
c
c Compute limits for Actual level and flow
c
c Forecast Mode:
c     Compute limits for discretionary level and flow
c     compute plan level before discretionary level
c
c Check for any forced flows, if not then check to see if any incremental
c deviations were assigned, otherwise compute limits for discretion
c
      onoff=1
      if(iflo.gt.0)then
          limdis=iflo
          symd='F '
          if(onoff.eq.1)then
              if(ontverbosity.eq.9)then
                  write(812,126)
                  write(812,121)limdis,symd
              endif
          endif
          go to 15
      endif
      if(qdev.ne.0)then
          limdis=limplan+qdev
          symd='FD'
          if(onoff.eq.1)then
              if(ontverbosity.eq.9)then
                  write(812,126)
                  write(812,122)limdis,symd,limplan,qdev
              endif
          endif
          go to 15
      endif
      if(isetflow.eq.1)then
          if(itotdev.lt.0.and.dlev.lt.75.37)then
              limdis=limplan
              symd=symp
              if(onoff.eq.1)then
                  if(ontverbosity.eq.9)then
                      write(812,126)
                      write(812,123)limdis,symd
                  endif
              endif
              goto 15
          endif
      endif
      call limits(mq,iqt,ipn,nrc,nrcd,jsi(1),dlev,limdis,symd)
   15 jlimd=limdis-limplan
      sum1=itotdev+jlimd
      itotdev=sum1
      dif=sum1*100./3228.94
      call round1(dif,aout,1)
      astagd=aout*.01
      dlevc=plevc-astagd
c
c Calling the lake st. lawrence routine to check for the forebay
c elevation ( must be greater than the criteria, ie. alert or
c minimum level)
c
      planlev=plevc
      dislev=dlevc
      limd=limdis
      limp=limplan
c
c Comuting upstream (ie. Ogdensburg, Cardinal, Iroquois Dam, Morrisburg and
c Long Sault)
c
      if(kcrt.ne.0)then
          call lawren(nm1,nd1,ny1,limd,limp,itotdev,dislev,planlev,
     +        jdis,dlevx1,jj,isum,as,jx)
          if(jx.ne.0)then
              limdis=jdis
              jlimd=jj
              itotdev=isum
              astagd=as
              dlevc=dlevx1
              symd='A '
          endif
      endif
c
c DCompute = TRUE, then compute downstream levels at Montreal and others
c (ie. Jetty1, Varennes, Sorel, Lac St. Pierre, 3 Rivers & Batiscan)
c
      if( DCompute ) then
          q_stlouis=limdis+q_stlaw_ont_diff
          IQstl=q_stlouis*10
          call downstream(nm1,nd1,ny1,limdis,IQstl,IQdpmi,
     +        IQstfran,IQstmau,IQRich,RJetty1,RVarennes,RSorel,
     +        RPierre,R3Rivers,RBatiscan,Tidal)
      endif
c
c For Simulation Only: Total Deviations are eliminated &
c computed level set equal to Actual/Discretionary level.
c
      if(elim.eq.1)then
          jj=0
          do i=1,nlim
              if(nm1.eq.nlimMon(i).and.nd1.eq.nlimDay(i).and.
     +        ny1.eq.nlimyear(i))then
                  jj=1
C change by MMM 12Aug2002     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C                     exit
                  GOTO 4444
C end change 12Aug2002        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              endif
          end do
C change by MMM 12Aug2002     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 4444 CONTINUE
C end change 12Aug2002        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          if(jj.eq.1)then
              plevc=dlevc
              itotdev=0
          endif
      endif
      if(ontverbosity.ge.5)then
          write(753,104)nm1,nd1,ny1,isup,nwt1785,nwts,
     +        nnts,isi,icsi,iadj,iasi,krc,isad,nrc,limplan,symp,
     +        jlimp,astagp,plevc,limdis,symd,jlimd,itotdev,astagd,
     +        dlevc
          if(ontverbosity.eq.9)then
              write(812,125)
          endif          
      endif
      return
C     Apparently unused line commented by mmm 13Mar09 to hush compiler
C  102 format(i3,'/',i2.2,'/',i4,f6.2)
  104 format(i3,'/',i2.2,'/',i4,i4,i6,i5,i4,4i5,i6,i5,2i5,1x,
     +    a,i4,f5.2,f8.2,i5,1x,a,i4,i5,f6.2,2f7.2)
  105 format(6x,'Averaging of Change-in-Sup for Adjustments:')
  106 format(6x,'Rule Curve Computations:')
C     Apparently unused line commented by mmm 13Mar09 to hush compiler
C  110 format(i5,i3,' - ',f6.2)
  121 format(6x,'ACTUAL Flow limits:',/,
     + 6x,'Forced Flows : flow=',i4,1x,a)
  122 format(6x,'ACTUAL Flow limits:',/,
     + 6x,'Deviation Flows : flow=',i4,1x,a,' [Plan flow (',i4,
     + ') + deviations (',i4,')]')
  123 format(6x,'ACTUAL Flow limits:',/,
     + 6x,'Setting actual flow to plan flow (Tot Dev < 0 & ',
     + 'L limit): flow=',i4,1x,a)
  125 format(80('_'))
  126 format(6x,74('-'))
  201 format(i2,'/',i2.2,'/',i4,'   Qtr=',i3)
  202 format(9x,'Bypass')
  203 format(9x,'summing last 8 "ksi" and "x .111"',/,
     +       9x,8i4,'; ave=',f6.2)
  204 format(9x,'summing last 4 "ksi" and "x .222"',/,
     +       9x,4i4,'; ave=',f6.2)
  205 format(9x,'Rounded "iadj" =',i4)
  206 format(9x,'ipn <= 13 or ipn = 48 ====> bypass')
  210 format(9x,'level=',f6.2,' month=',i3,' lasi=',i4)
      end
      subroutine limits(mon,qt,tqt,nrc,nrco,si,xlev,iq,sym)
      implicit integer*4 (a-z)
      real*8 xlev
      integer*4 isea(12,31),iwt(12,31)
      integer*4 isup,iflo,qdev,prevstlou
      character*2 sym,symax,symin,symm,symp
      character ztem(5)*100
      INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +        SupVerbosity, MLVerbosity, OntVerbosity
      COMMON/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
      common/q1/isea,iwt,isup,iflo,qdev
      common/r1/iprnt,kcrt,iforce,isetflow
      common/r2/prevstlou
      common/pd1/onoff
      common/pd2/ztem
c
c Computing Max & Min Preproject limits
c
      call maxp(tqt,si,qmaxp)
      call minp(tqt,si,qminp,symp)
c
c Comparing the rest of the limits
c
      call maxl(mon,qt,xlev,qmaxl)
      call maxi(tqt,qmaxi)
      call minm(tqt,qminm)
c
c--------------------------------------------------------------------      
c Checking the "J" limit or keeping the "RC" limit
c "imax" will store this value. "nrco" is last week's flow
c
      il=abs(nrc-nrco)
      if(il.gt.57)then
          jlim=nrco+57
          if(nrc.lt.nrco)then
              jlim=nrco-57
          endif
          imax=jlim
          symm='J '
      else
          imax=nrc
          symm='RC'
      endif
c--------------------------------------------------------------------      
c Getting the smallest of the maximum limits : L, P and I
c "jmax" will store this value.
c
      jmax=qmaxl
      symax='L '
c
c Using last week's Lake St. Louis flow as a check
c
      if(tqt.le.14.or.(prevstlou.gt.977.and.tqt.le.28))then
          if(qmaxp.ne.0)then
              if(qmaxp.lt.jmax)then
                  jmax=qmaxp
                  symax='P '
              endif
          endif
      endif
      if(qmaxi.ne.0)then
          if(qmaxi.lt.jmax)then
              jmax=qmaxi
              symax='I '
          endif
      endif
c--------------------------------------------------------------------      
c Getting the Biggest of the minimum limits : M and P*
c "imin" will store this value.
c
      imin=qminm
      symin='M '
      if(imin.lt.qminp)then
          imin=qminp
          symin=symp
      endif
c--------------------------------------------------------------------      
c Checking to see if the MAXimum < MINimum limit
c COMPARING MAXimum vs MINimum: 
c    If the Maximum Limit ("L", "I", or "P") is less than the Minimum 
c    Limit ("P" or "M"), then the Maximum Limit is the governing limit 
c    unless the Minimum Limit is the "M". The Minimum "M" Limit is the 
c    absolute minimum.
c
      if(jmax.lt.imin)then
          iq=jmax
          sym=symax
          if(symin.eq.'M ')then
              iq=imin
              sym=symin
          endif
      else
c
          if(imax.gt.jmax)then
c          
c Checking to see if either "J" or "RC" is GREATER than the MAXimum
c
              iq=jmax
              sym=symax
          elseif(imax.lt.imin)then
c          
c Checking to see if either "J" or "RC" is LESS than the MINimum
c
              iq=imin
              sym=symin
          else
c          
c Checking to see if either "J" or "RC" is in the MIDDLE
c
              iq=imax
              sym=symm
          endif
      endif
c--------------------------------------------------------------------      
      if(ontverbosity.eq.9)then
          write(812,115)
          if(onoff.eq.0)then
c
c This set is generated for the Plan Level
c
              write(812,107)xlev
          else
c
c This set is generated for the Actual Level
c
              write(812,108)xlev
          endif
          write(812,110)nrco,nrc,jlim,qmaxl,qmaxp,qminm,qminp,qmaxi
          write(812,111)iq,sym
          do i=1,5
              write(812,100)ztem(i)
          end do
          if(onoff.eq.1)then
              write(812,115)
          endif
      endif
      return
  100 format(a)
  107 format(6x,'PLAN Flow limits:  Plan Level =',f6.2)
  108 format(6x,'ACTUAL Flow limits:  Actual Level =',f6.2)
  110 format(6x,'Limits:  lastWk  RC    J    L    P    M    P*   I',/
     + 15x,i4,1x,8i5)
  111 format(6x,'SELECTED flow =',i5,1x,a)
  115 format(6x,74('-'))
      end
      subroutine lawren(nm1,nd1,ny1,limdis,limplan,idev,dlev,plev,
     +    jdis,dlevx,jj,isum,as,jx)
c
c Check Subroutine "relation" for the latest date of the coefficients.
c The Kingston level = Ontario level - 0.03 metres
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +        SupVerbosity, MLVerbosity, OntVerbosity
C     Apparently unused variables removed by mmm 13Mar09 to hush compiler
C     real*8 crt1(5),crt2(5),al(5),montreal_lev
      real*8 crt1(5),crt2(5),al(5)
      character month(12)*3,monts(12)*3
      logical DCompute
      common/DComp/DCompute
      COMMON/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
      common/c1/month,monts
      common/d2/xdislev,AveOntLev
      common/r1/iprnt,kcrt,iforce,isetflow
c
c Alert and Minimum depths for Ogdensburg, Cardinal, Iroquois Dam,
c Morrisburg and Long Sault Dam, respectively.
c
      data crt1/73.70,73.24,73.00,72.60,72.39/
      data crt2/74.02,73.26,73.00,72.44,71.93/
      jx=0
      idl=(xdislev+.009)*100.
   10 idle=(dlev+.009)*100.
      xd=(idle+idl)*.5
      call round1(xd,vd,1)
      AveOntLev=vd*.01
      ave1=AveOntLev-0.03
      flow=limdis*10
      call relation(ave1,flow,al)
      if(iforce.eq.1)goto 20
      do i=1,5
          if(kcrt.eq.1)then
              if(al(i).lt.crt1(i))goto 15
          else
              if(al(i).lt.crt2(i))goto 15
          endif
          go to 17
   15     limdis=limdis-1
          jj=limdis-limplan
          sum=idev+jj
          isum=sum
          dd=sum*100./3228.94
          call round1(dd,aout,1)
          as=aout*.01
          dlev=plev-as
          jdis=limdis
          dlevx=dlev
          jx=1
          goto 10
   17 end do
   20 if(ontverbosity.ge.5)then
          write(759,101)nm1,nd1,ny1,limdis,AveOntLev,al
      endif
      return
  101 format(1x,i2,'/',i2.2,'/',i4,i5,6f7.2)
      end
      Subroutine Downstream(nm,nd,ny,limdis,IQstl,IQdpmi, 
     +    IQstfran,IQstmau,IQRich,RJetty1,RVarennes,RSorel, 
     +    RPierre,R3Rivers,RBatiscan,Tidal)
      implicit integer*4 (i-n)
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      Integer*4 IQstl,IQdpmi,IQstfran,IQstmau,IQrich
      Double Precision RJetty1,RVarennes,RSorel,RPierre,R3Rivers,
     +    RBatiscan,Tidal,LevJetty1,LevVarennes,LevSorel,
     +    LevLacStPierre,Lev3Rivers,LevBatiscan
      integer   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +          SupVerbosity, MLVerbosity, OntVerbosity
      common/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
c
c New Downstream Equations, all flows are in m3s.
c
c Jetty 1:
      xx = (RJetty1 * (0.001757*IQstl + 0.000684*IQdpmi + 
     +   0.001161*IQstfran+0.000483*IQstmau))** 0.6587 + 0.9392*Tidal
      call round1(xx,x1,-2)
      LevJetty1=x1
c     
c Varennes:
      xx = (RVarennes * (0.001438*IQstl + 0.001377*IQdpmi + 
     +   0.001442*IQstfran+0.000698*IQstmau))** 0.6373 + 1.0578*Tidal
      call round1(xx,x1,-2)
      LevVarennes=x1
c
c Sorel:
      xx = (RSorel * (0.001075*IQstl + 0.001126*IQdpmi + 
     +   0.001854*IQstfran+0.000882*IQstmau))** 0.6331 + 1.2770*Tidal
      call round1(xx,x1,-2)
      LevSorel=x1
c
c Lac St. Pierre:
      xx = (RPierre * (0.000807*IQstl + 0.001199*IQdpmi + 
     +   0.001954*IQstfran+0.000976*IQstmau))** 0.6529 + 1.4772*Tidal
      call round1(xx,x1,-2)
      LevLacStPierre=x1
c
c Trois Rivieres:
      xx = (R3Rivers * (0.000584*IQstl + 0.000690*IQdpmi + 
     +   0.000957*IQRich+0.001197*IQstfran+0.000787*IQstmau))**0.7042 + 
     +   1.5895*Tidal
      call round1(xx,x1,-2)
      Lev3Rivers=x1
c
c Batiscan:
      xx = (RBatiscan * (0.000422*IQstl + 0.000553*IQdpmi +
     +   0.000903*IQRich+0.001023*IQstfran+0.000682*IQstmau))**0.6941 + 
     +   1.8303*Tidal
      call round1(xx,x1,-2)
      LevBatiscan=x1
c      
      if(ontverbosity.ge.5)then
          write(760,110)nm,nd,ny,limdis,LevJetty1,LevVarennes,
     +        LevSorel,LevLacStPierre,Lev3Rivers,LevBatiscan
          if(ontverbosity.eq.9)then
c             write(812,115)
              write(812,101)
              write(812,102)IQstl,IQdpmi,IQstfran,IQstmau,IQRich
              write(812,103)
              write(812,104)RJetty1,RVarennes,RSorel,RPierre,R3Rivers,
     +            RBatiscan,Tidal
              write(812,105)
              write(812,106)LevJetty1,LevVarennes,LevSorel,
     +            LevLacStPierre,Lev3Rivers,LevBatiscan
          endif
      endif
      return
  101 format(6x,'Montreal Computations:',/,
     +    11x,' Qstl    Qdpmi   Qstfran   Qstmau   Qrich   (cms)',/,
     +    10x,5('--------',1x))
  102 format(10x,5(i6,3x))
  103 format(/,10x,'RJetty1   RVaren   RSorel  RPierre  R3Rivers  ',
     +    'RBatis   Tidal',/,
     +    9x,7('--------',1x))
  104 format(9x,7(f7.3,2x))
  105 format(/,10x,'LJetty1   LVaren   LSorel  LPierre  L3Rivers  ',
     +    'LBatis (Metres)',/,
     +    9x,6('--------',1x))
  106 format(10x,6(f6.2,3x))
  110 format(1x,i2,'/',i2.2,'/',i4,i5,6(f6.2,1x))
C  115 format(6x,74('-'))
      end
      subroutine relation(kingston,q,dnlev)
c
c These Coefficients were developed in November 24, 1992
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      real*8 dnlev(1),coef(3,5),kingston
      data coef/63.2800,2.0925,0.4103,
     +          19.4908,2.3981,0.4169,
     +          23.6495,2.2886,0.4158,
     +          23.9537,2.2450,0.3999,
     +          22.9896,2.2381,0.3870/
      do i=1,5
          a0=coef(1,i)
          a1=coef(2,i)
          a2=1./coef(3,i)
          xxn=kingston-(q/(a0*(kingston-62.4)**a1))**a2
          call round1(xxn,xout,-2)
          dnlev(i)=xout
      end do
      return
      end
      subroutine rulecur(xi,mo,lasi,qrule)
c
c This computes the Basic Rule Curve discharge using
c     xi    = level, metric
c     mo    = month
c     lasi  = last week's adjusted supply indicator
c     qrule = new basic rule curve flow
c
      implicit integer*4 (a-z)
      real*8 xi,xo,xx,q,qo,con
      INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +        SupVerbosity, MLVerbosity, OntVerbosity
      COMMON/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
      xo=(xi-70)*100.
      call round1(xo,xx,1)
      if(xx.ge.447.or.(mo.gt.1.and.mo.lt.8))then
c
c Computed level greater than or equal to 74.31 m. or months between
c February and July inclusive
c
          con=20*.3048*.3048
          if(ontverbosity.ge.9)then
              write(812,211)xx,con
          endif
      else
c
c Computed level less than 74.31 m. between August and January inclusive
c
          con=20*.3048*.3048*5
          if(ontverbosity.ge.9)then
              write(812,212)xx,con
          endif
      endif
      if(lasi.gt.-45)then
          q=680+(xx-447)*con+lasi
          if(ontverbosity.ge.9)then
              write(812,213)q
          endif
      elseif(lasi.ge.-51.and.lasi.le.-45)then
          q=680+(xx-447)*con+(lasi+45)*3-45
          if(ontverbosity.ge.9)then
              write(812,214)q
          endif
      else
          q=680+(xx-447)*con+(lasi+51)*4-63
          if(ontverbosity.ge.9)then
              write(812,215)q
          endif
      endif
      call round1(q,qo,1)
      qrule=qo
      if(ontverbosity.ge.9)then
          write(812,216)qrule
      endif
      return
  211 format(9x,'xx=',f6.0,' Con=',f9.6,'  [20*.3048*.3048]')
  212 format(9x,'xx=',f6.0,' Con=',f9.6,'  [20*.3048*.3048*5]')
  213 format(9x,'q =',f6.1,'  [680+(xx-447)*con+lasi]')
  214 format(9x,'q =',f6.1,'  [680+(xx-447)*con+(lasi+45)*3-45]')
  215 format(9x,'q =',f6.1,'  [680+(xx-447)*con+(lasi+51)*4-63]')
  216 format(9x,'rc=',i6)
      end
      subroutine quarter(mon,day,mq,qt)
c
c Converting the days to quarter months
c
      implicit integer*4 (a-z)
      dimension i1(3),i2(3),ii(3)
      data i1/11,18,26/
      data i2/10,17,24/
      if(mon.eq.2)then
          do i = 1,3
              ii(i)=i2(i)
          end do
      else
          do i = 1,3
              ii(i)=i1(i)
          end do
      endif
      mq=mon
      qt=4
      if(day.le.3)then
          mq=mon-1
          if(mq.eq.0)mq=12
      elseif(day.le.ii(1))then
          qt=1
      elseif(day.le.ii(2))then
          qt=2
      elseif(day.le.ii(3))then
          qt=3
      endif
      return
      end
      subroutine maxl(mon,qt,xlev,qmaxl)
c
c Maximun flow "L" Limit
c This is based on quarter-months.
c
      implicit integer*4 (a-z)
      INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +        SupVerbosity, MLVerbosity, OntVerbosity
      real*8 xlev,qx,qo
      character ztem(5)*100
      common/pd2/ztem
      COMMON/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
      if(mon.le.3)then
          if(mon.eq.1)then
              if(xlev.lt.74.43)then
                  if(xlev.lt.74.22)then
                      qx=595
                  else
                      qx=(xlev-74.43+4.6725)*133.33333
                  endif
              else
                  qx=623
              endif
              goto 15
          elseif(mon.eq.2.and.qt.eq.1)then
              if(xlev.lt.74.42)then
                  if(xlev.lt.74.34)then
                      if(xlev.lt.74.22)then
                          qx=595
                      else
                          qx=(xlev-74.43+4.6725)*133.33333
                      endif
                  else
                      qx=(xlev-74.54+0.87143)*910.
                  endif
              else
                  qx=680
              endif
              goto 15
          elseif(mon.eq.3.or.(mon.eq.2.and.qt.eq.4))then
              if(xlev.lt.74.92)then
                  if(xlev.lt.74.42)goto 10
                  qx=(xlev-74.92+3.63671)*218.
              else
                  qx=793
              endif
              goto 15
          elseif(mon.eq.2.and.qt.ne.1)then
              if(xlev.lt.74.66)then
                  if(xlev.lt.74.42)goto 10
                  qx=(xlev-74.92+3.63671)*218.
              else
                  qx=736
              endif
              goto 15
          endif
   10     if(xlev.lt.74.34)then
              if(xlev.lt.74.22)then
                  qx=595
              else
                  qx=(xlev-74.43+4.6725)*133.33333
              endif
          else
              qx=(xlev-74.54+0.87143)*910.
          endif
      else
          if(xlev.lt.75.13)then
              if(xlev.lt.74.70)then
                  if(xlev.lt.74.54)then
                      if(xlev.lt.74.34)then
                          if(xlev.lt.74.22)then
                              qx=595
                          else
                              qx=(xlev-74.43+4.6725)*133.33333
                          endif
                      else
                          qx=(xlev-74.54+0.87143)*910.
                      endif
                  else
                      qx=(xlev-74.70+3.18095)*262.5
                  endif
              else
                  qx=(xlev-75.13)*100.+878.
              endif
          else
              qx=878
          endif
      endif
   15 call round1(qx,qo,1)
      qmaxl=qo
      if(ontverbosity.eq.9)then
          write(ztem(1),110)xlev,qmaxl
      endif
      return
  110 format(9x,'Max L : Level=',f6.2,', flow=',i4)
      end
      subroutine maxp(qt,si,qmaxp)
c
c Maximum flow "P" limit
c this is based on quarter month
c
      implicit integer*4 (a-z)
      dimension ipp(48)
      character ztem(5)*100
      INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +        SupVerbosity, MLVerbosity, OntVerbosity
      common/pd2/ztem
      COMMON/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
      data ipp/4*0,9*702,716,728,733,739,745,750,
     +         753,2*756,3*759,756,753,750,20*0/
      if(ipp(qt).ne.0)then
          qmaxp=ipp(qt)+si
      else
          qmaxp=0
      endif
      if(ontverbosity.eq.9)then
          write(ztem(2),110)ipp(qt),si,qmaxp
      endif
      return
  110 format(9x,'Max P : Constant=',i4,', SI=',i4,', PP flow=',i4)
      end
      subroutine maxi(qt,qmaxi)
c
c Maximum flow "I" limit, this is set thru "inidata.txt"
c this is based on quarter month
c
      implicit integer*4 (a-z)
      character ztem(5)*100
      INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +        SupVerbosity, MLVerbosity, OntVerbosity
      integer*4  q_stlaw_ont_diff,q_stlouis,q_desprairies
      common/rr1/q_stlaw_ont_diff,q_stlouis,q_desprairies
      common/pd2/ztem
      common/maxmin/maxilimit0,minplimit0
      COMMON/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
      ichk=0
      if(qt.ge.47)then
          if(maxilimit0.gt.0)then
              qmaxi=maxilimit0
              ichk=1
          else
c
c Total flow at Lake St. Louis = 793 cms x 10,
c Lake Ontario flow = 793 - (StLaw-Ont Diff)
c
              qmaxi=793-q_stlaw_ont_diff
          endif
      else
          qmaxi=0
      endif
      if(ontverbosity.eq.9)then
          if(ichk.eq.0)then
              write(ztem(3),110)q_stlaw_ont_diff,qmaxi
          else
              write(ztem(3),111)qmaxi
          endif
      endif
      return
  110 format(9x,'Max I : StLaw-Ont Diff=',i4,', flow=',i4)
  111 format(9x,'Max I : Forced flow=',i4)
      end
      subroutine minm(qt,qminm)
c
c Minimum flow "M" limit
c this is based on quarter month
c
      implicit integer*4 (a-z)
      dimension limit(48)
      character ztem(5)*100
      INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +        SupVerbosity, MLVerbosity, OntVerbosity
      common/pd2/ztem
      COMMON/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
      data limit/4*595,4*586,4*578,8*532,4*538,16*547,4*561,4*595/
      qminm=limit(qt)
      if(ontverbosity.eq.9)then
          write(ztem(4),110)qminm
      endif
      return
  110 format(9x,'Min M : flow=',i4)
      end
      subroutine minp(qt,si,qminp,symp)
c
c Minimum flow "P*" limit (this is based on quarter month)
c In this subroutine, The "P*" is distinquished from little "p" by:
c     P* = supply indicator plus PP constant  (for lake Ontario)
c     p  = 637 minus the (Ont-St.Law diff)/7  (for St. Law River)
c
      implicit integer*4 (a-z)
      character ztem(5)*100,symp*2
      INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +        SupVerbosity, MLVerbosity, OntVerbosity
      common/pd1/onoff
      common/pd2/ztem
      common/maxmin/maxilimit0,minplimit0
      COMMON/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
      real*8 x1,x2
      integer*4 ipa(48)
      integer*4 q_stlaw_ont_diff,q_stlouis,q_desprairies
      common/rr1/q_stlaw_ont_diff,q_stlouis,q_desprairies
      data ipa/18*0,643,657,671,685,694,700,705,711,714,716,719,722,
     +         4*725,722,714,705,700,694,688,683,680,677,3*674,0,0/
      ichk=0
      symp='P*'
      if(ipa(qt).ne.0)then
          qmp=ipa(qt)+si
          qminp=qmp
          if(minplimit0.gt.0)then
              if(qmp.gt.minplimit0)then
                  qminp=minplimit0
              endif
              ichk=1
          else
              x1=637.-q_stlaw_ont_diff/6.
              call round1(x1,x2,1)
              vflo=x2
              if(qminp.gt.vflo)then
                  qminp=vflo
              endif
          endif
      else
          qmp=0
          vflo=0
          qminp=0
      endif
      if(ontverbosity.eq.9)then
          if(ichk.eq.0)then
              write(ztem(5),110)ipa(qt),si,qmp,vflo
          else
              write(ztem(5),111)qminp
          endif
      endif
      return
  110 format(9x,'Min P*: Constant=',i4,', SI=',i4, ', PP flow=',i4,
     + ', StLaw-Ont Diff flow=',i4)
  111 format(9x,'Min P*: Forced flow=',i4)
      end
      subroutine round1(val1,result1,k)
c
c (Corps' way of rounding)
c Rounding up if ODD, rounding down if EVEN - MUST be within range 
c of -6 <= k <= 6.
c     val1   - value
c     result - rounded value
c     k      - decimal place to be rounded
c
      implicit integer*4 (a-z)
      real*8 val1,result1,xval,xx,tens,zx
      xval=abs(val1)
      if(k.gt.0)then
          n=k-1
      else
          n=k
      endif
      tens=10.**n
      xx=xval/tens
      ix=xx
      zx=xx-ix
c
      if(ABS(zx-.5).le.0.00000001)then
          ires=mod(ix,2)
          if(ires.gt.0)then
              add=1
          else
              add=0
          endif
      elseif(zx.lt.0.5)then
          add=0
      else
          add=1
      endif
      result1=(ix+add)*tens
      if(val1.lt.0.)result1=-result1
      return
      end
      subroutine date10(mon,iday,nyr,kday,imon,nday,let)
c
c Finding Calendar Dates:
c     mon  = month
c     iday = day
c     nyr  = year
c     kday = array to contain the days
c     imon = array to contain the months
c     nday = number of day back or ahead
c     let  = 'B' or 'b', Backward
c          = 'F' or 'f', Forward
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      integer*4 kday(1),imon(1),idat(12)
      character*(*) let
      data idat/31,28,31,30,31,30,31,31,30,31,30,31/
      ii=0
      if((let.eq.'B'.or.let.eq.'b').and.mon.eq.3)call leap1(nyr,ii)
      if((let.eq.'F'.or.let.eq.'f').and.mon.eq.2)call leap1(nyr,ii)
      idat(2)=28
      if(ii.eq.1)idat(2)=29
      k=mon
c
c Backwards
c
      if(let.eq.'B'.or.let.eq.'b')then
          id=iday+1
          do i=nday,1,-1
    5         id=id-1
              if(id.eq.0)then
                  k=k-1
                  if(k.eq.0)k=12
                  id=idat(k)+1
                  goto 5
              else
                  kday(i)=id
                  imon(i)=k
              endif
          end do
      else
c
c Forwards
c
          id=iday-1
          do i=1,nday
   10         id=id+1
              if(id.gt.idat(k))then
                  k=k+1
                  if(k.gt.12)k=1
                  id=0
                  goto 10
              else
                  kday(i)=id
                  imon(i)=k
              endif
          end do
      endif
      return
      end
      subroutine days10(monb,idb,yrb,mone,ide,yre,ndays)
c
c Computing the total number of days between two dates
c     idb,monb,yrb - beginning day, month and year
c     ide,mone,yre - ending day,    month and year
c     ndays - total number of days in-between the two dates inclusive
c
c Note: Does NOT include the starting day
c     examples: ( 1/ 1/95 to 1/5/95 => ndays=4)
c     examples: (12/31/94 to 1/5/95 => ndays=5)
c
      implicit integer*4 (a-z)
      real*8 cum
      integer*4 mtnn(12)
      common/b1/mtnn
      if(monb.eq.mone.and.yrb.eq.yre)then
          ndays=ide-idb
      else
          if((yrb.lt.yre).or.(yrb.eq.yre.and.monb.le.mone))then
              m1=monb
              y1=yrb
              m2=mone
              y2=yre
              cum=-idb
              rem=ide
              sign1=1
          else
              m1=mone
              y1=yre
              m2=monb
              y2=yrb
              cum=-ide
              rem=idb
              sign1=-1
          endif
          flg=1
          do while(flg.eq.1)
              mtn=mtnn(m1)
              if(m1.eq.2)then
                  call leap1(y1,ii)
                  if(ii.eq.1)mtn=29
              endif
              cum=cum+mtn
              m1=m1+1
              if(m1.gt.12)then
                  m1=1
                  y1=y1+1
              endif
              if(m1.eq.m2.and.y1.eq.y2)then
                  cum=cum+rem
                  goto 10
              endif
          end do
   10     ndays=cum*sign1
      endif
      return
      end
      subroutine leap1(nyr,ii)
c
c This subroutine checks if the year is a leap year
c     ii = 0, not a leap year
c        = 1, leap year
c
c Note: divisible by 400 ==> if yes, it is a leap year
c     divisible by 4 but NOT by 100 ==> if yes, it is a leap year
c     divisible by 4 AND by 100 ==> if yes, it is Not a leap year
c                                   (ie. 1700)
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      i400=mod(nyr,400)
      i4  =mod(nyr,4)
      i100=mod(nyr,100)
      ii=0
      if(i400.eq.0.or.(i4.eq.0.and.i100.ne.0))ii=1
      return
      end
      block data comdata
      implicit integer*4 (i-n)
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      character month(12)*3,monts(12)*3
      integer*4 isea(12,31),iwt(12,31)
      integer*4 isup,iflo,qdev
      integer*4 mtnn(12)
      common/b1/mtnn
      common/c1/month,monts
      common/q1/isea,iwt,isup,iflo,qdev
      data mtnn/31,28,31,30,31,30,31,31,30,31,30,31/
      data month/'JAN','FEB','MAR','APR','MAY','JUN',
     +           'JUL','AUG','SEP','OCT','NOV','DEC'/
      data monts/'Jan','Feb','Mar','Apr','May','Jun',
     +           'Jul','Aug','Sep','Oct','Nov','Dec'/
c
c "isea" is next week's Seasonal Adjustments
c
      data (isea(1,i), i=1,31)/4*0,27*-17/
      data (isea(2,i), i=1,29)/29*-17/
      data (isea(3,i), i=1,31)/27*-17,4*-23/
      data (isea(4,i), i=1,30)/4*-23,7*-28,8*-34,7*-40,4*-45/
      data (isea(5,i), i=1,31)/4*-45,7*-51,20*-57/
      data (isea(6,i), i=1,30)/26*-57,4*-51/
      data (isea(7,i), i=1,31)/4*-51,7*-45,8*-40,8*-34,4*-28/
      data (isea(8,i), i=1,31)/4*-28,7*-23,8*-17,8*-11,4*-6/
      data (isea(9,i), i=1,30)/4*-6,7*0,19*6/
      data (isea(10,i),i=1,31)/19*6,12*11/
      data (isea(11,i),i=1,30)/4*11,15*17,11*23/
      data (isea(12,i),i=1,31)/11*23,16*17,4*0/
c
c Obtained from metric graph (by Reg Young) - 1985
c
      data (iwt(1,i),i=1,31)/18*663,3*662,3*661,7*660/
      data (iwt(2,i),i=1,29)/29*660/
      data (iwt(3,i),i=1,31)/11*660,661,662,663,664,665,2*666,667,2*668,
     +                       669,2*670,671,672,673,675,677,679,680/
      data (iwt(4,i),i=1,30)/682,684,686,687,689,690,691,693,694,696,
     +                       697,699,700,701,702,703,705,706,707,708,
     +                       709,710,711,712,713,714,715,2*716,717/
      data (iwt(5,i),i=1,31)/718,2*719,720,721,723,724,725,726,727,728,
     +                       3*729,2*730,2*731,732,2*733,734,2*735,736,
     +                       736,3*737,2*738/
      data (iwt(6,i),i=1,30)/3*739,3*740,3*741,10*742,2*743,3*744,
     +                       6*745/
      data (iwt(7,i),i=1,31)/4*745,2*744,3*743,9*742,3*741,2*740,2*739,
     +                       2*738,737,736,2*735/
      data (iwt(8,i),i=1,31)/734,4*733,3*732,3*731,730,729,728,727,726,
     +                       2*725,724,723,2*722,721,720,2*719,718,
     +                       2*717,716,715/
      data (iwt(9,i),i=1,30)/2*714,713,2*712,711,2*710,709,2*708,707,
     +                       706,705,704,703,2*702,701,2*700,699,2*698,
     +                       2*697,696,695,694,693/
      data (iwt(10,i),i=1,31)/692,2*691,3*690,3*689,2*688,687,686,2*685,
     +                        684,2*683,682,681,2*680,679,678,3*677,
     +                        3*676,675/
      data (iwt(11,i),i=1,30)/675,3*674,2*673,3*672,2*671,3*670,2*669,
     +                        3*668,4*667,7*666/
      data (iwt(12,i),i=1,31)/18*666,3*665,3*664,7*663/
      end
