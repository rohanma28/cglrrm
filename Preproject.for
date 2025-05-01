      Subroutine Pre_project(Onnts, LouisOntFlowDiff, StLouisRetard,
     +    DesPrFlow, MilleFlow, OntMeanLev, LakeOntarioFlow,
     +    Qdpmi, Qstfran, Qstmau, QRich, RJetty1, RVarennes, 
     +    RSorel, RPierre, R3Rivers, RBatiscan, Tidal)
c
c     Weekly Preproject Computations
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h)
      implicit real*8 (o-z)
C
c     Inputs:
C     Onnts- The sum of NiagaraConverted, SupplyON, and FlowFromWelland.
C            This includes all Lake Ontario water balance terms, except
C            outflow.  
C     Apparently unused variables removed by mmm 13Mar09 to hush compiler
C     DOUBLE PRECISION Onnts, DesPrFlow, MilleFlow, OntDeviation
      DOUBLE PRECISION Onnts, DesPrFlow, MilleFlow
C     DOUBLE PRECISION StLouisRetard, StLouisFlow, OntForceFlow  
      DOUBLE PRECISION StLouisRetard
      DOUBLE PRECISION LouisOntFlowDiff
c
c     OntMeanLev      - Average Lake Ontario level for period
C     LakeOntarioFlow - The Lake Ontario outflow for the routing period.
      DOUBLE PRECISION OntMeanLev, LakeOntarioFlow
C
      integer*4 q_stlaw_ont_diff,q_stlouis,q_desprairies
      integer*4 q_milleiles,ice_retard
C     Apparently unused variables removed by mmm 13Mar09 to hush compiler
C     integer*4 erieout,Locsup,ice(12),md(8),jd(8)
      integer*4                 ice(12),md(8),jd(8)
      integer*4 mon,day,yr
      integer*4 isup,iflo,qdev
      integer*4 isea(12,31),iwt(12,31)
c
      INTEGER  GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +         SupVerbosity, MLVerbosity, OntVerbosity
C     Some compilers squawk about common blocks that don't total a multiple of
C     8 bytes, so use ipad to satisfy them when necessary.
C     Apparently unused variables removed by mmm 13Mar09 to hush compiler
C     integer*4 ipd, ipe
      integer*4      ipe
      logical DCompute
      common/DComp/DCompute
      COMMON/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
c
      common/ck1/alevmon,alevtot,totv
      common/d1/mon,day,yr,md,jd
      common/d2/xdislev,AveOntLev
      common/q1/isea,iwt,isup,iflo,qdev
      common/r1/iprnt,kcrt,iforce,isetflow
      common/rr1/q_stlaw_ont_diff,q_stlouis,q_desprairies
      common/rr2/ice_retard,q_milleiles
      common/t1/ipe,nm,nd,ny,plev,dlev,iws,iws1,nrcp,nrcd,lasi
      common/precond1/beglev
c
      DATA ice/-20,-20,-11,9*0/
c
c Setting the proper input data for model
c
      xx=Onnts*.1
      call round1(xx,x1,1)
      isup=x1
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
c StLouisRetard, DesPrFlow & MilleFlow are already in cms x 10 units
c
      xx=StLouisRetard
      call round1(xx,x1,1)
      ice_retard=x1
c      
      xx=DesPrFlow
      call round1(xx,x1,1)
      q_desprairies=x1
c      
      xx=MilleFlow
      call round1(xx,x1,1)
      q_milleiles=x1
c
c
c Computing the quarter-month
      if(ontverbosity.eq.9)then
          call quarter(nm,nd,mmm,mqt)
          mqtr=(mmm-1)*4+mqt
          write(812,105)nm,nd,ny,mqtr
      endif
c      
      call flowcomp(beglev,isup,ny,ice(nm),endlev,iflow,idif,dif)
      call lawren(nm,nd,ny,iflow,idum,idum,endlev,dum,idum,dum,idum,
     +    idum,dum,idum)
      if(ontverbosity.ge.5)then
          write(753,103)nm,nd,ny,beglev,isup,iflow,idif,dif,endlev
          if(ontverbosity.eq.9)then
              write(812,109)
          endif
      endif
      
c
c DCompute = TRUE, then compute downstream levels at Montreal and others
c (ie. Jetty1, Varennes, Sorel, Lac St. Pierre, 3 Rivers & Batiscan)
c
      if( DCompute ) then
          q_stlouis= iflow + q_stlaw_ont_diff 
          IQstl=q_stlouis*10
          call downstream(nm,nd,ny,iflow,IQstl,IQdpmi,IQstfran,
     +        IQstmau,IQRich,RJetty1,RVarennes,RSorel,RPierre,
     +        R3Rivers,RBatiscan,Tidal)
      endif
      if(ontverbosity.eq.9)then
          write(812,110)
      endif
c      
      OntMeanLev=endlev
      LakeOntarioFlow=iflow
      beglev=beglev+dif
c
c     Getting next week's date
c
      call date10(nm,nd,ny,jd,md,8,'f')
      if(nm.eq.12.and.md(8).eq.1)ny=ny+1
      nm=md(8)
      nd=jd(8)
c
      xdislev=endlev
c
      return
C     Apparently unused line commented by mmm 13Mar09 to hush compiler
C  101 format(1x,i2,'/',i2.2,'/',i4,i5,9f7.2)
C     Apparently unused line commented by mmm 13Mar09 to hush compiler
C  102 format(1x,i2,'/',i2.2,'/',i4,i5,8f7.2,'   N/A')
  103 format(5x,i2,'/',i2.2,'/',i4,f7.2,3i7,2f7.2)
  105 format(i2,'/',i2.2,'/',i4,'   Qtr=',i3)
  109 format(6x,69('-'))
  110 format(1x,74('-'))
      end
      subroutine flowcomp(wlev,isup,iyear,ik,plev,iq,idif,dif)
c
c     Calculate crustal movement to adjust 1903 flow to present
c
      implicit integer*4 (i-n)
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      INTEGER  GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +         SupVerbosity, MLVerbosity, OntVerbosity
      COMMON/RCPARI/GenVerbosity, TSDVerbosity, ChkVerbosity, 
     +              SupVerbosity, MLVerbosity, OntVerbosity
      adj=.0017*(iyear-1903)
      out1=57.719*(wlev-69.485)**1.5
      call round1(out1,outs,1)
      iq=outs
c
      if(ontverbosity.eq.9)then
          write(812,110)
          write(812,101)wlev,isup,adj,iq,ik
      endif
      do i=1,4
          idif=isup-iq
          dif=float(idif)/3228.94
          call round1(dif,di,-3)
          wlev1=wlev+di
          both=(wlev+wlev1)*.5
          call round1(both,alev,-3)
          out1=57.719*(alev-69.485)**1.5
          call round1(out1,outs,1)
          iq=outs+ik
          if(ontverbosity.eq.9)then
              if(i.eq.1)then
                  write(812,102)i,idif,di,both,outs,iq
              else
                  write(812,103)i,idif,di,both,outs,iq
              endif
          endif
      end do
      dif=di
      xx=alev+adj
      call round1(xx,plev,-2)
      if(ontverbosity.eq.9)then
          write(812,104)alev,plev
      endif
      return
  101 format(10x,'wlev, isup, adj, iq, ice:',/,
     + 12x,f6.2,i5,f5.2,2i5)      
  102 format(10x,'no, idif,di,both,outs,iq:',/,
     + 13x,i2,' -',i4,2f7.3,f6.1,i5)
  103 format(13x,i2,' -',i4,2f7.3,f6.1,i5)
  104 format(10x,'alev, plev:',2f6.2)
  110 format(6x,'Preproject Level and Flow Computation:')
      end

