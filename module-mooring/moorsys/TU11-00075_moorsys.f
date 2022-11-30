      SUBROUTINE MOORSYS(Iniflg,IHIS,time,Npj,Xpj,Rpj,RKpj)
C
      IMPLICIT  REAL*8  (A-H,O-Z)
C
      DIMENSION Xpj(6),Rpj(6),RKpj(6,6)
C
C     Iniflg = 0   :   Data in
C            = 1   :   Rpj(6)  Cal.
C            = 2   :   Rpj(6)  and RKpj(6,6) Cal.
C
C     IHIS = 0     :  Histerysis not Consider
C          = 1     :  Histerysis Consider  for Fender
C
      IF(Iniflg.ge.1)  go to 10
c                                                     ** STEP-1 **
      CALL INPUTP
C
      RETURN
C                                                     ** STEP-2 **
   10 CALL  ELEMNT( Iniflg ,IHIS , TIME , Npj ,Xpj )
C
      CALL  FINTOT( Npj , Rpj )
C
      IF(Iniflg.le.1)   return
c
      CALL  ASSEMP( Npj ,RKpj )
C
      RETURN
      END
      SUBROUTINE INPUTP
C
      IMPLICIT  REAL*8  (A-H,O-Z)
C
      COMMON / C00 / KR,KW,KM,IPAGE,IERR,NEP(20),NITLE(20),LTITL(20)
      COMMON / C01 / NBANE,NPONT,Npfend(50),Npfs(50),Npfe(50),Npcain(50)
     *              ,Npcs(50),Npce(50),Pitle(50),Xpi(50),Ypi(50),Zpi(50)
      COMMON / B00 /NBTYPE(400),NBIND(400),NOBN(400),KBMX(400),KBMI(400)
      COMMON / B01 / BDLT(400,24),BFOS(400,24),CDLT(400,24),CFOS(400,24)
      COMMON / B02 / NOFN(999),NFNB(999),GL
      COMMON / B03 / NOCN(999),NCNB(999)
      COMMON / B04 / Xfi(999),Yfi(999),Zfi(999),Thif(999)       
      COMMON / B05 / DLFN0(999),DLFN(999),FFN(999),RKFN(999),TFN(999,6)
      COMMON /B06/DLCN(999),FCN(999),FCNH(999),FICN(6,999),RKCN(6,6,999)
     *              ,TCN(999,6),FICN0(6,999),FCNMIN(999),FCNMAX(999)
      COMMON / B07 / UKFN(999),DKFN(999),UKCN(999),DKCN(999),DDFA(999)
      COMMON / B08 / RKBNI(400),UKBI(400,24),DKBI(400,24)
      COMMON / B09 / Cth0(999),Xci(999),Yci(999),Zci(999),Thic(999)
     *              ,Cdep(999),Cw(999),Cl(999),Chd(999)        
C
C     ----- 22/09/04 add for gfortran compiler by okubo
      CHARACTER(8) POINT
C     ---------------------------------------------------
      DIMENSION  POINT(24)  ,AA(3,3),V(3,3),W(3,3)
C
      DATA POINT /' POINT-1',' POINT-2',' POINT-3',' POINT-4',
     *            ' POINT-5',' POINT-6',' POINT-7',' POINT-8',
     *            ' POINT-9','POINT-10','POINT-11','POINT-12',
     *            'POINT-13','POINT-14','PONIT-15','POINT-16',
     *            'POINT-17','POINT-18','PONIT-19','POINT-20',
     *            'POINT-21','POINT-22','PONIT-23','POINT-24'/
C
C                              ** DATA INPUT **
      CALL LISTUP(KR,KW)
C
      READ(KR,100)  NITLE
      READ(KR,110)  NBANE,NPONT,GL,NEP
C
      IERR=0
C
      DO 10 I=1,400
       NBTYPE(I)=0
       NBIND(I)=0
   10 CONTINUE
C
      DO  12  I=1,999
      DLFN0(I)=0.0
      DLFN(I)=0.0
       FFN(I)=0.0
      UKFN(I)=0.0
      DKFN(I)=0.0
      DLCN(I)=0.0
       FCN(I)=0.0
      UKCN(I)=0.0
      DKCN(I)=0.0
      DDFA(I)=0.0
      FCNMIN(I)=-100000.0
      FCNMAX(I)= 100000.0
   12 CONTINUE
C
      IF(NBANE.EQ.0)  GO TO 25
C
      DO 20 I=1,NBANE
C
      READ(KR,200)  NOBN(I),KBMX(I),KBMI(I),RKBNI(I)
      J=NOBN(I)
      IF(NBIND(J).EQ.0) GO TO 14
      IERR=IERR+1
      WRITE(KW,209) J
   14 NBIND(J)=I
      KBM=KBMX(I)
      READ(KR,201) (BDLT(I,K),K=1,KBM)
      READ(KR,201) (BFOS(I,K),K=1,KBM)
      KBM1=KBM-1
      DO  16  K=1,KBM1
   16 UKBI(I,K)=(BFOS(I,K+1)-BFOS(I,K))/(BDLT(I,K+1)-BDLT(I,K))
C
      KBM=KBMI(I)
      IF(KBM.EQ.0)  GO  TO  20
      READ(KR,201)  (CDLT(I,K),K=1,KBM)
      READ(KR,201)  (CFOS(I,K),K=1,KBM)
      KBM1=KBM-1
      DO  18  K=1,KBM1
   18 DKBI(I,K)=(CFOS(I,K+1)-CFOS(I,K))/(CDLT(I,K+1)-CDLT(I,K))
C
   20 CONTINUE
C
   25 NFE=0
      NCE=0
C
      DO 40 IP=1,NPONT            
C
      READ(KR,120)  NO,Pitle(ip),Npfend(ip),Npcain(ip),Xpi(ip),Ypi(ip)
     *                                                ,Zpi(ip)
c
      if(Npfend(ip).eq.0)  go to 44
C
      Npfe(ip)=NFE+Npfend(ip)
      NFE=Npfe(ip)
      Npfs(ip)=Npfe(ip)-Npfend(ip)+1
C
      DO 42 i=Npfs(ip),Npfe(ip)
C
      READ(KR,130) NOFN(I),NFNB(I),Xfi(i),Yfi(i),Zfi(i),Thif(i)
C
C7   *****  TRANS-MATRIX, TFN
C
       X1=Xfi(I) 
       Y1=Yfi(I)
       Z1=Zfi(I)
       X2=X1+DCOS(Thif(I)*3.1415927/180.)
       Y2=Y1+DSIN(Thif(I)*3.1415927/180.)
       Z2=Z1       
C
      CALL  TRANSF(X1,Y1,Z1,X2,Y2,Z2,AA,V,W)
C
       J=1
        DO  76  K=1,3
         TFN(i,K) =AA(J,K)
         K2=K+3
         TFN(i,K2 ) =W(J,K)
   76   CONTINUE
C
      IF(NEP(1).LE.1)  GO  TO 41  
      WRITE(KW,202) X1,Y1,Z1,X2,Y2,Z2,(AA(1,K),K=1,3),(W(1,K),K=1,3)
C
   41 IF(NBIND(NFNB(I)).NE.0) GO TO 42
      IERR=IERR+1
      WRITE(KW,215)  NOFN(I),NFNB(I)
   42 CONTINUE
C
   44 if(Npcain(ip).eq.0)  go to 40
C
      Npce(ip)=NCE+Npcain(ip)
      NCE=Npce(ip)
      Npcs(ip)=Npce(ip)-Npcain(ip)+1
C
      DO 46 i=Npcs(ip),Npce(ip)
      READ(KR,140) NOCN(I),Cth0(I),Xci(i),Yci(i),Zci(i),Thic(i)
     *                    ,Cdep(i),Cw(i),Cl(i)
C
      CALL CASE3(Cdep(I),Cw(I),Cl(I),Cth0(I),Chd(I),ISPEC)
C
      IF(ISPEC.LE.2)  GO TO 46
      IERR=IERR+1
      WRITE(KW,270)  NO,NOCN(I),ISPEC,Cth0(I),Cdep(I),Cw(I),Cl(I),Chd(I)
   46 CONTINUE
C   
   40 continue
C
      IF(IERR.EQ.0)  GO TO 49
      WRITE(KW,216)  IERR
      STOP  'DATA ERROR'
C
C     OUTPUT OF INPUT DATA
C
   49 CALL WTITLE
C
      WRITE(KW,250)   NBANE,Npont,GL
C
      IF(NBANE.EQ.0)  GO TO 55
C
      WRITE(KW,229)
C
      DO 50 I=1,NBANE
      KBM=KBMX(I)
       WRITE(KW,230) I,NOBN(I),KBM,(POINT(K),K=1,KBM)
       WRITE(KW,231) (BDLT(I,K),K=1,KBM)
       WRITE(KW,232) (BFOS(I,K),K=1,KBM)
      KBM=KBMI(I)
      IF(KBM.EQ.0)  GO  TO  50
      WRITE(KW,240)  KBM,(POINT(K),K=1,KBM)
      WRITE(KW,231)  (CDLT(I,K),K=1,KBM)
      WRITE(KW,232)  (CFOS(I,K),K=1,KBM)
      WRITE(KW,241)  RKBNI(I)
   50 CONTINUE
C
   55 DO 60 IP=1,NPONT               
C
      WRITE(KW,260) IP,Pitle(ip),Npfend(IP),Npcain(IP),Xpi(IP),Ypi(ip)
     *                ,Zpi(ip),Npfs(ip),Npfe(ip),Npcs(ip),Npce(ip)
c
       if(Npfend(ip).eq.0)  go to 64
c
       WRITE(KW,233)
      DO 62 I=Npfs(ip),Npfe(ip)
       WRITE(KW,234) NOFN(I),NFNB(I),Xfi(i),Yfi(i),Zfi(i),Thif(i)
   62 CONTINUE
C
   64 if(Npcain(ip).eq.0)  go to 60
C
      WRITE(KW,235)
      DO 66 I=Npcs(ip),Npce(ip)
       WRITE(KW,238) NOCN(I),Cth0(i),Xci(i),Yci(i),Zci(i),Thic(i)
     *                              ,Cdep(i),Cw(i),Cl(i),Chd(I)
   66 CONTINUE
C
   60 CONTINUE
      RETURN
C
  100 FORMAT(20A4)
  110 FORMAT(2I5,F10.2,40X,20I1)
  120 FORMAT(I2,A8,2I5,3F10.2)
  130 FORMAT(2I5,4F10.2)
  140 FORMAT(I2,F8.2,7F10.2)
  200 FORMAT(3I5,F10.0)
  201 FORMAT(8F10.0)
  202 FORMAT(1H ,12F10.2)
  209 FORMAT(1H0,'** ERROR  BANE NO.',I3,' APPEARED AGAIN **')
  215 FORMAT(1H0,'** ERROR  FENDER OR CHAIN NO.',I3,
     *   '  HAS  UNDIFINED BANE NO.=',I3)
  216 FORMAT(1H0,I5,'ERROR WAS FOUND ,STOP CALCULATION')
  229 FORMAT(1H0,'*** BANE PROPERTY ***')
  230 FORMAT(1H ,'(I)  NOBN(I)  KBMX(I)  (UP)'
     *                 /1H ,I2,I7,I9,4X,8(2X,A8)/(1H ,22X,8(2X,A8)))
  231 FORMAT(1H ,13X,'DELTA(m) ',8F10.2/(1H ,22X,8F10.2))
  232 FORMAT(1H ,13X,'FORCE(kN)',8F10.2/(1H ,22X,8F10.2))
  233 FORMAT(/1H0,'*** FENDER INFORMATION ***' /1H , 'FENDER NO.',
     *'      NFNB','    Xfi(m)    Yfi(m)    Zfi(m)  Thif(deg)')
  234 FORMAT(1H ,I9,I10,4F10.2) 
  235 FORMAT(/1H0,'*** CHAIN INFORMATION ***' /1H , 'CHAIN  NO.',
     *'   T0(kN) ','    Xci(m)    Yci(m)    Zci(m)  Thic(deg)'
     *,'   Dep(m)   w(kN/m)     L(m)    Hd(m)')            
  238 FORMAT(1H ,I9,f10.2,4F10.2,F10.2,F10.5,2F10.2)               
  240 FORMAT(1H0,14X,'KBMI(I)  (DOWN)'/1H ,I18,4X,8(2X,A8)/(1H ,22X
     *          ,8(2X,A8)))
  241 FORMAT(1H0,10X,'KB(I)=',F10.2,'(kN/m)')
  250 FORMAT(1H0,'  NBANE=',I3,'  NPONT=',I3,'   GL=',F8.2,'(m)')
  260 FORMAT(/1H0,'PONTOON NO.=',I2,'  NAME=',A8,'  NUMBER OF FENDERS='
     *  ,I3,'   NUMBER OF CHAINS=',I3/1H0,4X,'(X,Y,Z)=',3F10.2
     * /1H0,10X,'Npfs,Npfe=',2I5,'   Npcs,Npce=',2I5)
  270 FORMAT(1H ,'CHain Data Error   Npj=',I3,'  CHain NO.=',I3,  
     *'   ISPEC=',I2  /' (Cth0,Cdep,Cw,Cl,Chd)',5F10.2)
      END
      SUBROUTINE  ELEMNT(Iniflg,IHIS,TIME,Npj,X)    
C
      IMPLICIT  REAL*8  (A-H,O-Z)
      COMMON / C00 / KR,KW,KM,IPAGE,IERR,NEP(20),NITLE(20),LTITL(20)
      COMMON / C01 / NBANE,NPONT,Npfend(50),Npfs(50),Npfe(50),Npcain(50)
     *              ,Npcs(50),Npce(50),Pitle(50),Xpi(50),Ypi(50),Zpi(50)
      COMMON / B00 /NBTYPE(400),NBIND(400),NOBN(400),KBMX(400),KBMI(400)
      COMMON / B01 / BDLT(400,24),BFOS(400,24),CDLT(400,24),CFOS(400,24)
      COMMON / B02 / NOFN(999),NFNB(999),GL
      COMMON / B03 / NOCN(999),NCNB(999)
      COMMON / B05 / DLFN0(999),DLFN(999),FFN(999),RKFN(999),TFN(999,6)
      COMMON /B06/DLCN(999),FCN(999),FCNH(999),FICN(6,999),RKCN(6,6,999)
     *              ,TCN(999,6),FICN0(6,999),FCNMIN(999),FCNMAX(999)
      COMMON / B07 / UKFN(999),DKFN(999),UKCN(999),DKCN(999),DDFA(999)
C
      DIMENSION X(6)
C
C  **   Iniflg =1  F CAL.                 **
C              =2  F CAL. AND K CAL.      **
C
C                IKK=1  MAKE K
      IKK=0
      IF(Iniflg.EQ.2) IKK=1
      
C
C    ***** FENDER  ***
C
      IF(Npfend(Npj).EQ.0) GO TO 222
C
      IF(NEP(12).GT.0)  WRITE(KW,930)
C
      DO 10 NF=Npfs(Npj),Npfe(Npj)
C
C      IF(INFE(NF) .EQ. 0) GO TO 111
C
       D=0.0
       DO 11  I=1,6
        D=D+TFN(NF,I)*X(I)
   11  CONTINUE
C
       DA= D      +DLFN0(NF)
       NN=NFNB(NF)
C
C    ***   FORCE AND   K  OF  EACH FENDER
C
      CALL  BFEND('FENDER  ',NF,NN,DLFN(NF),DDFA(NF),UKFN(NF),DKFN(NF)
     *                      ,  DA    ,  FA    , UKA    , DKA ,IHIS,Npj)
C
      FFN(NF)= FA
      DDFA(NF)=FA
      DLFN(NF)=DA
      UKFN(NF)=UKA
      DKFN(NF)=DKA
      RKFN(NF)=(UKA+DKA)/2.0
C     GO TO 10
C
C 111  FFN(NF)=0.0
C      RKFN(NF)=0.0
C     DLFN(NF)=0.0
C
   10 CONTINUE
C
      IF(NEP(2).EQ.0)   GO  TO  222
      WRITE(KW,900)   TIME,Npj,(NOFN(I),I=Npfs(Npj),Npfe(Npj))
      WRITE(KW,901)  (DLFN(I),I=Npfs(Npj),Npfe(Npj))
      WRITE(KW,910)  ( FFN(I),I=Npfs(Npj),Npfe(Npj))
      IF(NEP(2).LE.1) GO TO 222
      WRITE(KW,912)  (UKFN(I),I=Npfs(Npj),Npfe(Npj))
      WRITE(KW,914)  (DKFN(I),I=Npfs(Npj),Npfe(Npj))
      WRITE(KW,920)  (RKFN(I),I=Npfs(Npj),Npfe(Npj))
C
C    *****  CHAIN   *****
C
  222 IF(Npcain(Npj).EQ. 0) GO TO 233
C
      DO  30  NC=Npcs(Npj),Npce(Npj)
C
      CALL BCHAIN( Npj , NC,  IKK, X  )
C
   30 CONTINUE
C
      
      IF(NEP(3).EQ.0)   GO  TO  233
      WRITE(KW,800)   TIME,Npj,(NOCN(I),I=Npcs(Npj),Npce(Npj))
      WRITE(KW,801)  (DLCN(I),I=Npcs(Npj),Npce(Npj))
      WRITE(KW,810)  ( FCN(I),I=Npcs(Npj),Npce(Npj))
      WRITE(KW,811)  (FCNH(I),I=Npcs(Npj),Npce(Npj))
C
  233 RETURN
C
  900 FORMAT(1H0,4X,'FENDER DISP. AND FORCE    TIME=',F8.2,4X,'Npj=',I2
     *      /1H ,'FENDER NO.', 10I10/(1H ,10X,10I10))
  901 FORMAT(1H ,2X,'DLFN(I)=',10F10.4/(1H ,10X,10F10.4))
  910 FORMAT(1H ,2X,' FFN(I)=',10F10.1/(1H ,10X,10F10.1))
  912 FORMAT(1H ,2X,'UKFN(I)=',10F10.0/(1H ,10X,10F10.0))
  914 FORMAT(1H ,2X,'DKFN(I)=',10F10.0/(1H ,10X,10F10.0))
  920 FORMAT(1H ,2X,'RKFN(I)=',10F10.0/(1H ,10X,10F10.0))
  930 FORMAT(1H0,'SUB.ELEMNT CALL SUB.BANE',19X,'DB        DA     FB
     * FA    UKB    DKB    UKA    DKA    FAD    FUI    FDI    UKI    DKI
     * IND')
  800 FORMAT(1H0,4X,' CHAIN DISP. AND FORCE    TIME=',F8.2,4X,'Npj=',I2
     *      /1H ,' CHAIN NO.', 10I10/(1H ,10X,10I10))
  801 FORMAT(1H ,2X,'DLCN(I)=',10F10.4,1X/(1H ,10X,10F10.4,1X))
  810 FORMAT(1H ,2X,' FCN(I)=',10F10.2,1X/(1H ,10X,10F10.2,1X))
  811 FORMAT(1H ,2X,'FCNH(I)=',10F10.2,1X/(1H ,10X,10F10.2,1X))
      END
      SUBROUTINE BFEND(ANAME,IN,NN,DB,FB,UKB,DKB, DA,FA,UKA,DKA,IHIS
     *                      ,Npj)
C
      IMPLICIT  REAL*8  (A-H,O-Z)
      character*(4) aname
C
      COMMON / C00 / KR,KW,KM,IPAGE,IERR,NEP(20),NITLE(20),LTITL(20)
      COMMON / B00 /NBTYPE(400),NBIND(400),NOBN(400),KBMX(400),KBMI(400)
      COMMON / B01 / BDLT(400,24),BFOS(400,24),CDLT(400,24),CFOS(400,24)
      COMMON / B02 / NOFN(999),NFNB(999),GL
      COMMON / B08 / RKBNI(400),UKBI(400,24),DKBI(400,24)
C
      IB=NBIND(NN)
      IUMAX=KBMX(IB)-1
      IDMAX=KBMI(IB)-1
C         
      IF(DA.LT.BDLT(IB,1))  GO TO 33
      IF(IDMAX.EQ.-1)    GO  TO  100
      IF(IHIS.EQ.0)      GO  TO  100
C
C **   HISTERISIS  **
C
      DD=DA-DB
      IF(DD.LT.0.0)    GO  TO  50
C
C ***  ( UP )  ***
C
      FAD= FB + UKB*DD
C
      CALL FBANE( IB,IUMAX,IDMAX,DA, FAD, FUI, FDI , UKI, DKI, IND )
C
      GO  TO  ( 11, 12, 11, 14 )  , IND
C
   11 FA= FUI
      UKA=UKI
      DKA=RKBNI(IB)
      GO  TO  20
   12 FA= FAD
      UKA=RKBNI(IB)
      DKA=RKBNI(IB)
      GO  TO  20
   14 FA= FUI
      UKA=UKI
      DKA=UKA
      GO  TO  20
C
C ***  ( DOWN )  ***
C
   50 FAD= FB + DKB*DD
C
      CALL FBANE( IB,IUMAX,IDMAX,DA, FAD, FUI, FDI , UKI, DKI, IND )
C
      GO  TO  ( 11, 12, 23, 14 )  , IND
C
   23 FA= FDI
      UKA=RKBNI(IB)
      DKA=DKI
      GO  TO  20
C
   25 UKA=UKB
      DKA=DKB
      GO  TO  20
C
C **  NON-HISTERISIS  **
C
  100 CONTINUE
C
      DO  30  I=1,IUMAX
      IF(DA.GE.BDLT(IB,I).AND.DA.LT.BDLT(IB,I+1))  GO  TO  32
   30 CONTINUE
C
      I=IUMAX
C
   32 UKA=UKBI(IB,I)
      DKA= UKA
      FA= BFOS(IB,I) + UKA*( DA-BDLT(IB,I) )
      IND=0
      GO TO 20
C                         (  hanareru )
   33 UKA=0.0
      DKA=0.0
      FA=0.0
      IND=0
C
   20 IF( DA.GT.BDLT(IB,IUMAX+1))
     *        WRITE(KW,900)  ANAME,Npj,NOFN(IN),NN, DA
C
C  **  CHECK  WRITE  **
C
      IF(NEP(12).LE.0)  RETURN
      WRITE(KW,910) ANAME,Npj,NOFN(IN),DB,DA,FB,FA,UKB,DKB,UKA,DKA,
     *                  FAD,FUI,FDI,UKI,DKI,IND
C
      RETURN
  900 FORMAT(1H ,5X,'** WARNING  DISP. IS OUT OF RANGE. ',A8,'Npj=',I3,
     *'  FENDER NO.=',I3,  '  BANE NO.=',I3,'  DELTA=',F10.4)
  910 FORMAT(1H ,'SUB.BFEND',1X,A8,2I3,2F10.4,11F7.0,I2   )
      END
      SUBROUTINE FBANE( IB, IUMAX, IDMAX, DA,FAD, FUI,FDI,UKI,DKI,IND )
C
      IMPLICIT  REAL*8  (A-H,O-Z)
C
      COMMON / B01 / BDLT(400,24),BFOS(400,24),CDLT(400,24),CFOS(400,24)
      COMMON / B08 / RKBNI(400),UKBI(400,24),DKBI(400,24)
C
      DO  10  I=1,IUMAX
      IF(DA.GE.BDLT(IB,I).AND.DA.LT.BDLT(IB,I+1))  GO  TO  12
   10 CONTINUE
C
      I=IUMAX
C
   12 UKI=UKBI(IB,I)
      FUI= BFOS(IB,I) + UKI*( DA-BDLT(IB,I) )
C
      DO  20  I=1,IDMAX
      IF(DA.GE.CDLT(IB,I).AND.DA.LT.CDLT(IB,I+1))  GO  TO  22
   20 CONTINUE
C
      I=IDMAX
C
   22 DKI=DKBI(IB,I)
      FDI= CFOS(IB,I) + DKI*( DA-CDLT(IB,I) )
C
      IND=2
      IF(DABS(FUI-FDI).LE.0.001)   GO  TO  24
      IF(FAD.GE.FUI-0.001)   IND=1
      IF(FAD.LE.FDI+0.001)   IND=3
      RETURN
C
   24 IND=4
      RETURN
      END
      SUBROUTINE BCHAIN( Npj , NC, IKK, XXX  )
C --- [TU11-00035] commented out
C      USE NWTC_Library
C --- [TU11-00035] commented out

      IMPLICIT  REAL*8  (A-H,O-Z)
C
      COMMON / C00 / KR,KW,KM,IPAGE,IERR,NEP(20),NITLE(20),LTITL(20)
      COMMON / C01 / NBANE,NPONT,Npfend(50),Npfs(50),Npfe(50),Npcain(50)
     *              ,Npcs(50),Npce(50),Pitle(50),Xpi(50),Ypi(50),Zpi(50)
      COMMON /B06/DLCN(999),FCN(999),FCNH(999),FICN(6,999),RKCN(6,6,999)
     *              ,TCN(999,6),FICN0(6,999),FCNMIN(999),FCNMAX(999)
      COMMON / B09 / Cth0(999),Xci(999),Yci(999),Zci(999),Thic(999)
     *              ,Cdep(999),Cw(999),Cl(999),Chd(999)        
      COMMON / B10 /  W,SS,VF,HF,VD,HD,     DS      ,STL
     *               ,X(1000),Y(1000),TEN(1000),FAI(1000)
C
      DIMENSION  RK(6,6),F(6),BANE(4),DIR(3) ,AA(3,3),V(3,3),WW(3,3)
     *          ,XXX(6)

C --- [TU11-00032]
      REAL XX(6), TransMat(3,3)
C --- [TU11-00032]
C
C       IKK =1  MAKE K
C
      PAI=3.1415927
      TENR=PAI/180.
      TEND=180./PAI
C
      XI=Xci(NC)     
      YI=Yci(NC)          
      ZI=Zci(NC)     
      XJ=XI+DCOS(Thic(NC)*TENR)*Chd(NC)
      YJ=YI+DSIN(Thic(NC)*TENR)*Chd(NC)
      ZJ=ZI-Cdep(NC)  
C
C --- [TU11-00032] --- change transformation matrix
C      X1=XI+XXX(1)-XXX(6)*YI+XXX(5)*ZI
C      Y1=YI+XXX(2)+XXX(6)*XI-XXX(4)*ZI
C      Z1=ZI+XXX(3)+XXX(4)*YI-XXX(5)*XI

      XX(4)=XXX(4)
      XX(5)=XXX(5)
      XX(6)=XXX(6)
      
C --- [TU11-00073] --- change for large rotation.
C     CALL SmllRotTrans 
C    &     ( 'platform displacement',XX(4),XX(5),XX(6),TransMat )
      TransMat(1,1) =  COS(XX(6))*COS(XX(5))
      TransMat(1,2) =  SIN(XX(6))*COS(XX(5))
      TransMat(1,3) = -SIN(XX(5))
      TransMat(2,1) = -SIN(XX(6))*COS(XX(4))
     &                +COS(XX(6))*SIN(XX(5))*SIN(XX(4))
      TransMat(2,2) =  COS(XX(6))*COS(XX(4))
     &                +SIN(XX(6))*SIN(XX(5))*SIN(XX(4))
      TransMat(2,3) =  COS(XX(5))*SIN(XX(4))
      TransMat(3,1) =  SIN(XX(6))*SIN(XX(4))
     &                +COS(XX(6))*SIN(XX(5))*COS(XX(4))
      TransMat(3,2) = -COS(XX(6))*SIN(XX(4))
     &                +SIN(XX(6))*SIN(XX(5))*COS(XX(4))
      TransMat(3,3) =  COS(XX(5))*COS(XX(4))
C --- [TU11-00073] --- change for large rotation.
      X1=XXX(1)+TransMat(1,1)*XI+TransMat(2,1)*YI+TransMat(3,1)*ZI
      Y1=XXX(2)+TransMat(1,2)*XI+TransMat(2,2)*YI+TransMat(3,2)*ZI
      Z1=XXX(3)+TransMat(1,3)*XI+TransMat(2,3)*YI+TransMat(3,3)*ZI

C --- [TU11-00032] --- change transformation matrix
C
      SSS=DSQRT((X1-XJ)**2+(Y1-YJ)**2+(Z1-ZJ)**2)
C
      SS=Cl(NC)
       W=Cw(NC)
C
      VD=Z1-ZJ
      HD=DSQRT((X1-XJ)**2+(Y1-YJ)**2)
C                                     CHANGE  1992/2/8
C     IF(ILARGE.EQ.1)  GO TO 5
C
C--- start of Change on 2011/01/29 [TU11-00028]
C     X1=XI
C     Y1=YI
C     Z1=ZI
C--- end of Change on 2011/01/29 [TU11-00028]
C
    5 NBUN=50
C
      CALL CASE2(NBUN,ISPEC)
C
      GO TO (300,300,300,302,303)  ,ISPEC
C
  302 WRITE(KW,307) NC
      GO TO 306
  303 WRITE(KW,304)  NC,SS,SSS
  306 FAI(1)=PAI/4.0
      TEN(1)=FCNMAX(NC)
      HF    =FCNMAX(NC)
      GO TO 1000
C
C 300 IF(ISPEC.EQ.1)  WRITE(KW,309) NC
  300 CONTINUE                              
      IF(ISPEC.EQ.2.AND.NEP(11).GE.1) WRITE(KW,310)
      IF(ISPEC.EQ.3) WRITE(KW,305)
C
      IF(NEP(11).EQ.0)  GO TO 1000
      DFAI=FAI(1)*TEND
      WRITE(KW,312) NC,W,SS,HD,VD,HF,VF,TEN(1),DFAI,ISPEC
      IF(NEP(11).LE.3)  GO TO 1000
      WRITE(KW,315)
C
  304 FORMAT(1H0,6X,'*** WARNING, CABLE LENGTH IS TOO SHORT ***'/
     *       1H0,10X,'NC=',I3,'  (SS,SSS)',2F10.2)
  305 FORMAT(1H0,6X,'*** GIVEN CONDITIONS ARE NOT SATISFIED ( THE TOP ',
     *    'OF THE CABLE EXISTS BETWEEN THE TWO CABLE ENDS )  ***')
  307 FORMAT(1H0,6X,'*** WARNING, CATENARY PARAMETER CANNOT BE OBTAINED
     *  ', ' ( ARCCOSH VALUE CANNOT BE COMPUTED )  ***   NC=',I3)
  309 FORMAT(1H0,6X,'***  WARNING  CABLE',I3,' IS PIN PIN  ***')
  310 FORMAT(1H0,6X,'***  CABLE END  REACH THE BOTTOM  ***')
  312 FORMAT(1H0,3X,'NC=',I2,' (W,SS,HD,VD,HF,VF,T1,DFAI)',8F10.2,5X
     *  ,'ISPEC=',I2)
  315 FORMAT(1H0,7X,'NO.',10X,'S-CO.',10X,'X-CO.',10X,'Y-CO.',6X,
     *       'FAI(DEG.)',8X,'T.FORCE'/)
  320 FORMAT(1H ,5X,I5,5F15.3)
C
      NBUN1=NBUN+1
      DO 50 I=1,NBUN1
      S=DS*FLOAT(I-1)
      DFAI =FAI(I)*TEND
   50 WRITE(KW,320) I,S,X(I),Y(I),DFAI ,TEN(I)
C
 1000 FAI1=FAI(1)
      T1=TEN(1)
      T=HF/W
      SI=DSQRT(VD**2+2.0*VD*T)
C
C     ALFA=DATAN2(YJ-Y1,XJ-X1)
      ALFA=Thic(NC)*TENR            
      ALP =PAI-ALFA
C
      DIR(1)= DCOS(FAI1)*DCOS(ALFA)
      DIR(2)= DCOS(FAI1)*DSIN(ALFA)
      DIR(3)=-DSIN(FAI1)
C
      X2=X1+DIR(1)
      Y2=Y1+DIR(2)
      Z2=Z1+DIR(3)
C
      CALL TRANSF(X1,Y1,Z1,X2,Y2,Z2,AA,V,WW)
C
      J=1
      DO 76 K=1,3
      TCN(NC,K)=AA(J,K)
      K2=K+3
   76 TCN(NC,K2)=WW(J,K)
C
      DO 10 I=1,3
   10 F(I)=T1*DIR(I)
C
C --- start of change on 2011/01/29 [TU11-00029]
C     F(4)=-F(2)*Z1+F(3)*Y1 (original)
C     F(5)=-F(3)*X1+F(1)*Z1 (original)
C     F(6)=-F(1)*Y1+F(2)*X1 (original)
      F(4)=-F(2)*(Z1-XXX(3))+F(3)*(Y1-XXX(2))
      F(5)=-F(3)*(X1-XXX(1))+F(1)*(Z1-XXX(3))
      F(6)=-F(1)*(Y1-XXX(2))+F(2)*(X1-XXX(1))
C --- end of change on 2011/01/29 [TU11-00029]
C
      DO 12 I=1,6
      FICN(I,NC)= -F(I)
   12 FICN0(I,NC)=-F(I)
C
      DLCNDD=0.0
      DO 13 I=1,6
      F(I)=XXX(I)
   13 DLCNDD=DLCNDD+TCN(NC,I)*F(I)
C
      DLCN(NC)=DLCNDD
C
      FCNH(NC)=HF
      FCN(NC)= T1
C     HDCN(NC)=HF
C
      IF(IKK.EQ.0)  RETURN
C
      IF(T1.EQ.FCNMIN(NC))  RETURN
      IF(T1.EQ.FCNMAX(NC))  RETURN
C
      CALL CATENA(SS,SI,T,W,VD,HD,ALP,XI,YI,ZI,RK,ISPEC,BANE)
C
      DO 20 I=1,6
      DO 20 J=1,6
   20 RKCN(I,J,NC)=RK(I,J)
C
      IF(NEP(11).LE.2)  RETURN
C
      WRITE(KW,22)  DIR,BANE
   22 FORMAT(1H0,20X,'DIR(I)=',3F8.3,'   BANE(I)=',4F8.2)
C
      CALL MATOUT('(K) ',RK,KW)
C
      RETURN
      END
      SUBROUTINE  CASE2(NBUN,ISPEC)
C
      IMPLICIT  REAL*8  (A-H,O-Z)
      COMMON / B10 /  W,SS,VF,HF,VD,HD,     DS      ,STL
     *               ,X(1000),Y(1000),TEN(1000),FAI(1000)
C
      DIMENSION  DUMM(2,100)
C
      IF(SS.GT.DSQRT(HD**2+VD**2))  GO TO 11
      ISPEC=5
      RETURN
C
   11 CONTINUE
C
C     ***  CALCULATION OF CATENARY PARAMETER  ***
C
      IF(SS.LE.VD+HD)  GO TO 12
      NBUN=1
      ISPEC=2
      DS=VD
      HF=0.0001
      VF=W*VD
      X(1)=0.0
      Y(1)=0.0
      X(2)=0.0
      Y(2)=-VD
      TEN(1)=W*VD
      FAI(1)=90.0
      TEN(2)=0.0
      FAI(2)=0.0
      RETURN
C
   12 CPARA=1000.
   30 CONTINUE
      FUNCT=CPARA*DSINH(HD/CPARA)-SS
      IF(FUNCT.LE.0.0)  GO TO 31
      CPARA=CPARA+1000.
      GO TO 30
C
   31 CONTINUE
      DC=CPARA/2.0
      CPARA=CPARA-DC
      I=1
   10 CONTINUE
      FUNCT=CPARA*DSINH(HD/CPARA)-SS
      DC=DC/2.0
      IF(FUNCT.LE.0.)  CPARA=CPARA-DC
      IF(FUNCT.GT.0.)  CPARA=CPARA+DC
      IF(I.GE.15)  GO TO 20
      I=I+1
      GO TO 10
C
   20 CONTINUE
C
      SI=DSQRT( VD**2+2.*VD*CPARA)
      IF(SI.GT.SS)  GO TO 38
C
      CALL  CASE21(SS,VD,HD,CPARA)
C
      SI=DSQRT( VD**2+2.*VD*CPARA)
      AA=VD/CPARA+1.0
      XC=CPARA*DLOG(AA+DSQRT(AA**2-1.0))
      XB=XC
      HHD=XC
      SSS=SI
      ISPEC=2
      SAMP=1.0
      GO TO 39
C
   38 CALL  CASE22(SS,VD,HD,CPARA)
C
      BB=SS/(2.0*CPARA*DSINH(HD/(2.0*CPARA)))
      IF(BB.GE.1.0)  GO TO 37
      ISPEC=4
      RETURN
   37 XB=CPARA*DLOG(BB+DSQRT(BB**2-0.9999999))-HD/2.0
      IF(XB.LT.0.0)  ISPEC=3
      IF(XB.GE.0.0)  ISPEC=1
      HHD=XB
      SSS=SS
      SAMP=-1.0
   39 AA=HHD/CPARA
      YB=CPARA*(DCOSH(AA)-1.0)
      SB=CPARA*DSINH(AA)
      HF=W*CPARA
      DS=SSS/FLOAT(NBUN)
      H=HF
      N=1
   55 CONTINUE
      S=SB-DS*FLOAT(N-1)*SAMP
      AA=S/CPARA
      X0=CPARA*DLOG(AA+DSQRT(AA**2+1.0))
      X(N)=X0-XB
      Y(N)=CPARA*(DCOSH(X0/CPARA)-1.0)-YB
      IF(ISPEC.EQ.2) X(N)=-X(N)
      V=W*S
      TEN(N)=DSQRT(V**2+H**2)
      FAI(N)=DATAN(V/H)
      IF(N.GE.NBUN+1)  GO TO 54
      N=N+1
      GO TO 55
C
   54 VF=DSQRT(TEN(1)**2-HF**2)
C
      IF(ISPEC.EQ.2)  RETURN
C
      NBUN1=NBUN+1
      DO 60 I=1,NBUN1
      Y(I)=-Y(I)
      DUMM(1,I)=TEN(I)
      DUMM(2,I)=FAI(I)
   60 CONTINUE
C
      DO 65 I=1,NBUN1
      J=NBUN1+1-I
      TEN(I)=DUMM(1,J)
   65 FAI(I)=DUMM(2,J)
C
      VF=TEN(1)*DSIN(FAI(1))
      RETURN
      END
      SUBROUTINE  CASE21(SS,VD,HD,CPARA)
      IMPLICIT  REAL*8  (A-H,O-Z)
C
C     ***  CALCULATION OF CATENARY PARAMETER  ***
C
      DC=CPARA/2.0
      CPARA=CPARA-DC
      I=1
   10 CONTINUE
      AA=VD/CPARA+1.0
      XC=CPARA*DLOG(AA+DSQRT(AA**2-1.0))
      FUNCT=SS-DSQRT(VD**2+2.*VD*CPARA)+XC-HD
      DC=DC/2.0
      IF(FUNCT.LE.0.)  CPARA=CPARA+DC
      IF(FUNCT.GT.0.)  CPARA=CPARA-DC
      IF(I.GE.20)  GO TO 20
      I=I+1
      GO TO 10
C
   20 RETURN
      END
      SUBROUTINE  CASE22(SS,VD,HD,CPARA)
C
      IMPLICIT  REAL*8  (A-H,O-Z)
C
C     ***  CALCULATION OF CATENARY PARAMETER  ***
C
      A=HD/2.0
      B=-DSQRT(SS**2-VD**2)/2.0
      CPARA=1000.
   30 CONTINUE
      FUNCT=CPARA*DSINH( A/CPARA)+B
      IF(FUNCT.LE.0.0)  GO TO 31
      CPARA=CPARA+1000.
      GO TO 30
C
   31 CONTINUE
      DC=CPARA/2.0
      CPARA=CPARA-DC
      I=1
   10 CONTINUE
      FUNCT=CPARA*DSINH( A/CPARA)+B
      DC=DC/2.0
      IF(FUNCT.LE.0.)  CPARA=CPARA-DC
      IF(FUNCT.GT.0.)  CPARA=CPARA+DC
      IF(I.GE.20)  GO TO 20
      I=I+1
      GO TO 10
C
   20 RETURN
      END
      SUBROUTINE  CASE3(VD,W,SS,HF, HD ,ISPEC)
C
      IMPLICIT  REAL*8  (A-H,O-Z)
C
      ISPEC=1
C
      IF(HF.LE.0.001)  THEN
      ISPEC=2
      HD=SS-VD
      IF(HD.LE.0.0)  ISPEC=5
      RETURN
      ENDIF
C
      CPARA=HF/W
      SI=DSQRT( VD**2+2.*VD*CPARA)
      IF(SI.GT.SS)  GO TO 10
C                                
      ISPEC=2
      HD=CPARA*ARHYCS(VD/CPARA+1.)+(SS-SI)
      RETURN
C
   10 ISPEC=1
      HD=CPARA*ARHYCS((SS**2-VD**2)/2./CPARA**2)
C
      RETURN
      END
      SUBROUTINE CATENA(SS,SI,T,W,VD,HD,ALP,XI,YI,ZI,RK,ISPEC,BANE)
C
      IMPLICIT  REAL*8  (A-H,O-Z)
C
      COMMON / C00 / KR,KW,KM,IPAGE,IERR,NEP(20),NITLE(20),LTITL(20)
C
      DIMENSION  RK(6,6),V(3,3),BANE(4)
C
      IF(ISPEC.LE.3)  GO TO 40
      CSS=100000.
      CHS=100000.
      CHH=100000.
      GO TO 80
C
   40 IF(SI.LE.SS)  GO TO 50
      IF(SI.GT.SS)  GO TO 60
   50 XC=T*ARHYCS(VD/T+1.0)
      DS=0.0
                       GO TO 70
C  60 XC=T*ARHYCS((SS**2-VD**2)/2./T**2)
C  60 XC=HD
   60 AAA=(SS**2-VD**2)/2./T**2 + 1.0
      IF(AAA.GT.1.0)  XC=T*ARHYCS(AAA)
      IF(NEP(11).GE.1)  WRITE(KW,65) SS,SI,T,AAA,XC,HD,VD
      IF(AAA.LE.1.0)  XC=HD
   65 FORMAT(1H0,'(CATENA)   (SS,SI,T,AAA,XC,HD,VD)',7F10.3)
      XX=VD/(2.*T*DSINH(XC/2./T))
      DS=T*ARHYSN(XX)-XC/2.
   70 CONTINUE
C
      D=XC/T
      X=DS/T
      THITA0=  DATAN2(DSINH(X),1.0D0)*180./3.1415927
      CSS=W*DSINH(D)
      CSH=W*(DCOSH(X)*(DCOSH(D)-1.0)+DSINH(X)*DSINH(D))
      CHH=W*(D*DCOSH(X)*DCOSH(X+D)-DSINH(D))
      Q=D*DSINH(D)-2.*(DCOSH(D)-1.0)
      CSS=CSS/Q
      CSH=CSH/Q
      CHH=CHH/Q
   80 CHS=CSH
C
      BANE(1)=CSS
      BANE(2)=CHH
      BANE(3)=CSH
      BANE(4)=DSQRT(CSS**2+CHH**2)
C
      COSA=DCOS(ALP)
      SINA=DSIN(ALP)
C
      RK(1,1)=CSS*COSA**2
      RK(1,2)=-CSS*SINA*COSA
      RK(1,3)= -CSH*COSA
      RK(1,4)= RK(1,3)*YI-RK(1,2)*ZI
      RK(1,5)=-RK(1,3)*XI+RK(1,1)*ZI
      RK(1,6)=-RK(1,1)*YI+RK(1,2)*XI
      RK(2,1)=RK(1,2)
      RK(2,2)=CSS*SINA**2
      RK(2,3)=CSH*SINA
      RK(2,4)= RK(2,3)*YI-RK(2,2)*ZI
      RK(2,5)=-RK(2,3)*XI+RK(2,1)*ZI
      RK(2,6)=-RK(2,1)*YI+RK(2,2)*XI
      RK(3,1)=RK(1,3)
      RK(3,2)=RK(2,3)
      RK(3,3)=CHH
      RK(3,4)= RK(3,3)*YI-RK(3,2)*ZI
      RK(3,5)=-RK(3,3)*XI+RK(3,1)*ZI
      RK(3,6)=-RK(3,1)*YI+RK(3,2)*XI
C
      V(1,1)= 0.0
      V(1,2)= YI
      V(1,3)= ZI
      V(2,1)= XI
      V(2,2)= 0.0
      V(2,3)= ZI
      V(3,1)= XI
      V(3,2)= YI
      V(3,3)=0.0
C
       DO 72 I=1,3
      II=I+3
      I1=I+1
       IF(I1.GT.3) I1=I1-3
      I2=I1+1
       IF(I2.GT.3) I2=I2-3
       DO 73 J=1,6
      RK(II,J)=-RK(I1,J)*V(I,I2)+RK(I2,J)*V(I,I1)
   73 CONTINUE
   72 CONTINUE
      RETURN
      END
      SUBROUTINE  TRANSF(X1,Y1,Z1,X2,Y2,Z2,  AA,V,W)
C
      IMPLICIT  REAL*8  (A-H,O-Z)
      DIMENSION       AA(3,3),V(3,3),W(3,3)
C
       AX=X2-X1
       AY=Y2-Y1
       AZ=Z2-Z1
      BXY=DSQRT(AX**2 + AY**2)
C
      IF(BXY.GT.0.00001)  GO TO 10
C
      DO 3 I=1,3
      DO 3 J=1,3
    3 AA(I,J)=0.0
C
C     AA(3,3)=-AZ/DABS(AZ)    
      AA(1,3)=-AZ/DABS(AZ)
C
      GO TO 20
C
   10 PHY=DATAN2(AY,AX)
      THE=-DATAN2(AZ,BXY)
      CP=DCOS(PHY)
      CT=DCOS(THE)
      SP=DSIN(PHY)
      ST=DSIN(THE)
      AA(1,1)=CT*CP
      AA(1,2)=CT*SP
      AA(1,3)= -ST
      AA(2,1)= -SP
      AA(2,2)=  CP
      AA(2,3)= 0.0
      AA(3,1)=ST*CP
      AA(3,2)=ST*SP
      AA(3,3)= CT
C
   20 V(1,1)= 0.0
      V(1,2)= Z1
      V(1,3)=-Y1
      V(2,1)=-Z1
      V(2,2)= 0.0
      V(2,3)=X1
      V(3,1)=Y1
      V(3,2)=-X1
      V(3,3)=0.0
C
       DO 72 I=1,3
        DO 73 J=1,3
         CDM=0.0
          DO 74  K=1,3
          CDM=CDM+AA(I,K)*V(K,J)
   74      CONTINUE
          W(I,J)=CDM
   73   CONTINUE
   72  CONTINUE
      RETURN
      END
      SUBROUTINE  FINTOT(Npj,FI)
C
      IMPLICIT  REAL*8  (A-H,O-Z)
C
      COMMON / C00 / KR,KW,KM,IPAGE,IERR,NEP(20),NITLE(20),LTITL(20)
      COMMON / C01 / NBANE,NPONT,Npfend(50),Npfs(50),Npfe(50),Npcain(50)
     *              ,Npcs(50),Npce(50),Pitle(50),Xpi(50),Ypi(50),Zpi(50)
      COMMON / B05 / DLFN0(999),DLFN(999),FFN(999),RKFN(999),TFN(999,6)
      COMMON /B06/DLCN(999),FCN(999),FCNH(999),FICN(6,999),RKCN(6,6,999)
     *              ,TCN(999,6),FICN0(6,999),FCNMIN(999),FCNMAX(999)
      DIMENSION          FI(6)
C
      DO 10 I=1,6
      FI(I)=0.0
   10 CONTINUE
C
      IF(Npfend(Npj).EQ.0) GO TO 1
C
C     TOTAL FORCE ABOUT GRAVITY CENTER BY FENDER
C
      DO 20 J=Npfs(Npj),Npfe(Npj)
      DO 30 I=1,6
      FI(I)=FI(I)+FFN(J)*TFN(J,I )
   30 CONTINUE
   20 CONTINUE
C
    1 IF(Npcain(Npj).EQ.0) GO TO 2
C
C     TOTAL FORCE ABOUT GRAVITY CENTER BY CHAIN
C
      DO 40 J=Npcs(Npj),Npce(Npj)
      DO 50 I=1,6
      FI(I)=FI(I)+FICN(I,J)
   50 CONTINUE
   40 CONTINUE
C
C     TOTAL FORCE ABOUT GRAVITY CENTER BY MLINE
C
    2 RETURN                         
C
      END
      SUBROUTINE  ASSEMP(Npj,RK2)
C
      IMPLICIT  REAL*8  (A-H,O-Z)
C
      COMMON / C00 / KR,KW,KM,IPAGE,IERR,NEP(20),NITLE(20),LTITL(20)
      COMMON / C01 / NBANE,NPONT,Npfend(50),Npfs(50),Npfe(50),Npcain(50)
     *              ,Npcs(50),Npce(50),Pitle(50),Xpi(50),Ypi(50),Zpi(50)
      COMMON / B05 / DLFN0(999),DLFN(999),FFN(999),RKFN(999),TFN(999,6)
      COMMON /B06/DLCN(999),FCN(999),FCNH(999),FICN(6,999),RKCN(6,6,999)
     *              ,TCN(999,6),FICN0(6,999),FCNMIN(999),FCNMAX(999)
      DIMENSION RK2(6,6)
C
      DO  10  I=1,6
      DO  20  J=1,6
      RK2(I,J)=0.0
   20 CONTINUE
   10 CONTINUE
C
C     MATRIX BY FENDER
C
      IF(Npfend(Npj).EQ.0) GO TO 1
C
      DO 11 K=Npfs(Npj),Npfe(Npj)
      DO 12 I=1,6
      DO 13 J=1,6
      RK2(I,J)=RK2(I,J)+RKFN(K)*TFN(K,I)*TFN(K,J)
   13 CONTINUE
   12 CONTINUE
   11 CONTINUE
C
C     MATRIX BY CHAIN
C
    1 IF(Npcain(Npj).EQ.0) GO  TO  50
C
      DO 21 K=Npcs(Npj),Npce(Npj)
      DO 22 I=1,6
      DO 23 J=1,6
      RK2(I,J)=RK2(I,J)+RKCN(I,J,K)
   23 CONTINUE
   22 CONTINUE
   21 CONTINUE
C
   50 IF(NEP(14).LE.0)  RETURN
      WRITE(KW,900)     
  900 FORMAT(1H0,'SUB.ASSEM ')
      CALL  MATOUT('(K2)',RK2,KW)
      RETURN
      END
      SUBROUTINE   WTITLE
C
      COMMON / C00 / KR,KW,KM,IPAGE,IERR,NEP(20),NITLE(20),LTITL(20)
      WRITE(KW,100)    NITLE, IPAGE
       IPAGE=IPAGE+1
  100 FORMAT(1H1,20A4,30X,'PAGE',I3)
      RETURN
      END
      SUBROUTINE  MATOUT (NAME, A, KW)
C
      IMPLICIT  REAL*8  (A-H,O-Z)
      character*(4) name
      DIMENSION  A(6,6)
C
      WRITE(KW,100)    NAME
      WRITE(KW,101)   (A(1,J),J=1,6)
      WRITE(KW,102)   (A(2,J),J=1,6)
      WRITE(KW,103)   (A(3,J),J=1,6)
      WRITE(KW,104)   (A(4,J),J=1,6)
      WRITE(KW,105)   (A(5,J),J=1,6)
      WRITE(KW,106)   (A(6,J),J=1,6)
      RETURN
  100 FORMAT(1H0,10X,'**  MATRIX  NAME = ',A4,'  **'/1H0,20X,'SURGE(X)
     *     SWAY (Y)       HEAVE(Z)       ROLL (FAI)     PITCH(THI)     Y
     *AW (PSAI)' )
  101 FORMAT(1H ,5X,'SURGE(X)  ',6F15.2)
  102 FORMAT(1H ,5X,'SWAY (Y)  ',6F15.2)
  103 FORMAT(1H ,5X,'HEAVE(Z)  ',6F15.2)
  104 FORMAT(1H ,5X,'ROLL (FAI)',6F15.2)
  105 FORMAT(1H ,5X,'PITCH(THI)',6F15.2)
  106 FORMAT(1H ,5X,'YAW (PSAI)',6F15.2)
      END
      SUBROUTINE   LISTUP(L,N)
C     ****************************************************************  DATA015
C     **  THIS IS THE 'INPUT DATA PRINT OUT ROUTINE' FOR ANY JOB.   **  DATA020
C     **            INPUT -L --- FT05F001   CARD                    **  DATA025
C     **            OUTPUT-N --- FT01F001   PRINT                   **  DATA030
C     ****************************************************************  DATA040
C                                                 CALLED BY  INPUT      DATA045
      DIMENSION ATYPE(20)
              IPAGE=0                                                   DATA125
              IDATA=0                                                   DATA130
              LINE =0                                                   DATA135
C
10000 CONTINUE                                                          DATA140
              IPAGE=IPAGE+1                                             DATA145
          WRITE (N,101)  IPAGE
      IF(IPAGE.EQ.1)  WRITE(N,100)  L
  100 FORMAT(1H ,30X,'FILE NO.=',I3)
  101     FORMAT(1H1,T90,'PAGE',I6)                                     DATA155
          ASSIGN 1000  TO M                                             DATA160
          GO TO 20000                                                   DATA165
 1000 CONTINUE                                                          DATA170
          READ  (L,102,END=99999)  ATYPE                                DATA175
  102     FORMAT(20A4)                                                  DATA177
              LINE =LINE +1                                             DATA180
              IDATA=IDATA+1                                             DATA185
      WRITE(N,103)  ATYPE,IDATA
  103     FORMAT(1H ,20A4,T86,'DATA NO.',I6)                            DATA195
          IF(LINE.LT.50) GO TO  1000                                    DATA210
              LINE=0                                                    DATA215
          ASSIGN 10000  TO M                                            DATA220
20000 CONTINUE                                                          DATA225
          WRITE  (N,105)
  105     FORMAT ( 1H0,'*---.----1----.----2----.----3----.----4',
     *                 '----.----5----.----6----.----7----.----8' /)
          GO TO M,(1000,10000,90000)                                    DATA245
99999 CONTINUE                                                          DATA250
          ASSIGN 90000  TO M                                            DATA255
          GO TO 20000                                                   DATA257
90000 CONTINUE                                                          DATA258
          WRITE  (N,109)  IDATA
  109     FORMAT (1H ,T86,'TOTAL   ',I6)                                DATA265
          REWIND  L
C        --------                                                       DATA278
          RETURN                                                        DATA280
                              END
      FUNCTION ARHYSN(X)
      REAL*8 ARHYSN,X
      ARHYSN=DLOG(X+DSQRT(X**2+1.0))
      RETURN
      END
      FUNCTION ARHYCS(X)
      REAL*8 ARHYCS,X
      IF(X.GE.1.0) ARHYCS=DLOG(X+DSQRT(X**2-1.0))
      IF(X.LT.1.0) STOP 'ERROR IN FUNCTION ARHYCS X.LT.1.0'
      RETURN
      END
