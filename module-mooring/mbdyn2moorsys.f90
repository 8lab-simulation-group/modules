subroutine input_moorsys(ifname, ifnamelen, ofname, ofnamelen, *)
    IMPLICIT    NONE

    INTEGER(4)                      ifnamelen, ofnamelen
    CHARACTER(ifnamelen)            ifname
    CHARACTER(ofnamelen)            ofname

    CHARACTER(1024)                 Iptfile
    CHARACTER(1024)                 Outfile  

    REAL*8                          Xpj(6), Rpj(6), RKpj(6,6), DTIME

    REAL*8                          DLCN, FCN, FCNH
    INTEGER                         KR, KW, KM, IPAGE, IERR, NEP, NITLE, LTITL, Iniflg, IHIS, Npj
    INTEGER                         NBANE,NPONT,Npfend,Npfs,Npfe,Npcain

    COMMON / C00 / KR,KW,KM,IPAGE,IERR,NEP(20),NITLE(20),LTITL(20)
    COMMON / C01 / NBANE,NPONT,Npfend(50),Npfs(50),Npfe(50),Npcain(50)
    COMMON / B06 / DLCN(999),FCN(999),FCNH(999)

    KR     = 60
    KW     = 70
    IPAGE  = 1
    Iniflg = 0
    IHIS   = 0

    Iptfile = ifname
    Outfile = ofname

    OPEN(KR,file=Iptfile)
    OPEN(KW,file=Outfile)

    DTIME = 0.0
    CALL MOORSYS ( Iniflg, IHIS, DTIME, Npj, Xpj, Rpj, RKpj ) 

    print *, 'input completed'

end subroutine input_moorsys

subroutine call_moorsys(DTIME, Xpj, Rpj)
    IMPLICIT    NONE

    REAL*8                          Xpj(6), Rpj(6), RKpj(6,6), DTIME

    REAL*8                          DLCN, FCN, FCNH
    INTEGER                         KR, KW, KM, IPAGE, IERR, NEP, NITLE, LTITL, Iniflg, IHIS, Npj
    INTEGER                         NBANE,NPONT,Npfend,Npfs,Npfe,Npcain

    COMMON / C00 / KR,KW,KM,IPAGE,IERR,NEP(20),NITLE(20),LTITL(20)
    COMMON / C01 / NBANE,NPONT,Npfend(50),Npfs(50),Npfe(50),Npcain(50)
    COMMON / B06 / DLCN(999),FCN(999),FCNH(999)

    Iniflg = 1
    IHIS   = 0
    Npj = 1
    Rpj = 0.0
    CALL MOORSYS ( Iniflg, IHIS, DTIME, Npj, Xpj, Rpj, RKpj ) 

end subroutine call_moorsys
