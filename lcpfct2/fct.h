c     /FCT_SCRH/ Holds scratch arrays for use by LCPFCT and CNVFCT
          Real     SCRH(0:NPT),     SCR1(0:NPT),     DIFF(0:NPT)
          Real     FLXH(0:NPT),     FABS(0:NPT),     FSGN(0:NPT)
          Real     TERM(0:NPT),     TERP(0:NPT),     LNRHOT(0:NPT)
          Real     LORHOT(0:NPT),   RHOT(0:NPT),     RHOTD(0:NPT)
      COMMON /FCT_SCRH/SCRH,SCR1,DIFF,FLXH,FABS,FSGN,
     &                       TERM, TERP, LNRHOT, LORHOT, RHOT, RHOTD
c$OMP THREADPRIVATE(/FCT_SCRH/)

c     /FCT_GRID/ Holds geometry, grid, area and volume information
          Real     LO(0:NPT),       LN(0:NPT),       AH (0:NPT) 
          Real     RLN(0:NPT),      LH (0:NPT),      RLH(0:NPT) 
          Real     ROH(0:NPT),      RNH(0:NPT),      ADUGTH(0:NPT)
      COMMON /FCT_GRID/LO,LN,AH,RLN,LH,RLH,ROH,RNH,ADUGTH
c$OMP THREADPRIVATE(/FCT_GRID/)

c     /FCT_VELO/ Holds velocity-dependent flux coefficients
          Real     HADUDTH(0:NPT),  NULH(0:NPT),     MULH(0:NPT)
          Real     EPSH(0:NPT),     VDTODR(0:NPT)   
      COMMON /FCT_VELO/HADUDTH,NULH,MULH,EPSH,VDTODR
c$OMP THREADPRIVATE( /FCT_VELO/)

c     /FCT_MISC/ Holds the source array and diffusion coefficient
          Real     SOURCE(0:NPT),   DIFF1
       COMMON /FCT_MISC/SOURCE,DIFF1
c$OMP THREADPRIVATE( /FCT_MISC/)

c     /FCT_NDEX/ Holds a scalar list of special cell information . . .
          Real     SCALARS(NINDMAX)
          Integer  INDEX(NINDMAX), NIND
c      COMMON /FCT_NDEX/NIND,INDEX,SCALARS
       COMMON /FCT_NDEX/scalars,NIND,INDEX
c$OMP THREADPRIVATE( /FCT_NDEX/)
