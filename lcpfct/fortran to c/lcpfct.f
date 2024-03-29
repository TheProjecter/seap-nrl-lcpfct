C=======================================================================
Cc$mta inline
      Subroutine LCPFCT ( RHOO, RHON, I1, IN,
     &                    SRHO1, VRHO1, SRHON, VRHON, PBC )

C-----------------------------------------------------------------------
c
c     Originated: J.P. Boris         Code 4400, NRL          Feb 1987
c     Modified:  Laboratory for Computational Physics & Fluid Dynamics
c     Contact:    J.P. Boris, J.H. Gardner, A.M. Landsberg, or E.S. Oran
c
c     Description:  This routine solves generalized continuity equations
c     of the form  dRHO/dt = -div (RHO*V) + SOURCES in the user's choice
c     of Cartesian, cylindrical, or spherical coordinate systems.  A
c     facility is included to allow definition of other coordinates.
c     The grid can be Eulerian, sliding rezone, or Lagrangian and can
c     be arbitrarily spaced.  The algorithm is a low-phase-error FCT
c     algorithm, vectorized and optimized for a combination of speed and
c     flexibility.  A complete description appears in the NRL Memorandum
c     Report (1992), "LCPFCT - A Flux-Corrected Transport Algorithm For
c     Solving Generalized Continuity Equations".
c
c     Arguments:
c     RHOO   Real Array        grid point densities at start of step   I
c     RHON   Real Array        grid point densities at end of step     O
c     I1     Integer           first grid point of integration         I
c     IN     Integer           last grid point of intergration         I
c     SRHO1  Real Array        boundary guard cell factor at cell I1+1 I
c     VRHO1  Real Array        boundary value added to guard cell I1-1 I
c     SRHON  Real Array        boundary guard cell factor at cell IN+1 I
c     VRHON  Real Array        boundary value added to guard cell IN+1 I
c     PBC    Logical           periodic boundaries if PBC = .true.     I
c
c        In this routine the last interface at RADHN(INP) is the outer
c     boundary of the last cell indexed IN.  The first interface at
c     RADHN(I1) is the outer boundary of the integration domain before
c     the first cell indexed I1.
c
c     Language and Limitations:  LCPFCT is a package of FORTRAN 77 sub-
c     routines written in single precision (64 bits CRAY). The parameter
c     NPT is used to establish the internal FCT array dimensions at the
c     maximum size expected.  Thus NPT = 2002 means that continuity equa-
c     tions for systems up to 200 cells long in one direction can be
c     integrated.  Underflows can occur when the function being trans-
c     ported has a region of zeroes.  The calculations misconserve by
c     one or two bits per cycle.  Relative phase and amplitude errors
c     (for smooth functions) are typically a few percent for character-
c     istic lengths of 1 - 2 cells (wavelengths of order 10 cells).  The
c     jump conditions for shocks are generally accurate to better than 1
c     percent.  Common blocks are used to transmit all data between the
c     subroutines in the LCPFCT package.
c
c     Auxiliary Subroutines:  CNVFCT, CONSERVE, COPYGRID, MAKEGRID,
c     NEW_GRID, RESIDIFF, SET_GRID, SOURCES, VELOCITY, ZERODIFF, and
c     ZEROFLUX.  The detailed documentation report provided (or the
c     listing below) explains the definitions and use of the arguments
c     to these other subroutines making up the LCPFCT package.  These
c     routines are not called from LCPFCT itself but are controlled by
c     calls from the user.  Subroutines MAKEGRID, VELOCITY and SOURCES
c     in this package must first be called to set the grid geometry,
c     velocity-dependent flux and diffusion coefficients, and external
c     source arrays used by LCPFCT.  The other subroutines may be called
c     to perform other functions such as to modify boundary conditions,
c     to perform special grid operations, or compute conservation sums.
c
C-----------------------------------------------------------------------

          Implicit  NONE
          Integer   NPT, I1, IN, I1P, INP, I
          Real      BIGNUM, SRHO1, VRHO1, SRHON, VRHON, RHO1M, RHONP
          Real      RHOT1M, RHOTNP, RHOTD1M, RHOTDNP
          Logical   PBC
          Parameter ( NPT = 2002 )
          Parameter ( BIGNUM = 1.0E38 )
c     BIGNUM = Machine Dependent Largest Number - Set By The User!!!!

          Real     RHOO(NPT),     RHON(NPT)

c     /FCT_SCRH/ Holds scratch arrays for use by LCPFCT and CNVFCT
          Real     SCRH(NPT),     SCR1(NPT),     DIFF(NPT)
          Real     FLXH(NPT),     FABS(NPT),     FSGN(NPT)
          Real     TERM(NPT),     TERP(NPT),     LNRHOT(NPT)
          Real     LORHOT(NPT),   RHOT(NPT),     RHOTD(NPT)
      COMMON /FCT_SCRH/SCRH,SCR1,DIFF,FLXH,FABS,FSGN,
     &                       TERM, TERP, LNRHOT, LORHOT, RHOT, RHOTD
c$OMP THREADPRIVATE(/FCT_SCRH/)

c     /FCT_GRID/ Holds geometry, grid, area and volume information
          Real     LO(NPT),       LN(NPT),       AH (NPT)
          Real     RLN(NPT),      LH (NPT),      RLH(NPT)
          Real     ROH(NPT),      RNH(NPT),      ADUGTH(NPT)
      COMMON /FCT_GRID/LO,LN,AH,RLN,LH,RLH,ROH,RNH,ADUGTH
c$OMP THREADPRIVATE(/FCT_GRID/)

c     /FCT_VELO/ Holds velocity-dependent flux coefficients
          Real     HADUDTH(NPT),  NULH(NPT),     MULH(NPT)
          Real     EPSH(NPT),     VDTODR(NPT)
      COMMON /FCT_VELO/HADUDTH,NULH,MULH,EPSH,VDTODR
c$OMP THREADPRIVATE( /FCT_VELO/)

c     /FCT_MISC/ Holds the source array and diffusion coefficient
          Real     SOURCE(NPT),   DIFF1
       COMMON /FCT_MISC/SOURCE,DIFF1
c$OMP THREADPRIVATE( /FCT_MISC/)

C-----------------------------------------------------------------------
          I1P = I1 + 1
          INP = IN + 1

c     Calculate the convective and diffusive fluxes . . .
C-----------------------------------------------------------------------
          If ( PBC ) Then
             RHO1M = RHOO(IN)
             RHONP = RHOO(I1)
          Else
             RHO1M = SRHO1*RHOO(I1) + VRHO1
             RHONP = SRHON*RHOO(IN) + VRHON
          End If

          DIFF(I1) = NULH(I1) * ( RHOO(I1) - RHO1M )
          FLXH(I1) = HADUDTH(I1) * ( RHOO(I1) + RHO1M )

          Do 1 I = I1P, IN
             FLXH(I) = HADUDTH(I) * ( RHOO(I) + RHOO(I-1) )
 1                   DIFF(I) = NULH(I) * ( RHOO(I) - RHOO(I-1) )

          DIFF(INP) = NULH(INP) * ( RHONP - RHOO(IN) )
          FLXH(INP) = HADUDTH(INP) * ( RHONP + RHOO(IN) )

c     Calculate LORHOT, the transported mass elements, and LNRHOT, the
c     transported & diffused mass elements  . . .
C-----------------------------------------------------------------------
          Do 2 I = I1, IN
             LORHOT(I) = LO(I)*RHOO(I) + SOURCE(I) + (FLXH(I)-FLXH(I+1))
             LNRHOT(I) = LORHOT(I) + (DIFF(I+1) - DIFF(I))
             RHOT(I)  = LORHOT(I)*RLN(I)
 2                   RHOTD(I) = LNRHOT(I)*RLN(I)

c     Evaluate the boundary conditions for RHOT and RHOTD . . .
C-----------------------------------------------------------------------
          If ( PBC ) Then
             RHOT1M = RHOT(IN)
             RHOTNP = RHOT(I1)
             RHOTD1M = RHOTD(IN)
             RHOTDNP = RHOTD(I1)
          Else
             RHOT1M = SRHO1*RHOT(I1) + VRHO1
             RHOTNP = SRHON*RHOT(IN) + VRHON
             RHOTD1M = SRHO1*RHOTD(I1) + VRHO1
             RHOTDNP = SRHON*RHOTD(IN) + VRHON
          End If

c     Calculate the transported antiduffusive fluxes and transported
c     and diffused density differences . . .
C-----------------------------------------------------------------------
          FLXH(I1) = MULH(I1) * ( RHOT(I1) - RHOT1M )
          DIFF(I1) = RHOTD(I1) - RHOTD1M
          FABS(I1) = ABS ( FLXH(I1) )
          FSGN(I1) = SIGN ( DIFF1, DIFF(I1) )

          Do 3 I = I1P, IN
             FLXH(I) = MULH(I) * ( RHOT(I) - RHOT(I-1) )
 3                   DIFF(I) = RHOTD(I) - RHOTD(I-1)

          FLXH(INP) = MULH(INP) * ( RHOTNP - RHOT(IN) )
          DIFF(INP) = RHOTDNP - RHOTD(IN)

c     Calculate the magnitude & sign of the antidiffusive flux followed
c     by the flux-limiting changes on the right and left . . .
C-----------------------------------------------------------------------
          Do 4 I = I1, IN
             FABS(I+1) = ABS ( FLXH(I+1) )
             FSGN(I+1) = SIGN ( DIFF1, DIFF(I+1) )
             TERM(I+1) = FSGN(I+1)*LN(I)*DIFF(I)
 4                   TERP(I) = FSGN(I)*LN(I)*DIFF(I+1)

          If ( PBC ) Then
             TERP(INP) = TERP(I1)
             TERM(I1) = TERM(INP)
          Else
             TERP(INP) = BIGNUM
             TERM(I1) = BIGNUM
          End If

c     Correct the transported fluxes completely and then calculate the
c     new Flux-Corrected Transport densities . . .
C-----------------------------------------------------------------------
          FLXH(I1) = FSGN(I1) * AMAX1 ( 0.0,
     &                  AMIN1 ( TERM(I1), FABS(I1), TERP(I1) ) )

          Do 5 I = I1, IN
             FLXH(I+1) = FSGN(I+1) * AMAX1 ( 0.0,
     &                AMIN1 ( TERM(I+1), FABS(I+1), TERP(I+1) ) )
             RHON(I) = RLN(I) * ( LNRHOT(I) + (FLXH(I) - FLXH(I+1)) )
 5                   SOURCE(I) = 0.0

      Return
      End

C=======================================================================

Cc$mta inline
      Subroutine MAKEGRID ( RADHO, RADHN, I1, INP, ALPHA )

C-----------------------------------------------------------------------
c
c     Description:  This Subroutine initializes geometry variables and
c     coefficients. It should be called first to initialize the grid.
c     The grid must be defined for all of the grid interfaces from I1 to
c     INP.  Subsequent calls to VELOCITY and LCPFCT can work on only
c     portions of the grid, however, to perform restricted integrations
c     on separate line segments.
c
c     Arguments:
c     RADHO    Real Array(INP)    old cell interface positions         I
c     RADHN    Real Array(INP)    new cell interface positions         I
c     I1       Integer            first cell interface                 I
c     INP      Integer            last cell interface                  I
c     ALPHA    Integer            = 1 for cartesian geometry           I
c                                 = 2 for cylindrical geometry         I
c                                 = 3 for spherical geometry           I
c                                 = 4 general geometry (user supplied) I
c
C-----------------------------------------------------------------------

          Implicit  NONE
          Integer   NPT, I1, I1P, I, IN, INP, ALPHA
          Parameter ( NPT = 2002 )

          Real     RADHO(INP), RADHN(INP), PI, FTPI

c     /FCT_SCRH/ Holds scratch arrays for use by LCPFCT and CNVFCT
          Real     SCRH(NPT),     SCR1(NPT),     DIFF(NPT)
          Real     FLXH(NPT),     FABS(NPT),     FSGN(NPT)
          Real     TERM(NPT),     TERP(NPT),     LNRHOT(NPT)
          Real     LORHOT(NPT),   RHOT(NPT),     RHOTD(NPT)
       COMMON /FCT_SCRH/SCRH,SCR1,DIFF,FLXH,FABS,FSGN,
     &                       TERM, TERP, LNRHOT, LORHOT, RHOT, RHOTD
c$OMP THREADPRIVATE(/FCT_SCRH/)

c     /FCT_GRID/ Holds geometry, grid, area and volume information
          Real     LO(NPT),       LN(NPT),       AH (NPT)
          Real     RLN(NPT),      LH (NPT),      RLH(NPT)
          Real     ROH(NPT),      RNH(NPT),      ADUGTH(NPT)
       COMMON /FCT_GRID/LO,LN,AH,RLN,LH,RLH,ROH,RNH,ADUGTH
c$OMP THREADPRIVATE( /FCT_GRID/)

          DATA     PI, FTPI /3.1415927, 4.1887902/

C-----------------------------------------------------------------------
          I1P = I1 + 1
          IN = INP - 1

c  Store the old and new grid interface locations from input and then
c  update the new and average interface and grid coefficients . . .
C-----------------------------------------------------------------------
          Do 1 I = I1, INP
             ROH(I) = RADHO(I)
 1                   RNH(I) = RADHN(I)

c  Select the choice of coordinate systems . . .
C-----------------------------------------------------------------------
          Go To (100, 200, 300, 400), ALPHA

c  Cartesian coordinates . . .
C-----------------------------------------------------------------------
 100            AH(INP) = 1.0
          Do 101 I = I1, IN
             AH(I) = 1.0
             LO(I) = ROH(I+1) - ROH(I)
 101                 LN(I) = RNH(I+1) - RNH(I)
          Go To 500

c  Cylindrical Coordinates: RADIAL . . .
C-----------------------------------------------------------------------
 200           DIFF(I1) = RNH(I1)*RNH(I1)
          SCRH(I1) = ROH(I1)*ROH(I1)
          AH(INP) = PI*(ROH(INP) + RNH(INP))
          DO 201 I = I1, IN
             AH(I) = PI*(ROH(I) + RNH(I))
             SCRH(I+1) = ROH(I+1)*ROH(I+1)
             LO(I) = PI*(SCRH(I+1) - SCRH(I))
             DIFF(I+1) = RNH(I+1)*RNH(I+1)
 201                 LN(I) = PI*(DIFF(I+1) - DIFF(I))
          Go To 500

c  Spherical Coordinates: RADIAL . . .
C-----------------------------------------------------------------------
 300           SCR1(I1) = ROH(I1)*ROH(I1)*ROH(I1)
          DIFF(I1) = RNH(I1)*RNH(I1)*RNH(I1)
          SCRH(INP) = (ROH(INP) + RNH(INP))*ROH(INP)
          AH(INP) = FTPI*(SCRH(INP) + RNH(INP)*RNH(INP))
          DO 301 I = I1, IN
             SCR1(I+1) = ROH(I+1)*ROH(I+1)*ROH(I+1)
             DIFF(I+1) = RNH(I+1)*RNH(I+1)*RNH(I+1)
             SCRH(I) = (ROH(I) + RNH(I))*ROH(I)
             AH(I) = FTPI*(SCRH(I) + RNH(I)*RNH(I))
             LO(I) = FTPI*(SCR1(I+1) - SCR1(I))
 301                 LN(I) = FTPI*(DIFF(I+1) - DIFF(I))
          Go To 500

c  Special Coordinates: Areas and Volumes are User Supplied . . .
C-----------------------------------------------------------------------
 400           Continue

c  Additional system independent geometric variables . . .
C-----------------------------------------------------------------------
 500                Do 501 I = I1, IN
 501                           RLN(I) = 1.0/LN(I)
          LH(I1)  = LN(I1)
          RLH(I1) = RLN(I1)
          Do 502 I = I1P, IN
             LH(I) =  0.5*(LN(I) + LN(I-1))
 502                 RLH(I) = 0.5*(RLN(I) + RLN(I-1))
          LH(INP)  = LN(IN)
          RLH(INP) = RLN(IN)
          Do 503 I = I1, INP
 503                 ADUGTH(I) = AH(I)*(RNH(I) - ROH(I))

      Return
      End

C=======================================================================

Cc$mta inline
      Subroutine VELOCITY ( UH, I1, INP, DT )

C-----------------------------------------------------------------------
c
c     Description:   This subroutine calculates all velocity-dependent
c     coefficients for the LCPFCT and CNVFCT routines. This routine
c     must be called before either LCPFCT or CNVFCT is called.  MAKEGRID
c     must be called earlier to set grid and geometry data used here.
c
c     Arguments:
c     UH     Real Array(NPT)   flow velocity at cell interfaces        I
c     I1     Integer           first cell interface of integration     I
c     INP    Integer           last cell interface = N + 1             I
c     DT     Real              stepsize for the time integration       I
c
C-----------------------------------------------------------------------

          Implicit NONE
          Integer  NPT, I1, I1P, I, IN, INP
          Parameter ( NPT = 2002 )

          Real     UH(INP), DT, RDT, DTH, DT2, DT4, ONE3RD, ONE6TH

c     /FCT_SCRH/ Holds scratch arrays for use by LCPFCT and CNVFCT
          Real     SCRH(NPT),     SCR1(NPT),     DIFF(NPT)
          Real     FLXH(NPT),     FABS(NPT),     FSGN(NPT)
          Real     TERM(NPT),     TERP(NPT),     LNRHOT(NPT)
          Real     LORHOT(NPT),   RHOT(NPT),     RHOTD(NPT)
       COMMON /FCT_SCRH/SCRH,SCR1,DIFF,FLXH,FABS,FSGN,
     &                       TERM, TERP, LNRHOT, LORHOT, RHOT, RHOTD
c$OMP THREADPRIVATE( /FCT_SCRH/)

c     /FCT_GRID/ Holds geometry, grid, area and volume information
          Real     LO(NPT),       LN(NPT),       AH (NPT)
          Real     RLN(NPT),      LH (NPT),      RLH(NPT)
          Real     ROH(NPT),      RNH(NPT),      ADUGTH(NPT)
       COMMON /FCT_GRID/LO,LN,AH,RLN,LH,RLH,ROH,RNH,ADUGTH
c$OMP THREADPRIVATE( /FCT_GRID/)

c     /FCT_VELO/ Holds velocity-dependent flux coefficients
          Real     HADUDTH(NPT),  NULH(NPT),     MULH(NPT)
          Real     EPSH(NPT),     VDTODR(NPT)
       COMMON /FCT_VELO/HADUDTH,NULH,MULH,EPSH,VDTODR
c$OMP THREADPRIVATE(/FCT_VELO/)

C-----------------------------------------------------------------------
          I1P = I1 + 1
          IN = INP - 1

c     Calculate 0.5*Interface Area * Velocity Difference * DT (HADUDTH).
c     Next calculate the interface epsilon (EPSH = V*DT/DX).  Then find
c     the diffusion (NULH) and antidiffusion (MULH) coefficients.  The
c     variation with epsilon gives fourth-order accurate phases when the
c     grid is uniform, the velocity constant, and SCRH is set to zero.
c     With SCRH nonzero (as below) slightly better results are obtained
c     in some of the tests.  Optimal performance, of course, depends on
c     on the application.
C-----------------------------------------------------------------------
          RDT = 1.0/DT
          DTH = 0.5*DT
          ONE6TH = 1.0/6.0
          ONE3RD = 1.0/3.0
          Do 1 I = I1, INP
             HADUDTH(I) = DT*AH(I)*UH(I) - ADUGTH(I)
             EPSH(I) = HADUDTH(I)*RLH(I)
             SCRH(I) = AMIN1 ( ONE6TH, ABS(EPSH(I)) )
             SCRH(I) = ONE3RD*SCRH(I)**2
             HADUDTH(I) = 0.5*HADUDTH(I)
             NULH(I) =  ONE6TH + ONE3RD*(EPSH(I) + SCRH(I))*
     &                                  (EPSH(I) - SCRH(I))
             MULH(I) =  0.25 - 0.5*NULH(I)
             NULH(I) = LH(I)*(NULH(I) + SCRH(I))
             MULH(I) = LH(I)*(MULH(I) + SCRH(I))
 1                   DIFF(I) = UH(I) - RDT*(RNH(I) - ROH(I))

c     Now calculate VDTODR for CNVFCT . . .
C-----------------------------------------------------------------------
          DT2 = 2.0*DT
          DT4 = 4.0*DT
          VDTODR(I1) = DT2*DIFF(I1)/(RNH(I1P)-RNH(I1) +
     &                               ROH(I1P)-ROH(I1))

          Do 2 I = I1P, IN
 2                   VDTODR(I) = DT4*DIFF(I)/(RNH(I+1)-RNH(I-1) +
     &                                ROH(I+1)-ROH(I-1))

          VDTODR(INP) = DT2*DIFF(INP)/(RNH(INP)-RNH(IN) +
     &                                 ROH(INP)-ROH(IN))

       Return
       End

C=======================================================================

Cc$mta inline
      Subroutine SOURCES ( I1, IN, DT, MODE, C, D, D1, DN)

C-----------------------------------------------------------------------
c
c     Description:   This Subroutine accumulates different source terms.
c
c     Arguments:
c     I1     Integer           first cell to be integrated             I
c     IN     Integer           last cell  to be integrated             I
c     DT     Real              stepsize for the time integration       I
c     MODE   Integer           = 1   computes + DIV (D)                I
c                              = 2   computes + C*GRAD (D)             I
c                              = 3   adds + D to the sources           I
c                              = 4   + DIV (D) from interface data     I
c                              = 5   + C*GRAD (D) from interface data  I
c                              = 6   + C for list of scalar indices    I
c     C      Real Array(NPT)   Array of source variables               I
c     D      Real Array(NPT)   Array of source variables               I
c     D1     Real              first boundary value of D               I
c     DN     Real              last  boundary value of D               I
c
C-----------------------------------------------------------------------

          Implicit NONE
          Integer  NPT, NINDMAX, MODE, IS, I, I1, IN, I1P, INP
          Parameter ( NPT = 2002, NINDMAX = 150 )

          Real     C(NPT), D(NPT), DT, DTH, DTQ, D1, DN

c     /FCT_NDEX/ Holds a scalar list of special cell information . . .
          Real     SCALARS(NINDMAX)
          Integer  INDEX(NINDMAX), NIND
c      COMMON /FCT_NDEX/NIND,INDEX,SCALARS
       COMMON /FCT_NDEX/scalars,NIND,INDEX
c$OMP THREADPRIVATE( /FCT_NDEX/)

c     /FCT_SCRH/ Holds scratch arrays for use by LCPFCT and CNVFCT
          Real     SCRH(NPT),     SCR1(NPT),     DIFF(NPT)
          Real     FLXH(NPT),     FABS(NPT),     FSGN(NPT)
          Real     TERM(NPT),     TERP(NPT),     LNRHOT(NPT)
          Real     LORHOT(NPT),   RHOT(NPT),     RHOTD(NPT)
       COMMON /FCT_SCRH/SCRH,SCR1,DIFF,FLXH,FABS,FSGN,
     &                       TERM, TERP, LNRHOT, LORHOT, RHOT, RHOTD
c$OMP THREADPRIVATE( /FCT_SCRH/)

c     /FCT_GRID/ Holds geometry, grid, area and volume information
          Real     LO(NPT),       LN(NPT),       AH (NPT)
          Real     RLN(NPT),      LH (NPT),      RLH(NPT)
          Real     ROH(NPT),      RNH(NPT),      ADUGTH(NPT)
      COMMON /FCT_GRID/LO,LN,AH,RLN,LH,RLH,ROH,RNH,ADUGTH
c$OMP THREADPRIVATE( /FCT_GRID/)

c     /FCT_MISC/ Holds the source array and diffusion coefficient
          Real     SOURCE(NPT),   DIFF1
      COMMON /FCT_MISC/SOURCE,DIFF1
c$OMP THREADPRIVATE(/FCT_MISC/)

C-----------------------------------------------------------------------
         I1P = I1 + 1
         INP = IN + 1
         DTH = 0.5*DT
         DTQ = 0.25*DT
         Go To  ( 101, 202, 303, 404, 505, 606 ), MODE

c  + DIV(D) is computed conservatively and added to SOURCE . . .
C-----------------------------------------------------------------------
 101         SCRH(I1) = DT*AH(I1)*D1
         SCRH(INP) = DT*AH(INP)*DN
         Do 1 I = IN, I1P, -1
            SCRH(I) = DTH*AH(I)*(D(I) + D(I-1))
 1                 SOURCE(I) = SOURCE(I) + (SCRH(I+1) - SCRH(I))
         SOURCE(I1) = SOURCE(I1) + (SCRH(I1P) - SCRH(I1))
      Return

c  + C*GRAD(D) is computed efficiently and added to the SOURCE . . .
C-----------------------------------------------------------------------
 202      SCRH(I1) = DTH*D1
         SCRH(INP) = DTH*DN
         Do 2 I = IN, I1P, -1
            SCRH(I) = DTQ*(D(I)+D(I-1))
            DIFF(I) = SCRH(I+1) - SCRH(I)
 2                 SOURCE(I) = SOURCE(I)
     &                + C(I)*(AH(I+1)+AH(I))*DIFF(I)
         SOURCE(I1) = SOURCE(I1) + C(I1)*(AH(I1P)+AH(I1))*
     &                 (SCRH(I1P)-SCRH(I1))
      Return

c  + D is added to SOURCE in an explicit formulation . . .
C-----------------------------------------------------------------------
 303      Do 3 I = I1, IN
 3                  SOURCE(I) = SOURCE(I) + DT*LO(I)*D(I)
      Return

c  + DIV(D) is computed conservatively from interface data . . .
C-----------------------------------------------------------------------
 404      SCRH(INP) = DT*AH(INP)*DN
         SCRH( I1) = DT*AH( I1)*D1
         Do 4 I = IN, I1P, -1
            SCRH(I)   = DT*AH(I)*D(I)
 4                 SOURCE(I) = SOURCE(I)+SCRH(I+1)-SCRH(I)
         SOURCE(I1) = SOURCE(I1) + SCRH(I1P) - SCRH(I1)
      Return

c  + C*GRAD(D) is computed using interface data . . .
C-----------------------------------------------------------------------
 505      SCRH( I1) = DTH*D1
         SCRH(INP) = DTH*DN
         Do 5 I = IN, I1P, -1
            SCRH(I) = DTH*D(I)
            DIFF(I) = SCRH(I+1) - SCRH(I)
 5                 SOURCE(I) = SOURCE(I)
     &                + C(I)*(AH(I+1)+AH(I))*DIFF(I)
         SOURCE(I1) = SOURCE(I1) + C(I1)*(AH(I1P)+AH(I1))*
     &                     (SCRH(I1P)-SCRH(I1))
      Return

c  + C for source terms only at a list of indices . . .
C-----------------------------------------------------------------------
 606      Do 6 IS = 1, NIND
            I = INDEX(IS)
 6                 SOURCE(I) = SOURCE(I) + SCALARS(IS)

      Return
      End

C=======================================================================

Cc$mta inline
      Subroutine CNVFCT ( RHOO, RHON, I1, IN,
     &                    SRHO1, VRHO1, SRHON, VRHON, PBC )

C-----------------------------------------------------------------------
c
c     Originated: J.P. Boris         Code 4400, NRL          Feb 1987
c     Modified:  Laboratory for Computational Physics & Fluid Dynamics
c     Contact:    J.P. Boris, J.H. Gardner, A.M. Landsberg, or E.S. Oran
c
c     Description:  This routine solves an advective continuity equation
c     of the form  dRHO/dt = -V*grad(RHO) + SOURCES in the user's choice
c     of Cartesian, cylindrical, or spherical coordinate systems.  A
c     facility is included to allow definition of other coordinates.
c     The grid can be Eulerian, sliding rezone, or Lagrangian and can
c     be arbitrarily spaced.  The algorithm is a low-phase-error FCT
c     algorithm, vectorized and optimized for a combination of speed and
c     flexibility.  A complete description appears in the NRL Memorandum
c     Report (1992), "LCPFCT - A Flux-Corrected Transport Algorithm For
c     Solving Generalized Continuity Equations".
c
c     Arguments:
c     RHOO   Real Array        grid point densities at start of step   I
c     RHON   Real Array        grid point densities at end of step     O
c     I1     Integer           first grid point of integration         I
c     IN     Integer           last grid point of intergration         I
c     SRHO1  Real Array        boundary guard cell factor at cell I1+1 I
c     VRHO1  Real Array        boundary value added to guard cell I1-1 I
c     SRHON  Real Array        boundary guard cell factor at cell IN+1 I
c     VRHON  Real Array        boundary value added to guard cell IN+1 I
c     PBC    Logical           periodic boundaries if PBC = .true.     I
c
c        In this routine the last interface at RADHN(INP) is the outer
c     boundary of the last cell indexed IN.  The first interface at
c     RADHN(I1) is the outer boundary of the integration domain before
c     the first cell indexed I1.  The description of CNVFCT and the
c     roles played by the auxiliary library routines is the same for
c     LCPFCT given above.
c
C-----------------------------------------------------------------------

          Implicit  NONE
          Integer   NPT, I1, IN, I1P, INP, I
          Real      BIGNUM, SRHO1, VRHO1, SRHON, VRHON, RHO1M, RHONP
          Real      RHOT1M, RHOTNP, RHOTD1M, RHOTDNP
          Logical   PBC
          Parameter ( NPT = 2002 )
          Parameter ( BIGNUM = 1.0E38 )
c     BIGNUM = Machine Dependent Largest Number - Set By The User!!!!

          Real     RHOO(NPT),     RHON(NPT)

c     /FCT_SCRH/ Holds scratch arrays for use by LCPFCT and CNVFCT
          Real     SCRH(NPT),     SCR1(NPT),     DIFF(NPT)
          Real     FLXH(NPT),     FABS(NPT),     FSGN(NPT)
          Real     TERM(NPT),     TERP(NPT),     LNRHOT(NPT)
          Real     LORHOT(NPT),   RHOT(NPT),     RHOTD(NPT)
      COMMON /FCT_SCRH/SCRH,SCR1,DIFF,FLXH,FABS,FSGN,
     &                       TERM, TERP, LNRHOT, LORHOT, RHOT, RHOTD
c$OMP THREADPRIVATE(/FCT_SCRH/)

c     /FCT_GRID/ Holds geometry, grid, area and volume information
          Real     LO(NPT),       LN(NPT),       AH (NPT)
          Real     RLN(NPT),      LH (NPT),      RLH(NPT)
          Real     ROH(NPT),      RNH(NPT),      ADUGTH(NPT)
      COMMON /FCT_GRID/LO,LN,AH,RLN,LH,RLH,ROH,RNH,ADUGTH
c$OMP THREADPRIVATE( /FCT_GRID/)

c     /FCT_VELO/ Holds velocity-dependent flux coefficients
          Real     HADUDTH(NPT),  NULH(NPT),     MULH(NPT)
          Real     EPSH(NPT),     VDTODR(NPT)
      COMMON /FCT_VELO/HADUDTH,NULH,MULH,EPSH,VDTODR
c$OMP THREADPRIVATE(/FCT_VELO/)

c     /FCT_MISC/ Holds the source array and diffusion coefficient
          Real     SOURCE(NPT),   DIFF1
      COMMON /FCT_MISC/SOURCE,DIFF1
c$OMP THREADPRIVATE( /FCT_MISC/)

C-----------------------------------------------------------------------
          I1P = I1 + 1
          INP = IN + 1

c     Calculate the convective and diffusive fluxes . . .
C-----------------------------------------------------------------------
          If ( PBC ) Then
             RHO1M = RHOO(IN)
             RHONP = RHOO(I1)
          Else
             RHO1M = SRHO1*RHOO(I1) + VRHO1
             RHONP = SRHON*RHOO(IN) + VRHON
          End If

          DIFF(I1) = NULH(I1) * ( RHOO(I1) - RHO1M )
          FLXH(I1) = VDTODR(I1) * ( RHOO(I1) - RHO1M )

          Do 1 I = I1P, IN
             DIFF(I) = ( RHOO(I) - RHOO(I-1) )
             FLXH(I) = VDTODR(I) * DIFF(I)
 1                   DIFF(I) = NULH(I) * DIFF(I)

          DIFF(INP) = NULH(INP) * ( RHONP - RHOO(IN) )
          FLXH(INP) = VDTODR(INP) * ( RHONP - RHOO(IN) )

c     Calculate LORHOT, the transported mass elements, and LNRHOT, the
c     transported & diffused mass elements  . . .
C-----------------------------------------------------------------------
          Do 2 I = I1, IN
             LORHOT(I) = LN(I) * (RHOO(I) - 0.5*(FLXH(I+1) + FLXH(I)))
     &                   + SOURCE(I)
             LNRHOT(I) = LORHOT(I) + (DIFF(I+1) - DIFF(I))
             RHOT(I)  = LORHOT(I)*RLN(I)
 2                   RHOTD(I) = LNRHOT(I)*RLN(I)

c     Evaluate the boundary conditions for RHOT and RHOTD . . .
C-----------------------------------------------------------------------
          If ( PBC ) Then
             RHOT1M = RHOT(IN)
             RHOTNP = RHOT(I1)
             RHOTD1M = RHOTD(IN)
             RHOTDNP = RHOTD(I1)
          Else
             RHOT1M = SRHO1*RHOT(I1) + VRHO1
             RHOTNP = SRHON*RHOT(IN) + VRHON
             RHOTD1M = SRHO1*RHOTD(I1) + VRHO1
             RHOTDNP = SRHON*RHOTD(IN) + VRHON
          End If

c     Calculate the transported antiduffusive fluxes and transported
c     and diffused density differences . . .
C-----------------------------------------------------------------------
          FLXH(I1) = MULH(I1) * ( RHOT(I1) - RHOT1M )
          DIFF(I1) = RHOTD(I1) - RHOTD1M
          FABS(I1) = ABS ( FLXH(I1) )
          FSGN(I1) = SIGN ( DIFF1, DIFF(I1) )

          Do 3 I = I1P, IN
             FLXH(I) = MULH(I) * ( RHOT(I) - RHOT(I-1) )
 3                   DIFF(I) = RHOTD(I) - RHOTD(I-1)

          FLXH(INP) = MULH(INP) * ( RHOTNP - RHOT(IN) )
          DIFF(INP) = RHOTDNP - RHOTD(IN)

c     Calculate the magnitude & sign of the antidiffusive flux followed
c     by the flux-limiting changes on the right and left . . .
C-----------------------------------------------------------------------
          Do 4 I = I1, IN
             FABS(I+1) = ABS ( FLXH(I+1) )
             FSGN(I+1) = SIGN ( DIFF1, DIFF(I+1) )
             TERM(I+1) = FSGN(I+1)*LN(I)*DIFF(I)
 4                   TERP(I) = FSGN(I)*LN(I)*DIFF(I+1)

          If ( PBC ) Then
             TERP(INP) = TERP(I1)
             TERM(I1) = TERM(INP)
          Else
             TERP(INP) = BIGNUM
             TERM(I1) = BIGNUM
          End If

c     Correct the transported fluxes completely and then calculate the
c     new Flux-Corrected Transport densities . . .
C-----------------------------------------------------------------------
          FLXH(I1) = FSGN(I1) * AMAX1 ( 0.0,
     &                  AMIN1 ( TERM(I1), FABS(I1), TERP(I1) ) )

          Do 5 I = I1, IN
             FLXH(I+1) = FSGN(I+1) * AMAX1 ( 0.0,
     &                AMIN1 ( TERM(I+1), FABS(I+1), TERP(I+1) ) )
             RHON(I) = RLN(I) * ( LNRHOT(I) + (FLXH(I) - FLXH(I+1)) )
 5                   SOURCE(I) = 0.0

      Return
      End

C=======================================================================

Cc$mta inline
      Subroutine CONSERVE ( RHO, I1, IN, CSUM )

C-----------------------------------------------------------------------
c
c     Description:   This routine computes the ostensibly conserved sum.
c     Beware your boundary conditions and note that only one continuity
c     equation is summed for each call to this subroutine.
c
c     Arguments:
c     RHO    Real Array(NPT)   cell values for physical variable 'RHO' I
c     I1     Integer           first cell to be integrated             I
c     IN     Integer           last cell  to be integrated             I
c     CSUM   Real              value of the conservation sum of rho    O
c
C-----------------------------------------------------------------------

          Implicit NONE
          Integer  NPT, I, I1, IN
          Parameter ( NPT = 2002 )
          Real     CSUM, RHO(NPT)

c     /FCT_GRID/ Holds geometry, grid, area and volume information
          Real     LO(NPT),       LN(NPT),       AH (NPT)
          Real     RLN(NPT),      LH (NPT),      RLH(NPT)
          Real     ROH(NPT),      RNH(NPT),      ADUGTH(NPT)
      COMMON /FCT_GRID/LO,LN,AH,RLN,LH,RLH,ROH,RNH,ADUGTH
c$OMP THREADPRIVATE(/FCT_GRID/)

c  Compute the ostensibly conserved total mass (BEWARE B.C.) . . .
C-----------------------------------------------------------------------
          CSUM = 0.0
          Do 80 I = I1, IN
 80                    CSUM = CSUM + LN(I)*RHO(I)

      Return
      End

C=======================================================================

Cc$mta inline
      Subroutine COPYGRID ( MODE, I1, IN )

C-----------------------------------------------------------------------
c
c     Description:   This Subroutine makes a complete copy of the grid
c     variables defined by the most recent call to MAKEGRID from cell
c     I1 to IN including the boundary values at interface IN+1 when the
c     argument MODE = 1.  When MODE = 2, these grid variables are reset
c     from common block OLD_GRID.  This routine is used where the same
c     grid is needed repeatedly after some of the values have been over-
c     written, for example, by a grid which moves between the halfstep
c     and the whole step.
c
c     Argument:
c     I1       Integer         first cell index                        I
c     IN       Integer         last cell index                         I
c     MODE     Integer         = 1 grid variables copied into OLD_GRID I
c                              = 2 grid restored from OLD_GRID common  I
c
C-----------------------------------------------------------------------

          Implicit  NONE
          Integer   NPT, I, MODE, I1, IN
          Parameter ( NPT = 2002 )

c     /OLD_GRID/ Holds geometry, grid, area and volume information
          Real     LOP(NPT),      LNP(NPT),      AHP(NPT)
          Real     RLNP(NPT),     RLHP(NPT),     LHP(NPT)
          Real     ROHP(NPT),     RNHP(NPT),     ADUGTHP(NPT)
      COMMON /OLD_GRID/LOP,LNP,AHP,RLNP,LHP,RLHP,
     &                       ROHP,     RNHP,      ADUGTHP
c$OMP THREADPRIVATE( /OLD_GRID/)

c     /FCT_GRID/ Holds geometry, grid, area and volume information
          Real     LO(NPT),       LN(NPT),       AH (NPT)
          Real     RLN(NPT),      LH (NPT),      RLH(NPT)
          Real     ROH(NPT),      RNH(NPT),      ADUGTH(NPT)
      COMMON /FCT_GRID/LO,LN,AH,RLN,LH,RLH,ROH,RNH,ADUGTH
c$OMP THREADPRIVATE( /FCT_GRID/)

C-----------------------------------------------------------------------
          If ( MODE .eq. 1 ) Then
             Do 101 I = I1, IN
                LOP(I)  = LO(I)
                LNP(I)  = LN(I)
 101                       RLNP(I) = RLN(I)
             Do 102 I = I1, IN+1
                AHP(I)  = AH(I)
                LHP(I)  = LH(I)
                RLHP(I) = RLH(I)
                ROHP(I) = ROH(I)
                RNHP(I) = RNH(I)
 102                       ADUGTHP(I) = ADUGTH(I)
          Else If ( MODE .eq. 2 ) Then
             Do 201 I = I1, IN
                LO(I)  = LOP(I)
                LN(I)  = LNP(I)
 201                       RLN(I) = RLNP(I)
             Do 202 I = I1, IN+1
                AH(I)  = AHP(I)
                LH(I)  = LHP(I)
                RLH(I) = RLHP(I)
                ROH(I) = ROHP(I)
                RNH(I) = RNHP(I)
 202                       ADUGTH(I) = ADUGTHP(I)
          Else
             Write ( 6, 1001 ) MODE
          End If
 1001     Format ( ' COPYGRID Error! MODE =', I3, ' (not 1 or 2!)' )

      Return
      End

C=======================================================================

      Block Data FCTBLK

C-----------------------------------------------------------------------
          Implicit  NONE
          Integer   NPT
          Parameter ( NPT = 2002 )

c     /FCT_MISC/ Holds the source array and diffusion coefficient
          Real     SOURCE(NPT),   DIFF1
      COMMON /FCT_MISC/SOURCE,DIFF1
c$OMP THREADPRIVATE(/FCT_MISC/)

          Data SOURCE / NPT*0.0 /, DIFF1 / 0.999 /

      End

C=======================================================================

Cc$mta inline
      Subroutine NEW_GRID ( RADHN, I1, INP, ALPHA )

C-----------------------------------------------------------------------
c
c     Description:  This Subroutine initializes geometry variables and
c     coefficients when the most recent call to MAKEGRID used the same
c     set of values RADHO and only the new interface locations RADHN are
c     different.  NEW_GRID is computationally more efficienty than the
c     complete grid procedure in MAKEGRID because several formulae do
c     not need to be recomputed.  The grid should generally be defined
c     for the entire number of grid interfaces from 1 to INP,  however
c     subsets of the entire grid may be reinitialized with care.
c
c     Arguments:
c     RADHN    Real Array(INP)    new cell interface positions         I
c     I1       Integer            first interface index                I
c     INP      Integer            last interface index                 I
c     ALPHA    Integer            = 1 for cartesian geometry           I
c                                 = 2 for cylindrical geometry         I
c                                 = 3 for spherical geometry           I
c                                 = 4 general geometry (user supplied) I
c
C-----------------------------------------------------------------------
          Implicit  NONE
          Integer   NPT, I1, I1P, I, IN, INP, ALPHA
          Parameter ( NPT = 2002 )

          Real     RADHN(INP), PI, FTPI

c     /FCT_SCRH/ Holds scratch arrays for use by LCPFCT and CNVFCT
          Real     SCRH(NPT),     SCR1(NPT),     DIFF(NPT)
          Real     FLXH(NPT),     FABS(NPT),     FSGN(NPT)
          Real     TERM(NPT),     TERP(NPT),     LNRHOT(NPT)
          Real     LORHOT(NPT),   RHOT(NPT),     RHOTD(NPT)
      COMMON /FCT_SCRH/SCRH,SCR1,DIFF,FLXH,FABS,FSGN,
     &                       TERM, TERP, LNRHOT, LORHOT, RHOT, RHOTD
c$OMP THREADPRIVATE( /FCT_SCRH/)

c     /FCT_GRID/ Holds geometry, grid, area and volume information
          Real     LO(NPT),       LN(NPT),       AH (NPT)
          Real     RLN(NPT),      LH (NPT),      RLH(NPT)
          Real     ROH(NPT),      RNH(NPT),      ADUGTH(NPT)
      COMMON /FCT_GRID/LO,LN,AH,RLN,LH,RLH,ROH,RNH,ADUGTH
c$OMP THREADPRIVATE( /FCT_GRID/)

          DATA     PI, FTPI /3.1415927, 4.1887902/

C-----------------------------------------------------------------------
          I1P = I1 + 1
          IN = INP - 1

c  Store the old and new grid interface locations from input and then
c  update the new and average interface and grid coefficients . . .
C-----------------------------------------------------------------------
          Do 1 I = I1, INP
 1                   RNH(I) = RADHN(I)

c  Select the choice of coordinate systems . . .
C-----------------------------------------------------------------------
          Go To (100, 200, 300, 400), ALPHA

c  Cartesian coordinates . . .
C-----------------------------------------------------------------------
 100            AH(INP) = 1.0
          Do 101 I = I1, IN
 101                 LN(I) = RNH(I+1) - RNH(I)
          Go To 500

c  Cylindrical Coordinates: RADIAL . . .
C-----------------------------------------------------------------------
 200           DIFF(I1) = RNH(I1)*RNH(I1)
          AH(INP) = PI*(ROH(INP) + RNH(INP))
          DO 201 I = I1, IN
             AH(I) = PI*(ROH(I) + RNH(I))
             DIFF(I+1) = RNH(I+1)*RNH(I+1)
 201                 LN(I) = PI*(DIFF(I+1) - DIFF(I))
          Go To 500

c  Spherical Coordinates: RADIAL . . .
C-----------------------------------------------------------------------
 300           DIFF(I1) = RNH(I1)*RNH(I1)*RNH(I1)
          SCRH(INP) = (ROH(INP) + RNH(INP))*ROH(INP)
          AH(INP) = FTPI*(SCRH(INP) + RNH(INP)*RNH(INP))
          DO 301 I = I1, IN
             DIFF(I+1) = RNH(I+1)*RNH(I+1)*RNH(I+1)
             SCRH(I) = (ROH(I) + RNH(I))*ROH(I)
             AH(I) = FTPI*(SCRH(I) + RNH(I)*RNH(I))
 301                 LN(I) = FTPI*(DIFF(I+1) - DIFF(I))
          Go To 500

c  Special Coordinates: Areas and Volumes are User Supplied . . .
C-----------------------------------------------------------------------
 400           Continue

c  Additional system independent geometric variables . . .
C-----------------------------------------------------------------------
 500                Do 501 I = I1, IN
 501                           RLN(I) = 1.0/LN(I)
          LH(I1)  = LN(I1)
          RLH(I1) = RLN(I1)
          Do 502 I = I1P, IN
             LH(I) =  0.5*(LN(I) + LN(I-1))
 502                 RLH(I) = 0.5*(RLN(I) + RLN(I-1))
          LH(INP)  = LN(IN)
          RLH(INP) = RLN(IN)
          Do 503 I = I1, INP
 503                 ADUGTH(I) = AH(I)*(RNH(I) - ROH(I))

      Return
      End

C=======================================================================

Cc$mta inline
      Subroutine RESIDIFF ( DIFFA )

C-----------------------------------------------------------------------
c
c     Description:   Allows the user to give FCT some residual numerical
c     diffusion by making the anti-diffusion coefficient smaller.
c
c     Arguments:
c     DIFFA  Real    Replacement residual diffusion coefficient        I
c                    Defaults to 0.999 but could be as high as 1.0000
c
C-----------------------------------------------------------------------
          Implicit NONE
          Integer  NPT
          Real     DIFFA
          Parameter ( NPT = 2002 )

c     /FCT_MISC/ Holds the source array and diffusion coefficient
          Real     SOURCE(NPT),   DIFF1
      COMMON /FCT_MISC/SOURCE,DIFF1
c$OMP THREADPRIVATE( /FCT_MISC/)

          DIFF1 = DIFFA

      Return
      End

C=======================================================================

Cc$mta inline
      Subroutine SET_GRID ( RADR, I1, IN )

C-----------------------------------------------------------------------
c
c     Description:  This subroutine includes the radial factor in the
c     cell volume for polar coordinates. It must be preceeded by a call
c     to MAKE_GRID with ALPHA = 1 to establish the angular dependence of
c     the cell volumes and areas and a call to COPY_GRID to save this
c     angular dependence.  The angular coordinate is measured in radians
c     (0 to 2 pi) in cylindrical coordinates and cos theta (1 to -1) in
c     spherical coordinates.  SET_GRID is called inside the loop over
c     radius in a multidimensional model to append the appropriate
c     radial factors when integrating in the angular direction.
c
c     Arguments:
c     RADR     Real               radius of cell center                I
c     I1       Integer            first cell index                     I
c     IN       Integer            last cell index                      I
c
C-----------------------------------------------------------------------
          Implicit  NONE
          Integer   NPT, I1, I1P, I, IN, INP
          Real      RADR
          Parameter ( NPT = 2002 )

c     /OLD_GRID/ Holds geometry, grid, area and volume information
          Real     LOP(NPT),      LNP(NPT),      AHP(NPT)
          Real     RLNP(NPT),     RLHP(NPT),     LHP(NPT)
          Real     ROHP(NPT),     RNHP(NPT),     ADUGTHP(NPT)
      COMMON /OLD_GRID/LOP,LNP,AHP,RLNP,LHP,RLHP,
     &                       ROHP,     RNHP,      ADUGTHP
c$OMP THREADPRIVATE(/OLD_GRID/)

c     /FCT_GRID/ Holds geometry, grid, area and volume information
          Real     LO(NPT),       LN(NPT),       AH (NPT)
          Real     RLN(NPT),      LH (NPT),      RLH(NPT)
          Real     ROH(NPT),      RNH(NPT),      ADUGTH(NPT)
      COMMON /FCT_GRID/LO,LN,AH,RLN,LH,RLH,ROH,RNH,ADUGTH
c$OMP THREADPRIVATE(/FCT_GRID/)

C-----------------------------------------------------------------------
          I1P = I1 + 1
          INP = IN + 1

C  Multiply each volume element by the local radius
          DO 100 I = I1, IN
          LN(I) = LNP(I)*RADR
 100           LO(I) = LOP(I)*RADR

c  Additional system independent geometric variables . . .
C-----------------------------------------------------------------------
 500                Do 501 I = I1, IN
 501                           RLN(I) = 1.0/LN(I)
          LH(I1)  = LN(I1)
          RLH(I1) = RLN(I1)
          Do 502 I = I1P, IN
             LH(I) =  0.5*(LN(I) + LN(I-1))
 502                 RLH(I) = 0.5*(RLN(I) + RLN(I-1))
          LH(INP)  = LN(IN)
          RLH(INP) = RLN(IN)

      Return
      End

C=======================================================================

Cc$mta inline
      Subroutine ZERODIFF ( IND )

C-----------------------------------------------------------------------
c
c     Description:   This Subroutine sets the FCT diffusion and anti-
c     diffusion parameters to zero at the specified cell interface to
c     inhibit unwanted diffusion across the interface.  This routine is
c     used for inflow and outflow boundary conditions.  If argument IND
c     is positive, the coefficients at that particular interface are
c     reset.  If IND is negative, the list of NIND indices in INDEX are
c     used to reset that many interface coefficients.
c
c     Argument:
c     IND    Integer              index of interface to be reset       I
c
C-----------------------------------------------------------------------
          Implicit  NONE
          Integer   NPT, NINDMAX, IND, IS, I
          Parameter ( NPT = 2002, NINDMAX = 150 )

c     /FCT_NDEX/ Holds a scalar list of special cell information . . .
          Real     SCALARS(NINDMAX)
          Integer  INDEX(NINDMAX), NIND
c     COMMON /FCT_NDEX/NIND,INDEX,SCALARS
       COMMON /FCT_NDEX/scalars,NIND,INDEX
c$OMP THREADPRIVATE( /FCT_NDEX/)

c     /FCT_VELO/ Holds velocity-dependent flux coefficients
          Real     HADUDTH(NPT),  NULH(NPT),     MULH(NPT)
          Real     EPSH(NPT),     VDTODR(NPT)
      COMMON /FCT_VELO/HADUDTH,NULH,MULH,EPSH,VDTODR
c$OMP THREADPRIVATE(/FCT_VELO/)

C-----------------------------------------------------------------------
          If ( IND .gt. 0 ) Then
              NULH(IND) = 0.0
              MULH(IND) = 0.0
          Else If ( IND .le. 0 ) Then
             If ( NIND.lt.1 .or. NIND.gt.NINDMAX .or. IND.eq.0 ) Then
                Write ( 6,* ) ' ZERODIFF Error! IND, NIND =', IND, NIND
                Stop
             End If
             Do IS = 1, NIND
                I = INDEX(IS)
                NULH(I) = 0.0
                MULH(I) = 0.0
             End Do
          End If

      Return
      End

C=======================================================================

Cc$mta inline
      Subroutine ZEROFLUX ( IND )

C-----------------------------------------------------------------------
c
c     Description:   This Subroutine sets all the velocity dependent FCT
c     parameters to zero at the specified cell interface to inhibit
c     transport fluxes AND diffusion of material across the interface.
c     This routine is needed in solid wall boundary conditions.  If IND
c     is positive, the coefficients at that particular interface are
c     reset.  If IND is negative, the list of NIND indices in INDEX are
c     used to reset that many interface coefficients.
c
c     Argument:
c     IND    Integer              index of interface to be reset       I
c
C-----------------------------------------------------------------------
          Implicit  NONE
          Integer   NPT, NINDMAX, IND, IS, I
          Parameter ( NPT = 2002, NINDMAX = 150 )

c     /FCT_NDEX/ Holds a scalar list of special cell information . . .
          Real     SCALARS(NINDMAX)
          Integer  INDEX(NINDMAX), NIND
c     COMMON /FCT_NDEX/NIND,INDEX,SCALARS
       COMMON /FCT_NDEX/scalars,NIND,INDEX
c$OMP THREADPRIVATE(/FCT_NDEX/)

c     /FCT_VELO/ Holds velocity-dependent flux coefficients
          Real     HADUDTH(NPT),  NULH(NPT),     MULH(NPT)
          Real     EPSH(NPT),     VDTODR(NPT)
      COMMON /FCT_VELO/HADUDTH,NULH,MULH,EPSH,VDTODR
c$OMP THREADPRIVATE( /FCT_VELO/ )

C-----------------------------------------------------------------------
          If ( IND .gt. 0 ) Then
             HADUDTH(IND) = 0.0
             NULH(IND) = 0.0
             MULH(IND) = 0.0
          Else If ( IND .le. 0 ) Then
             If ( NIND.lt.1 .or. NIND.gt.NINDMAX .or. IND.eq.0 ) Then
                Write ( 6,* ) ' ZEROFLUX Error! IND, NIND =', IND, NIND
                Stop
             End If
             Do IS = 1, NIND
                I = INDEX(IS)
                HADUDTH(I) = 0.0
                NULH(I) = 0.0
                MULH(I) = 0.0
             End Do
          End If

      Return
      End

C=======================================================================