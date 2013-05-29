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