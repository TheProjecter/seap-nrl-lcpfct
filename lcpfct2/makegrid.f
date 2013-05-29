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
          Integer   I1, I1P, I, IN, INP, ALPHA

          Real     RADHO(0:INP), RADHN(0:INP), PI, FTPI

          include 'prm.h'
          include 'fct.h'

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
