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
          Integer  MODE, IS, I, I1, IN, I1P, INP

          include 'prm.h'
          include 'fct.h'

          Real     C(0:NPT), D(0:NPT), DT, DTH, DTQ, D1, DN

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
