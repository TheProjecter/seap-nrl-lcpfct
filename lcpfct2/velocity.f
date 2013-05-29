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
          Integer  I1, I1P, I, IN, INP

          Real     UH(0:INP), DT, RDT, DTH, DT2, DT4, ONE3RD, ONE6TH

	  include 'prm.h'
	  include 'fct.h'


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
