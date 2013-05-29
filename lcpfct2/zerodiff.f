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
          Integer   IND, IS, I

          include 'prm.h'
	  include 'fct.h'

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
