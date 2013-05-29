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
          Integer   IND, IS, I

          include 'prm.h'
          include 'fct.h'

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
