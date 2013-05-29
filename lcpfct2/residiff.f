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
          Real     DIFFA

          include 'prm.h'
	  include 'fct.h'

          DIFF1 = DIFFA

      Return
      End
