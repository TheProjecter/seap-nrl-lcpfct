C234567
      subroutine fctblk_init
C  This subroutine was added to the openmp version of lcpfct to
C  facilitate initialization of the threat private fct_misc common
C  block.  Inserted by Rick Roberts on or about 9/17/2002.  The
C  authors of the lcpfct code are aware of its inclusion.
      implicit none
      integer npt, i
      Parameter (npt = 2002)
c     /FCT_MISC/ Holds the source array and diffusion coefficient
          Real     SOURCE(0:NPT),   DIFF1
       COMMON /FCT_MISC/SOURCE,DIFF1
c$OMP THREADPRIVATE( /FCT_MISC/)

      do i = 1,npt
         source(i)=0.0
      end do
      diff1 = 0.0
      return
      end
