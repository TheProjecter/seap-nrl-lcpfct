C=======================================================================

      Subroutine GASDYNP( K1, KN, BC1, BCN, DT, nx ,ny, j)

C-----------------------------------------------------------------------
c
c  This routine integrates the gasdynamic equations using the momentum
c  component RVRN as the direction of integration and the momentum RVTN 
c  as the transverse direction.  In 2D models the two directions of 
c  integration are chosen by exchanging RVRN and RVTN in Common.

c  K1  . . .   Index of the integration's first cell       
c  KN  . . .   Index of the integration's last cell        
c  BC1 . . .   Indicates boundary condition on integration at K1
c  BCN . . .   Indicates boundary condition on integration at K1
c  DT  . . .   Timestep for the integrations of this step
c
C-----------------------------------------------------------------------

         Implicit  NONE
         Integer    K1, K1P, BC1, BCN, K, KN, KNP, IT,ndim,idim,K1M,nx,ny,j,ii
         include 'prm.h'
         Logical   PBC
         Real      SBC1,  SRV1,  SBCN,     SRVN,  VRHO1, VRHON
         Real      VRVR1, VRVRN, VRVT1,    VRVTN, VERG1, VERGN
         Real      MPINT(0:NPT), VEL(0:NPT),   UNIT(0:NPT),  ZERO(0:NPT)
         Real      RHOO(0:NPT),  rvvo(0:NPT,3),  ERGO(0:NPT)
         Real      VINT(0:NPT),  PRE(0:NPT),   MPVINT(0:NPT)
         Real      DTSUB,      DT,         RELAX

         Real      RHO_IN,     PRE_IN,     VEL_IN,     GAMMA0
         Real      RHOAMB,     PREAMB,     VELAMB,     GAMMAM
         Real      RHON(0:NPT),  RVVN(0:NPT,3),   ERGN(0:NPT), RHONI(0:NPT)
         real      rhobc(2),prebc(2),sbci(3),sbcj(3),vrvvi(3),vrvvj(3)
         real      bc(4,0:npt),bcx(4,0:npt,0:npt),bcv(2,0:npt,0:npt,2),bce(2,0:npt,0:npt),err
	 integer   ibc(0:npt,0:npt)
       COMMON /ARRAYS/RHON,RVVN,ERGN,RELAX,rhobc,prebc,ndim
c$OMP THREADPRIVATE(/ARRAYS/)
C      TASKCOMMON /ARRAYS/RHON,RVRN,RVTN,ERGN,RELAX
      COMMON  /INVARIANTS/ 
     &                        RHO_IN, PRE_IN, VEL_IN, GAMMA0,
     &                        RHOAMB, PREAMB, VELAMB, GAMMAM
         Data      UNIT / NPT1*1.0 /,       ZERO / NPT1*0.0 /
 	common /tau_prf/ profiler
        common /fct_bc / bc,bcx,bcv,bce,ibc
 	integer profiler(2,10),isdim(3,0:2)
	data isdim /1,1,1,-1,1,1,1,1,1/

c  Prepare for the time integration. Index K is either I or J depending
c  on the definitions of RVRN and RVTN. Copies of the physical variable 
c  are needed to recover values for the whole step integration . . .
C-----------------------------------------------------------------------
      KNP = KN + 1
      K1P = K1 + 1
      K1M = K1 - 1
      PBC = .false.
      If ( BC1.eq.3 .OR. BCN.eq.3 ) PBC = .true.
      Do 50 K = K1, KN
         RHOO(K) = RHON(K)
        rvvo(k,1)=rvvn(k,1)
	rvvo(k,2)=rvvn(k,2)
c	rvvo(k,3)=rvvn(k,3)
 50     ERGO(K) = ERGN(K)

c  Integrate first the half step then the whole step . . .
C-----------------------------------------------------------------------
      Do 500 IT = 1, 2
  
         DTSUB = 0.5*DT*FLOAT(IT)

         Do 100 K = K1, KN
            rhoni(k)=1.0/rhon(k)
            vel  (k)=rvvn(k,1)*rhoni(k)
           pre(k)=gammam*(ergn(k)-0.5*(rvvn(k,1)**2+rvvn(k,2)**2)*rhoni(k))
c           pre(k)=gammam*(ergn(k)-0.5*(rvvn(k,1)**2+rvvn(k,2)**2+rvvn(k,3)**2)*rhoni(k))
 100   continue      


c  Calculate the interface velocities and pressures as weighted values
c  of the cell-centered values computed just above . . .
C-----------------------------------------------------------------------
         Do 200 K = K1+1, KN
            VINT(K)  =  0.5*(VEL(K) + VEL(K-1))
            MPINT(K) = -0.5*(PRE(K) + PRE(K-1))
 200        MPVINT(K) = -0.5*(PRE(K)*VEL(K) + PRE(K-1)*VEL(K-1))

c  The unweighted interface averages can be computed as follows . . .


c  Call the FCT utility routines and set the boundary conditions. Other
c  boundary conditions could be added for inflow, outflow, etc . . .
c  BC1, BCN = 1  => ideal solid wall or axis boundary condition 
c  BC1, BCN = 2  => an extrapolative outflow boundary condition 
c  BC1, BCN = 3  => periodic boundary conditions . . .
c  BC1, BCN = 4  => specified boundary values (e.g. shock tube problem)
C-----------------------------------------------------------------------
         Go To ( 310, 320, 330, 340 ), BC1
 310        VINT  (K1) = 0.0
            MPINT (K1) = - PRE(K1)
            MPVINT(K1) = 0.0
            Go To 350
 320        VINT(K1)   = VEL(K1)*(1.0 - RELAX) 
!           MPINT(K1)  = - PRE(K1)*(1.0 - RELAX) - RELAX*PRE_IN
            MPINT(K1)  = - PRE(K1)*(1.0 - RELAX) - RELAX*prebc(1)
            MPVINT(K1) = MPINT(K1)*VINT(K1)
            Go To 350
 330               MPVINT(K1) = 1.0/( RHON(K1) + RHON(KN) )
            VINT(K1)   = (VEL(K1)*RHON(KN)+VEL(KN)*RHON(K1)) *MPVINT(K1)
            MPINT(K1)  = -(PRE(K1)*RHON(KN)+PRE(KN)*RHON(K1))*MPVINT(K1)
            MPVINT(K1) = -(PRE(K1)*VEL(K1)*RHON(KN) 
     &                + PRE(KN)*VEL(KN)*RHON(K1))*MPVINT(K1)
            Go To 350
 340               VINT(K1)   = VEL_IN
            MPINT(K1)  = - PRE_IN
            MPVINT(K1) = - PRE_IN*VEL_IN

 350            Go To ( 410, 420, 430, 440 ), BCN
 410        VINT  (KNP) = 0.0
            MPINT (KNP) = - PRE(KN)
            MPVINT(KNP) = 0.0
            Go To 450
 420        VINT(KNP)   = VEL(KN)*(1.0 - RELAX)
            MPINT(KNP)  = - PRE(KN)*(1.0 - RELAX) - RELAX*prebc(2)
            MPVINT(KNP) = MPINT(KNP)*VINT(KNP)
            Go To 450
 430               VINT(KNP)   = VINT(K1)
            MPINT(KNP)  = MPINT(K1)
            MPVINT(KNP) = MPVINT(K1)
            Go To 450
 440               VINT(KNP)   = VELAMB
            MPINT(KNP)  = - PREAMB
            MPVINT(KNP) = - PREAMB*VELAMB
 450            Continue

c  The velocity dependent FCT coefficients are set and the boundary 
c  condition calculations are completed.  Here the periodic boundary 
c  conditions require no action as (S)lope and (V)alue boundary value
c  specifiers are ignored in LCPFCT when PBC = .true.
C-----------------------------------------------------------------------
         Call VELOCITY ( VINT, K1, KNP, DTSUB )

         Go To ( 510, 520, 550, 540 ), BC1
 510            Call ZEROFLUX ( K1 )
            SBC1  = 1.0
            SRV1  = -1.0
            sbci(:)= 1.0
            sbci(1)=-1.0
            VRHO1 = 0.0
            VRVR1 = 0.0
            VRVT1 = 0.0
            vrvvi(:)=0.0
            VERG1 = 0.0
            Go To 550
 520               Call ZERODIFF ( K1 )
            SBC1  = 1.0 - RELAX
            SRV1  = 1.0 - RELAX
            sbci(:)=1.0-relax
!           VRHO1 = RELAX*RHO_IN
            VRHO1 = RELAX*rhobc(1)
            VRVR1 = 0.0
            VRVT1 = 0.0
            vrvvi(:)=0.0
!           VERG1 = RELAX*PRE_IN/GAMMAM
            VERG1 = RELAX*prebc(1)/GAMMAM
            Go To 550
 540               SBC1  = 0.0
            SRV1  = 0.0
            VRHO1 = RHO_IN
            VRVR1 = RHO_IN*VEL_IN
            VRVT1 = 0.0
            VERG1 = PRE_IN/GAMMAM + 0.5*RHO_IN*VEL_IN**2

 550            Go To ( 610, 620, 650, 640 ), BCN
 610                   Call ZEROFLUX ( KNP )
            SBCN = 1.0
            SRVN = -1.0
            sbcj(:)= 1.0
            sbcj(1)=-1.0
            VRHON = 0.0
            VRVRN = 0.0
            VRVTN = 0.0
            vrvvj(:)=-0.0
            VERGN = 0.0
            Go To 650
 620               Call ZERODIFF ( KNP )
            SBCN  = 1.0 - RELAX
            SRVN  = 1.0 - RELAX
            sbcj(:)=1.0-relax
            VRHON = RELAX*rhobc(2)
            VRVRN = 0.0
            VRVTN = 0.0
            vrvvj(:)=0.0
            VERGN = RELAX*prebc(2)/GAMMAM
            Go To 650
 640               SBCN  = 0.0
            SRVN  = 0.0
            VRHON = RHOAMB
            VRVRN = RHOAMB*VELAMB
            VRVTN = 0.0
            VERGN = PREAMB/GAMMAM + 0.5*RHOAMB*VELAMB**2
 650            Continue

c  Integrate the continuity equations using LCPFCT . . .  
C-----------------------------------------------------------------------
 	call tau_profile_start(profiler(1,6))
C	bc(1,k1m)= sbc1
C	bc(2,k1m)=vrho1
C	bc(1,knp)= sbcn
C	bc(2,knp)=vrhon

	bc(1,k1m:knp)=bcx(1,k1m:knp,j)
	bc(2,k1m:knp)=bcx(2,k1m:knp,j)
	bc(3,k1m:knp)=bcx(3,k1m:knp,j)
	bc(4,k1m:knp)=bcx(4,k1m:knp,j)

C	bc(1,0:nx+1)=bcx(1,0:nx+1,j)
C	bc(2,0:nx+1)=bcx(2,0:nx+1,j)
C	bc(3,0:nx+1)=bcx(3,0:nx+1,j)
C	bc(4,0:nx+1)=bcx(4,0:nx+1,j)
	rhoo(k1m)=rhoo(k1)
	rhoo(knp)=rhoo(kn)
         Call LCPFCTP( RHOO, RHON, K1,KN, SBC1,VRHO1, SBCN,VRHON, PBC, nx, ny)
 	call tau_profile_stop (profiler(1,6))

         call SOURCES( K1,KN, DTSUB, 5, UNIT, MPINT, 
     &                                        MPINT(K1),  MPINT(KNP)  )    

 	call tau_profile_start(profiler(1,6))
	do idim=1,ndim
C	bc(1,k1m)= sbci(idim)
C	bc(2,k1m)=vrvvi(idim)
C	bc(1,knp)= sbcj(idim)
C	bc(2,knp)=vrvvj(idim)
	bc(1,k1m:knp)=bcv(1,k1m:knp,j,idim)*isdim(idim,ibc(k1m:knp,j))
	bc(2,k1m:knp)=bcv(2,k1m:knp,j,idim)*isdim(idim,ibc(k1m:knp,j))
	rvvo(k1m,idim)=rvvo(k1,idim)
	rvvo(knp,idim)=rvvo(kn,idim)
	call lcpfctp(rvvo(0,idim),rvvn(0,idim),k1,kn,sbci(idim),vrvvi(idim),sbcj(idim),vrvvj(idim),pbc,nx,ny)
	end do
 	call tau_profile_stop (profiler(1,6))

         call SOURCES( K1,KN, DTSUB, 4, UNIT, MPVINT, 
     &                                        MPVINT(K1), MPVINT(KNP) )    

 	call tau_profile_start(profiler(1,6))
C	bc(1,k1m)= sbc1
C	bc(2,k1m)=verg1
C	bc(1,knp)= sbcn
C	bc(2,knp)=vergn
	bc(1,k1m:knp)=bce(1,k1m:knp,j)
	bc(2,k1m:knp)=bce(2,k1m:knp,j)
	ergo(k1m)=ergo(k1)
	ergo(knp)=ergo(kn)
         Call LCPFCTP( ERGO, ERGN, K1,KN, SBC1,VERG1, SBCN,VERGN, PBC,nx,ny )
 	call tau_profile_stop (profiler(1,6))

 500      Continue       ! End of halfstep-wholestep loop.
      Return
      End

C=======================================================================
