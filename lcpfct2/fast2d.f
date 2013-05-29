C=======================================================================

          Program FAST2D
CC$mta parallel off

C-----------------------------------------------------------------------
c
c  BURSTING DIAPHRAGM "MUZZLE FLASH" - LCPFCT TEST # 4       August 1992
c
c  The problem begins with 1000:1 pressure and 100:1 density ratios 
c  across a diaphragm inside a solid cylindrical barrel.  The ideal wall
c  of the barrel is 10 cells thick (1.0 cm) with its inner radius given
c  as 1.5 cm and its outer radius of 2.5 cm.  The run starts at time 
c  t = 0.0 when the diaphragm at interface J = 11 (inside the barrel) is
c  ruptured.  The flow then expands upward in a 1D manner, spilling out 
c  of the barrel in a 2D flow which eventually reaches the boundaries at
c  R = 4.0 cm and Z = 4.0 cm where a very simple extrapolative outflow
c  condition is expressed through the LCPFCT boundary conditions values.
c  The outflow condition used here includes a slow relaxation to ambient 
c  conditions far from the origin.
c
C-----------------------------------------------------------------------

         Implicit  NONE

	 include 'prm.h'

         Integer   I, J, IJ, NSIZE,nsave,nr0,nz0,nx,nxp,ialfx,ndim
         real      dr0,dz0,time0,timel
         Parameter ( NSIZE = 516, nsave=100000, nr0=129, nz0=129, time0=3.3612e-5)
	 parameter ( nx=2*nsize)
         parameter ( dr0=0.1,dz0=0.1 )
         Integer   NR, NRP, IALFR,     BC_AXIS, BC_WALL, BC_OUTF
         Integer   NZ, NZP, IALFZ,     LOUT, MAXSTP, IPRINT
         Parameter (NR = NSIZE, NZ = 2*NSIZE)
         Integer   ICIN, ICOUT, JCTOP,     JSTEP, ISTEP
	 integer icount,irate,imax,jcount
	 real twall,twall0,tcpu,ctime
         Real      DR,         DZ,         DT,         TIME, DX
         Real      COURANT,    DTNEW,      VTYPICAL,   RELAX
         Real      DTMAX,      VZMAX,      R(0:NPT),     Z(0:NPT), X(0:NPT)
         Real RHO(nx,NZ), RVR(nx,NZ), RVZ(nx,NZ), ERG(nx,NZ)
         Real RHOT(nz,nx), RVRT(nz,nx), RVZT(nz,nx), ERGT(nz,nx)
         Real val(4,nx,nz)
         Real*4 RHO4(nx,NZ), RVR4(nx,NZ), RVZ4(nx,NZ), ERG4(nx,NZ)
         Real      DIN(9)


         Real      RHO_IN,     PRE_IN,     VEL_IN,     GAMMA0
         Real      RHOAMB,     PREAMB,     VELAMB,     GAMMAM
         Real      RHON(0:NPT),  RVRN(0:NPT),  RVTN(0:NPT),  ERGN(0:NPT)
         real      rhobc(2), prebc(2),rvvn(0:npt,3)
         real      bc(4,0:npt)
	 character(4) astep
         logical lil,lir,lim
       COMMON /ARRAYS/RHON,rvvn,ERGN,RELAX,rhobc,prebc,ndim
c$OMP THREADPRIVATE(/ARRAYS/)
C      TASKCOMMON /ARRAYS/RHON,RVRN,RVTN,ERGN,RELAX
      COMMON  /INVARIANTS/ 
     &                         RHO_IN, PRE_IN, VEL_IN, GAMMA0,
     &                         RHOAMB, PREAMB, VELAMB, GAMMAM

 	common /tau_prf/ profiler
	common /fct_bc/ bc

  	integer profiler(2,10)
  	call tau_profile_init()
  	call tau_profile_set_node(0)
  	call tau_profile_timer(profiler(1,1),'fast2d')
  	call tau_profile_timer(profiler(1,2),'integx')
  	call tau_profile_timer(profiler(1,3),'integy')
  	call tau_profile_timer(profiler(1,5),'gasdyn')
  	call tau_profile_timer(profiler(1,6),'lcpfct')
  	call tau_profile_start(profiler(1,1))
	ndim=2

c  The 2D barrel explosion program control parameters are specified.
c  (Change here to run other cases) . . .
C-----------------------------------------------------------------------
c      NR      =  129  ! Number of cells in the first (radial R) direction
!     IALFR   =   2  ! Sets cylindrical coordinates in the R direction
      IALFX   =   1  ! Sets cartesian coordinates in the R direction
      DR      =  0.1 ! Cell size (e.g., cm) in the radial direction
      DX      =  0.1
c      NZ      =  129  ! Number of cells in the second (axial Z) direction
      IALFZ   =   1  ! Sets Cartesian coordinates in the Z direction
      DZ      =  0.1 ! Cell size (e.g., cm) in the axial direction
      LOUT    =   6  ! Logical unit number of printed output device
      BC_AXIS =   1  ! Cylindrical axis set as an impermeable wall
      BC_OUTF =   2  ! Outer boundaries set as extrapolative outflow
      BC_WALL =   1  ! Walls of the barrel set as an ideal solid wall
      MAXSTP  = 700  ! Maximum number of timesteps of length DT
c     MAXSTP  =  65  ! Maximum number of timesteps of length DT
c     MAXSTP  =  26  ! Maximum number of timesteps of length DT
c     IPRINT  = 100  ! Initial frequency of validation printout results
      IPRINT  =  25  ! Initial frequency of validation printout results
      COURANT =  0.4 ! Approximate maximum Cournat number allowed
      DTMAX = 2.0E-7 ! Maximum timestep allowed in the computation
      DT    = 1.0E-9 ! Initial (small guess) for starting timestep
      
       DR=DR0*NR0/NR
       DX=DR0*NR0/NR
       DZ=DZ0*NR0/NR
       TIMEL=TIME0*NR/NR0
!      courant=courant*dr0/dr

c  Initialize the test problem geometry, a cylindrical shell JCTOP cells 
c  high in Z (indexed by J) which extends from the left of cell ICIN to 
c  the right of cell ICOUT in X (indexed by I).  
C-----------------------------------------------------------------------
      ICIN    =  16  ! Number of the innermost radial cell in the barrel
      ICOUT   =  25  ! Number of the outermost radial cell in the barrel
      JCTOP   =  20  ! Number of the uppermost axial cell in the barrel
      icin  = icin *float(nr)/nr0
      icout = icout*float(nr)/nr0
      jctop = jctop*float(nr)/nr0
      GAMMA0 = 1.4   ! Gas constant (adiabatic index)
      RHOAMB = 0.00129       ! Initialization and relaxation BCN = 2
      PREAMB = 1.013E+6      ! Initialization and relaxation BCN = 2
      rhoamb = rhoamb*(dr0*dz0)/(dr*dz)
!     preamb = preamb*(dr0/dr)
      VELAMB = 0.0           ! Initialization and relaxation BCN = 2
      RHO_IN = 100.0*RHOAMB  ! Initialization and relaxation BC1 = 2
      PRE_IN = 1000.0*PREAMB ! Initialization and relaxation BC1 = 2
      VEL_IN = 0.0           ! Initialization and relaxation BC1 = 2
      RELAX =  0.002 ! Relaxation rate, used when BC1 or BCN = 2
      GAMMAM = GAMMA0 - 1.0

c  Determine the cell interface locations, here a uniform grid . . .
C-----------------------------------------------------------------------
      NRP = NR + 1
      nxp = nx + 1
      NZP = NZ + 1
      Do 100 I = 1, nxp
 100         X(I) = DX*FLOAT(I-1)
      Do 110 J = 1, NZP
 110         Z(J) = DZ*FLOAT(J-1)

c  Fill the arrays with air at STP and behind the diaphragm increase the
c  density by 100 to 1 and the pressure by 1000 to 1 . . .
C-----------------------------------------------------------------------
      Do 200 J = 1, NZ
      Do 200 I = 1, nx
         val(1,I,J) = RHOAMB
         val(2,I,J) = 0.0
         val(3,I,J) = 0.0
 200     val(4,I,J) = PREAMB/GAMMAM

      Do  J = 1, JCTOP/2
      Do  I = -(icin-1) + nr, ICIN - 1 + nr
         val(4,I,J) = PRE_IN/GAMMAM
         val(1,I,J) = RHO_IN
      end do
      end do    

c  Mark the unused cells inside the cylindrical 'barrel' so they will
c  show up distinctly compared to ambient values in the plots.  This
c  has no effect as the simulation does not access these values . . .
C-----------------------------------------------------------------------
C      Do  J = 1, JCTOP
C      Do  I = -icout + nr, -icin + nr
C         val(4,I,J) = 20.0*val(4,I,J)
C         val(1,I,J) = 20.0*val(1,I,J)
C      end do
C      Do  I = ICIN + nr, ICOUT + nr
C         val(4,I,J) = 20.0*val(4,I,J)
C         val(1,I,J) = 20.0*val(1,I,J)
C      end do
C      end do

c  Begin loop over the timesteps . . .
C-----------------------------------------------------------------------
      TIME = 0.0
C234567
		rho4=val(1,:,:)
C		rvr4=rvr
C		rvz4=rvz
C		erg4=erg
		write(astep,'(i4.4)') 0
		open(12,file='fast2d.'//astep//'.dat',status='unknown',form='unformatted')
!        	write(12)rho4,rvr4,rvz4,erg4
        	write(12)rho4
	        close(12)
		
	call bcfcn(icout,icin,jctop,nr,nx,nz,bc_wall,bc_outf,relax,gammam,rhoamb,preamb)

      Do 9999 ISTEP = 1, MAXSTP

c  Compute the next timestep based on a 'Courant' number COURANT . . .
C-----------------------------------------------------------------------
         VZMAX = 0.0
         Do 240 J = 1, NZ
         Do 240 I = 1, Nx
            VTYPICAL = val(4,I,J)/val(1,I,J) 
 240        VZMAX = AMAX1 ( VTYPICAL, VZMAX ) 
         VZMAX = SQRT ( VZMAX )
         DTNEW = COURANT*AMIN1(DR,DZ)/VZMAX
         DT = AMIN1 ( DTMAX, 1.25*DT, DTNEW )

c  The results are printed when required . . .
C-----------------------------------------------------------------------
         JSTEP = ISTEP - 1
         If ( MOD(JSTEP,25) .eq. 0 ) Write ( 6, 1005 ) JSTEP, TIME, DT
         If ( MOD(JSTEP,IPRINT) .eq. 0 ) Then
            Write ( LOUT, 1000 ) ISTEP, NR, NZ, TIME, DT
            If ( ISTEP .ge. 4*IPRINT ) IPRINT = 2*IPRINT
            Write ( LOUT, 1010 ) 
            Do 230 IJ = 1, NSIZE
               DIN(1) = val(1,1,IJ)/RHOAMB
               DIN(2) = 0.01*val(3,1,IJ)/val(1,1,IJ)
               DIN(3) = (GAMMAM/PREAMB)*(val(4,1,IJ) - 0.5*
     &                  (val(2,1,IJ)**2 + val(3,1,IJ)**2)/val(1,1,IJ))
               DIN(4) = val(1,IJ,NZ)/RHOAMB
               DIN(5) = 0.01*val(2,IJ,NZ)/val(1,IJ,NZ)
               DIN(6) = (GAMMAM/PREAMB)*(val(4,IJ,NZ) - 0.5*
     &                  (val(2,IJ,NZ)**2 + val(3,IJ,NZ)**2)/val(1,IJ,NZ))
               DIN(7) = val(1,NR,IJ)/RHOAMB
               DIN(8) = 0.01*val(3,NR,IJ)/val(1,NR,IJ)
               DIN(9) = (GAMMAM/PREAMB)*(val(4,NR,IJ) - 0.5*
     &                  (val(2,NR,IJ)**2 + val(3,NR,IJ)**2)/val(1,NR,IJ))
 230                     Write ( LOUT, 1011 ) IJ, ( DIN(I), I = 1, 9 ) 
            Write ( LOUT, 1011 )
         Endif

c  Integrate the fluid equations in the radial direction (indexed by I). 
c  The outer boundary condition at interface I = NR+1 is an extra- 
c  polation from the interior cell values with a slow relaxation to
c  the known distant ambient conditions . . .
C-----------------------------------------------------------------------
c$omp parallel  default (none)
c$omp+ shared  (R, NRP, IALFR, RHO, RVR, RVZ, ERG, JCTOP, ICIN, 
c$omp+          BC_AXIS, BC_WALL, DT, ICOUT, BC_OUTF, Z, NZP,
c$omp+          IALFZ)
c$omp+ private ( J, I)

         call fctblk_init

         Call RESIDIFF ( 0.999 )
         Call MAKEGRID ( X, X, 1, nxp, ialfx )
       rhobc(1)=rhoamb
       rhobc(2)=rhoamb
       prebc(1)=preamb
       prebc(2)=preamb


c  Pick up the data from the 2D arrays in the radial direction, setting
c  the temporary, compact 1D arrays for GASDYN . . .
C-----------------------------------------------------------------------
c$omp do
 	call tau_profile_start(profiler(1,2))
         Do 300 J = 1, NZ
            Do 400 I = 1, nx
               RHON(I  ) = val(1,I,J)
               rvvn(I,1) = val(2,I,J)
               rvvn(I,2) = val(3,I,J)
 400           ERGN(I  ) = val(4,I,J)


c  Integrate along the radials inside and outside the cylinder  . . .
C-----------------------------------------------------------------------
 	call tau_profile_start(profiler(1,5))
            If ( J .le. JCTOP ) Then
               Call GASDYNP(            1, -(icout+1)+nr, BC_OUTF, BC_WALL, DT, nx,nz,j)
               Call GASDYNP( -(icin-1)+nr, ICIN-1    +nr, BC_WALL, BC_WALL, DT, nx,nz,j)
               Call GASDYNP(   ICOUT+1+nr,            nx, BC_WALL, BC_OUTF, DT, nx,nz,j)

c  Integrate along the radials (indexing in I) above the cylinder
c  which reach from the axis to the outer boundary . ..
C-----------------------------------------------------------------------
            Else
               Call GASDYNP( 1, nx, bc_outf, BC_OUTF, DT , nx,nz,j)
            End If
C	call gasdynp(1,nx,bc_outf,bc_outf,dt,nx,nz,j)
 	call tau_profile_stop (profiler(1,5))

c  Put the data back into the 2D arrays in the radial direction . . .    
C-----------------------------------------------------------------------
            Do 500 I = 1, nx
               val(1,I,J) = RHON(I  )
               val(2,I,J) = rvvn(I,1)
               val(3,I,J) = rvvn(I,2)
 500           val(4,I,J) = ERGN(I  )
 300                         Continue     ! End loop integrating the NZ rows.
 	call tau_profile_stop (profiler(1,2))
c$omp end do

c  Integrate along the axials (indexing in J) which reach from the 
c  lower active J cell (1 or JCTOP+1) to the upper boundary.  The
c  upper boundary condition at interface J = NZ+1 (BCN = 2 ) is an 
c  extrapolation from the interior cell values with a slow relaxation 
c  to the known distant ambient conditions . . .
C-----------------------------------------------------------------------
         Call MAKEGRID ( Z, Z, 1, NZP, IALFZ )
       rhobc(2)=rhoamb
       prebc(2)=preamb
c$omp do 
 	call tau_profile_start(profiler(1,3))
         Do 600 I = 1, nx
            lil=i.ge.-icout+nr.and.i.le.-icin +nr
            lir=i.ge. icin +nr.and.i.le. icout+nr
            lim=i.gt.-icin +nr.and.i.lt. Icin +nr
c  Pick up the data from the 2D arrays in the axial direction, setting
c  the temporary, compact 1D arrays for GASDYN . . .
C-----------------------------------------------------------------------
	do j=1,nz
		rhon(j  )=val(1,i,j)
		rvvn(j,2)=val(2,i,j)
		rvvn(j,1)=val(3,i,j)
		ergn(j  )=val(4,i,j)
	end do

c  Integrate along the axials either from the lower solid boundary at 
c  interface J = 1 or from the top of the barrel at J = 21 for cells
c  with I = ICIN to ICOUT . . .
C-----------------------------------------------------------------------
 	call tau_profile_start(profiler(1,5))
            If ( lil.or.lir) Then
               Call GASDYN ( JCTOP+1, NZ, BC_WALL, BC_OUTF, DT)
            else if (lim) then
	       rhobc(1)=rho_in
	       prebc(1)=pre_in
               Call GASDYN (       1, NZ, BC_WALL, BC_OUTF, DT )
!              Call GASDYN (       1, NZ, BC_OUTF, BC_OUTF, DT )
            Else
	       rhobc(1)=rhoamb
	       prebc(1)=preamb
               Call GASDYN (       1, NZ, BC_OUTF, BC_OUTF, DT )
            End If 
 	call tau_profile_stop (profiler(1,5))

c  Put the data back into the 2D arrays in the axial direction . . .    
C-----------------------------------------------------------------------
	do j=1,nz
		val(1,i,j)=rhon(j  )
		val(2,i,j)=rvvn(j,2)
		val(3,i,j)=rvvn(j,1)
		val(4,i,j)=ergn(j  )
	end do
 600                         Continue     ! End loop integrating the NR columns.
 	call tau_profile_stop (profiler(1,3))
c$omp end do
C$omp end parallel

         TIME = TIME + DT

	if (mod(istep,nsave).eq.0) then
		rho4=val(1,:,:)
C		rvr4=rvr
C		rvz4=rvz
C		erg4=erg
		write(astep,'(i4.4)') istep
		open(12,file='fast2d.'//astep//'.dat',status='unknown',form='unformatted')
!        	write(12)rho4,rvr4,rvz4,erg4
        	write(12)rho4
	        close(12)

	end if
c	if (time.ge.timel) goto 1111

 9999     Continue        ! End of the timestep loop.
 1111	continue
  	call tau_profile_stop(profiler(1,2))

		rho4=val(1,:,:)
C		rvr4=rvr
C		rvz4=rvz
C		erg4=erg
		open(12,file='fast2d.dat',status='unknown',form='unformatted')
!        	write(12)rho4,rvr4,rvz4,erg4
        	write(12)rho4
	        close(12)

  	call tau_profile_stop(profiler(1,1))
  	call tau_profile_exit('Normal Termination')
        
 1000  Format ('1', /, '   LCPFCT Test # 4 - FAST2D Barrel Explosion:',
     1        '  Step =', I4, /, 5X, I3, ' x', I3, ' Uniform Grid.',
     2        '  Time=', 1PE12.4, ' and   DT =', E12.4 )
 1005   Format (' After step ', I5, '  TIME = ', 1PE12.4,
     &           ' and timestep DT = ', E12.4 )
 1010    Format (1X, /, ' Fluid variables on selected lines ',
     1           /, ' I, J', '    RHO axis VZ     PRE ',
     2                      '     RHO top  VR     PRE ',
     3                      '     RHO wall VZ     PRE ', / )
 1011     Format ( I3, 1X, 3(1X, F8.2, F8.2, F8.2) ) 

      Stop
      End

C=======================================================================
