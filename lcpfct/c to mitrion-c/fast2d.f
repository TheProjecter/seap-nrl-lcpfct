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

         Integer   NPT, I, J, IJ, NSIZE,nsave,nr0,nz0
         real      dr0,dz0
         Parameter ( NPT = 2002, NSIZE = 1025, nsave=1000, nr0=129, 
     &    		nz0=129)
         parameter ( dr0=0.1,dz0=0.1 )
         Integer   NR, NRP, IALFR,     BC_AXIS, BC_WALL, BC_OUTF
         Integer   NZ, NZP, IALFZ,     LOUT, MAXSTP, IPRINT
         Parameter (NR = NSIZE, NZ = NSIZE)
         Integer   ICIN, ICOUT, JCTOP,     JSTEP, ISTEP
        integer icount,irate,imax,jcount
        real twall,twall0,tcpu,ctime
         Real      DR,         DZ,         DT,         TIME
         Real      COURANT,    DTNEW,      VTYPICAL,   RELAX
         Real      DTMAX,      VZMAX,      R(NPT),     Z(NPT)
         Real RHO(NR,NZ), RVR(NR,NZ), RVZ(NR,NZ), ERG(NR,NZ)
         Real*4 RHO4(NR,NZ), RVR4(NR,NZ), RVZ4(NR,NZ), ERG4(NR,NZ)
         Real      DIN(9)

         Real      RHO_IN,     PRE_IN,     VEL_IN,     GAMMA0
         Real      RHOAMB,     PREAMB,     VELAMB,     GAMMAM
         Real      RHON(NPT),  RVRN(NPT),  RVTN(NPT),  ERGN(NPT)
        character(4) astep
       COMMON /ARRAYS/RHON,RVRN,RVTN,ERGN,RELAX
c$OMP THREADPRIVATE(/ARRAYS/)
C      TASKCOMMON /ARRAYS/RHON,RVRN,RVTN,ERGN,RELAX
      COMMON  /INVARIANTS/
     &                         RHO_IN, PRE_IN, VEL_IN, GAMMA0,
     &                         RHOAMB, PREAMB, VELAMB, GAMMAM


c       integer profiler(2,1)

c  The 2D barrel explosion program control parameters are specified.
c  (Change here to run other cases) . . .
C-----------------------------------------------------------------------
c      NR      =  129  ! Number of cells in the first (radial R) direction
      IALFR   =   2  ! Sets cylindrical coordinates in the R direction
      DR      =  0.1 ! Cell size (e.g., cm) in the radial direction
c      NZ      =  129  ! Number of cells in the second (axial Z) direction
      IALFZ   =   1  ! Sets Cartesian coordinates in the Z direction
      DZ      =  0.1 ! Cell size (e.g., cm) in the axial direction
      LOUT    =   6  ! Logical unit number of printed output device
      BC_AXIS =   1  ! Cylindrical axis set as an impermeable wall
      BC_OUTF =   2  ! Outer boundaries set as extrapolative outflow
      BC_WALL =   1  ! Walls of the barrel set as an ideal solid wall
c      MAXSTP  = 801  ! Maximum number of timesteps of length DT
      MAXSTP  =  26  ! Maximum number of timesteps of length DT
c     IPRINT  = 100  ! Initial frequency of validation printout results
      IPRINT  =  25  ! Initial frequency of validation printout results
      COURANT =  0.4 ! Approximate maximum Cournat number allowed
      DTMAX = 2.0E-7 ! Maximum timestep allowed in the computation
      DT    = 1.0E-9 ! Initial (small guess) for starting timestep

       DR=(DR0*NR0)/NR
       DZ=(DZ0*NZ0)/NZ
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
      jctop = jctop*float(nz)/nz0
      GAMMA0 = 1.4   ! Gas constant (adiabatic index)
      RHOAMB = 0.00129       ! Initialization and relaxation BCN = 2
      PREAMB = 1.013E+6      ! Initialization and relaxation BCN = 2
      rhoamb = rhoamb*(dr0*dz0)/(dr*dz)
!      preamb = preamb*(dr0/dr)
      VELAMB = 0.0           ! Initialization and relaxation BCN = 2
      RHO_IN = 100.0*RHOAMB  ! Initialization and relaxation BC1 = 2
      PRE_IN = 1000.0*PREAMB ! Initialization and relaxation BC1 = 2
      VEL_IN = 0.0           ! Initialization and relaxation BC1 = 2
      RELAX =  0.002 ! Relaxation rate, used when BC1 or BCN = 2
      GAMMAM = GAMMA0 - 1.0

c  Determine the cell interface locations, here a uniform grid . . .
C-----------------------------------------------------------------------
      NRP = NR + 1
      NZP = NZ + 1
      Do 100 I = 1, NRP
 100         R(I) = DR*FLOAT(I-1)
      Do 110 J = 1, NZP
 110         Z(J) = DZ*FLOAT(J-1)

c  Fill the arrays with air at STP and behind the diaphragm increase the
c  density by 100 to 1 and the pressure by 1000 to 1 . . .
C-----------------------------------------------------------------------
      Do 200 J = 1, NZ
      Do 200 I = 1, NR
         RHO(I,J) = RHOAMB
         RVR(I,J) = 0.0
         RVZ(I,J) = 0.0
 200         ERG(I,J) = PREAMB/GAMMAM
  
      Do 210 J = 1, JCTOP/2
      Do 210 I = 1, ICIN - 1
         ERG(I,J) = PRE_IN/GAMMAM
 210         RHO(I,J) = RHO_IN

c  Mark the unused cells inside the cylindrical 'barrel' so they will
c  show up distinctly compared to ambient values in the plots.  This
c  has no effect as the simulation does not access these values . . .
C-----------------------------------------------------------------------
      Do 220 J = 1, JCTOP
      Do 220 I = ICIN, ICOUT
         ERG(I,J) = 20.0*ERG(I,J)
 220         RHO(I,J) = 20.0*RHO(I,J)

c  Begin loop over the timesteps . . .
C-----------------------------------------------------------------------
      TIME = 0.0
C234567
c       call tau_profile_init()
c       call tau_profile_set_node(0)
c       call tau_profile_timer(profiler(1,1),'loop')
c       call tau_profile_start(profiler(1,1))
      Do 9999 ISTEP = 1, MAXSTP

c  Compute the next timestep based on a 'Courant' number COURANT . . .
C-----------------------------------------------------------------------
         VZMAX = 0.0
         Do 240 J = 1, NZ
         Do 240 I = 1, NR
            VTYPICAL = ERG(I,J)/RHO(I,J)
 240               VZMAX = AMAX1 ( VTYPICAL, VZMAX )
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
               DIN(1) = RHO(1,IJ)/RHOAMB
               DIN(2) = 0.01*RVZ(1,IJ)/RHO(1,IJ)
               DIN(3) = (GAMMAM/PREAMB)*(ERG(1,IJ) - 0.5*
     &                  (RVR(1,IJ)**2 + RVZ(1,IJ)**2)/RHO(1,IJ))
               DIN(4) = RHO(IJ,NZ)/RHOAMB
               DIN(5) = 0.01*RVR(IJ,NZ)/RHO(IJ,NZ)
               DIN(6) = (GAMMAM/PREAMB)*(ERG(IJ,NZ) - 0.5*
     &                  (RVR(IJ,NZ)**2 + RVZ(IJ,NZ)**2)/RHO(IJ,NZ))
               DIN(7) = RHO(NR,IJ)/RHOAMB
               DIN(8) = 0.01*RVZ(NR,IJ)/RHO(NR,IJ)
               DIN(9) = (GAMMAM/PREAMB)*(ERG(NR,IJ) - 0.5*
     &                  (RVR(NR,IJ)**2 + RVZ(NR,IJ)**2)/RHO(NR,IJ))
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
         Call MAKEGRID ( R, R, 1, NRP, IALFR )

c  Pick up the data from the 2D arrays in the radial direction, setting
c  the temporary, compact 1D arrays for GASDYN . . .
C-----------------------------------------------------------------------
c$omp do
         Do 300 J = 1, NZ
            Do 400 I = 1, NR
               RHON(I) = RHO(I,J)
               RVRN(I) = RVR(I,J)
               RVTN(I) = RVZ(I,J)
 400                     ERGN(I) = ERG(I,J)

c  Integrate along the radials inside and outside the cylinder  . . .
C-----------------------------------------------------------------------
            If ( J .le. JCTOP ) Then
               Call GASDYN ( 1, ICIN-1, BC_AXIS, BC_WALL, DT )
               Call GASDYN ( ICOUT+1, NR, BC_WALL, BC_OUTF, DT)

c  Integrate along the radials (indexing in I) above the cylinder
c  which reach from the axis to the outer boundary . ..
C-----------------------------------------------------------------------
            Else
               Call GASDYN ( 1, NR, BC_AXIS, BC_OUTF, DT )
            End If

c  Put the data back into the 2D arrays in the radial direction . . .
C-----------------------------------------------------------------------
            Do 500 I = 1, NR
               RHO(I,J) = RHON(I)
               RVR(I,J) = RVRN(I)
               RVZ(I,J) = RVTN(I)
 500                     ERG(I,J) = ERGN(I)
 300                         Continue     ! End loop integrating the NZ rows.
c$omp end do

c  Integrate along the axials (indexing in J) which reach from the
c  lower active J cell (1 or JCTOP+1) to the upper boundary.  The
c  upper boundary condition at interface J = NZ+1 (BCN = 2 ) is an
c  extrapolation from the interior cell values with a slow relaxation
c  to the known distant ambient conditions . . .
C-----------------------------------------------------------------------
         Call MAKEGRID ( Z, Z, 1, NZP, IALFZ )
c$omp do
         Do 600 I = 1, NR

c  Pick up the data from the 2D arrays in the axial direction, setting
c  the temporary, compact 1D arrays for GASDYN . . .
C-----------------------------------------------------------------------
            Do 700 J = 1, NZ
               RHON(J) = RHO(I,J)
               RVTN(J) = RVR(I,J)
               RVRN(J) = RVZ(I,J)
 700                     ERGN(J) = ERG(I,J)

c  Integrate along the axials either from the lower solid boundary at
c  interface J = 1 or from the top of the barrel at J = 21 for cells
c  with I = ICIN to ICOUT . . .
C-----------------------------------------------------------------------
            If ( I.ge.ICIN .and. I.le.ICOUT ) Then
               Call GASDYN ( JCTOP+1, NZ, BC_WALL, BC_OUTF, DT)
            Else
               Call GASDYN ( 1, NZ, BC_WALL, BC_OUTF, DT )
            End If

c  Put the data back into the 2D arrays in the axial direction . . .
C-----------------------------------------------------------------------
            Do 800 J = 1, NZ
               RHO(I,J) = RHON(J)
               RVR(I,J) = RVTN(J)
               RVZ(I,J) = RVRN(J)
 800                     ERG(I,J) = ERGN(J)
 600                         Continue     ! End loop integrating the NR columns.
c$omp end do
C$omp end parallel

         TIME = TIME + DT

       if (mod(istep,nsave).eq.0) then
       do j=1,nz
       do i=1,nr
               rho4(i,j)=rho(i,j)
               rvr4(i,j)=rvr(i,j)
               rvz4(i,j)=rvz(i,j)
               erg4(i,j)=erg(i,j)
       end do
       end do
               write(astep,'(i4.4)') istep
               open(12,file='fast2d.'//astep//'.dat',
     &         		status='unknown',form='unformatted')
               write(12)rho4,rvr4,rvz4,erg4
               close(12)
       end if

 9999     Continue        ! End of the timestep loop.
c       call tau_profile_stop(profiler(1,1))
c       call tau_profile_exit('Normal Termination')

       do j=1,nz
       do i=1,nr
	       rho4(i,j)=rho(i,j)
	       rvr4(i,j)=rvr(i,j)
	       rvz4(i,j)=rvz(i,j)
	       erg4(i,j)=erg(i,j)
       end do
       end do

       open(12,file='fast2d.dat',status='unknown',form='unformatted')
        write(12)rho4,rvr4,rvz4,erg4
        close(12)

c  Added section to write a .general file
c-----------------------------------------------------------------------
       open(12,file='fast2d.general',status='unknown',form='formatted')
        write(12,*)'file = C:\\fct_test\\fast2d.mxu.dat'
        write(12,*)'grid = ',nsize,' x',nsize
        write(12,*)'format = lsb ieee'
        write(12,*)'interleaving = record'
        write(12,*)'field = rho, velocity, erg'
        write(12,*)'structure = scalar, 2-vector, scalar'
        write(12,*)'type = float, float, float'
        write(12,*)'dependency = positions, positions, positions'
        write(12,*)'positions = regular, regular, 0, 1, 0, 1'
        write(12,*)'end'
        close(12)

 1000  Format ('1', /, '   LCPFCT Test # 4 - FAST2D Barrel Explosion:',
     1        '  Step =', I4, /, 5X, I3, ' x ', I3, ' Uniform Grid.',
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