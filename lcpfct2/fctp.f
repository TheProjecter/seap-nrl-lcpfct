C=======================================================================
Cc$mta inline
      Subroutine LCPFCTP( RHOO, RHON, I1, IN, 
     &                    SRHO1, VRHO1, SRHON, VRHON, PBC, nx, ny)

C-----------------------------------------------------------------------
c 
c     Originated: J.P. Boris         Code 4400, NRL          Feb 1987
c     Modified:  Laboratory for Computational Physics & Fluid Dynamics
c     Contact:    J.P. Boris, J.H. Gardner, A.M. Landsberg, or E.S. Oran
c
c     Description:  This routine solves generalized continuity equations
c     of the form  dRHO/dt = -div (RHO*V) + SOURCES in the user's choice
c     of Cartesian, cylindrical, or spherical coordinate systems.  A
c     facility is included to allow definition of other coordinates.
c     The grid can be Eulerian, sliding rezone, or Lagrangian and can 
c     be arbitrarily spaced.  The algorithm is a low-phase-error FCT 
c     algorithm, vectorized and optimized for a combination of speed and 
c     flexibility.  A complete description appears in the NRL Memorandum
c     Report (1992), "LCPFCT - A Flux-Corrected Transport Algorithm For 
c     Solving Generalized Continuity Equations".
c
c     Arguments:    
c     RHOO   Real Array        grid point densities at start of step   I
c     RHON   Real Array        grid point densities at end of step     O
c     I1     Integer           first grid point of integration         I
c     IN     Integer           last grid point of intergration         I
c     SRHO1  Real Array        boundary guard cell factor at cell I1+1 I
c     VRHO1  Real Array        boundary value added to guard cell I1-1 I
c     SRHON  Real Array        boundary guard cell factor at cell IN+1 I 
c     VRHON  Real Array        boundary value added to guard cell IN+1 I
c     PBC    Logical           periodic boundaries if PBC = .true.     I
c
c        In this routine the last interface at RADHN(INP) is the outer
c     boundary of the last cell indexed IN.  The first interface at
c     RADHN(I1) is the outer boundary of the integration domain before 
c     the first cell indexed I1.
c
c     Language and Limitations:  LCPFCT is a package of FORTRAN 77 sub-
c     routines written in single precision (64 bits CRAY). The parameter
c     NPT is used to establish the internal FCT array dimensions at the
c     maximum size expected.  Thus NPT = 2002 means that continuity equa-
c     tions for systems up to 200 cells long in one direction can be 
c     integrated.  Underflows can occur when the function being trans-
c     ported has a region of zeroes.  The calculations misconserve by
c     one or two bits per cycle.  Relative phase and amplitude errors
c     (for smooth functions) are typically a few percent for character-
c     istic lengths of 1 - 2 cells (wavelengths of order 10 cells).  The 
c     jump conditions for shocks are generally accurate to better than 1 
c     percent.  Common blocks are used to transmit all data between the 
c     subroutines in the LCPFCT package.
c
c     Auxiliary Subroutines:  CNVFCT, CONSERVE, COPYGRID, MAKEGRID, 
c     NEW_GRID, RESIDIFF, SET_GRID, SOURCES, VELOCITY, ZERODIFF, and 
c     ZEROFLUX.  The detailed documentation report provided (or the
c     listing below) explains the definitions and use of the arguments
c     to these other subroutines making up the LCPFCT package.  These
c     routines are not called from LCPFCT itself but are controlled by
c     calls from the user.  Subroutines MAKEGRID, VELOCITY and SOURCES 
c     in this package must first be called to set the grid geometry, 
c     velocity-dependent flux and diffusion coefficients, and external 
c     source arrays used by LCPFCT.  The other subroutines may be called 
c     to perform other functions such as to modify boundary conditions,
c     to perform special grid operations, or compute conservation sums.
c
C-----------------------------------------------------------------------
          Implicit  NONE
          Integer   I1, IN, I1P, INP, I, I1Q, INM, I1M, nx, ny
          Real      SRHO1, VRHO1, SRHON, VRHON, RHO1M, RHONP
          Real      RHOT1M, RHOTNP, RHOTD1M, RHOTDNP 
          Logical   PBC

          include 'prm.h'
          include 'fct.h'
          real bc(4,0:npt),bcx(4,0:npt,0:npt),bcv(2,0:npt,0:npt,2)
          common /fct_bc/ bc,bcx,bcv

          Real     RHOO(0:NPT),     RHON(0:NPT)

c         Real DIFFP(NPT),FLXHP(NPT),DIFFQ(NPT),FLXHQ(NPT)
          Real flxhi1,diffi1,flxhi1p,diffi1p,flxhi1q,diffi1q
          Real flxhi,diffi,flxhip,diffip,flxhin,diffin,flxhinp,diffinp
          Real flxhpi1,diffpi1,flxhpi1p,diffpi1p,flxhpi,diffpi,diffpim
          Real flxhpinm,diffpinm,flxhpinp,diffpinp,flxhpin,diffpin
          Real flxhqi1,flxhqim,flxhqil,flxhqi,flxhqinn,flxhqinm,flxhqin,flxhqinp
          Real lorhoti1,lnrhoti1,lorhoti1p,lnrhoti1p,lorhoti,lnrhoti
          Real lnrhotin,lnrhotim,lorhotin,lnrhotil,lnrhotinn,lnrhotinm
          Real rhoti1,rhotdi1,rhoti1p,rhotdi1p,rhotin,rhotdin,rhotim,rhotdim
          Real rhotinp,rhotdinp,rhoti,rhotdi,rhotinm,rhotdinm
          Real fabsi1,fsgni1,fabsi1p,fsgni1p,fabsi,fsgni,fabsim,fsgnim
          Real fabsin,fsgnin,fabsinm,fsgninm,fabsinp,fsgninp
          Real termi1,terpi1,termi1p,terpi1p,termi,terpi,termim,terpim
          Real terminp,terpinp,terminm,terpinm,termin,terpin
       

C-----------------------------------------------------------------------
          I1M = I1 - 1
          I1P = I1 + 1
	  I1Q = I1 + 2
          INP = IN + 1
	  INM = IN - 1

c     Calculate the convective and diffusive fluxes . . .
C-----------------------------------------------------------------------
C          If ( PBC ) Then
C             RHO1M = RHOO(IN)
C             RHONP = RHOO(I1)
C          Else
C             RHO1M = SRHO1*RHOO(I1) + VRHO1 
C             RHONP = SRHON*RHOO(IN) + VRHON
C          End If

C	flxhi1 =hadudth(i1 )*(rhoo(i1 )+rho1m    )
C	diffi1 =nulh   (i1 )*(rhoo(i1 )-rho1m    )
C	flxhi1p=hadudth(i1p)*(rhoo(i1p)+rhoo(i1 ))
C	diffi1p=nulh   (i1p)*(rhoo(i1p)-rhoo(i1 ))
C	flxhi1q=hadudth(i1q)*(rhoo(i1q)+rhoo(i1p))
C	diffi1q=nulh   (i1q)*(rhoo(i1q)-rhoo(i1p))

	flxhi1 =hadudth(i1 )*(bc(1,i1 )*rhoo(i1 )+bc(2,i1 )+bc(1,i1m)*rhoo(i1m)+bc(2,i1m))
	diffi1 =nulh   (i1 )*(bc(1,i1 )*rhoo(i1 )+bc(2,i1 )-bc(1,i1m)*rhoo(i1m)-bc(2,i1m))
	flxhi1p=hadudth(i1p)*(bc(1,i1p)*rhoo(i1p)+bc(2,i1p)+bc(1,i1 )*rhoo(i1 )+bc(2,i1 ))
	diffi1p=nulh   (i1p)*(bc(1,i1p)*rhoo(i1p)+bc(2,i1p)-bc(1,i1 )*rhoo(i1 )-bc(2,i1 ))
	flxhi1q=hadudth(i1q)*(bc(1,i1q)*rhoo(i1q)+bc(2,i1q)+bc(1,i1p)*rhoo(i1p)+bc(2,i1p))
	diffi1q=nulh   (i1q)*(bc(1,i1q)*rhoo(i1q)+bc(2,i1q)-bc(1,i1p)*rhoo(i1p)-bc(2,i1p))

        flxhi  =flxhi1q
        diffi  =diffi1q

	lorhoti1  =lo   (i1 )*rhoo (i1 )+source(i1 )+(flxhi1 -flxhi1p)
	lnrhoti1  =lorhoti1                         +(diffi1p-diffi1 )
	rhoti1    =lorhoti1   *rln (i1 )
	rhotdi1   =lnrhoti1   *rln (i1 )
	lorhoti1p =lo   (i1p)*rhoo(i1p)+source(i1p)+(flxhi1p-flxhi1q)
	lnrhoti1p =lorhoti1p                       +(diffi1q-diffi1p)
	rhoti1p   =lorhoti1p *rln (i1p)
	rhotdi1p  =lnrhoti1p *rln (i1p)
        lnrhotil  =lnrhoti1
        lnrhotim  =lnrhoti1p
        rhotim    =rhoti1p
        rhotdim   =rhotdi1p

C          If ( PBC ) Then
C             RHOT1M  = RHOTIN
C             RHOTD1M = RHOTDIN
C          Else
C             RHOT1M  = SRHO1*RHOTI1  + VRHO1 
C             RHOTD1M = SRHO1*RHOTDI1 + VRHO1 
C          End If
C
C	flxhpi1   =mulh(i1 )*(rhoti1  -rhot1m )
C	diffpi1   =          (rhotdi1 -rhotd1m)
C	fabsi1    =abs(       flxhpi1 )
C	fsgni1    =sign(diff1,diffpi1 )
C	flxhpi1p  =mulh(i1p)*(rhoti1p -rhoti1 )
C	diffpi1p  =          (rhotdi1p-rhotdi1)
C        diffpim   =diffpi1p

	rhot1m =rhoti1
	rhotd1m=rhotdi1
	flxhpi1   =mulh(i1 )*(bc(1,i1 )*rhoti1 +bc(2,i1p)-bc(1,i1m)*rhot1m -bc(2,i1m))
	diffpi1   =          (bc(1,i1 )*rhotdi1+bc(2,i1p)-bc(1,i1m)*rhotd1m-bc(2,i1m))
	fabsi1    =abs(       flxhpi1 )
	fsgni1    =sign(diff1,diffpi1 )
	flxhpi1p  =mulh(i1p)*(bc(1,i1p)*rhoti1p +bc(2,i1p)-bc(1,i1 )*rhoti1 -bc(2,i1 ))
	diffpi1p  =          (bc(1,i1p)*rhotdi1p+bc(2,i1p)-bc(1,i1 )*rhotdi1-bc(2,i1 ))
        diffpim   =diffpi1p

 	fabsi1p=abs(       flxhpi1p)
 	fsgni1p=sign(diff1,diffpi1p)
 	termi1p=fsgni1p*ln(i1 )*diffpi1 +bc(3,i1p)*bignum
 	terpi1 =fsgni1 *ln(i1 )*diffpi1p+bc(4,i1 )*bignum
	fabsim =fabsi1p
	fsgnim =fsgni1p
	termim =termi1p
	termi1 =bc(3,i1)*bignum

C          If ( PBC ) Then
C             TERMI1 = TERMINP
C          Else
C             TERMI1 = BIGNUM
C          End If

	flxhqi1=fsgni1 *amax1(0.0,amin1(termi1,fabsi1,terpi1))
        flxhqil=flxhqi1

	do i=i1q,inm
		flxhip  =hadudth(i+1)*(bc(1,i+1)*rhoo(i+1)+bc(1,i  )*rhoo(i  )+bc(2,i  ))
		diffip  =nulh   (i+1)*(bc(1,i+1)*rhoo(i+1)-bc(1,i  )*rhoo(i  )-bc(2,i  ))
		lorhoti =lo   (i)*rhoo(i)+source(i)+(flxhi    -flxhip   )
		lnrhoti =lorhoti                   +(diffip   -diffi    )
		rhoti   =lorhoti *rln (i)
		rhotdi  =lnrhoti *rln (i)
		flxhpi   =mulh(i)*(bc(1,i  )*rhoti +bc(2,i  )-bc(1,i-1)*rhotim -bc(2,i-1))
		diffpi   =        (bc(1,i  )*rhotdi+bc(2,i  )-bc(1,i-1)*rhotdim-bc(2,i-1))
 		fabsi  =abs(       flxhpi  )
 		fsgni  =sign(diff1,diffpi  )
 		termi  =fsgni *ln(i-1)*diffpim+bc(3,i  )*bignum
 		terpim =fsgnim*ln(i-1)*diffpi +bc(4,i-1)*bignum
		flxhqim  =fsgnim *amax1(0.0,amin1(termim,fabsim,terpim))
		rhon  (i-2)=rln (i-2)*(lnrhotil+(flxhqil-flxhqim))
		source(i-2)=0.0
                flxhi   =flxhip
                diffi   =diffip
                diffpim =diffpi
                flxhqil =flxhqim
                lnrhotil=lnrhotim
                lnrhotim=lnrhoti
                rhotim  =rhoti
                rhotdim =rhotdi
                fabsim  =fabsi
                fsgnim  =fsgni
		termim  =termi
	end do

        flxhin=flxhi
        diffin=diffi
C	flxhinp=nulh   (inp)*(rhonp+rhoo(in))
C	diffinp=hadudth(inp)*(rhonp-rhoo(in))

	flxhinp=nulh   (inp)*(bc(1,inp)*rhoo(inp)+bc(2,inp)+bc(1,in)*rhoo(in)+bc(2,in))
	diffinp=hadudth(inp)*(bc(1,inp)*rhoo(inp)+bc(2,inp)-bc(1,in)*rhoo(in)-bc(2,in))

	lorhotin =lo   (in)*rhoo(in)+source(in)+(flxhin -flxhinp)
	lnrhotin =lorhotin                     +(diffinp-diffin )
	rhotin =lorhotin *rln (in)
	rhotdin=lnrhotin *rln (in)

C          If ( PBC ) Then
C             RHOTNP  = RHOTI1
C             RHOTDNP = RHOTDI1
C          Else
C             RHOTNP  = SRHON*RHOTIN  + VRHON
C             RHOTDNP = SRHON*RHOTDIN + VRHON
C          End If
C
C        rhotinm =rhotim
C        rhotdinm=rhotdim
C        diffpinm=diffpim
C	flxhpin =mulh(in )*(rhotin -rhotinm )
C	diffpin =          (rhotdin-rhotdinm)
C	flxhpinp=mulh(inp)*(rhotnp -rhotin  )
C	diffpinp=          (rhotdnp-rhotdin )

        rhotinm =rhotim
        rhotdinm=rhotdim
        diffpinm=diffpim
	rhotnp=rhotin
	rhotdnp=rhotdin
	flxhpin =mulh(in )*(bc(1,in )*rhotin +bc(2,in )-bc(1,inm)*rhotinm -bc(2,inm))
	diffpin =          (bc(1,in )*rhotdin+bc(2,in )-bc(1,inm)*rhotdinm-bc(2,inm))
	flxhpinp=mulh(inp)*(bc(1,inp)*rhotnp +bc(2,inp)-bc(1,in )*rhotin  -bc(2,in ))
	diffpinp=          (bc(1,inp)*rhotdnp+bc(2,inp)-bc(1,in )*rhotdin -bc(2,in ))

	fabsinm=fabsim
	fsgninm=fsgnim
	terminm=termim
	
 	fabsin =abs(       flxhpin )
 	fsgnin =sign(diff1,diffpin )
 	termin =fsgnin *ln(inm)*diffpinm+bc(3,in )*bignum
 	terpinm=fsgninm*ln(inm)*diffpin +bc(4,inm)*bignum
 	fabsinp=abs(       flxhpinp)
 	fsgninp=sign(diff1,diffpinp)
 	terminp=fsgninp*ln(in )*diffpin +bc(3,inp)*bignum
 	terpin =fsgnin *ln(in )*diffpinp+bc(4,in )*bignum

C          If ( PBC ) Then
C             TERPINP = TERPI1
C          Else
C             TERPINP = BIGNUM
C          End If

        lnrhotinn=lnrhotil
        lnrhotinm=lnrhotim
        flxhqinn=flxhqil
	flxhqinm=fsgninm *amax1(0.0,amin1(terminm,fabsinm,terpinm))
	rhon  (in-2)=rln (in-2)*(lnrhotinn+(flxhqinn-flxhqinm))
	source(in-2)=0.0
	flxhqin =fsgnin  *amax1(0.0,amin1(termin ,fabsin ,terpin ))
	rhon  (inm )=rln (inm )*(lnrhotinm+(flxhqinm-flxhqin ))
	source(inm )=0.0
	flxhqinp=fsgninp *amax1(0.0,amin1(terminp,fabsinp,terpinp))
	rhon  (in  )=rln (in  )*(lnrhotin +(flxhqin -flxhqinp))
	source(in  )=0.0

      Return
      End
