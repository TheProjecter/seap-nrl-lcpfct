/*
 * converted almost verbatim from lcpfct.f
 * 		interger -> int
 * 		real -> float
 * 		logical -> boolean
 * 		var(size) -> var[size]
 * 		parameter -> #define
 * 		common -> struct
 * all variable names are converted into lower case letters (except constants),
 * but underscores and other stuff are kept the same so whoever wrote this originally
 * can still read this converted script
 *
 * from here on out, I will religiously copy all comments in code
 * for better comparability between original code and my translated
 * version ... comments in original code use //
 *
 * 					-- mabel xu, jun 25 2007
 */

/* #include <stdio.h>		/* for printing debugging messages */
#include <stdbool.h>	/* support on UNIX for bool vars */
#include <math.h>		/* support for fabs() */

//--------------------------------------------------------------------
//
//     Originated: J.P. Boris         Code 4400, NRL          Feb 1987
//     Modified:  Laboratory for Computational Physics & Fluid Dynamics
//     Contact:    J.P. Boris, J.H. Gardner, A.M. Landsberg, or E.S. Oran
//
//     Description:  This routine solves generalized continuity equations
//     of the form  dRHO/dt = -div (RHO*V) + SOURCES in the user's choice
//     of Cartesian, cylindrical, or spherical coordinate systems.  A
//     facility is included to allow definition of other coordinates.
//     The grid can be Eulerian, sliding rezone, or Lagrangian and can
//     be arbitrarily spaced.  The algorithm is a low-phase-error FCT
//     algorithm, vectorized and optimized for a combination of speed and
//     flexibility.  A complete description appears in the NRL Memorandum
//     Report (1992), "LCPFCT - A Flux-Corrected Transport Algorithm For
//     Solving Generalized Continuity Equations".
//
//     Arguments:
//     RHOO   Real Array        grid point densities at start of step   I
//     RHON   Real Array        grid point densities at end of step     O
//     I1     Integer           first grid point of integration         I
//     IN     Integer           last grid point of intergration         I
//     SRHO1  Real Array        boundary guard cell factor at cell I1+1 I
//     VRHO1  Real Array        boundary value added to guard cell I1-1 I
//     SRHON  Real Array        boundary guard cell factor at cell IN+1 I
//     VRHON  Real Array        boundary value added to guard cell IN+1 I
//     PBC    Logical           periodic boundaries if PBC = .true.     I
//
//       In this routine the last interface at RADHN(INP) is the outer
//     boundary of the last cell indexed IN.  The first interface at
//     RADHN(I1) is the outer boundary of the integration domain before
//     the first cell indexed I1.
//
//     Language and Limitations:  LCPFCT is a package of FORTRAN 77 sub-
//     routines written in single precision (64 bits CRAY). The parameter
//     NPT is used to establish the internal FCT array dimensions at the
//     maximum size expected.  Thus NPT = 2002 means that continuity equa-
//     tions for systems up to 200 cells long in one direction can be
//     integrated.  Underflows can occur when the function being trans-
//     ported has a region of zeroes.  The calculations misconserve by
//     one or two bits per cycle.  Relative phase and amplitude errors
//     (for smooth functions) are typically a few percent for character-
//     istic lengths of 1 - 2 cells (wavelengths of order 10 cells).  The
//     jump conditions for shocks are generally accurate to better than 1
//     percent.  Common blocks are used to transmit all data between the
//     subroutines in the LCPFCT package.
//
//     Auxiliary Subroutines:  CNVFCT, CONSERVE, COPYGRID, MAKEGRID,
//     NEW_GRID, RESIDIFF, SET_GRID, SOURCES, VELOCITY, ZERODIFF, and
//     ZEROFLUX.  The detailed documentation report provided (or the
//     listing below) explains the definitions and use of the arguments
//     to these other subroutines making up the LCPFCT package.  These
//     routines are not called from LCPFCT itself but are controlled by
//     calls from the user.  Subroutines MAKEGRID, VELOCITY and SOURCES
//     in this package must first be called to set the grid geometry,
//     velocity-dependent flux and diffusion coefficients, and external
//     source arrays used by LCPFCT.  The other subroutines may be called
//     to perform other functions such as to modify boundary conditions,
//     to perform special grid operations, or compute conservation sums.
//
//--------------------------------------------------------------------

/* for compatability between c and fortran: */
/* c arrays start at 0, fortran arrays start at 1 */
#define I1 i1 - 1
#define IN in - 1

/* here be extra functions that mimic intrinsic functions in fortran 77 */
float float_abs(float a)	/* regular abs() deals with ints, fabs() deals with doubles */
{
	return (float)fabs((double)a);
}
float amax1 (float a, float b)		/* this particular "amax1" only supports 2 variables! woot! */
{
	if (a >= b)
		return a;
	else
		return b;
}
float amin1 (float a, float b, float c)		/* this particular "amin1" only supports 3 variables! woot! */
{
	if (a <= b && a <= c)
		return a;
	else if (b <= a && b <= c)
		return b;
	else
		return c;
}
float sign (float a, float b)	/* does whatever sign does ... I don't pretend to understand this at all */
{
	if (b >= 0)
		return float_abs(a);
	else
		return -float_abs(a);
}

#define NPT 2002
#define BIGNUM 1e38f
//	BIGNUM = Machine Dependent Largest Number - Set By The User!!!!

//	/FCT_SCRH/ Holds scratch arrays for use by LCPFCT and CNVFCT
struct
{
	float scrh[NPT];
	float scr1[NPT];
	float diff[NPT];
	float flxh[NPT];
	float fabs[NPT];
	float fsgn[NPT];
	float term[NPT];
	float terp[NPT];
	float lnrhot[NPT];
	float lorhot[NPT];
	float rhot[NPT];
	float rhotd[NPT];
} fct_scrh__;
//	$OMP THREADPRIVATE(/FCT_SCRH/)

//	/FCT_GRID/ Holds geometry, grid, area and volume information
struct
{
	float lo[NPT];
	float ln[NPT];
	float ah[NPT];
	float rln[NPT];
	float lh[NPT];
	float rlh[NPT];
	float roh[NPT];
	float rnh[NPT];
	float adugth[NPT];
} fct_grid__;
//	$OMP THREADPRIVATE(/FCT_GRID/)

//	/FCT_VELO/ Holds velocity-dependent flux coefficients
struct
{
	float hadudth[NPT];
	float nulh[NPT];
	float mulh[NPT];
	float epsh[NPT];
	float vdtodr[NPT];
} fct_velo__;
//	$OMP THREADPRIVATE( /FCT_VELO/)

//	/FCT_MISC/ Holds the source array and diffusion coefficient
struct
{
	float source[NPT];
	float diff1_;
} fct_misc__;
//	$OMP THREADPRIVATE( /FCT_MISC/)

void lcpfct_ (float *rhoo, float *rhon, int* i1, int* in,
				float* srho1, float* vrho1, float* srhon,
				float* vrhon, bool* pbc)
{
	/* passed vars and constants went here originally */

	static int i1p, inp, i;
	static float rho1m, rhonp, rhot1m, rhotnp, rhotd1m, rhotdnp;

	/* COMMONs originally went here (now structs) */

//-----------------------------------------------------------------------
	i1p = *I1 + 1;
	inp = *IN + 1;

//	Calculate the convective and diffusive fluxes . . .
//-----------------------------------------------------------------------
	if ( *pbc )
	{
		rho1m = rhoo[*IN];
		rhonp = rhoo[*I1];
	}
	else
	{
		rho1m = *srho1 * rhoo[*I1] + *vrho1;
		rhonp = *srhon * rhoo[*IN] + *vrhon;
	}

	fct_scrh__.diff[*I1] = fct_velo__.nulh[*I1] * ( rhoo[*I1] - rho1m );
	fct_scrh__.flxh[*I1] = fct_velo__.hadudth[*I1] * ( rhoo[*I1] + rho1m );

	for (i = i1p; i <= *IN; i++)
	{
		fct_scrh__.flxh[i] = fct_velo__.hadudth[i] * ( rhoo[i] + rhoo[i-1] );
		fct_scrh__.diff[i] = fct_velo__.nulh[i] * ( rhoo[i] - rhoo[i-1] );
	}

	fct_scrh__.diff[inp] = fct_velo__.nulh[inp] * ( rhonp - rhoo[*IN] );
	fct_scrh__.flxh[inp] = fct_velo__.hadudth[inp] * ( rhonp + rhoo[*IN] );

//	Calculate LORHOT, the transported mass elements, and LNRHOT, the
//	transported & diffused mass elements . . .
//-----------------------------------------------------------------------
	for (i = *I1; i <= *IN; i++)
	{
		fct_scrh__.lorhot[i] = fct_grid__.lo[i] * rhoo[i] + fct_misc__.source[i] + (fct_scrh__.flxh[i] - fct_scrh__.flxh[i+1]);
		fct_scrh__.lnrhot[i] = fct_scrh__.lorhot[i] + (fct_scrh__.diff[i+1] - fct_scrh__.diff[i]);
		fct_scrh__.rhot[i] = fct_scrh__.lorhot[i] * fct_grid__.rln[i];
		fct_scrh__.rhotd[i] = fct_scrh__.lnrhot[i] * fct_grid__.rln[i];
	}

//	Evaluate the boundary conditions for RHOT and RHOTD . . .
//-----------------------------------------------------------------------
	if ( *pbc )
	{
		rhot1m = fct_scrh__.rhot[*IN];
		rhotnp = fct_scrh__.rhot[*I1];
		rhotd1m = fct_scrh__.rhotd[*IN];
		rhotdnp = fct_scrh__.rhotd[*I1];
	}
	else
	{
		rhot1m = *srho1 * fct_scrh__.rhot[*I1] + *vrho1;
		rhotnp = *srhon * fct_scrh__.rhot[*IN] + *vrhon;
		rhotd1m = *srho1 * fct_scrh__.rhotd[*I1] + *vrho1;
		rhotdnp = *srhon * fct_scrh__.rhotd[*IN] + *vrhon;
	}

//	Calculate the transported antiduffusive fluxes and transported
//	and diffused density differences . . .
//-----------------------------------------------------------------------
	fct_scrh__.flxh[*I1] = fct_velo__.mulh[*I1] * ( fct_scrh__.rhot[*I1] - rhot1m );
	fct_scrh__.diff[*I1] = fct_scrh__.rhotd[*I1] - rhotd1m;
	fct_scrh__.fabs[*I1] = float_abs ( fct_scrh__.flxh[*I1] );
	fct_scrh__.fsgn[*I1] = sign ( fct_misc__.diff1_, fct_scrh__.diff[*I1] );

	for (i = i1p; i <= *IN; i++)
	{
		fct_scrh__.flxh[i] = fct_velo__.mulh[i] * ( fct_scrh__.rhot[i] - fct_scrh__.rhot[i-1] );
		fct_scrh__.diff[i] = fct_scrh__.rhotd[i] - fct_scrh__.rhotd[i-1];
	}

	fct_scrh__.flxh[inp] = fct_velo__.mulh[inp] * ( rhotnp - fct_scrh__.rhot[*IN] );
	fct_scrh__.diff[inp] = rhotdnp - fct_scrh__.rhotd[*IN];

//	Calculate the magnitude & sign of the antidiffusive flux followed
//	by the flux-limiting changes on the right and left . . .
//-----------------------------------------------------------------------
	for (i = *I1; i <= *IN; i++)
	{
		fct_scrh__.fabs[i+1] = float_abs ( fct_scrh__.flxh[i+1] );
		fct_scrh__.fsgn[i+1] = sign ( fct_misc__.diff1_, fct_scrh__.diff[i+1] );
		fct_scrh__.term[i+1] = fct_scrh__.fsgn[i+1] * fct_grid__.ln[i] * fct_scrh__.diff[i];
		fct_scrh__.terp[i] = fct_scrh__.fsgn[i] * fct_grid__.ln[i] * fct_scrh__.diff[i+1];
	}

	if ( *pbc )
	{
		fct_scrh__.terp[inp] = fct_scrh__.terp[*I1];
		fct_scrh__.term[*I1] = fct_scrh__.term[inp];
	}
	else
	{
		fct_scrh__.terp[inp] = BIGNUM;
		fct_scrh__.term[*I1] = BIGNUM;
	}

//	Correct the transported fluxes completely and then calculate the
//	new Flux-Corrected Transport densities . . .
//-----------------------------------------------------------------------
	fct_scrh__.flxh[*I1] = fct_scrh__.fsgn[*I1] * amax1 ( 0.0, amin1 ( fct_scrh__.term[*I1], fct_scrh__.fabs[*I1], fct_scrh__.terp[*I1] ) );

	for (i = *I1; i <= *IN; i++)
	{
		fct_scrh__.flxh[i+1] = fct_scrh__.fsgn[i+1] * amax1 ( 0.0, amin1 ( fct_scrh__.term[i+1], fct_scrh__.fabs[i+1], fct_scrh__.terp[i+1] ) );
		rhon[i] = fct_grid__.rln[i] * ( fct_scrh__.lnrhot[i] + (fct_scrh__.flxh[i] - fct_scrh__.flxh[i+1]) );
		fct_misc__.source[i] = 0.0;
	}
}
