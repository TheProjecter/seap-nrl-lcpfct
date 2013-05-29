/* helper c methods */
float float_abs(float a)	/* regular abs() deals with ints, fabs() deals with doubles */
{
	return (float)fabs((double)a);
}

float amax1(float a, float b)		/* this particular "amax1" only supports 2 variables */
{
	if (a >= b)
		return a;
	else
		return b;
}

float amin1(float a, float b, float c)		/* this particular "amin1" only supports 3 variables */
{
	if (a <= b && a <= c)
		return a;
	else if (b <= a && b <= c)
		return b;
	else
		return c;
}

float sign(float a, float b)	/* adds sign of second number to the first number */
{
	if (b >= 0)
		return float_abs(a);
	else
		return -float_abs(a);
}

/* for compatability between c and fortran: */
/* c arrays start at 0, fortran arrays start at 1 */
#define I1 i1 - 1
#define IN in - 1

#include <stdbool.h>	// handle booleans
#include <math.h>		// support for fabs()

/* constants used in the program */
#define NPT 2002
#define BIGNUM 1e38f
//	BIGNUM = Machine Dependent Largest Number - Set By The User!!!!

/* only some of the structs are included here since not all are needed anymore */
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

struct
{
	float bc[4][NPT + 1];
	float bcx[4][NPT + 1][NPT + 1];
	float bcv[2][NPT + 1][NPT + 1][2];
} fct_bc__;

/* main lcpfct single loop algorithm */
void lcpfctp_ (float *rhoo, float *rhon, int* i1, int* in,
				float *srho1, float *vrho1, float *srhon, float *vrhon, bool pbc, int nx, int ny)
{
	static int i1p, i1q, i1m;
	static int inp, inq, inm;
	static float flxhi, flxhi1, flxhi1p, flxhi1q, flxhip, flxhpi, flxhpi1, flxhin, flxhinp, flxhpi1p, flxhpim, flxhpin, flxhpinm, flxhpinp, flxhqi1, flxhqil, flxhqim, flxhqin, flxhqinm, flxhqinn, flxhqinp;
	static float diffi, diffi1, diffi1p, diffi1q, diffip, diffpi, diffpi1, diffin, diffinp, diffpi1p, diffpim, diffpin, diffpinm, diffpinp, diffqi1, diffqil, diffqim, diffqin, diffqinm, diffqinn, diffqinp;
	static float lorhoti, lorhoti1, lorhoti1p, lorhotil, lorhotim, lorhotin, lorhotinm, lorhotinn;
	static float lnrhoti, lnrhoti1, lnrhoti1p, lnrhotil, lnrhotim, lnrhotin, lnrhotinm, lnrhotinn;
	static float rhoti, rhoti1, rhoti1p, rhotim, rhot1m, rhotin, rhotinm, rhotnp;
	static float rhotdi, rhotdi1, rhotdi1p, rhotdim, rhotd1m, rhotdin, rhotdinm, rhotdnp;
	static float fabsi, fabsi1, fabsi1p, fabsim, fabsin, fabsinm, fabsinp;
	static float fsgni, fsgni1, fsgni1p, fsgnim, fsgnin, fsgninm, fsgninp;
	static float termi, termi1, termi1p, termim, termin, terminm, terminp;
	static float terpi, terpi1, terpi1p, terpim, terpin, terpinm, terpinp;
	// the first dimension of bc[][] is adjusted with -1 as usual with arrays from FORTRAN to c, but the second dimension already counts from 0 to NPT inclusive for its indices, so a +1 is used to "unadjust" ...
	// why not just don't adjust instead of adjusting and then "unadjusting"?  this way the I1 and IN #define substitutions can still be used, especially since i1m i1n i1p inp and inp all depend on i1 and in

	i1m = *I1 - 1;
	i1p = *I1 + 1;
	i1q = *I1 + 2;
	inp = *IN + 1;
	inm = *IN - 1;

	flxhi1 = fct_velo__.hadudth[*I1] * (fct_bc__.bc[1-1][*I1+1] * rhoo[*I1] + fct_bc__.bc[2-1][*I1+1] + fct_bc__.bc[1-1][i1m+1] * rhoo[i1m] + fct_bc__.bc[2-1][i1m+1]);
	diffi1 = fct_velo__.nulh[*I1] * (fct_bc__.bc[1-1][*I1+1] * rhoo[*I1] + fct_bc__.bc[2-1][*I1+1] - fct_bc__.bc[1-1][i1m+1] * rhoo[i1m] - fct_bc__.bc[2-1][i1m+1]);
	flxhi1p = fct_velo__.hadudth[i1p] * (fct_bc__.bc[1-1][i1p+1] * rhoo[i1p] + fct_bc__.bc[2-1][i1p+1] + fct_bc__.bc[1-1][*I1+1] * rhoo[*I1] + fct_bc__.bc[2-1][*I1+1]);
	diffi1p = fct_velo__.nulh[i1p] * (fct_bc__.bc[1-1][i1p+1] * rhoo[i1p] + fct_bc__.bc[2-1][i1p+1] - fct_bc__.bc[1-1][*I1+1] * rhoo[*I1] - fct_bc__.bc[2-1][*I1+1]);
	flxhi1q = fct_velo__.hadudth[i1q] * (fct_bc__.bc[1-1][i1q+1] * rhoo[i1q] + fct_bc__.bc[2-1][i1q+1] + fct_bc__.bc[1-1][i1p+1] * rhoo[i1p] + fct_bc__.bc[2-1][i1p+1]);
	diffi1q = fct_velo__.nulh[i1q] * (fct_bc__.bc[1-1][i1q+1] * rhoo[i1q] + fct_bc__.bc[2-1][i1q+1] - fct_bc__.bc[1-1][i1p+1] * rhoo[i1p] - fct_bc__.bc[2-1][i1p+1]);

	flxhi = flxhi1q;
	diffi = diffi1q;

	lorhoti1 = fct_grid__.lo[*I1] * rhoo[*I1] + fct_misc__.source[*I1] + (flxhi1 - flxhi1p);
	lnrhoti1 = lorhotil + (diffi1p - diffi1);
	rhoti1 = lorhoti1 * fct_grid__.rln[*I1];
	rhotdi1 = lnrhoti1 * fct_grid__.rln[*I1];
	lorhoti1p = fct_grid__.lo[i1p] * rhoo[*I1] + fct_misc__.source[i1p] + (flxhi1p - flxhi1q);
	lnrhoti1p = lorhoti1p;
	rhoti1p = lorhoti1p * fct_grid__.rln[i1p];
	rhotdi1p = lnrhoti1p * fct_grid__.rln[i1p];
	lnrhotil = lnrhoti1;
	lnrhotim = lnrhoti1p;
	rhotim = rhoti1p;
	rhotdim = rhotdi1p;

	rhot1m = rhoti1;
	rhotd1m = rhotdi1;
	flxhpi1 = fct_velo__.mulh[*I1] * (fct_bc__.bc[1-1][*I1+1] * rhoti1 + fct_bc__.bc[2-1][i1p+1] - fct_bc__.bc[1-1][i1m+1] * rhot1m - fct_bc__.bc[2-1][i1m+1]);
	diffpi1 =                        (fct_bc__.bc[1-1][*I1+1] * rhotdi1 + fct_bc__.bc[2-1][i1p+1] - fct_bc__.bc[1-1][i1m+1] * rhotd1m - fct_bc__.bc[2-1][i1m+1]);
	fabsi1 = float_abs(flxhpi1);
	fsgni1 = sign(fct_misc__.diff1_, diffpi1);
	flxhpi1p = fct_velo__.mulh[i1p] * (fct_bc__.bc[1-1][i1p+1] * rhoti1p + fct_bc__.bc[2-1][i1p+1] - fct_bc__.bc[1-1][*I1+1] * rhoti1 - fct_bc__.bc[2-1][*I1+1]);
	diffpi1p =                        (fct_bc__.bc[1-1][i1p+1] * rhotdi1p + fct_bc__.bc[2-1][i1p+1] - fct_bc__.bc[1-1][*I1+1] * rhotdi1 - fct_bc__.bc[2-1][*I1+1]);
	diffpim = diffpi1p;

	fabsi1p = float_abs(flxhpi1p);
	fsgni1p = sign(fct_misc__.diff1_, diffpi1p);
	termi1p = fsgni1p * fct_grid__.ln[*I1] * diffpi1 + fct_bc__.bc[3-1][i1p+1] * BIGNUM;
	terpi1 = fsgni1 * fct_grid__.ln[*I1] * diffpi1p + fct_bc__.bc[4-1][*I1+1] * BIGNUM;
	fabsim = fabsi1p;
	fsgnim = fsgni1p;
	termim = termi1p;
	termi1 = fct_bc__.bc[3-1][*I1+1] * BIGNUM;

	flxhqi1 = fsgni1 * amax1(0.0, amin1(termi1, fabsi1, terpi1));
	flxhqil = flxhqi1;

	int i;	// compiler complains if counter declared in for statement in c (but ok in c++?!)
	for (i = i1q; i <= inm; i++)
	{
		flxhip = fct_velo__.hadudth[i+1] * (fct_bc__.bc[1-1][i+1+1] * rhoo[i+1] + fct_bc__.bc[1-1][i+1] * rhoo[i] + fct_bc__.bc[2-1][i+1]);
		diffip = fct_velo__.nulh[i+1] * (fct_bc__.bc[1-1][i+1+1] * rhoo[i+1] - fct_bc__.bc[1-1][i+1] * rhoo[i] - fct_bc__.bc[2-1][i+1]);
		lorhoti = fct_grid__.lo[i] * rhoo[i] + fct_misc__.source[i] + (flxhi - flxhip);
		lnrhoti = lorhoti + (diffip - diffi);
		rhoti = lorhoti * fct_grid__.rln[i];
		rhotdi = lnrhoti * fct_grid__.rln[i];
		flxhpi = fct_velo__.mulh[i] * (fct_bc__.bc[1-1][i+1] * rhoti + fct_bc__.bc[2-1][i+1] - fct_bc__.bc[1-1][i-1+1] * rhotim - fct_bc__.bc[2-1][i-1+1]);
		diffpi =                      (fct_bc__.bc[1-1][i+1] * rhotdi + fct_bc__.bc[2-1][i+1] - fct_bc__.bc[1-1][i-1+1] * rhotdim - fct_bc__.bc[2-1][i-1+1]);
		fabsi = float_abs(flxhpi);
		fsgni = sign(fct_misc__.diff1_, diffpi);
		termi = fsgni * fct_grid__.ln[i-1] * diffpim + fct_bc__.bc[3-1][i+1] * BIGNUM;
		terpim = fsgnim * fct_grid__.ln[i-1] * diffpi + fct_bc__.bc[4-1][i-1+1] * BIGNUM;
		flxhqim = fsgnim * amax1(0.0, amin1(termim, fabsim, terpim));
		rhon[i-2] = fct_grid__.rln[i-2] * (lnrhotil + (flxhqil - flxhqim));
		fct_misc__.source[i-2] = 0.0;
		flxhi = flxhip;
		diffi = diffip;
		diffpim = diffpi;
		flxhqil = flxhqim;
		lnrhotil = lnrhotim;
		lnrhotim = lnrhoti;
		rhotim = rhoti;
		rhotdim = rhotdi;
		fabsim = fabsi;
		fsgnim = fsgni;
		termim = termi;
	}

	flxhin = flxhi;
	diffin = diffi;

	flxhinp = fct_velo__.nulh[inp] * (fct_bc__.bc[1-1][inp+1] * rhoo[inp] + fct_bc__.bc[2-1][inp+1] + fct_bc__.bc[1-1][*IN+1] * rhoo[*IN] + fct_bc__.bc[2-1][*IN+1]);
	diffinp = fct_velo__.hadudth[inp] * (fct_bc__.bc[1-1][inp+1] * rhoo[inp] + fct_bc__.bc[2-1][inp+1] - fct_bc__.bc[1-1][*IN+1] * rhoo[*IN] - fct_bc__.bc[2-1][*IN+1]);

	lorhotin = fct_grid__.lo[*IN] * rhoo[*IN] + fct_misc__.source[*IN] + (flxhin - flxhinp);
	lnrhotin = lorhotin + (diffinp - diffin);
	rhotin = lorhotin * fct_grid__.rln[*IN];
	rhotdin = lnrhotin * fct_grid__.rln[*IN];

	rhotinm = rhotim;
	rhotdinm = rhotdim;
	diffpinm = diffpim;
	rhotnp = rhotin;
	rhotdnp = rhotdin;
	flxhpin = fct_velo__.mulh[*IN] * (fct_bc__.bc[1-1][*IN+1] * rhotin + fct_bc__.bc[2-1][*IN+1] - fct_bc__.bc[1-1][inm+1] * rhotinm - fct_bc__.bc[2-1][inm+1]);
	diffpin =                        (fct_bc__.bc[1-1][*IN+1] * rhotdin + fct_bc__.bc[2-1][*IN+1] - fct_bc__.bc[1-1][inm+1] * rhotdinm - fct_bc__.bc[2-1][inm+1]);
	flxhpinp = fct_velo__.mulh[inp] * (fct_bc__.bc[1-1][inp+1] * rhotnp + fct_bc__.bc[2-1][inp+1] - fct_bc__.bc[1-1][*IN+1] * rhotin - fct_bc__.bc[2-1][*IN+1]);
	diffpinp =                        (fct_bc__.bc[1-1][inp+1] * rhotdnp + fct_bc__.bc[2-1][inp+1] - fct_bc__.bc[1-1][*IN+1] * rhotdin - fct_bc__.bc[2-1][*IN+1]);

	fabsinm = fabsim;
	fsgninm = fsgnim;
	terminm = termim;

	fabsin = float_abs(flxhpin);
	fsgnin = sign(fct_misc__.diff1_, diffpin);
	termin = fsgnin * fct_grid__.ln[inm] * diffpinm + fct_bc__.bc[3-1][*IN+1] * BIGNUM;
	terpinm = fsgninm * fct_grid__.ln[inm] * diffpin + fct_bc__.bc[4-1][inm+1] * BIGNUM;
	fabsinp = float_abs(flxhpinp);
	fsgninp = sign(fct_misc__.diff1_, diffpinp);
	terminp = fsgninp * fct_grid__.ln[*IN] * diffpin + fct_bc__.bc[3-1][inp+1] * BIGNUM;
	terpin = fsgnin * fct_grid__.ln[*IN] * diffpinp + fct_bc__.bc[4-1][*IN+1] * BIGNUM;

	lnrhotinn = lnrhotil;
	lnrhotinm = lnrhotim;
	flxhqinn = flxhqil;
	flxhqinm = fsgninm * amax1(0.0, amin1(terminm, fabsinm, terpinm));
	rhon[*IN-2] = fct_grid__.rln[*IN-2] * (lnrhotinn + (flxhqinn - flxhqinm));
	fct_misc__.source[*IN-2] = 0.0;
	flxhqin = fsgnin * amax1(0.0, amin1(termin, fabsin, terpin));
	rhon[inm] = fct_grid__.rln[inm] * (lnrhotinm + (flxhqinm - flxhqin));
	fct_misc__.source[inm] = 0.0;
	flxhqinp = fsgninp * amax1(0.0, amin1(terminp, fabsinp, terpinp));
	rhon[*IN] = fct_grid__.rln[*IN] * (lnrhotin + (flxhqin - flxhqinp));
	fct_misc__.source[*IN] = 0.0;

	return;
}
