#include <stdio.h>		/* for printing debugging messages */
#include <stdbool.h>	/* support on UNIX for bool vars */
#include <math.h>		/* support for fabs() */

/* all the extern struct stuff from the original c program put here instead ... since they need to be processed */
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

/* number of places to move decmial point by */
#define move 1.0e12

void main (float *rhoo, float *rhon, float* i1, float* in,
			float* srho1, float* vrho1, float* srhon,
			float* vrhon, bool* pbc)
{
	/* declare int structs */
	//	/FCT_SCRH/ Holds scratch arrays for use by LCPFCT and CNVFCT
	struct
	{
		int scrh[NPT];
		int scr1[NPT];
		int diff[NPT];
		int flxh[NPT];
		int fabs[NPT];
		int fsgn[NPT];
		int term[NPT];
		int terp[NPT];
		int lnrhot[NPT];
		int lorhot[NPT];
		int rhot[NPT];
		int rhotd[NPT];
	} scrh;
	//	$OMP THREADPRIVATE(/FCT_SCRH/)

	//	/FCT_GRID/ Holds geometry, grid, area and volume information
	struct
	{
		int lo[NPT];
		int ln[NPT];
		int ah[NPT];
		int rln[NPT];
		int lh[NPT];
		int rlh[NPT];
		int roh[NPT];
		int rnh[NPT];
		int adugth[NPT];
	} grid;
	//	$OMP THREADPRIVATE(/FCT_GRID/)

	//	/FCT_VELO/ Holds velocity-dependent flux coefficients
	struct
	{
		int hadudth[NPT];
		int nulh[NPT];
		int mulh[NPT];
		int epsh[NPT];
		int vdtodr[NPT];
	} velo;
	//	$OMP THREADPRIVATE( /FCT_VELO/)

	//	/FCT_MISC/ Holds the source array and diffusion coefficient
	struct
	{
		int source[NPT];
		int diff1_;
	} misc;
	//	$OMP THREADPRIVATE( /FCT_MISC/)
	
	/* arguments passed to lcpfct */
	int irhoo[rhoo.length];
	int irhon[rhon.length];
	int isrho1, ivrho1, isrhon, ivrhon;
	
	/* convert floats to ints */
	for (int i = 0; i < NPT; i++)
	{
		scrh.scrh[i] = fct_scrh__.scrh[i] * move;
		scrh.scr1[i] = fct_scrh__.scr1[i] * move;
		scrh.diff[i] = fct_scrh__.diff[i] * move;
		scrh.flxh[i] = fct_scrh__.flxh[i] * move;
		scrh.fabs[i] = fct_scrh__.flxh[i] * move;
		scrh.fsgn[i] = fct_scrh__.fsgn[i] * move;
		scrh.term[i] = fct_scrh__.term[i] * move;
		scrh.terp[i] = fct_scrh__.terp[i] * move;
		scrh.lnrhot[i] = fct_scrh__.lnrhot[i] * move;
		scrh.lorhot[i] = fct_scrh__.lorhot[i] * deimcalMove;
		scrh.rhot[i] = fct_scrh__.rhot[i] * move;
		scrh.rhotd[i] = fct_scrh__.rhotd[i] * move;
		
		grid.lo[i] = fct_grid__.io[i] * move;
		grid.ln[i] = fct_grid__.ln[i] * move;
		grid.ah[i] = fct_grid__.ah[i] * move;
		grid.rln[i] = fct_grid__.rln[i] * move;
		grid.lh[i] = fct_grid__.lh[i] * move;
		grid.rlh[i] = fct_grid__.rlh[i] * move;
		grid.roh[i] = fct_grid__.roh[i] * move;
		grid.rnh[i] = fct_grid__.rnh[i] * move;
		grid.adugth[i] = fct_grid__.adugth[i] * move;
		
		velo.hadudth[i] = fct_velo__.hadudth[i] * move;
		velo.nulh[i] = fct_velo__.nulh[i] * move;
		velo.mulh[i] = fct_velo__.mulh[i] * move;
		velo.epsh[i] = fct_velo__.epsh[i] * move;
		velo.vdtodr[i] = fct_velo__.vdtodr[i] * move;
		
		misc.source[i] = fct_misc__.source[i] * move;
		misc.diff1_ = fct_misc__.diff1_ * move;
	}
	
	for (int i = 0; i < rhoo.length; i++)
	{
		irhoo[i] = rhoo[i] * move;
	}
	for (int i = 0; i < rhon.length; i++)
	{
		irhon[i] = rhon[i] * move;
	}
	isrho1 = srho1 * move;
	ivrho1 = vrho1 * move;
	isrhon = srhon * move;
	ivrhon = vrhon * move;
	
	/* plan of attack:
	 * 	method 1 - convert parameters just like structs ... create new vars, initialize, then set equal
	 *		requires thinking :(
	 * 		needs extra variables
	 *		would be something like new_num = *old_num, then pass int& new_num
	 *		arrays would be like new_array = old_array, then pass &new_array
	 * 	method 2 - convert parameters "in real time" within the function call
	 * 		bad because arrays are tricky things and simply typing is not going to work
	 * 		and besides, you cannot cast from a pointer of one type to a pointer of another type
	 *			except that you can cast a pointer to void and back ... ^_^
	 */
	 
	/* call lcpfct fpga function thing with "corrected" values */
	int fpga;
	err_e e;
	fpga = fpga_open("/dev/ufp0", 0_RDWRIO_SYNC, &e);
	fpga(int &rhoo, int &rhon, int& i1, int& in,
			int& isrho1, int& ivrho1, int& isrhon,
			int& ivrhon, int& pbc);
	
	/* convert ints to floats */
	for (int i = 0; i < irhoo.length; i++)
	{
		rhoo[i] = irhoo[i] / move;
	}
	for (int i = 0; i < irhon.length; i++)
	{
		rhon[i] = irhon[i] / move;
	}
	srho1 = isrho1 / move;
	vrho1 = ivrho1 / move;
	srhon = isrhon / move;
	vrhon = ivrhon / move;
	
	for (int i = 0; i < NPT; i++)
	{
		fct_scrh__.scrh[i] = scrh.scrh[i] / move;
		fct_scrh__.scr1[i] = scrh.scr1[i] / move;
		fct_scrh__.diff[i] = scrh.diff[i] / move;
		fct_scrh__.flxh[i] = scrh.flxh[i] / move;
		fct_scrh__.fabs[i] = scrh.flxh[i] / move;
		fct_scrh__.fsgn[i] = scrh.fsgn[i] / move;
		fct_scrh__.term[i] = scrh.term[i] / move;
		fct_scrh__.terp[i] = scrh.terp[i] / move;
		fct_scrh__.lnrhot[i] = scrh.lnrhot[i] / move;
		fct_scrh__.lorhot[i] = scrh.lorhot[i] / deimcalMove;
		fct_scrh__.rhot[i] = scrh.rhot[i] / move;
		fct_scrh__.rhotd[i] = scrh.rhotd[i] / move;
		
		fct_grid__.lo[i] = grid.io[i] / move;
		fct_grid__.ln[i] = grid.ln[i] / move;
		fct_grid__.ah[i] = grid.ah[i] / move;
		fct_grid__.rln[i] = grid.rln[i] / move;
		fct_grid__.lh[i] = grid.lh[i] / move;
		fct_grid__.rlh[i] = grid.rlh[i] / move;
		fct_grid__.roh[i] = grid.roh[i] / move;
		fct_grid__.rnh[i] = grid.rnh[i] / move;
		fct_grid__.adugth[i] = grid.adugth[i] / move;
		
		fct_velo__.hadudth[i] = velo.hadudth[i] / move;
		fct_velo__.nulh[i] = velo.nulh[i] / move;
		fct_velo__.mulh[i] = velo.mulh[i] / move;
		fct_velo__.epsh[i] = velo.epsh[i] / move;
		fct_velo__.vdtodr[i] = velo.vdtodr[i] / move;
		
		fct_misc__.source[i] = misc.source[i] / move;
		fct_misc__.diff1_ = misc.diff1_ / move;
	}
}