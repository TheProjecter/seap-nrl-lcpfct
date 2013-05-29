#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fcntl.h>
#include <argp.h>
#include <einlib.h>  
#include <unistd.h>

#define NPT = 2002
#define MEM_SIZE     (8*0x100000) //in bytes
#define MEM_DISTANCE (8*0x100000)
#define MEM_OFFSET   (0x40 * 0x100000)
#define TYPE_VAL  0x0UL
#define TYPE_ADDR 0x1UL

typedef unsigned float u_64;
u_64* rams[4];  
void *ftr_mem;

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

int main(void)
{
	int i, fp_id, num_bytes, pagesize, buf_size, memali;
	const char *loadfile = "top.bin.ufp";
	err_e err, e;
	void * bp;
	u_64 val;
	long long data[2][1024];
	char line_in[256];
	FILE* file;
	char* filename;  
	int k;

	printf("Start program\n");

	// Set FPGA ready
	fp_id = fpga_open("/dev/ufp0", O_RDWR|O_SYNC, &e);
	num_bytes = fpga_load(fp_id, loadfile, &err);
	fpga_reset(fp_id, &e);
	fpga_start(fp_id, &e);

	// Set RAMS ready
	rams[0] = fpga_memmap(fp_id, MEM_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, MEM_OFFSET + 0*MEM_DISTANCE, &e);
	rams[1] = fpga_memmap(fp_id, MEM_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, MEM_OFFSET + 1*MEM_DISTANCE, &e);
	rams[2] = fpga_memmap(fp_id, MEM_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, MEM_OFFSET + 2*MEM_DISTANCE, &e);
	rams[3] = fpga_memmap(fp_id, MEM_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, MEM_OFFSET + 3*MEM_DISTANCE, &e);

	// Convert structs to arrays
	u_64 scrh[], grid[], velo[], misc[];

	for (int i = 0; i < NPT; i++)
	{
		scrh[i + 0 * NPT] = fct_scrh__.scrh[i];
		scrh[i + 1 * NPT] = fct_scrh__.scr1[i];
		scrh[i + 2 * NPT] = fct_scrh__.diff[i];
		scrh[i + 3 * NPT] = fct_scrh__.flxh[i];
		scrh[i + 4 * NPT] = fct_scrh__.fabs[i];
		scrh[i + 5 * NPT] = fct_scrh__.fsgn[i];
		scrh[i + 6 * NPT] = fct_scrh__.term[i];
		scrh[i + 7 * NPT] = fct_scrh__.terp[i];
		scrh[i + 8 * NPT] = fct_scrh__.lnorhot[i];
		scrh[i + 9 * NPT] = fct_scrh__.lorhot[i];
		scrh[i + 10 * NPT] = fct_scrh__.rhot[i];
		scrh[i + 11 * NPT] = fct_scrh__.rhotd[i];

		grid[i + 0 * NPT] = fct_grid__.lo[i];
		grid[i + 1 * NPT] = fct_grid__.ln[i];
		grid[i + 2 * NPT] = fct_grid__.ah[i];
		grid[i + 3 * NPT] = fct_grid__.rln[i];
		grid[i + 4 * NPT] = fct_grid__.lh[i];
		grid[i + 5 * NPT] = fct_grid__.rlh[i];
		grid[i + 6 * NPT] = fct_grid__.roh[i];
		grid[i + 7 * NPT] = fct_grid__.rnh[i];
		grid[i + 8 * NPT] = fct_grid__.adugth[i];

		velo[i + 0 * NPT] = fct_velo__.hadudth[i];
		velo[i + 1 * NPT] = fct_velo_.nulh[i];
		velo[i + 2 * NPT] = fct_velo__.mulh[i];
		velo[i + 3 * NPT] = fct_velo__.epsh[i];
		velo[i + 4 * NPT] = fct_velo__.vdtodr[i];

		misc[i + 0 * NPT] = fct_misc__.source[i];
		misc[i + 1 * NPT] = fct_misc__.diff1_[i];
	}
	
	// Initalize vectors and input data on rams
	u_64 word0 = (u_64) scrh[0];
	u_64 word1 = (u_64) grid[0];
	u_64 word2 = (u_64) velo[0];
	u_64 word3 = (u_64) misc[0];
	u_64 *word_0 = rams[0];
	u_64 *word_1 = rams[1];
	u_64 *word_2 = rams[2];
	u_64 *word_3 = rams[3];
	
	for( i = 0; i < 12 * NPT; i++)
    {
		*word_0 ++= (u_64) scrh[i];
		*word_1 ++= (u_64) grid[i];
		*word_2 ++= (u_64) velo[i];
		*word_3 ++= (u_64) misc[i];
	}

	// Set memory for FPGA use
	pagesize = getpagesize();
	buf_size = ((128*1024*1024 + (pagesize - 1)) / pagesize) * pagesize;
	memali   = posix_memalign(&bp, pagesize, buf_size);
	fpga_register_ftrmem(fp_id, bp, buf_size, &e);
	ftr_mem = bp;

	// Setup Mitrion run
	void * addr = ftr_mem;
	fpga_wrt_appif_val(fp_id, 0x00000000000000001UL, (0x01*sizeof(u_64)) + 0, TYPE_VAL, &e);
	fpga_wrt_appif_val(fp_id, 88, ((0x40+1)*sizeof(u_64)), TYPE_VAL, &e);
	fpga_wrt_appif_val(fp_id, (u_64)addr, (0x40*sizeof(u_64)), TYPE_ADDR, &e);
	fpga_wrt_appif_val(fp_id, 0x000000000000000000UL, (0x02 * sizeof(u_64)), TYPE_VAL, &e);
	fpga_wrt_appif_val(fp_id, 0x000000000000000000UL, (0x01 * sizeof(u_64)), TYPE_VAL, &e);

	// Loop and wait until FPGA calculation is done
	do
	{
		fpga_rd_appif_val(fp_id, (void*)&val, (0x02 * sizeof(u_64)), &e);
	} while((val&1) == 0);

	// put data from rams back into structs
	*word_0 = rams[0];
	*word_1 = rams[1];
	*word_2 = rams[2];
	*word_3 = rams[3];
	word0 = *word_0;
	word1 = *word_1;
	word2 = *word_2;
	word3 = *word_3;
	for(i = 0;i < 12 * NPT; i++)
	{
		word0 = *word_0; word_0++;
		word1 = *word_1; word_1++;
		word2 = *word_2; word_2++;
		word3 = *word_3; word_3++;
	}

	// convert arrays back into structs
	for (int i = 0; i < NPT; i++)
	{
		fct_scrh__.scrh[i] = scrh[i + 0 * NPT];
		fct_scrh__.scr1[i] = scrh[i + 1 * NPT];
		fct_scrh__.diff[i] = scrh[i + 2 * NPT];
		fct_scrh__.flxh[i] = scrh[i + 3 * NPT];
		fct_scrh__.fabs[i] = scrh[i + 4 * NPT];
		fct_scrh__.fsgn[i] = scrh[i + 5 * NPT];
		fct_scrh__.term[i] = scrh[i + 6 * NPT];
		fct_scrh__.terp[i] = scrh[i + 7 * NPT];
		fct_scrh__.lnorhot[i] = scrh[i + 8 * NPT];
		fct_scrh__.lorhot[i] = scrh[i + 9 * NPT];
		fct_scrh__.rhot[i] = scrh[i + 10 * NPT];
		fct_scrh__.rhotd[i] = scrh[i + 11 * NPT];

		fct_grid__.lo[i] = grid[i + 0 * NPT];
		fct_grid__.ln[i] = grid[i + 1 * NPT];
		fct_grid__.ah[i] = grid[i + 2 * NPT];
		fct_grid__.rln[i] = grid[i + 3 * NPT];
		fct_grid__.lh[i] = grid[i + 4 * NPT];
		fct_grid__.rlh[i] = grid[i + 5 * NPT];
		fct_grid__.roh[i] = grid[i + 6 * NPT];
		fct_grid__.rnh[i] = grid[i + 7 * NPT];
		fct_grid__.adugth[i] = grid[i + 8 * NPT];

		fct_velo__.hadudth[i] = velo[i + 0 * NPT];
		fct_velo_.nulh[i] = velo[i + 1 * NPT];
		fct_velo__.mulh[i] = velo[i + 2 * NPT];
		fct_velo__.epsh[i] = velo[i + 3 * NPT];
		fct_velo__.vdtodr[i] = velo[i + 4 * NPT];

		fct_misc__.source[i] = misc[i + 0 * NPT];
		fct_misc__.diff1_[i] = misc(i + 1 * NPT);
	}
	
	printf("End of program ");
	return 0;
}
