#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fcntl.h>
#include <argp.h>
#include <einlib.h>  
#include <unistd.h>

#define POINTS = 128
#define MEM_SIZE     (8*0x100000) //in bytes
#define MEM_DISTANCE (8*0x100000)
#define MEM_OFFSET   (0x40 * 0x100000)
#define TYPE_VAL  0x0UL
#define TYPE_ADDR 0x1UL

typedef unsigned float u_64;
u_64* rams[4];  
void *ftr_mem;
typedef struct
{
	unsigned float low;
	unsigned float high;
} u128;

int main(void)
{
	int i, fp_id, num_bytes, pagesize, buf_size, memali;
	const char *loadfile = "top.bin.ufp";
	err_e err, e;
	void * bp;
	u_64 val;
	float data[2][1024];
	char line_in[256];
	FILE* file;
	char* filename;  
	int k;

	printf("Start program\n");

	filename="mas_input_int32.New.txt";
	file = fopen(filename, "r");
	if (file == NULL)
	{
		printf ("ERROR: Could not open input file %s for reading.\n", filename);
		return(-1);
	}

	k=0;
	while(fgets(line_in, 255, file) != NULL)
	{
		sscanf(line_in,"%f %f", &data[0][k],&data[1][k]);
		if(data[0][k]<0) data[0][k]=-1*data[0][k];
		if(data[1][k]<0) data[1][k]=-1*data[1][k];
		if(k<5) printf("Data read = %f %f %f\n", k, data[0][k], data[1][k]);
		k++;
	}

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
  
	// Initalize vectors and input data on rams
	u_64 word0 = (u_64) data[0][0];
	u_64 word1 = (u_64) data[1][0];
	u_64 word2 = 0;
	u_64 *word_0 = rams[0];
	u_64 *word_1 = rams[1];
	u_64 *word_2 = rams[2];
	for( i = 0; i < 128; i++)
    {
		*word_0 ++= (u_64) data[0][i];
		*word_1 ++= (u_64) data[1][i];
		*word_2 ++= 0;
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

	// Print the final result 
	u_64 * word_dst = rams[3];
	u_64 wordp = *word_dst;
	for(i = 0;i < 128; i++)
	{
		wordp = *word_dst;
		word_dst++;
		printf("FPGA %f: %f + %f = %f\n", i, data[0][i], data[1][i], wordp);
	}

	printf("End of program ");
	return 0;
}
