#include <assert.h>
#include <getopt.h>
#include <iostream>
#include <stdio.h>
#include "allvars.h"

using namespace std;

enum { GAS, HALO, DISK, BULGE, STARS, BNDRY };


typedef struct 
{
    float  Pos[3];
    float  Vel[3];
    float  Mass;
    int    Type;

    float  Rho, U, Temp, Ne;
} particle_data_t;

int load_snapshot(char *fname, int files, bool only_show_header);

void help()
{
    cerr << "Usage: gadget_header <gadget-file>" << endl;
    exit(2);
}

int main(int argc, char **argv)
{
    int only_show_header = false;

    if (argc < 2) help();

    while (1)
    {
        int option_index = 0;
        int c = getopt_long(argc, argv, "h", NULL, &option_index);

        if (c == -1) break;

        switch (c)
        {
            case 'h': only_show_header = true; break;
        }
    }

    if (optind >= argc) help();

    load_snapshot(argv[optind], 1, only_show_header);

    return 0;
}

int load_snapshot(char *fname, int files, bool only_show_header)
{
    FILE *fd;
    char   buf[256];
    uint32_t dummy;
    char   id[5];
    int    i,k,ntot_withmasses;
    int    n,pc,pc_new,pc_sph;

    io_header iohdr;
    unsigned int block0, block1;

    int Ngas = 0;
    int NumPart = 0;

    particle_data_t P;

    int SnapFormat = 0;
    
#define SKIP0 \
    do { \
        fread(&block0, sizeof(dummy), 1, fd); \
        if (SnapFormat == 0) { \
            fread(&id, 4, 1, fd); \
            int _i; \
            for (_i=0; _i < 4; _i++) \
                if (id[_i] != "HEAD"[_i]) break; \
            SnapFormat = 1; \
            if (_i == 4) SnapFormat = 2; \
            printf("Autodetected GADGET version %i\n", SnapFormat); \
            fseek(fd, -4L, SEEK_CUR); \
        } \
        if (SnapFormat == 2) { \
            fread(&id, 4, 1, fd); \
            id[4] = 0; \
            printf("Found block '%s'\n", id); \
        } \
    } while (0)
#define SKIP1 do { fread(&block1, sizeof(dummy), 1, fd); } while (0)
#define TEST_BLOCK_LEN if (block0 != block1) { cerr << "Block size mismatch: " << block0 << " " << block1 << endl; }
//#define SKIP do { fseek(fd, sizeof(dummy), SEEK_CUR); } while (0)

    for(i=0, pc=1; i<files; i++, pc=pc_new)
    {

        if(files>1)
	        snprintf(buf, 256, "%s.%d",fname,i);
        else
	        snprintf(buf, 256, "%s",fname);

        if(!(fd=fopen(buf,"r")))
	    {
	        printf("ERROR: Can`t open file %s\n",buf);
	        exit(0);
	    }

        SKIP0;
        fread(&iohdr, sizeof(iohdr), 1, fd);
        SKIP1;
        TEST_BLOCK_LEN;

        cout << "----------------------------------------------------------------" << endl;
        cout << "Gadget Header (" << buf << ")" << endl;
        cout << "----------------------------------------------------------------" << endl;
        cout << "Particle Count: "
             << "Gas("   << iohdr.npartTotal[GAS  ] << ") "
             << "Halo("  << iohdr.npartTotal[HALO ] << ") "
             << "Disk("  << iohdr.npartTotal[DISK ] << ") " << endl
             << "                " 
             << "Bulge(" << iohdr.npartTotal[BULGE] << ") "
             << "Stars(" << iohdr.npartTotal[STARS] << ") "
             << "Bndry(" << iohdr.npartTotal[BNDRY] << ") " << endl;
        cout << "Time:        " << iohdr.time        << endl;
        cout << "Redshift:    " << iohdr.redshift    << endl;
        cout << "Omega0:      " << iohdr.Omega0      << endl;
        cout << "OmegaLambda: " << iohdr.OmegaLambda << endl;
        cout << "HubbleParam: " << iohdr.HubbleParam << endl;
        cout << "BoxSize:     " << iohdr.BoxSize     << endl;
        cout << "----------------------------------------------------------------" << endl << endl;

        if (only_show_header) continue;

        cout << "Verifying..." << endl;

        if(files==1)
	    {
	        for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	            NumPart+= iohdr.npart[k];
	        Ngas= iohdr.npart[0];
	    }
        else
	    {
	        for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	            NumPart+= iohdr.npartTotal[k];
	        Ngas= iohdr.npartTotal[0];
	    }

        for(k=0, ntot_withmasses=0; k<5; k++)
	    {
	        if(iohdr.mass[k]==0)
	            ntot_withmasses+= iohdr.npart[k];
	    }
	    //printf("N: %d\n",ntot_withmasses);

        SKIP0;
        for(k=0,pc_new=pc;k<6;k++)
	    {
	        for(n=0;n<iohdr.npart[k];n++)
	        {
	            fread(&P.Pos, sizeof(float), 3, fd);
                int sum = (P.Pos[0] == 0) + (P.Pos[1] == 0) + (P.Pos[2] == 0);
                if (sum > 1)
                {
                    cerr << "Suspicious position: [" << pc_new << " | " << k << "] " 
                         << P.Pos[0] << " " << P.Pos[1] << " " << P.Pos[2] << endl;
                }
	            pc_new++;
	        }
	    }
        SKIP1;
        TEST_BLOCK_LEN;
    
        SKIP0;
        for(k=0,pc_new=pc;k<6;k++)
	    {
	        for(n=0;n<iohdr.npart[k];n++)
	        {
	            fread(&P.Vel, sizeof(float), 3, fd);
                int sum = (P.Vel[0] == 0) + (P.Vel[1] == 0) + (P.Vel[2] == 0);
                if (sum > 1)
                {
                    cerr << "Suspicious velocity: [" << pc_new << " | " << k << "] "
                         << P.Vel[0] << " " << P.Vel[1] << " " << P.Vel[2] << endl;
                }
	            pc_new++;
	        }
	    }
        SKIP1;
        TEST_BLOCK_LEN;
    

        SKIP0;
        for(k=0,pc_new=pc;k<6;k++)
	    {
	        for(n=0;n<iohdr.npart[k];n++)
	        {
                unsigned int Id;
	            fread(&Id, sizeof(int), 1, fd);
	            pc_new++;
	        }
	    }
        SKIP1;
        TEST_BLOCK_LEN;


        if(ntot_withmasses>0)
	        SKIP0;

        for(k=0, pc_new=pc; k<6; k++)
	    {
            for(n=0;n<iohdr.npart[k];n++)
            {
                if(iohdr.mass[k]==0)
                {
                    fread(&P.Mass, sizeof(float), 1, fd);
                    if (P.Mass <= 0)
                    {
                        cerr << "Suspicious Mass: [" << pc_new << " | " << k << "] " << P.Mass << endl;
                    }
                }
                pc_new++;
            }
	    }
        if(ntot_withmasses>0)
        {
	        SKIP1;
            TEST_BLOCK_LEN;
        }

        if(iohdr.npart[0]>0)
	    {
	        SKIP0;
	        for(n=0, pc_sph=pc; n<iohdr.npart[0];n++)
	        {
	            fread(&P.U, sizeof(float), 1, fd);
                pc_sph++;
	        }
	        SKIP1;
            TEST_BLOCK_LEN;

	        SKIP0;
	        for(n=0, pc_sph=pc; n<iohdr.npart[0];n++)
	        {
	            fread(&P.Rho, sizeof(float), 1, fd);
                pc_sph++;
	        }
	        SKIP1;
            TEST_BLOCK_LEN;

	        if(iohdr.flag_cooling)
	        {
	            SKIP0;
	            for(n=0, pc_sph=pc; n<iohdr.npart[0];n++)
		        {
		            fread(&P.Ne, sizeof(float), 1, fd);
                    pc_sph++;
		        }
	            SKIP1;
                TEST_BLOCK_LEN;
	        }
	    }

        fclose(fd);
    }


    return 0;
}
