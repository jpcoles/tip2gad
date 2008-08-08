#include <assert.h>
#include <iostream>
#include <stdio.h>
#include "allvars.h"

using namespace std;

enum { GAS, HALO, DISK, BULGE, STARS, BNDRY };

int main(int argc, char **argv)
{
    FILE *in;
    io_header iohdr;
    int block0, block1;

    if (argc < 2)
    {
        cerr << "Usage: gadget_header <gadget-file>" << endl;
        exit(2);
    }

    in = fopen(argv[1], "r");

    if (in == NULL)
    {
        cerr << "ERROR: Can't open file " << argv[1] << "." << endl;
        exit(1);
    }

    if (fread(&block0, sizeof(block0), 1, in) < 1)
    {
        cerr << "ERROR: Unexpected EOF." << endl;
        exit(1);
    }

    if (fread(&iohdr, sizeof(iohdr), 1, in) < 1)
    {
        cerr << "ERROR: Unexpected EOF." << endl;
        exit(1);
    }

    if (fread(&block1, sizeof(block1), 1, in) < 1)
    {
        cerr << "ERROR: Unexpected EOF." << endl;
        exit(1);
    }

    fclose(in);

    assert(block0 == block1);

    cout << "----------------------------------------------------------------" << endl;
    cout << "Gadget Header" << endl;
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

    return 0;
}

