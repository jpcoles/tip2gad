#include <assert.h>
#include <getopt.h>
#include <iostream>
#include <stdio.h>
#include "ftipsy.hpp"

using namespace std;

int load_snapshot(char *fname, bool only_show_header);

void help()
{
    cerr << "Usage: verify_tipsy <tipsy-file>" << endl;
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

    load_snapshot(argv[optind], only_show_header);

    return 0;
}

int load_snapshot(char *fname, bool only_show_header)
{
    ifTipsy in;
    TipsyHeader h;

    in.open(fname, "standard");
    if(!in.is_open())
    {
        printf("ERROR: Can`t open file %s\n", fname);
        exit(0);
    }

    in >> h;

    cout << "----------------------------------------------------------------" << endl;
    cout << "Tipsy Header (" << fname << ")" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << "Particle Count: "
         << "Total(" << h.h_nBodies << ") "
         << "Gas("   << h.h_nGas << ") "
         << "Dark("  << h.h_nDark << ") "
         << "Star("  << h.h_nStar << ") " << endl;
    cout << "Time:        " << h.h_time        << endl;
    cout << "Redshift:    " << (1.0F / h.h_time - 1.0F)    << endl;
    cout << "----------------------------------------------------------------" << endl << endl;

    if (only_show_header) return 0;

#if 0
    for (int i=0; i < h.h_nGas; i++)
    {

    }
#endif

    return 0;
}
