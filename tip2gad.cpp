#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <getopt.h>
#include <assert.h>

#include "ftipsy.hpp"
#include "allvars.h"

using namespace std;

#define VERSION "0.8"

#define DEFAULT_MEM 128

#define DIR_SEP '/'

#define IO_HEAD (IO_NBLOCKS)

enum { GAS, HALO, DISK, BULGE, STARS, BNDRY };

int verbosity = 1;
#define VL(__l) if (verbosity >= __l)


#define READ_GAS   1
#define READ_DARK  2
#define READ_STARS 4
#define WRITE_DATA 8
#define READ_DONE  16

class Config
{
public:
    float Omega0;
    float OmegaLambda;
    float BoxSize;
    float HubbleParam;
    float mass[6];
    float time;
    unsigned int SnapFormat;
    string TipsyFile;
    string OutputDir;
    float RhoCritical;

    Config() :
        Omega0(0), OmegaLambda(0),
        BoxSize(0), HubbleParam(0),
        time(0), SnapFormat(1),
        TipsyFile(), OutputDir(),
        RhoCritical(1.0F)
    {
        memset(mass, 0, sizeof(mass));
    }
};


class FileWriter
{
protected:
    ofstream of;

public:

    FileWriter(string filename) : of(filename.c_str(), ios::binary) { }
    virtual ~FileWriter() { close(); }

    bool is_open() 
    { return of.is_open(); }

    void close() 
    { if (is_open()) of.close(); }

    virtual inline void write(unsigned int val, size_t &offs)
    { of.write((const char *)&(val), sizeof(val)); }

    virtual inline void write(float val, size_t &offs)
    { of.write((const char *)&(val), sizeof(val)); }

    virtual inline void write(double val, size_t &offs) 
    { write((float)val, offs); }

    virtual inline void write(const char *val, size_t val_size, size_t &offs)
    { of.write(val, val_size); }
};

class SingleFileWriter : public FileWriter
{
public:

    SingleFileWriter(string filename) : FileWriter(filename) { }
    ~SingleFileWriter() { };

    inline void write(unsigned int val, size_t &offs)
    { of.seekp(offs); FileWriter::write(val, offs); offs = of.tellp(); }

    inline void write(float val, size_t &offs)
    { of.seekp(offs); FileWriter::write(val, offs); offs = of.tellp(); }

    inline void write(const char *val, size_t val_size, size_t &offs)
    { of.seekp(offs); FileWriter::write(val, val_size, offs); offs = of.tellp(); }
};

class Block 
{
public:
    const char *ext;
    char id[4];
    size_t offs;
    unsigned int size;

    FileWriter *fw;

    Block(const char *bext, const char *bid) : ext(bext), offs(0), size(0)
    { memcpy(id, bid, 4); }

    bool is_open() 
    { return fw->is_open(); }

    void close() 
    { if (is_open()) fw->close(); }

    virtual inline void write(unsigned int val)
    { fw->write((const char *)&(val), sizeof(val), offs); }

    virtual inline void write(float val)
    { fw->write((const char *)&(val), sizeof(val), offs); } //cerr << id << " " << offs << endl;}

    virtual inline void write(double val) 
    { write((float)val); }

    virtual inline void write(const char *val, size_t val_size)
    { fw->write(val, val_size, offs); }

    //~Block() { delete fw; }
};

Block block[IO_NBLOCKS+1] = 
{ 
    Block(".pos",  "POS "),
    Block(".vel",  "VEL "),
    Block(".id",   "ID  "),
    Block(".mass", "MASS"),
    Block(".u",    "U   "),
    Block(".rho",  "RHO "),
    Block(".hsml", "HSML"),
    Block(".pot",  "POT "),
    Block(".acce", "ACCE"),
    Block(".endt", "ENDT"),
    Block(".tstp", "TSTP"),
    Block(".head", "HEAD") 
};

void create_gadget_header(TipsyHeader &h, Config &cfg, struct io_header &iohdr);
int read_config_file(char *file, Config &cfg);
double convertToDouble(const std::string& s, bool failIfLeftoverChars = true);
int convertToInt(const std::string& s, bool failIfLeftoverChars = true);
void help() __attribute__ ((noreturn));

//#define SCALE_POS(p)  (((p) + 0.5F)*(1.0e3F*iohdr.HubbleParam)*fact_scale*(1.0F+iohdr.redshift))
#define SCALE_POS(p0, p1)  \
    do { \
        const float _s = (1.0e3F*iohdr.HubbleParam)*fact_scale*(1.0F+iohdr.redshift); \
        p1[0] = (p0[0] + 0.5F) * _s; \
        p1[1] = (p0[1] + 0.5F) * _s; \
        p1[2] = (p0[2] + 0.5F) * _s; \
    } while (0)
        //cerr << _s << " " << fact_scale << " " << p1[0] << " " << p1[1] << " " << p1[2] << endl;\

//#define SCALE_VEL(v)  ((v) / sqrt(iohdr.time) * fact_vel)
#define SCALE_VEL(v0, v1)  \
    do { \
        v1[0] = (v0[0]) / sqrt(iohdr.time) * fact_vel; \
        v1[1] = (v0[1]) / sqrt(iohdr.time) * fact_vel; \
        v1[2] = (v0[2]) / sqrt(iohdr.time) * fact_vel; \
    } while (0)
#define SCALE_MASS(m) ((m) * cfg.RhoCritical * iohdr.HubbleParam * fact_mass)
#define SCALE_TEMP(t) ((t) * 1.5F * kel2kms2)
#define SCALE_RHO(r)  ((r) * rho_z / 1.0e10F)

void create_gadget_header(TipsyHeader &h, Config &cfg, struct io_header &iohdr)
{
    iohdr.npart[GAS  ] = iohdr.npartTotal[GAS  ] = h.h_nSph;     // Gas
    iohdr.npart[HALO ] = iohdr.npartTotal[HALO ] = h.h_nDark;    // Halo
    iohdr.npart[DISK ] = iohdr.npartTotal[DISK ] = 0;            // Disk
    iohdr.npart[BULGE] = iohdr.npartTotal[BULGE] = 0;            // Bulge
    iohdr.npart[STARS] = iohdr.npartTotal[STARS] = h.h_nStar;    // Stars
    iohdr.npart[BNDRY] = iohdr.npartTotal[BNDRY] = 0;            // Bndry

    iohdr.mass[GAS  ] = cfg.mass[GAS  ];
    iohdr.mass[HALO ] = cfg.mass[HALO ];
    iohdr.mass[DISK ] = cfg.mass[DISK ];
    iohdr.mass[BULGE] = cfg.mass[BULGE];
    iohdr.mass[STARS] = cfg.mass[STARS];
    iohdr.mass[BNDRY] = cfg.mass[BNDRY];

    iohdr.time             = (cfg.time > 0) ? cfg.time : h.h_time;
    iohdr.redshift         = (iohdr.time == 0) ? iohdr.time : (1./iohdr.time) - 1.;

    iohdr.num_files        = 1;

    iohdr.BoxSize          = cfg.BoxSize;
    iohdr.Omega0           = cfg.Omega0;
    iohdr.OmegaLambda      = cfg.OmegaLambda;
    iohdr.HubbleParam      = cfg.HubbleParam;

    iohdr.npartTotalHighWord[GAS  ] = (((long)iohdr.npart[GAS  ]) >> 32) & 0x00000000ffffffffL;
    iohdr.npartTotalHighWord[HALO ] = (((long)iohdr.npart[HALO ]) >> 32) & 0x00000000ffffffffL;
    iohdr.npartTotalHighWord[DISK ] = (((long)iohdr.npart[DISK ]) >> 32) & 0x00000000ffffffffL;
    iohdr.npartTotalHighWord[BULGE] = (((long)iohdr.npart[BULGE]) >> 32) & 0x00000000ffffffffL;
    iohdr.npartTotalHighWord[STARS] = (((long)iohdr.npart[STARS]) >> 32) & 0x00000000ffffffffL;
    iohdr.npartTotalHighWord[BNDRY] = (((long)iohdr.npart[BNDRY]) >> 32) & 0x00000000ffffffffL;;

    iohdr.flag_entropy_instead_u = 0;
    iohdr.flag_sfr               = 0;     // Unused
    iohdr.flag_feedback          = 0;     // Unused
    iohdr.flag_cooling           = 0;     // Unused
    iohdr.flag_stellarage        = 0;     // Unused
    iohdr.flag_metals            = 0;     // Unused
}

void show_header(struct io_header &iohdr)
{
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
}


double convertToDouble(const std::string& s, bool failIfLeftoverChars)
{
    std::istringstream i(s);
    double x;
    char c;
    if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
    {
        cerr << "ERROR: Problems converting " << s << " to a number." << endl;
        exit(1);
    }
        
    return x;
}

int convertToInt(const std::string& s, bool failIfLeftoverChars)
{
    std::istringstream i(s);
    int x;
    char c;
    if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
    {
        cerr << "ERROR: Problems converting " << s << " to a number." << endl;
        exit(1);
    }
        
    return x;
}

int read_config_file(char *file, Config &cfg)
{
    ifstream in;
    string key, val;
    string x;

    in.open(file, ios::in);

    while (in.good())
    {
        key.erase();
        val.erase();
        in >> key >> val;

        if (key.empty()) continue;

        if (key.at(0) == '#') continue;

        if (val.empty())
        {
            cerr << "WARNING: Missing value for " << key << "." << endl;
            continue;
        }

        if (key == "SnapFormat")   cfg.SnapFormat  = convertToInt(val);
        else if (key == "Omega0")       cfg.Omega0      = convertToDouble(val);
        else if (key == "OmegaLambda")  cfg.OmegaLambda = convertToDouble(val);
        else if (key == "BoxSize")      cfg.BoxSize     = convertToDouble(val);
        else if (key == "HubbleParam")  cfg.HubbleParam = convertToDouble(val);
        else if (key == "Mass[0]")      cfg.mass[0]     = convertToDouble(val);
        else if (key == "Mass[1]")      cfg.mass[1]     = convertToDouble(val);
        else if (key == "Mass[2]")      cfg.mass[2]     = convertToDouble(val);
        else if (key == "Mass[3]")      cfg.mass[3]     = convertToDouble(val);
        else if (key == "Mass[4]")      cfg.mass[4]     = convertToDouble(val);
        else if (key == "Mass[5]")      cfg.mass[5]     = convertToDouble(val);
        else if (key == "TipsyFile")    cfg.TipsyFile   = val; 
        else if (key == "OutputDir")    cfg.OutputDir   = val; 
        else if (key == "RhoCritical")  cfg.RhoCritical = convertToDouble(val); 
        else 
        {
            cerr << "WARNING: Unrecognized key " << key << "." << endl;
            continue;
        }
    }

#define ERR(var) \
    if (cfg.var <= 0) { cerr << #var " missing or invalid in config file." << endl; exit(1); }

    ERR(Omega0);
    ERR(OmegaLambda);
    ERR(HubbleParam);
    ERR(BoxSize);

    return 0;
}


void open_block_files(Config &cfg, bool useOneFile)
{
    string basename;

    basename = cfg.TipsyFile.substr(cfg.TipsyFile.find_last_of(DIR_SEP) + 1);

    if (!cfg.OutputDir.empty())
        basename = cfg.OutputDir + DIR_SEP + basename;

    if (useOneFile)
    {
        SingleFileWriter *out;

        basename += ".gad";

        VL(1) cerr << "Creating a single file called " << basename << endl;

        out = new SingleFileWriter(basename);
        if (!out->is_open())
        {
            cerr << "Unable to open " << basename << ". Aborting." << endl;
            exit(1);
        }

        for (int i=0; i < IO_NBLOCKS+1; i++)
            block[i].fw = out;

        /*====================================================================
         * If we use one file then we need to compute the offsets within that 
         * file where each of the blocks begin. We also create the proper 
         * file writer here. 
         *
         * If we use multiple files then the offsets within each file is zero.
         *==================================================================*/
        block[IO_HEAD].offs = 0;
        block[IO_HEAD].size = sizeof(struct io_header);

        int offs0 = 0;
        if (cfg.SnapFormat == 2) offs0 = 4;

        int offs = sizeof(unsigned int) + offs0 + block[IO_HEAD].size + sizeof(unsigned int);
        for (int i=0; i < IO_NBLOCKS; i++)
        {
            block[i].offs = offs;
            offs +=
                sizeof(unsigned int)
                + offs0
                + block[i].size
                + sizeof(unsigned int);

            VL(2) cerr << block[i].id << " " << block[i].offs << " " << block[i].size << endl;
        }
#if 0
        block[0].offs = sizeof(unsigned int) + block[IO_HEAD].size + sizeof(unsigned int);
        for (int i=1; i < IO_NBLOCKS; i++) 
        {
            block[i].offs 
                = block[i-1].offs 
                + sizeof(unsigned int)
                + block[i-1].size 
                + sizeof(unsigned int);

        }
#endif
    }
    else
    {
        VL(1) cerr << "Creating multiple files called " << endl;

        for (int i=0; i < IO_NBLOCKS+1; i++)
        {
            string blockfile(basename + block[i].ext);

            VL(1) cerr << "\t" << blockfile;
            VL(2) cerr << " " << block[i].id << " " << block[i].offs << " " << block[i].size;
            VL(1) cerr << endl;

            block[i].fw = new FileWriter(blockfile);

            if (!block[i].fw->is_open())
            {
                cerr << "ERROR: Unable to open " << blockfile << ". Aborting." << endl;
                exit(1);
            }

            block[i].offs = 0;
        }
    }
}

size_t parse_memory_arg(const char *arg)
{
    string a(arg);
    double mem;

    mem = convertToDouble(a, 0);

    return (size_t)(mem * 1024 * 1024);
}

void help() 
{
    //cerr << "tip2gad [-f tipsy file] [-v] [-q] <config file>" << endl
    cerr << "Usage: tip2gad [OPTIONS] <config file>" << endl
         << endl
         << "OPTIONS:" << endl
         << "    -f tipsy_file      Override the tipsy file specified in the config file." << endl
         << "    -m <N>             Convert file in blocks of up to N megabytes (default " 
         <<                                                               DEFAULT_MEM << ")."<<endl
         << "                       If the amount of available memory is larger than the size"<<endl
         << "                       of the tipsy file this can significantly decrease the" << endl
         << "                       running time because the whole file can be converted" << endl
         << "                       in one read-convert-write operation." << endl
         << "    -v                 Increase verbosity." << endl
         << "    -q                 Decrease verbosity (make quiet)." << endl
         << "    --single           Output a whole gadget file (default)." << endl
         << "    --multi            Output the individual parts of the gadget file." << endl
         << "    --version          Print version information." << endl
         << "    --help             This help screen." << endl
         << endl
         << "The output file will have the same name as in the input file plus '.gad' if" << endl
         << "a single output file is requested. The --multi option produces files with" << endl
         << "the name of the input file plus an extension for each of the gadget blocks." << endl
         << endl
         << "The configuration file has the following format. Options with defaults" << endl
         << "can be left out. TipsyFile can be left out only if the file is specified" << endl
         << "on the command line. Lines with '#' as the first character are ignored." << endl
         << endl
         << "TipsyFile      <file>   -- The tipsy file to convert"  << endl
         << "OutputDir      <dir>    -- The directory to place the gadget file (default: .)" << endl
         << "SnapFormat     [n]      -- n must be 1 or 2 (default: 1)" << endl
         << "Omega0         <n>      -- Value of Omega-naught" << endl
         << "OmegaLambda    <n>      -- Value of Omega lambda" << endl
         << "BoxSize        <n>      -- Length of one side of the box in Mpc" << endl
         << "HubbleParam    <n>      -- The value of the Hubble constant" << endl
         << "RhoCritical    <n>      -- Critical density (default: 1)" << endl
         << "Mass[0]        <n>      -- The mass of particle type 0 (default: 0)" << endl
         << "Mass[1]        <n>      -- The mass of particle type 1 (default: 0)" << endl
         << "Mass[2]        <n>      -- The mass of particle type 2 (default: 0)" << endl
         << "Mass[3]        <n>      -- The mass of particle type 3 (default: 0)" << endl
         << "Mass[4]        <n>      -- The mass of particle type 4 (default: 0)" << endl
         << "Mass[5]        <n>      -- The mass of particle type 5 (default: 0)" << endl
         << endl
         << "A mass type with a value of 0 means that the masses of each particle of that" << endl
         << "type are not the same and will have to be stored individually in the output file." << endl
         << endl
         << "Send questions, comments, bug reports to Jonathan Coles <jonathan@physik.uzh.ch>" << endl
         << endl;
    
    exit(2);
}

void version()
{
    cerr << "tip2gad v" << VERSION << endl
         << "Tipsy to Gadget file converter." << endl
         << "Written by Jonathan Coles <jonathan@physik.uzh.ch>" << endl;

    exit(0);
}

int main(int argc, char **argv)
{
    ifTipsy in;
    TipsyHeader h;
    struct io_header iohdr;

    bool useOneFile = true;

    Config cfg;

    static struct option long_options[] = {
        {"single",  no_argument, 0, 0},
        {"multi",  no_argument, 0, 0},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 0},
        {0, 0, 0, 0}
    };

    char *tipsy_file = NULL;

    size_t use_mem = DEFAULT_MEM * 1024 * 1024;

    if (argc < 2) help();

    /*========================================================================
     * Process the command line flags
     *======================================================================*/
    while (1)
    {
        int option_index = 0;
        int c = getopt_long(argc, argv, "hf:vqm:",
                            long_options, &option_index);

        if (c == -1) break;

        switch (c)
        {
            case 0:
                if (!strcmp("single", long_options[option_index].name))
                    useOneFile = true;
                else if (!strcmp("multi", long_options[option_index].name))
                    useOneFile = false;
                else if (!strcmp("version", long_options[option_index].name))
                    version();
                break;

            case 'f': tipsy_file = optarg; break;
            case 'h': help(); break;
            case 'v': verbosity++; break;
            case 'q': verbosity = 0; break;
            case 'm': use_mem = parse_memory_arg(optarg); break;
            case '?': break;
        }
    }

    if (optind >= argc) help();

    size_t pbuffer_size = use_mem / ((3+3+1+1+1+1+1+1+1) * sizeof(float));

    if (pbuffer_size == 0)
    {
        cerr << "Memory size given is too small." << endl;
        exit(1);
    }

    VL(2) cout << "Using " << use_mem << " bytes to load particles." << endl;
    VL(2) cout << "Buffers are " << pbuffer_size << " elements long." << endl;

    if (read_config_file(argv[optind], cfg))
    {
        cerr << "ERROR: Can't read config file " << argv[optind] << ". Aborting." << endl;
        exit(1);
    }

    if (tipsy_file != NULL) cfg.TipsyFile = tipsy_file;

    if (cfg.TipsyFile.empty())
    {
        cerr << "No Tipsy file specified in config file or on command line." << endl;
        exit(1);
    }

    in.open(cfg.TipsyFile.c_str(), "standard");
    if (!in.is_open()) 
    {
        in.open(cfg.TipsyFile.c_str(), "native");
        if (!in.is_open()) 
        {
            cerr << "ERROR: Unable to open Tipsy binary " << cfg.TipsyFile << endl;
            exit(2);
        }
    }

    /*========================================================================
     * Read in the tipsy header and convert extract as much informaton as
     * we can into the gadget header format.
     *======================================================================*/
    in >> h; 
    create_gadget_header(h, cfg, iohdr);

    block[IO_HEAD].size     = sizeof(iohdr);
    block[IO_POS].size      = h.h_nBodies * 3 * sizeof(float);
    block[IO_VEL].size      = h.h_nBodies * 3 * sizeof(float);
    block[IO_ID].size       = h.h_nBodies * 1 * sizeof(unsigned int);
    block[IO_U].size        = h.h_nSph    * 1 * sizeof(float);
    block[IO_RHO].size      = h.h_nSph    * 1 * sizeof(float);
    block[IO_HSML].size     = h.h_nSph    * 1 * sizeof(float);
    block[IO_POT].size      = h.h_nBodies * 1 * sizeof(float);
    block[IO_ACCEL].size    = 0; //h.h_nBodies * 3 * sizeof(float);
    block[IO_DTENTR].size   = 0; //h.h_nSph    * 1 * sizeof(float);
    block[IO_TSTP].size     = 0; //h.h_nBodies * 1 * sizeof(float);

    block[IO_MASS].size = 0;
    if (iohdr.mass[GAS]   == 0) block[IO_MASS].size += h.h_nSph  * 1 * sizeof(float);
    if (iohdr.mass[HALO]  == 0) block[IO_MASS].size += h.h_nDark * 1 * sizeof(float);
    if (iohdr.mass[STARS] == 0) block[IO_MASS].size += h.h_nStar * 1 * sizeof(float);

    cerr << block[IO_MASS].size << endl;
    cerr << "****" << endl;

    /*========================================================================
     * Conversion values for scaling between units in tipsy files and units
     * in gadget files.
     *======================================================================*/
    float H0         = 100. * iohdr.HubbleParam;
    float kel2kms2   = BOLTZMANN/(0.6*PROTONMASS)*1.e-10;
    float rho_crit   = 1.9e-29 * pow(iohdr.HubbleParam,2);      /*! Units of g/cm^3 for Omega_m=1 */
    float rho_bar    = 0.019 * 1.9e-29 *pow(iohdr.HubbleParam,2); /*! Baryon density (indep of H) */
    float rho_msol   = 2.7755*10.*pow(iohdr.HubbleParam,2);       /*! Units of Msol/Mpc^3  */
    float volume     = pow((cfg.BoxSize/1.0e3/iohdr.HubbleParam),3);
    float fact_mass  = rho_msol * volume;
    float fact_scale = cfg.BoxSize/1.0e3/iohdr.HubbleParam/(1.+iohdr.redshift);
    float fact_vel   = (cfg.BoxSize/1.0e3/iohdr.HubbleParam * H0 / 2.894405) / (1+iohdr.redshift);
    float rho_z      = 1.e10*rho_msol/pow(iohdr.HubbleParam,2)/1.e9; 

    /*========================================================================
     * Open the output files. If we are not using a single output file then
     * a seperate file will be created for each gadget i/o block.
     *======================================================================*/
    open_block_files(cfg, useOneFile);

    /*========================================================================
     * Write each block size and optionally the block marker. 
     *======================================================================*/
    for (int i=0; i < IO_NBLOCKS+1; i++) 
    {
        if (block[i].size) 
        {
            /*========================================================================
             * If we are writing in SnapFormat >= 2 then each block has a marker
             * that is four bytes long.
             *======================================================================*/
            if (cfg.SnapFormat >= 2) 
            {
                block[i].write(block[i].size + 4);
                block[i].write(block[i].id, 4);
            }
            else
            {
                block[i].write(block[i].size);
            }
        }
    }

    VL(1) show_header(iohdr);


    /*========================================================================
     * Write all the particle information. 
     *
     * Gadget files expect gas particles first, then dark particles, then
     * stars.
     *======================================================================*/

    unsigned int two_percent = (int)(h.h_nBodies * 0.02);

    VL(1) cerr << "   ""                                                   ] 100%\r";
    VL(1) cerr << "0% [";

    //float mass = -1;
    bool sameMass = true;

    float *pos  = new float[3 * pbuffer_size], *pos_p   = pos;
    float *vel  = new float[3 * pbuffer_size], *vel_p   = vel;
    float *mass = new float[1 * pbuffer_size], *mass_p  = mass;
    float *u    = new float[1 * pbuffer_size], *u_p     = u;
    float *rho  = new float[1 * pbuffer_size], *rho_p   = rho;
    float *temp = new float[1 * pbuffer_size], *temp_p  = temp;
    float *hsml = new float[1 * pbuffer_size], *hsml_p  = hsml;
    float *pot  = new float[1 * pbuffer_size], *pot_p   = pot;
    unsigned int *id  = new unsigned int[1 * pbuffer_size], *id_p    = id;

    TipsyGasParticle  g;
    TipsyDarkParticle d;
    TipsyStarParticle s;

    int read_state = READ_GAS;

    for(size_t i=0, read=0; read_state != READ_DONE; )
    {
        if (i == h.h_nSph)                         read_state = READ_DARK;
        if (i == h.h_nSph + h.h_nDark)             read_state = READ_STARS;
        if (i == h.h_nSph + h.h_nDark + h.h_nStar) read_state = WRITE_DATA;

        if (read == pbuffer_size) read_state |= WRITE_DATA;

        switch (read_state)
        {
            case READ_GAS:
                in >> g;
                SCALE_POS(g.pos, pos_p);        pos_p  += 3;
                SCALE_VEL(g.vel, vel_p);        vel_p  += 3;
                *id_p     = i;                  id_p   += 1;
                *temp_p   = SCALE_TEMP(g.temp); temp_p += 1;
                *rho_p    = SCALE_RHO(g.rho);   rho_p  += 1;
                *hsml_p   = g.hsmooth;          hsml_p += 1;
                *pot_p    = g.phi;              pot_p  += 1;

                if (iohdr.mass[GAS] == 0) 
                { *mass_p = SCALE_MASS(g.mass); mass_p += 1; }

                break;

            case READ_DARK:
                in >> d;
                SCALE_POS(d.pos, pos_p);        pos_p  += 3;
                SCALE_VEL(d.vel, vel_p);        vel_p  += 3;
                *id_p     = i;                  id_p   += 1;
                *pot_p    = d.phi;              pot_p  += 1;

                if (iohdr.mass[HALO] == 0) 
                { *mass_p = SCALE_MASS(d.mass); mass_p += 1; }

                break;

            case READ_STARS:
                in >> s;
                SCALE_POS(s.pos, pos_p);        pos_p  += 3;
                SCALE_VEL(s.vel, vel_p);        vel_p  += 3;
                *id_p     = i;                  id_p   += 1;
                *pot_p    = s.phi;              pot_p  += 1;

                if (iohdr.mass[STARS] == 0) 
                { *mass_p = SCALE_MASS(s.mass); mass_p += 1; }

                break;

            case WRITE_DATA:
            case (READ_STARS | WRITE_DATA):
            case (READ_DARK  | WRITE_DATA):
            case (READ_GAS   | WRITE_DATA):
#define WRITE_GROUP(ptr, ptrptr, label, type) \
    do { \
        ptrdiff_t l=ptrptr-ptr; \
        if (l) block[label].write((char *)ptr, l*sizeof(type));\
        ptrptr = ptr;\
    } while(0)

                WRITE_GROUP(pos,  pos_p,  IO_POS,  float);
                WRITE_GROUP(vel,  vel_p,  IO_VEL,  float);
                WRITE_GROUP(id,   id_p,   IO_ID,   unsigned int);
                WRITE_GROUP(mass, mass_p, IO_MASS, float);
                WRITE_GROUP(u,    u_p,    IO_U,    float);
                WRITE_GROUP(rho,  rho_p,  IO_RHO,  float);
                WRITE_GROUP(hsml, hsml_p, IO_HSML, float);
                WRITE_GROUP(pot,  pot_p,  IO_POT,  float);

                read = 0;

                read_state &= ~WRITE_DATA;
                if (i == h.h_nBodies) read_state = READ_DONE;

                continue;

            case READ_DONE:
                break;

            default:
                cerr << "read_state = " << read_state << endl;
                assert(0);
                break;

        }

        VL(1) if ((i % two_percent) == 0) cerr << ".";

        i++; 
        read++;
    }

    
#if 0
    for(unsigned int i=0; i < h.h_nSph; i += read, id += read)
    { 
        for (unsigned int read=0; read < pbuffer_size && i < h.h_nSph; read++)
        {
            in >> g[read]; 
            SCALE_POS(g.pos[read]);
            SCALE_VEL(g.vel[read]);
        }

        for (unsigned int j=0; j < read; j++)
        {
            pos[j][0] = SCALE_POS(g.pos[0]);
            pos[j][1] = SCALE_POS(g.pos[1]);
            pos[j][2] = SCALE_POS(g.pos[2]);
        }
        block[IO_POS].write((char *)&(pos), sizeof(pos));



        //SCALE_POS(g.pos);
        block[IO_POS].write((char *)&(g.pos), sizeof(g.pos));

        //SCALE_VEL(g.vel);
        block[IO_VEL].write((char *)&(g.vel), sizeof(g.vel));

        block[IO_ID].write(id);

        if (iohdr.mass[GAS] <= 0) 
        {
            float scl_mass = SCALE_MASS(g.mass);
            if (mass < 0) mass = scl_mass;
            if (mass != scl_mass) sameMass = false;

            block[IO_MASS].write(SCALE_MASS(g.mass));
        }

        block[IO_U].write(SCALE_TEMP(g.temp));
        block[IO_RHO].write(SCALE_RHO(g.rho));
        block[IO_HSML].write(g.hsmooth);

        block[IO_POT].write(g.phi);

        VL(1) if ((id % two_percent) == 0) cerr << ".";
    }

#if 0
    if (sameMass)
    {
        block[IO_MASS].zero
    }
#endif

    for(unsigned int i=0; i < h.h_nDark; i++, id++) 
    { 
        TipsyDarkParticle d;

        in >> d; 
        SCALE_POS(d.pos);
        block[IO_POS].write((char *)&(d.pos), sizeof(d.pos));

        SCALE_VEL(d.vel);
        block[IO_VEL].write((char *)&(d.vel), sizeof(d.vel));

        block[IO_ID].write(id);

        if (iohdr.mass[HALO] == 0) block[IO_MASS].write(SCALE_MASS(d.mass));

        block[IO_POT].write(d.phi);

        VL(1) if ((id % two_percent) == 0) cerr << ".";
    }

    for(unsigned int i=0; i < h.h_nStar; i++, id++) 
    { 
        TipsyStarParticle s;

        in >> s; 
        SCALE_POS(s.pos);
        block[IO_POS].write((char *)&(s.pos), sizeof(s.pos));

        SCALE_VEL(s.vel);
        block[IO_VEL].write((char *)&(s.vel), sizeof(s.vel));

        block[IO_ID].write(id);

        if (iohdr.mass[STARS] == 0) block[IO_MASS].write(SCALE_MASS(s.mass));

        block[IO_POT].write(s.phi);

        VL(1) if ((id % two_percent) == 0) cerr << ".";
    }
#endif

    VL(1) cerr << endl;


    /*========================================================================
     * Write the gadget header.
     *======================================================================*/
    block[IO_HEAD].write((const char *)&iohdr, (int)sizeof(iohdr));

    /*========================================================================
     * Finish up with the size of the block again. This is the end marker
     * expected by FORTRAN.
     *======================================================================*/
    for (int i=0; i < IO_NBLOCKS+1; i++) 
    {
        if (block[i].size) 
        {
            if (cfg.SnapFormat >= 2) 
                block[i].write(block[i].size + 4);
            else
                block[i].write(block[i].size);
        }
    }

    VL(1) cerr << "Closing files";

    /*========================================================================
     * Close all files and release memory.
     *======================================================================*/
    for (int i=0; i < IO_NBLOCKS+1; i++) 
    {
        if (block[i].fw != NULL && block[i].is_open()) 
        {
            block[i].close();
            delete block[i].fw;
            block[i].fw = NULL;
            //delete block[i];
        }

        VL(1) cerr << ".";
    }

    VL(1) cerr << endl;

    in.close();

    delete pos;
    delete vel;
    delete id;
    delete mass;
    delete u;
    delete rho;
    delete temp;
    delete hsml;
    delete pot;

    return 0;
}

