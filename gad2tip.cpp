#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <getopt.h>
#include <assert.h>

#include "ftipsy.hpp"

using namespace std;

#define VERSION "0.8"

#define DEFAULT_MEM 128

#define GADGET_FORMAT_DETECT 0
#define GADGET_FORMAT_1      1
#define GADGET_FORMAT_2      2

#define  PROTONMASS        1.6726e-24
#define  BOLTZMANN         1.3806e-16

enum { GAS, HALO, DISK, BULGE, STARS, BNDRY };

int verbosity = 1;
#define VL(__l) if (verbosity >= __l)

#define MIN(a,b) (((a) < (b)) ? (a) : (b))


#define WRITE_GAS   1
#define WRITE_DARK  2
#define WRITE_STARS 4
#define READ_DATA   8
#define READ_DONE  16


class Config
{
public:
    unsigned int SnapFormat;
    float RhoCritical;
    bool sort;
    bool CosmoScaling;

    Config() :
        SnapFormat(0),
        RhoCritical(1.0F),
        sort(false),
        CosmoScaling(false)
    {
    }
};

class GadgetBaseParticle
{
public:
    float pos[3];
    float vel[3];
    uint32_t id;
    float mass;
    float pot;
    float acc[3];
    float dt;
};

class GadgetGasParticle : public GadgetBaseParticle
{
public:
    float u;
    float rho;
    float hsml;
    float endt;
};

class GadgetHaloParticle : public GadgetBaseParticle
{
};

class GadgetDiskParticle : public GadgetBaseParticle
{
};

class GadgetBulgeParticle : public GadgetBaseParticle
{
};

class GadgetStarParticle : public GadgetBaseParticle
{
};

class GadgetBoundaryParticle : public GadgetBaseParticle
{
};

class GadgetHeader
{
public:
    uint32_t npart[6];                        /*!< number of particles of each type in this file */
    double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
                                            stored in the mass-block of the snapshot file, otherwise they are omitted */
    double   time;                         /*!< time of snapshot file */
    double   redshift;                     /*!< redshift of snapshot file */
    uint32_t flag_sfr;                        /*!< flags whether the simulation was including star formation */
    uint32_t flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
    uint32_t npartTotal[6];          /*!< total number of particles of each type in this snapshot. This can be
                                            different from npart if one is dealing with a multi-file snapshot. */
    uint32_t flag_cooling;                    /*!< flags whether cooling was included  */
    uint32_t num_files;                       /*!< number of files in multi-file snapshot */
    double   BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
    double   Omega0;                       /*!< matter density in units of critical density */
    double   OmegaLambda;                  /*!< cosmological constant parameter */
    double   HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
    uint32_t flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
    uint32_t flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
    uint32_t npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
    uint32_t flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
    char     fill[60];	               /*!< fills to 256 Bytes */
};

class GadgetFile;

class BaseBlock
{
protected:

    GadgetFile *gf;

public:
    virtual void setGadgetFile(GadgetFile *_gf)
    {
        this->gf = _gf;
    }
    virtual ~BaseBlock() {}

    virtual bool empty() { return false; }
    virtual size_t getByteSize() = 0;
    virtual char *getID() = 0;
    virtual void setFileOffset(size_t offs) = 0;
    virtual void setBufferSize(int _bufsize) = 0;
    virtual void setNumberOfElements(int _N) = 0;
    virtual bool verify() = 0;
    virtual void fill(GadgetGasParticle *p) = 0;
    virtual void fill(GadgetHaloParticle *p) = 0;
    virtual void fill(GadgetDiskParticle *p) = 0;
    virtual void fill(GadgetBulgeParticle *p) = 0;
    virtual void fill(GadgetStarParticle *p) = 0;
    virtual void fill(GadgetBoundaryParticle *p) = 0;
};

template <class T, int Nelem>
class Block : public BaseBlock
{
    const char *ext;
    char id[5];
    uint32_t size;
    unsigned int have_read;
    size_t bufsize;
    size_t buf_offset;
    T *buffer;

    size_t block_offset;;
    size_t block_start_offset;

protected:
    uint32_t N;
    virtual void read_next(T *);

public:
    Block() {}
    Block(const char *bid, uint32_t _N)
        : have_read(0),
          bufsize(0),
          buf_offset(0),
          buffer(NULL),
          block_offset(0),
          N(_N)
    {
        strcpy(id, bid);
    }

    virtual ~Block() 
    {
        if (buffer != NULL) free(buffer);
    }

    virtual void setFileOffset(size_t offs)
    {
        block_start_offset = offs;
    }

    virtual void setBufferSize(int _bufsize)
    {
        const size_t bytes = Nelem * sizeof(T);
        assert(_bufsize >= bytes);
        bufsize = _bufsize / bytes;
        buf_offset = bufsize;
        buffer = (T *)realloc(buffer, bufsize * bytes);
        assert(buffer != NULL);
        VL(3) cout << "Set block " << id << " buffer size to " << (bufsize*bytes) << endl;
    }

    virtual void setNumberOfElements(int _N)
    {
        N = _N;
    }

    virtual size_t getByteSize()
    {
        return N * Nelem * sizeof(T);
    }

    virtual char *getID()
    {
        return id;
    }

    virtual void fill(GadgetGasParticle *p)  {}
    virtual void fill(GadgetHaloParticle *p)  {}
    virtual void fill(GadgetDiskParticle *p)  {}
    virtual void fill(GadgetBulgeParticle *p)  {}
    virtual void fill(GadgetStarParticle *p)  {}
    virtual void fill(GadgetBoundaryParticle *p) {}
    virtual bool verify();
};

class GadgetFile
{
private:
    //friend template<class T, Nelem> class Block<T,Nelem>;
    vector<BaseBlock*> blocks;

public:
    size_t buflen;
    int current_block_offset;
    int SnapFormat;
    ifstream *is;
    GadgetHeader h;

    uint32_t nBodies;

public:
    GadgetFile(char *filename, int snap_format=0, int _buflen=DEFAULT_MEM * 1024 * 1024) 
        : buflen(_buflen),
          current_block_offset(0), 
          SnapFormat(snap_format)
    {
        is = new ifstream(filename, ios::in | ios::binary);

        uint32_t len0, len1;
        char id[5];
        is->read((char *)&len0, sizeof(len0));
        if (SnapFormat == 0)
        {
            is->read(id,    sizeof(id)-1);
            id[4] = 0;
            SnapFormat = 1;
            if (!strcmp(id, "HEAD")) SnapFormat = 2;
            is->seekg(sizeof(len0));
            VL(1) cout << "Auto-detected GADGET file format " << SnapFormat << endl;
        }

        if (SnapFormat == 2)
        {
            is->read(id, sizeof(id)-1);
            id[4] = 0;
            if (strcmp(id, "HEAD")) throw 2;
        }

        is->read((char *)&h,  sizeof(h));
        is->read((char *)&len1, sizeof(len1));
        if (len0 != len1) throw 1;
        current_block_offset = len0 + 2*sizeof(uint32_t); // + (SnapFormat==2?4:1);

        nBodies = h.npart[GAS] 
                + h.npart[HALO] 
                + h.npart[STARS] 
                + h.npart[DISK] 
                + h.npart[BULGE];
    }

    virtual ~GadgetFile()
    {
        close();
    }

    virtual void close()
    {
        if (is != NULL && is->is_open()) is->close();
    }

    virtual inline bool is_open()
    {
        return is != NULL && is->is_open();
    }

    virtual bool addBlock(BaseBlock *b)
    {
        b->setGadgetFile(this);

        if (!b->empty())
        {
            b->setFileOffset(current_block_offset);
#if 1
            if (!b->verify())
            {
                cerr << "Block " << b->getID() << " does not appear to exist in file!" << endl;
                return false;
            }
#endif

            VL(2) cout << "Adding block " << b->getID() << endl;

            is->clear();
            current_block_offset += b->getByteSize() + 2*sizeof(uint32_t) + (SnapFormat==2?4:0);;

            blocks.push_back(b);
        }
        return true;
    }

    virtual inline ifstream *getStream() { return is; }

    virtual void setBufferSize(size_t size)
    {
        if (blocks.size() > 0)
        {
            size_t *block_sizes = new size_t[blocks.size()];
            size_t max_size = 0;
            int i=0;
            for (vector<BaseBlock*>::const_iterator iter = blocks.begin(); iter != blocks.end(); ++iter, i++)
            { 
                block_sizes[i] = (*iter)->getByteSize();
                cerr << block_sizes[i] << endl;
                if (block_sizes[i] > max_size) max_size = block_sizes[i];
            }

#if 0
            for (i=0; i < blocks.size(); i++)
            {
                block_sizes[i] /= max_size * blocks.size();
                cerr << block_sizes[i] << endl;
            }
#endif

            i=0;
            for (vector<BaseBlock*>::const_iterator iter = blocks.begin(); iter != blocks.end(); ++iter, i++)
            { 
                (*iter)->setBufferSize(size * block_sizes[i] / max_size / blocks.size());
            }

            delete block_sizes;
        }
    }

#define DO_FILL(p) \
    do {  \
        for (vector<BaseBlock*>::const_iterator iter = blocks.begin(); iter != blocks.end(); ++iter) \
        { (*iter)->fill(&(p)); } \
    } while(0)

    virtual GadgetFile& operator>>(GadgetGasParticle &p)      { DO_FILL(p); return *this; }
    virtual GadgetFile& operator>>(GadgetHaloParticle &p)     { DO_FILL(p); return *this; }
    virtual GadgetFile& operator>>(GadgetDiskParticle &p)     { DO_FILL(p); return *this; }
    virtual GadgetFile& operator>>(GadgetBulgeParticle &p)    { DO_FILL(p); return *this; }
    virtual GadgetFile& operator>>(GadgetStarParticle &p)     { DO_FILL(p); return *this; }
    virtual GadgetFile& operator>>(GadgetBoundaryParticle &p) { DO_FILL(p); return *this; }
};

template <class T, int Nelem> 
void Block<T,Nelem>::read_next(T *b)
{
    //cerr << have_read << endl;
    //if (have_read == getByteSize()) return 0;
    assert (bufsize != 0);
    if (buf_offset + Nelem > bufsize) /* > because the buffer many have been resized */
    {
        size_t left_to_read = N - have_read;
        size_t read = MIN(bufsize, left_to_read);
        gf->is->seekg(block_offset);
        gf->is->read((char *)buffer, read * Nelem * sizeof(T));
        block_offset += read * Nelem * sizeof(T);
        buf_offset = 0;
    }

    for (int i=0; i < Nelem; i++)
    {
        have_read++;
        b[i] = buffer[buf_offset++];
    }
}

template <class T, int Nelem> 
bool Block<T,Nelem>::verify()
{
    char bid[5];
    uint32_t len0, len1;
    VL(3) cout << "block_start_offset: " << block_start_offset << endl;
    gf->is->seekg(block_start_offset);
    gf->is->read((char *)&len0, sizeof(len0));
    VL(3) cout << "block length (1st marker): " << len0 << endl;
    if (gf->SnapFormat == 2)
    {
        gf->is->read((char *)bid, sizeof(id)-1);
        bid[4] = 0;
        cerr << "'" << bid << "'" << endl;
        cerr << "'" << id << "'" << endl;
        if (strcmp(id, bid)) return false;
    }
    block_offset = gf->is->tellg();
    gf->is->seekg(block_start_offset + len0 + sizeof(len0));
    gf->is->read((char *)&len1, sizeof(len1));
    VL(3) cout << "block length (2nd marker): " << len1 << endl;

    return len0 == len1;
}

class PosBlock : public Block<float, 3>
{
public:
    PosBlock() : Block<float,3>("POS ", 0) { }
    virtual void setGadgetFile(GadgetFile *_gf) 
    { 
        BaseBlock::setGadgetFile(_gf); 
        setNumberOfElements(_gf->nBodies);
    }
    virtual void fill(GadgetGasParticle *p)      { read_next(p->pos); }
    virtual void fill(GadgetHaloParticle *p)     { read_next(p->pos); }
    virtual void fill(GadgetDiskParticle *p)     { read_next(p->pos); }
    virtual void fill(GadgetBulgeParticle *p)    { read_next(p->pos); }
    virtual void fill(GadgetStarParticle *p)     { read_next(p->pos); }
    virtual void fill(GadgetBoundaryParticle *p) { read_next(p->pos); }
};

class VelBlock : public Block<float,3>
{
public:
    VelBlock() : Block<float,3>("VEL ", 0) { }
    virtual void setGadgetFile(GadgetFile *_gf) 
    { 
        BaseBlock::setGadgetFile(_gf); 
        setNumberOfElements(_gf->nBodies);
    }
    virtual void fill(GadgetGasParticle *p)      { read_next(p->vel); }
    virtual void fill(GadgetHaloParticle *p)     { read_next(p->vel); }
    virtual void fill(GadgetDiskParticle *p)     { read_next(p->vel); }
    virtual void fill(GadgetBulgeParticle *p)    { read_next(p->vel); }
    virtual void fill(GadgetStarParticle *p)     { read_next(p->vel); }
    virtual void fill(GadgetBoundaryParticle *p) { read_next(p->vel); }
};

class IDBlock : public Block<uint32_t,1>
{
public:
    IDBlock() : Block<uint32_t,1>("ID  ", 0) { }
    virtual void setGadgetFile(GadgetFile *_gf) 
    { 
        BaseBlock::setGadgetFile(_gf); 
        setNumberOfElements(_gf->nBodies);
    }
    virtual void fill(GadgetGasParticle *p)      { read_next(&(p->id)); }
    virtual void fill(GadgetHaloParticle *p)     { read_next(&(p->id)); }
    virtual void fill(GadgetDiskParticle *p)     { read_next(&(p->id)); }
    virtual void fill(GadgetBulgeParticle *p)    { read_next(&(p->id)); }
    virtual void fill(GadgetStarParticle *p)     { read_next(&(p->id)); }
    virtual void fill(GadgetBoundaryParticle *p) { read_next(&(p->id)); }
};

class MassBlock : public Block<float,1>
{
    int nPart;
public:
    MassBlock() : Block<float,1>("MASS", 0) { }
    virtual void setGadgetFile(GadgetFile *_gf) 
    { 
        BaseBlock::setGadgetFile(_gf); 
        nPart = 0;
        nPart += (gf->h.mass[GAS]   == 0) ? gf->h.npart[GAS]   : 0;
        nPart += (gf->h.mass[HALO]  == 0) ? gf->h.npart[HALO]  : 0;
        nPart += (gf->h.mass[DISK]  == 0) ? gf->h.npart[DISK]  : 0;
        nPart += (gf->h.mass[BULGE] == 0) ? gf->h.npart[BULGE] : 0;
        nPart += (gf->h.mass[STARS] == 0) ? gf->h.npart[STARS] : 0;
        nPart += (gf->h.mass[BNDRY] == 0) ? gf->h.npart[BNDRY] : 0;
        setNumberOfElements(nPart);
    }
    virtual void fill(GadgetGasParticle *p)      { if (gf->h.mass[GAS]   == 0) read_next(&(p->mass)); else p->mass = gf->h.mass[GAS];   }
    virtual void fill(GadgetHaloParticle *p)     { if (gf->h.mass[HALO]  == 0) read_next(&(p->mass)); else p->mass = gf->h.mass[HALO];  }
    virtual void fill(GadgetDiskParticle *p)     { if (gf->h.mass[DISK]  == 0) read_next(&(p->mass)); else p->mass = gf->h.mass[DISK];  }
    virtual void fill(GadgetBulgeParticle *p)    { if (gf->h.mass[BULGE] == 0) read_next(&(p->mass)); else p->mass = gf->h.mass[BULGE]; }
    virtual void fill(GadgetStarParticle *p)     { if (gf->h.mass[STARS] == 0) read_next(&(p->mass)); else p->mass = gf->h.mass[STARS]; }
    virtual void fill(GadgetBoundaryParticle *p) { if (gf->h.mass[BNDRY] == 0) read_next(&(p->mass)); else p->mass = gf->h.mass[BNDRY]; }
    virtual size_t getByteSize()
    {
        int n=0;
        n += (gf->h.mass[GAS]   == 0) ? gf->h.npart[GAS]   : 0;
        n += (gf->h.mass[HALO]  == 0) ? gf->h.npart[HALO]  : 0;
        n += (gf->h.mass[DISK]  == 0) ? gf->h.npart[DISK]  : 0;
        n += (gf->h.mass[BULGE] == 0) ? gf->h.npart[BULGE] : 0;
        n += (gf->h.mass[STARS] == 0) ? gf->h.npart[STARS] : 0;
        n += (gf->h.mass[BNDRY] == 0) ? gf->h.npart[BNDRY] : 0;
        return n * 1 * sizeof(float);
    }
    virtual bool empty() { return nPart == 0; }
};

class UBlock : public Block<float,1>
{
public:
    UBlock() : Block<float,1>("U   ", 0) { }
    virtual void setGadgetFile(GadgetFile *_gf) 
    { 
        BaseBlock::setGadgetFile(_gf); 
        setNumberOfElements(_gf->h.npart[GAS]);
    }
    virtual void fill(GadgetGasParticle *p)      { if (N) read_next(&(p->u)); }
    virtual bool empty() { return N == 0; }
};

class RhoBlock : public Block<float,1>
{
public:
    RhoBlock() : Block<float,1>("RHO ", 0) { }
    virtual void setGadgetFile(GadgetFile *_gf) 
    { 
        BaseBlock::setGadgetFile(_gf); 
        setNumberOfElements(_gf->h.npart[GAS]);
    }
    virtual void fill(GadgetGasParticle *p)      { if (N) read_next(&(p->rho)); }
    virtual bool empty() { return N == 0; }
};

class HsmlBlock : public Block<float,1>
{
public:
    HsmlBlock() : Block<float,1>("HSML", 0) { }
    virtual void setGadgetFile(GadgetFile *_gf) 
    { 
        BaseBlock::setGadgetFile(_gf); 
        setNumberOfElements(_gf->h.npart[GAS]);
    }
    virtual void fill(GadgetGasParticle *p)      { if (N) read_next(&(p->hsml)); }
    virtual bool empty() { return N == 0; }
};

class PotBlock : public Block<float,1>
{
public:
    PotBlock() : Block<float,1>("POT ", 0) { }
    virtual void setGadgetFile(GadgetFile *_gf) 
    { 
        BaseBlock::setGadgetFile(_gf); 
        setNumberOfElements(_gf->nBodies);
    }
    virtual void fill(GadgetGasParticle *p)      { read_next(&(p->pot)); }
    virtual void fill(GadgetHaloParticle *p)     { read_next(&(p->pot)); }
    virtual void fill(GadgetDiskParticle *p)     { read_next(&(p->pot)); }
    virtual void fill(GadgetBulgeParticle *p)    { read_next(&(p->pot)); }
    virtual void fill(GadgetStarParticle *p)     { read_next(&(p->pot)); }
    virtual void fill(GadgetBoundaryParticle *p) { read_next(&(p->pot)); }
};

class AcceBlock : public Block<float,3>
{
public:
    AcceBlock() : Block<float,3>("ACCE", 0)  { }
    virtual void setGadgetFile(GadgetFile *_gf) 
    { 
        BaseBlock::setGadgetFile(_gf); 
        setNumberOfElements(_gf->nBodies);
    }
    virtual void fill(GadgetGasParticle *p)      { read_next(p->acc); }
    virtual void fill(GadgetHaloParticle *p)     { read_next(p->acc); }
    virtual void fill(GadgetDiskParticle *p)     { read_next(p->acc); }
    virtual void fill(GadgetBulgeParticle *p)    { read_next(p->acc); }
    virtual void fill(GadgetStarParticle *p)     { read_next(p->acc); }
    virtual void fill(GadgetBoundaryParticle *p) { read_next(p->acc); }
};

class EndtBlock : public Block<float,1>
{
public:
    EndtBlock() : Block<float,1>("ENDT", 0) { }
    virtual void setGadgetFile(GadgetFile *_gf) 
    { 
        BaseBlock::setGadgetFile(_gf); 
        setNumberOfElements(_gf->h.npart[GAS]);
    }
    virtual void fill(GadgetGasParticle *p)      { read_next(&(p->endt)); }
};

class TstpBlock : public Block<float,1>
{
public:
    TstpBlock() : Block<float,1>("TSTP", 0) { }
    virtual void setGadgetFile(GadgetFile *_gf) 
    { 
        BaseBlock::setGadgetFile(_gf); 
        setNumberOfElements(_gf->nBodies);
    }
    virtual void fill(GadgetGasParticle *p)      { read_next(&(p->dt)); }
    virtual void fill(GadgetHaloParticle *p)     { read_next(&(p->dt)); }
    virtual void fill(GadgetDiskParticle *p)     { read_next(&(p->dt)); }
    virtual void fill(GadgetBulgeParticle *p)    { read_next(&(p->dt)); }
    virtual void fill(GadgetStarParticle *p)     { read_next(&(p->dt)); }
    virtual void fill(GadgetBoundaryParticle *p) { read_next(&(p->dt)); }
};

class Scaler
{
public:
    virtual ~Scaler() {}
    virtual void pos(float *a, float *b) { b[0] = a[0]; b[1] = a[1]; b[2] = a[2]; }
    virtual void vel(float *a, float *b) { b[0] = a[0]; b[1] = a[1]; b[2] = a[2]; }
    virtual float mass(float a) { return a; }
    virtual float temp(float a) { return a; }
    virtual float rho(float a)  { return a; }
};

class CosmoScaler : public Scaler
{
    float RhoCritical, HubbleParam, BoxSize, Redshift, Time;
    float fact_mass, fact_vel, kel2kms2, rho_z;

public:
    CosmoScaler(float _RhoCritical, float _HubbleParam, 
                float _BoxSize, float _Redshift, float _Time)
    {
        float rho_msol, z_p1, volume;

        RhoCritical = _RhoCritical;
        HubbleParam = _HubbleParam;
        BoxSize     = _BoxSize;
        Redshift    = _Redshift;
        Time        = _Time;

        kel2kms2   = BOLTZMANN/(0.6*PROTONMASS)*1.e-10;
        rho_msol   = 2.7755*10.*pow(HubbleParam,2);       /*! Units of Msol/Mpc^3  */
        rho_z      = 1.e10*rho_msol/pow(HubbleParam,2)/1.e9; 

        z_p1       = 1.0F + Redshift;
        volume     = pow((BoxSize/1.0e3/HubbleParam),3);
        fact_mass  = rho_msol * volume;
        fact_vel   = (BoxSize/10/2.894405)/z_p1;

        if (BoxSize < 1e-15)
            cerr << "WARNING: BoxSize is extremely small!" << endl;
        if (fact_vel < 1e-15)
            cerr << "WARNING: fact_vel is extremely small! " << "What is the BoxSize and (z+1)?" << endl;
        if (RhoCritical < 1e-15)
            cerr << "WARNING: RhoCritical is extremely small!" << endl;
        if (HubbleParam < 1e-15)
            cerr << "WARNING: HubbleParam is extremely small!" << endl;
        if (fact_mass < 1e-15)
            cerr << "WARNING: fact_mass is extremely small! What are the BoxSize and HubbleParam?" << endl;
        if (rho_z < 1e-15)
            cerr << "WARNING: rho_z is extremely small! What is the HubbleParam?" << endl;

    }

    virtual ~CosmoScaler() {}

    virtual void pos(float *a, float *b) 
    {
        b[0] = (a[0] / BoxSize) - 0.5F;
        b[1] = (a[1] / BoxSize) - 0.5F;
        b[2] = (a[2] / BoxSize) - 0.5F;
    }

    virtual void vel(float *a, float *b) 
    {
        b[0] = a[0] * sqrt(Time) / fact_vel;
        b[1] = a[1] * sqrt(Time) / fact_vel;
        b[2] = a[2] * sqrt(Time) / fact_vel; 
    }
                                                                  
    virtual float mass(float m) { return m / RhoCritical / HubbleParam / fact_mass; }
    virtual float temp(float t) { return t / 1.5F / kel2kms2; }
    virtual float rho(float r)  { return r / rho_z * 1.0e10F; }
};

/*============================================================================
 * Function prototypes
 *==========================================================================*/
void create_tipsy_header(TipsyHeader &h, Config &cfg, GadgetHeader &iohdr);
int read_gadget_header(GadgetHeader &iohdr);
int read_config_file(char *file, Config &cfg);
double convertToDouble(const std::string& s, bool failIfLeftoverChars = true);
int convertToInt(const std::string& s, bool failIfLeftoverChars = true);
void help() __attribute__ ((noreturn));

/*============================================================================
 * Global configuration
 *==========================================================================*/
Config cfg;
Scaler *scaler;

void show_headers(GadgetHeader &iohdr, TipsyHeader &h)
{
    cout << "----------------------------------------------------------------" << endl;
    cout << "Gadget Header" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << "Particle Count: "
         << "Gas("   << iohdr.npart[GAS  ] << ") "
         << "Halo("  << iohdr.npart[HALO ] << ") "
         << "Disk("  << iohdr.npart[DISK ] << ") " << endl
         << "                " 
         << "Bulge(" << iohdr.npart[BULGE] << ") "
         << "Stars(" << iohdr.npart[STARS] << ") "
         << "Bndry(" << iohdr.npart[BNDRY] << ") " << endl;
    cout << "Gas("   << iohdr.npartTotal[GAS  ] << ") "
         << "Halo("  << iohdr.npartTotal[HALO ] << ") "
         << "Disk("  << iohdr.npartTotal[DISK ] << ") " << endl
         << "                " 
         << "Bulge(" << iohdr.npartTotal[BULGE] << ") "
         << "Stars(" << iohdr.npartTotal[STARS] << ") "
         << "Bndry(" << iohdr.npartTotal[BNDRY] << ") " << endl;
    cout << "Particle Mass: "
         << "Gas("   << iohdr.mass[GAS  ] << ") "
         << "Halo("  << iohdr.mass[HALO ] << ") "
         << "Disk("  << iohdr.mass[DISK ] << ") " << endl
         << "                " 
         << "Bulge(" << iohdr.mass[BULGE] << ") "
         << "Stars(" << iohdr.mass[STARS] << ") "
         << "Bndry(" << iohdr.mass[BNDRY] << ") " << endl;
    cout << "Time:        " << iohdr.time        << endl;
    cout << "Redshift:    " << iohdr.redshift    << endl;
    cout << "Omega0:      " << iohdr.Omega0      << endl;
    cout << "OmegaLambda: " << iohdr.OmegaLambda << endl;
    cout << "HubbleParam: " << iohdr.HubbleParam << endl;
    cout << "BoxSize:     " << iohdr.BoxSize     << endl;
    cout << "----------------------------------------------------------------" << endl << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << "Tipsy Header" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << "Particle Count: "
         << "Total(" << h.h_nBodies << ") "
         << "Gas("   << h.h_nGas << ") "
         << "Dark("  << h.h_nDark << ") "
         << "Star("  << h.h_nStar << ") " << endl;
    cout << "Time:        " << h.h_time        << endl;
    cout << "Redshift:    " << (1.0F / h.h_time - 1.0F)    << endl;
    cout << "Dimensions:  " << h.h_nDims << endl;
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
    cerr << "Usage: gad2tip [OPTIONS] <gadget-file> <tipsy-output-file>" << endl
         << endl
         //<< "    <tipsy-file>       A standard or native TIPSY format simulation file.
         << "OPTIONS:" << endl
         << "    -1                 GADGET file is in version 1 format. (default: auto-detect)" << endl
         << "    -2                 GADGET file is in version 2 format." << endl
         << "    --sort             Arrange particles in tipsy file to be consistent with IDs in" << endl
         << "                       the GADGET file. (slow!)" << endl
         << "    --scale-cosmo      Rescale values based on cosmological parameters in GADGET file. (default)" << endl
         << "    --no-scale         Do not rescale." << endl
         << "    -m <N>             Convert file in blocks of up to N megabytes (default " 
         <<                                                               DEFAULT_MEM << ")."<<endl
         << "                       If the amount of available memory is larger than the size"<<endl
         << "                       of the tipsy file this can significantly decrease the" << endl
         << "                       running time because the whole file can be converted" << endl
         << "                       in one read-convert-write operation." << endl
         << "    -v                 Increase verbosity." << endl
         << "    -q                 Decrease verbosity (make quiet)." << endl
         << "    --version          Print version information." << endl
         << "    --help             This help screen." << endl
         << endl
         << "The output file will have the same name as in the input file plus '.gad' if" << endl
         << "a single output file is requested. The --multi option produces files with" << endl
         << "the name of the input file plus an extension for each of the gadget blocks." << endl
         << endl
         << "WARNING: IDs are not handled properly. Particles are transferred to the tipsy" << endl
         << "         file in the order they appear in the gadget file. Converting multiple" << endl
         << "         gadget files from the same simulation will not preserve particle IDs." << endl
         << "         Use --sort to get around this." << endl
         << endl
         << "Send questions, comments, bug reports to Jonathan Coles <jonathan@physik.uzh.ch>" << endl
         << endl;
    
    exit(2);
}

void version()
{
    cerr << "gad2tip v" << VERSION << endl
         << "Gadget to Tipsy file converter." << endl
         << "Written by Jonathan Coles <jonathan@physik.uzh.ch>" << endl;

    exit(0);
}

void create_tipsy_header(TipsyHeader &h, Config &_cfg, GadgetHeader &iohdr)
{
    h.h_time    = _cfg.CosmoScaling ? 1. / (iohdr.redshift + 1.) : iohdr.time;
    h.h_nGas    = iohdr.npart[GAS];
    h.h_nDark   = iohdr.npart[HALO];
    h.h_nStar   = iohdr.npart[STARS] + iohdr.npart[DISK] + iohdr.npart[BULGE];
    h.h_nBodies = h.h_nGas + h.h_nDark + h.h_nStar;
    h.h_nDims   = 3;
}

void scale_to_tipsy(GadgetGasParticle &Gg, TipsyGasParticle &Tg)
{
    scaler->pos(Gg.pos, Tg.pos);
    scaler->vel(Gg.vel, Tg.vel);
    Tg.mass    = scaler->mass(Gg.mass);
    Tg.phi     = Gg.pot;
    Tg.temp    = Gg.u;
    Tg.rho     = Gg.rho;
    Tg.hsmooth = Gg.hsml;
    Tg.metals  = 0;
}

void scale_to_tipsy(GadgetHaloParticle &Gh, TipsyDarkParticle &Td)
{
    scaler->pos(Gh.pos, Td.pos);
    scaler->vel(Gh.vel, Td.vel);
    Td.mass = scaler->mass(Gh.mass);
    Td.phi  = Gh.pot;
    Td.eps  = 0;
}

void scale_to_tipsy(GadgetStarParticle &Gs, TipsyStarParticle &Ts)
{
    scaler->pos(Gs.pos, Ts.pos);
    scaler->vel(Gs.vel, Ts.vel);
    Ts.mass    = scaler->mass(Gs.mass);
    Ts.phi     = Gs.pot;
    Ts.eps     = 0;
    Ts.tform   = 0;
    Ts.metals  = 0;
}

void scale_to_tipsy(GadgetDiskParticle &Gd, TipsyStarParticle &Ts)
{
    scaler->pos(Gd.pos, Ts.pos);
    scaler->vel(Gd.vel, Ts.vel);
    Ts.mass    = scaler->mass(Gd.mass);
    Ts.phi     = Gd.pot;
    Ts.eps     = 0;
    Ts.tform   = 0;
    Ts.metals  = 0;
}

void scale_to_tipsy(GadgetBulgeParticle &Gb, TipsyStarParticle &Ts)
{
    scaler->pos(Gb.pos, Ts.pos);
    scaler->vel(Gb.vel, Ts.vel);
    Ts.mass    = scaler->mass(Gb.mass);
    Ts.phi     = Gb.pot;
    Ts.eps     = 0;
    Ts.tform   = 0;
    Ts.metals  = 0;
}

void scale_to_tipsy(GadgetBoundaryParticle &Gy, TipsyStarParticle &Ts)
{
    scaler->pos(Gy.pos, Ts.pos);
    scaler->vel(Gy.vel, Ts.vel);
    Ts.mass    = scaler->mass(Gy.mass);
    Ts.phi     = Gy.pot;
    Ts.eps     = 0;
    Ts.tform   = 0;
    Ts.metals  = 0;
}

int main(int argc, char **argv)
{
    ofTipsy out;
    TipsyHeader h;
    char block_list[256] = {"PVIMURH"};

    static struct option long_options[] = {
        {"scale-cosmo",  no_argument, 0, 0},
        {"no-scale",  no_argument, 0, 0},
        {"sort",  no_argument, 0, 0},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 0},
        {"blocks", required_argument, 0, 'b'},
        {0, 0, 0, 0}
    };

    size_t use_mem = DEFAULT_MEM * 1024 * 1024;

    if (argc < 2) help();

    cfg.SnapFormat = GADGET_FORMAT_DETECT;

    /*========================================================================
     * Process the command line flags
     *======================================================================*/
    while (1)
    {
        int option_index = 0;
        int c = getopt_long(argc, argv, "hvqm:12b:",
                            long_options, &option_index);

        if (c == -1) break;

        switch (c)
        {
            case 0:
                if (!strcmp("scale-cosmo", long_options[option_index].name))
                    cfg.CosmoScaling = true;
                else if (!strcmp("no-scale", long_options[option_index].name))
                    cfg.CosmoScaling = false;
                else if (!strcmp("sort", long_options[option_index].name))
                    cfg.sort = true;
                else if (!strcmp("version", long_options[option_index].name))
                    version();
                break;

            case 'b': 
                strncpy(block_list, optarg, sizeof(block_list));
                block_list[sizeof(block_list)-1] = 0; 
                break;
            case '1': cfg.SnapFormat = GADGET_FORMAT_1; break;
            case '2': cfg.SnapFormat = GADGET_FORMAT_2; break;
            case 'h': help(); break;
            case 'v': verbosity++; break;
            case 'q': verbosity = 0; break;
            case 'm': use_mem = parse_memory_arg(optarg); break;
            case '?': help(); break;
        }
    }

    if (optind >= argc) help();

    VL(2) cout << "Using " << use_mem << " bytes to load particles." << endl;

    GadgetFile in(argv[optind], cfg.SnapFormat);
    if (!in.is_open())
    {
        cerr << "ERROR: Unable to open GADGET file " << argv[optind] << endl;
        exit(1);
    }
    optind++;

    out.open(argv[optind], "standard");
    if (!out.is_open()) 
    {
        out.open(argv[optind], "native");
        if (!out.is_open()) 
        {
            cerr << "ERROR: Unable to open Tipsy binary " << argv[optind] << endl;
            exit(1);
        }
    }
    optind++;

    for (int i=0; i < strlen(block_list); i++)
    {
        switch (block_list[i])
        {
            case 'P': if (!in.addBlock(new PosBlock())) exit(1);  break;
            case 'V': if (!in.addBlock(new VelBlock())) exit(1);  break;
            case 'I': if (!in.addBlock(new IDBlock())) exit(1);   break;
            case 'M': if (!in.addBlock(new MassBlock())) exit(1); break;
            case 'U': if (!in.addBlock(new UBlock())) exit(1);    break;
            case 'R': if (!in.addBlock(new RhoBlock())) exit(1);  break;
            case 'H': if (!in.addBlock(new HsmlBlock())) exit(1); break;
            case 'p': if (!in.addBlock(new PotBlock())) exit(1);  break;
            case 'A': if (!in.addBlock(new AcceBlock())) exit(1); break;
            case 'E': if (!in.addBlock(new EndtBlock())) exit(1); break;
            case 'T': if (!in.addBlock(new TstpBlock())) exit(1); break;
            default:
                cerr << "Block list code " << block_list[i] << " not recognized. Ignored." << endl;
                break;
        }
    }

    in.setBufferSize(use_mem);

    create_tipsy_header(h, cfg, in.h);

    if (cfg.CosmoScaling)
        scaler = new CosmoScaler(cfg.RhoCritical, in.h.HubbleParam, in.h.BoxSize, in.h.redshift, in.h.time);
    else
        scaler = new Scaler();


    VL(1) show_headers(in.h, h);

    out << h;

    /*========================================================================
     * Write all the particle information. 
     *
     * Gadget files expect gas particles first, then dark particles, then
     * stars.
     *======================================================================*/

    VL(1) cerr << "   ""                                                   ] 100%\r";
    VL(1) cerr << "0% [";

    unsigned int two_percent = (int)(h.h_nBodies * 0.02);
#define STATUS VL(1) if (id>0 && (id % two_percent) == 0) cerr << ".";

    tipsypos::offset_type i;

    GadgetGasParticle Gg;
    GadgetStarParticle Gs;
    GadgetHaloParticle Gh;
    GadgetDiskParticle Gd;
    GadgetBulgeParticle Gb;
    GadgetBoundaryParticle Gy;

    TipsyGasParticle  Tg;
    TipsyStarParticle Ts;
    TipsyDarkParticle Td;

    uint32_t id=0;
    for (uint32_t j=0; j < in.h.npart[GAS]; j++, id++)
    {
        in >> Gg; scale_to_tipsy(Gg, Tg); out << Tg;
        STATUS;
    }

    for (uint32_t j=0; j < in.h.npart[HALO]; j++, id++)
    {
        in >> Gh; scale_to_tipsy(Gh, Td); out << Td;
        STATUS;
    }

    /* Disk, Bulge, and Boundary particles are treated like stars */

    for (uint32_t j=0; j < in.h.npart[DISK]; j++, id++)
    {
        in >> Gd; scale_to_tipsy(Gd, Ts); out << Ts;
        STATUS;
    }

    for (uint32_t j=0; j < in.h.npart[BULGE]; j++, id++)
    {
        in >> Gb; scale_to_tipsy(Gb, Ts); out << Ts;
        STATUS;
    }

    for (uint32_t j=0; j < in.h.npart[STARS]; j++, id++)
    {
        in >> Gs; scale_to_tipsy(Gs, Ts); out << Ts;
        STATUS;
    }

    for (uint32_t j=0; j < in.h.npart[BNDRY]; j++, id++)
    {
        in >> Gy; scale_to_tipsy(Gy, Ts); out << Ts;
        STATUS;
    }

    VL(1) cerr << endl;

    VL(1) cerr << "Closing files" << endl;
    in.close();
    out.close();

    return 0;
}
