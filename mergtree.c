/* 25-9-2002 GM */
/* Multimass, multi-type version. Considers HiRes parts as Type 1*/
/* finds bigger DM HiRes group, writes its cm on stdout and writes
   its particels in separate file reminds fof idxes (idxfof.z###.dat) */
/* The part commented with --- is to select massive grps for output files */
/* this version also asks for an external index file; then finds
   the progenitor of the group given in the ext file at the current redshift.
   A group is a progenitor if al least PROGPERC of its mass is made of
   particles belonging to the input group. */
/* REM the input index file must contain GADGET label of groups, ONE
   group will be loaded*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define NTYP        10      /* Max number of different particle types */
#define PROGPERC    0.5     /* percentage of mass needed to be a progenitor */

#define UdM     1.e10       /* Unita' di misura "gadget naturali */
#define UdL     1.0         /* Msol, Kpc, km/s */
#define UdV     1.0
#define NMIN    100         /* default di gius e' 30 */
#define UdMA    5.97e10     /* Unita' di misura Anna  */
#define UdLA    20.0        /* Msol, Kpc, km/s */
#define UdVA    113.3

#define DO_NOTHING  0
#define MAKE_FOF    1
#define MAKE_GROUPS 2
#define MAKE_TREE   4

#define GAS     0
#define DM      1

struct io_header_1
{
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam; 
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;

struct particle_data 
{
    float  Pos[3];
    float  Vel[3];
    float  Mass;
    int    Type;

    float  Rho, U, Temp, Ne;
} *P;

struct group_data 
{
  float Pos[3];
  float Vel[3];
  float Mass;
  float HiResMass;
  int npart;
  float rad;
  int isProgenitor;
} *G;

typedef struct
{
    int gid;
    int npart;
    int *ps;
} gmap_t;

typedef struct
{
    int id;
    float Pos[3];
    float Vel[3];
    float Mass;
    float rad;
    int npart;
} group_t;

void EvalGroups0(int *idx,char *fname);
int EvalGroups1(int *idx);
void EvalProgenitors(char *fname, int *idx, int *sonidx, int npson, int Ngroups, gmap_t *gmap, int gmapLen);
int read_pmap(char *fname, 
              float fEpsf, 
              float fPeriodf1, float fPeriodf2, float fPeriodf3, 
              int nf, int *pmap);
int write_pmap(char *fname, 
               float fEpsf, 
               float fPeriodf1, float fPeriodf2, float fPeriodf3, 
               int nf, int *pmap);
int read_gmap(char *fname, gmap_t **gmap0, int *nGroups);
int read_groups(char *fname, group_t **groups0, int *nGroups0);
void fof_(float ,float ,float ,float , int , struct particle_data *, int* );
//int load_son_idx(int*);
int find_progenitors(int*,int*,int,int);
int load_snapshot(char *fname, int files);

void CenterGroups(int Ngroups,int* idx);
int allocate_memory(void);

int nMin = NMIN;
float MassMin = -1,MassMax = -1;
int     NumPart, Ngas;
int *Id;
double  Time, Redshift;

void help()
{
    fprintf(stderr, 
"Usage: mergtree --make-fof <snapshot>\n"
"Usage: mergtree [OPTIONS] -M MassMin MassMax <snapshot>\n"
"Usage: mergtree [OPTIONS] <snapshot> <index_file>\n"
"\n"
"    snapshot                   The GADGET snapshot file to read.\n"
"    index_file                 An optional index file. If this is specified\n"
"                               then the progenitors will be found.\n"
"    -M MassMin MassMax         Set the mass range.\n"
"\n"
"OPTIONS:\n"
"    -fof fof_file              Read the FOF data from fof_file instead of\n"
"                               calculating it.\n"
"\n"
    );

    exit(2);
}

/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int main(int argc, char **argv)
{
    char dir[256];
    char *snapshot_fname = NULL, 
         *fof_fname = NULL, 
         *groups_fname = NULL,
         *gmap_fname = NULL;
    int  files;
    float fPeriodf1,fPeriodf2,fPeriodf3,fEpsf;
    int i,nf;
    int *pmap; /* A map between a particle and the group it is in */
    int nGroups, nMapGroups, Ngroups;
    group_t *groups;
    gmap_t *gmap;

    int mode = DO_NOTHING;

    /*=======================================================================
     * Parse command line arguments
     *=====================================================================*/
    for (i=1; i < argc; i++)
    {
        if (!strcmp("--make-fof", argv[i]))
        {
            if (mode != DO_NOTHING) help();
            if (++i < argc) snapshot_fname = argv[i]; else help();
            if (++i < argc) fof_fname = argv[i]; else help();
            mode = MAKE_FOF;
        }
        else if (!strcmp("--make-groups", argv[i]))
        {
            if (mode != DO_NOTHING) help();
            if (++i < argc) snapshot_fname = argv[i]; else help();
            if (++i < argc) fof_fname      = argv[i]; else help();
            mode = MAKE_GROUPS;
        }
        else if (!strcmp("--make-tree", argv[i]))
        {
            if (mode != DO_NOTHING) help();
            if (++i < argc) groups_fname   = argv[i]; else help();
            if (++i < argc) snapshot_fname = argv[i]; else help();
            if (++i < argc) fof_fname      = argv[i]; else help();
            if (++i < argc) gmap_fname     = argv[i]; else help();
            mode = MAKE_TREE;
        }
        else if (!strcmp("-M", argv[i]))
        {
            if (++i < argc) MassMin = atof(argv[i]); else help();
            if (++i < argc) MassMax = atof(argv[i]); else help();
        }
        else if (!strcmp("-nmin", argv[i]))
        {
            if (++i < argc) nMin = atof(argv[i]); else help();
        }
    }

    /*=======================================================================
     * Check that the user gave us sensible parameters
     *=====================================================================*/
    if (mode == DO_NOTHING) help();

    /*=======================================================================
     * Load in the snapshot (GADGET format)
     *=====================================================================*/
    files=1;                               /* number of files per snapshot */
    fprintf(stderr, "Loading snapshot %s.\n", snapshot_fname);
    load_snapshot(snapshot_fname, files);

    printf("Npart: %d; Time: %lf; Redshift: %lf\n",NumPart,header1.time, header1.redshift);

    nf=NumPart;
    printf("First particle:  %f %f %f \n",P[1].Pos[0],P[1].Pos[1],P[1].Pos[2]);
    printf("Medium particle: %f %f %f \n",P[nf/2].Pos[0],P[nf/2].Pos[1],P[nf/2].Pos[2]);
    printf("Last particle:   %f %f %f \n",P[nf].Pos[0],P[nf].Pos[1],P[nf].Pos[2]);
  

    /*=======================================================================
     * Prepare to either create or load the FOF data
     *=====================================================================*/
    fEpsf =  header1.BoxSize/ pow( NumPart, 1./3.);
    printf("Distanza interparticellare media: %f\n",fEpsf);
    fEpsf *= 0.2; 
/*  fEpsf *= 0.15; */
/*  fEpsf *= 0.1; */
    fPeriodf1=fPeriodf2=fPeriodf3=(float) header1.BoxSize;
    pmap=malloc((NumPart+1)*sizeof(int));
    if(pmap==NULL) 
    { 
        printf("Errore di allocazione memoria indici.\n");
        exit(-1); 
    }

    /* LE UNITA' NATURALI GADGET LAVORANO IN 10^10 MSOL. SE SONO STATE 
    CAMBIATE, OCCORRE MODIFICARE QUI */

    if (MassMin != -1 && MassMax != -1) 
    { 
        fprintf(stderr, "Mass range: [%2.2e,%2.2e]\n", MassMin, MassMax);
        MassMin /= UdM; MassMax /= UdM; 
    }
    else
    {
        fprintf(stderr, "Mass range: Unbounded\n");
        MassMin = 0;
        MassMax = 1e33;
    }

    switch (mode)
    {
        case MAKE_FOF:
            fprintf(stderr, "Making FOFs.\n");
            fof_(fEpsf,fPeriodf1,fPeriodf2,fPeriodf3,NumPart,P,pmap);
            write_pmap(fof_fname, fEpsf,fPeriodf1,fPeriodf2,fPeriodf3,NumPart,pmap);
            break;

        case MAKE_GROUPS:
            fprintf(stderr, "Making Groups.\n");
            if (read_pmap(fof_fname, fEpsf,fPeriodf1,fPeriodf2,fPeriodf3,NumPart,pmap))
                exit(1);
            EvalGroups0(pmap,snapshot_fname);
            break;

        case MAKE_TREE:
            fprintf(stderr, "Making Progenitor Trees.\n");
            if (read_groups(groups_fname, &groups, &nGroups)) exit(1);
            if (read_gmap(gmap_fname, &gmap, &nMapGroups))    exit(1);

            /*=======================================================================
             * Sanity check that the group files match up.
             *=====================================================================*/
#if 0
            if (nGroups != nMapGroups)
            {
                fprintf(stderr, 
                    "Mismatched number of groups between group and map files (%i != %i).\n",
                    nGroups, nMapGroups);
                exit(1);
            }
#endif

#if 0
            int nPartInGroups=0, nPartInMap=0;
            for (i=0; i < nGroups; i++)
            {
                nPartInGroups += groups[i].npart;
                nPartInMap    += gmap[i].npart;
            }

            fprintf(stderr, "HERE 4\n");

            if (nGroups != nMapGroups)
            {
                fprintf(stderr, "Mismatched number of particles between group and map files.\n");
                exit(1);
            }
#endif

            /*=======================================================================
             * Now load the particle->group map.
             *=====================================================================*/
            if (read_pmap(fof_fname, fEpsf,fPeriodf1,fPeriodf2,fPeriodf3,NumPart,pmap))
                exit(1);


            //sonidx=malloc( (NumPart+1)*sizeof(int));
            //npson=load_son_idx(sonidx, index_fname);
            //if(npson <= 0) { printf("ERROR: group non find.\n"); exit(-1); }
            //printf("<<Son>> group with %d particles loaded.\n",npson);


            mkdir("GROUPS", -1);
            assert(!chdir("GROUPS"));
            Ngroups = EvalGroups1(pmap);
            chdir("..");

            gmap_t *gmap2 = (gmap_t *)calloc(Ngroups+1, sizeof(gmap_t));
            for (i=1; i < NumPart+1; i++)
                gmap2[pmap[i]].npart++;

            for (i=0; i < Ngroups+1; i++)
            {
                gmap2[i].ps = (int *)malloc(gmap2[i].npart * sizeof(int));
                assert(gmap2[i].ps != NULL);
                gmap2[i].npart = 0;
            }

            for (i=1; i < NumPart+1; i++)
            {
                gmap2[pmap[i]].ps[ gmap2[pmap[i]].npart++ ] = i;
            }
#if 1

            printf("Finding progenitors with mass %%: %f\n",PROGPERC);

            for (i=0; i < nGroups; i++)
            {

                //fprintf(stderr, "Group %i [Mass: %f]\n", i, groups[i].Mass);
#if 1
                if (MassMin != -1 && MassMax != -1)
                {
                    if (!(MassMin < groups[i].Mass && groups[i].Mass < MassMax))
                    {
                        fprintf(stderr, "Outside mass range.  Group %i/%i.\n", i+1, nGroups+1);;
                        continue;
                    }
                }

                snprintf(dir, 256, "%05i", groups[i].id);
                fprintf(stderr, "Making directory %s. Group %i/%i.\n", dir, i+1, nGroups+1);;
                mkdir(dir, -1);
                assert(!chdir(dir));

                int *ps = NULL;
                int npart = 0;
                int j;
                for (j=0; j < nMapGroups; j++)
                {
                    if (groups[i].id == gmap[j].gid)
                    {
                        assert(ps == NULL);
                        ps    = gmap[j].ps;
                        npart = gmap[j].npart;

                        assert(groups[i].npart == gmap[j].npart);

                        /* Intentionally not breaking here to check the assert()'ion if there
                         * is somehow a duplicate groups. */
                    }
                }

                assert(ps != NULL);

                EvalProgenitors(snapshot_fname, pmap, ps, npart, Ngroups, gmap2, Ngroups);
                
                chdir("..");
#endif
            }
#endif
            break;

        default:
            fprintf(stderr, "Unknown internal mode: %i\n", mode);
            exit(1);
    }

    return 0;
}

int read_groups(char *fname, group_t **groups0, int *nGroups0)
{
    FILE *in;
    int ret=0;
    int read;
    int nGroups = 0;
    int allocdGroups = 0;
    group_t *groups = NULL;

    in = fopen(fname, "r");
    if (in == NULL) ret = 1;

    while (!ret && !feof(in))
    {
        if (allocdGroups == nGroups)
        {
            if  (allocdGroups == 0) allocdGroups = 512; 
            else allocdGroups *= 2;

            groups = (group_t *)realloc(groups, allocdGroups * sizeof(group_t));
            if (groups == NULL) { ret = 2; break; }
        }

        group_t *g = &(groups[nGroups]);
        read = 
            fscanf(in,"%d %f %f %f %f %f %f %f %e %d\n",
                &(g->id),
                &(g->Pos[0]), &(g->Pos[1]), &(g->Pos[2]),
                &(g->Vel[0]), &(g->Vel[1]), &(g->Vel[2]),
                &(g->rad), &(g->Mass), &(g->npart));

        if (read == 0) continue; /* check for EOF at top of loop */

        if (ret != 10) { ret = 3; break; }

        g->Mass /= UdM;

        nGroups++;
    }

    fclose(in);

    switch (ret)
    {
        case 0:
            *groups0 = (group_t *)realloc(groups, nGroups * sizeof(group_t)); assert(*groups0 != NULL);
            *nGroups0 = nGroups;
            break;
        case 1:
            fprintf(stderr, "ERROR: Can't open %s.\n", fname);
            break;
        case 2:
            fprintf(stderr, "ERROR: Can't allocate memory for groups in %s.\n", fname);
            break;
        case 3:
            fprintf(stderr, "ERROR: Unexpected EOF in %s.\n", fname);
            break;
    }

    return ret;
}

int read_gmap(char *fname, gmap_t **gmap0, int *nGroups0)
{
    FILE *in;
    int nGroups = 0;
    int allocdGroups = 0;
    gmap_t *gmap = NULL;

    int nParticles;
    int gid;
    int *ps;

    int ret = 0;

    in = fopen(fname, "r");
    if (in == NULL) ret = 1;

    while (!ret && !feof(in))
    {
        if (allocdGroups == nGroups)
        {
            if  (allocdGroups == 0) allocdGroups = 512; 
            else allocdGroups *= 2;

            gmap = (gmap_t *)realloc(gmap, allocdGroups * sizeof(gmap_t)); 
            if (gmap == NULL) { ret = 2; break; }
        }

        if (fread(&nParticles, sizeof(int), 1, in) <= 0) { continue; }
        if (fread(&gid,        sizeof(int), 1, in) <= 0) { return 1; }

        //fprintf(stderr, "Trying to read %i particles. Group %i\n", nParticles, gid);

        ps = (int *)malloc(nParticles * sizeof(int)); assert(ps != NULL);
        ret=fread(ps, sizeof(int), nParticles, in);
        if (ret < nParticles) { ret = 3; break; }

        gmap[nGroups].gid   = gid;
        gmap[nGroups].npart = nParticles;
        gmap[nGroups].ps    = ps;

        nGroups++;
    }

    fclose(in);

    switch (ret)
    {
        case 0:
            *gmap0 = (gmap_t *)realloc(gmap, nGroups * sizeof(gmap_t)); assert(*gmap0 != NULL);
            *nGroups0 = nGroups;
            break;
        case 1:
            fprintf(stderr, "ERROR: Can't open %s.\n", fname);
            break;
        case 2:
            fprintf(stderr, "ERROR: Can't allocate memory for groups in %s.\n", fname);
            break;
        case 3:
            fprintf(stderr, "ERROR: Unexpected EOF in %s.\n", fname);
            break;
    }

    return ret;
}

int read_pmap(char *fname, 
              float fEpsf0, 
              float fPeriodf10, float fPeriodf20, float fPeriodf30, 
              int NumPart0, int *pmap)
{
    FILE *in;
    float fEpsf; 
    float fPeriodf1, fPeriodf2, fPeriodf3; 
    int NumPart;
    int ret = 0;

    in = fopen(fname, "r");
    if (in == NULL) ret = 1;

#define READ(var, type, count) if (fread(&(var), sizeof(type), count, in) < count) { ret = 2; break }
#define CHECK(var0, var1) if (var0 != var1) { ret = 3; break; }
    
    while (!ret)
    {
        READ(fEpsf, fEpsf, 1); CHECK(fEpsf, fEpsf0);
        READ(fPeriodf1, fPeriodf1, 1); CHECK(fPeriodf1, fPeriodf10);
        READ(fPeriodf2, fPeriodf2, 1); CHECK(fPeriodf2, fPeriodf20);
        READ(fPeriodf3, fPeriodf3, 1); CHECK(fPeriodf3, fPeriodf30);
        READ(NumPart,   NumPart,   1); CHECK(NumPart,   NumPart0);
        READ(pmap, int, NumPart);
    }
        
    fclose(in);

    switch (ret)
    {
        case 0: break; /* Normal case */
        case 1:
            fprintf(stderr, "ERROR: Cannot open %s\n", fname);
            break;
        case 2:
            fprintf(stderr, "ERROR: Unexpected EOF in %s\n", fname);
            break;
        case 3:
            fprintf(stderr, "ERROR: Parameters used to make %s are different from the snapshot.\n", fname);
            break;
    }

    return ret;
}

int write_pmap(char *fname, 
               float fEpsf, 
               float fPeriodf1, float fPeriodf2, float fPeriodf3, 
               int NumPart, int *pmap)
{
    FILE *out;

    out = fopen(fname, "w");
    if (out == NULL) return 1;
    
    fwrite(&fEpsf,     sizeof(fEpsf),     1, out);
    fwrite(&fPeriodf1, sizeof(fPeriodf1), 1, out);
    fwrite(&fPeriodf2, sizeof(fPeriodf2), 1, out);
    fwrite(&fPeriodf3, sizeof(fPeriodf3), 1, out);
    fwrite(&NumPart,   sizeof(NumPart),   1, out);
    fwrite(pmap,       sizeof(int), NumPart, out);

    fclose(out);

    return 0;
}

void EvalGroups0(int *idx,char *fname)
{
    float TotMass;
    float *hist,PMassMin,mbig;
    float x,y,z,vx,vy,vz,m;
    int i,j,Ngroups,NunderMin,Ncho,ibig;
    FILE *f,*f2;
    char flnm[256],fltmp[256],bfr[256];
    float *hhist;

    FILE *idx_out, *dm_out, *gas_out;

    printf("\n\nTrovo gruppi fof con estremi di massa dati.\n");
/* --- 
printf("Introdurre gli estremi di massa (in Msol):\n");
scanf("%f %f",&MassMin,&MassMax);
printf("Estremi in massa: %e %e\n",MassMin,MassMax);
   --- */

    hist=calloc((NumPart+1), sizeof(float));
    hhist=calloc((NumPart+1), sizeof(float));
    if(hist==NULL) { printf("Impossibile allocare memoria per l'istogramma delle masse.\n"); exit(-1); }
    if(hhist==NULL) { printf("Impossibile allocare memoria per l'istogramma delle masse.\n"); exit(-1); }

    PMassMin=1.e33;
    for(i=1;i<NumPart+1;i++) 
        if(PMassMin > P[i].Mass && P[i].Type != GAS) PMassMin=P[i].Mass;
    printf("Massa min. particelle: %e Msol\n",PMassMin*UdM);
/* ATTENIONE! i numeri dei gruppi partono da 1 non da 0, idem le
   particelle; idx parte da zero quindi la particella i ha indice idx[i-1] */
/* idx[] e' il vettore che associa ogni particella i al gruppo idx[i-1];
   se idx[i]=0 la particella non appartiene ai gruppi */
/* siccome e' indicizzato da idx[], anche hist[] parte da 1 */
    TotMass=0.0;
    Ngroups=0;
    for(i=1;i<NumPart+1;i++) 
    {
// if(idx[i-1]>0&&P[i].Mass<=PMassMin+0.001) {
        if(idx[i-1]>0) 
        {
            hist[idx[i-1]] += P[i].Mass;
            if(P[i].Type==DM) hhist[idx[i-1]] += P[i].Mass;
            TotMass += P[i].Mass;
            if(idx[i-1]>Ngroups) Ngroups=idx[i-1];
        }
    }

    printf("Massa totale in gruppi fof: %e Msol.\n",TotMass*UdM);
    printf("Numero gruppi fof: %d\n",Ngroups);
    NunderMin=0;
    for(i=1;i<Ngroups+1;i++) 
    {
        if(hist[i]/PMassMin<nMin) 
        {
            hist[i]=0;
            NunderMin++;
        }
    }
    printf("Numero gruppi sotto il minimo numero di particelle (%d, massa %e): %d\n", nMin,nMin*PMassMin*UdM,NunderMin);
    printf("Gruppi sopra soglia: %d\n",Ngroups-NunderMin);

    G=malloc( (Ngroups+1)*sizeof(struct group_data));
    if(G==NULL) { printf("Impossibile allocare la memoria per i dati dei gruppi.\n"); exit(-1); }

    for(i=1;i<Ngroups+1;i++) { G[i].Mass=hist[i]; G[i].HiResMass=hhist[i]; }
 
    CenterGroups(Ngroups,idx);

    sprintf(bfr,"Groups.z%4.2f.dat",header1.redshift);
    f=fopen(bfr,"w");
    for(i=1;i<Ngroups+1;i++) 
    {
        if(G[i].Mass>0.0) 
        {
            fprintf(f,"%d %f %f %f %f %f %f %f %e %d\n",
                i,
                G[i].Pos[0],G[i].Pos[1],G[i].Pos[2],
                G[i].Vel[0],G[i].Vel[1],G[i].Vel[2],
                G[i].rad,G[i].Mass*UdM,G[i].npart);
        }
    }
    fclose(f);

    mbig=0.0; ibig=0;
    for(i=1; i<Ngroups+1; i++)
        if (G[i].HiResMass>mbig) { mbig=G[i].HiResMass; ibig=i;}

    printf("Bigger HiRes group (%d) mass %e\n",ibig,mbig);
    printf("Coordinates %f %f %f\n",G[ibig].Pos[0],G[ibig].Pos[1],G[ibig].Pos[2]);

    f2=fopen("d.inp","w"); /* questo file serve a DensitySphere */
    fprintf(f2,"%s\n",fname);
    fprintf(f2,"%f %f %f\n",G[ibig].Pos[0],G[ibig].Pos[1],G[ibig].Pos[2]);
    fprintf(f2,"120.0\n");   /* radius for dicot. search */
    fprintf(f2,"0.3 0.7\n"); /* omega0,omegaL */
    fprintf(f2,"178.\n");    /* DELTA (at z=0 for omega0=1; usually 178 */
    fprintf(f2,"0.05\n");    /* error on DELTA  in (0.0,1.0)*/
    fclose(f2);

    Ncho=0;

    for(i=1;i<Ngroups+1;i++) 
        if( (MassMin < hist[i] && hist[i] < MassMax) || i==ibig) Ncho++;
    printf("Numero gruppi scelti: %d\n",Ncho);


    printf("Output file: grp####.dat\n");
    sprintf(bfr,"idxfof.z%4.2f.dat",header1.redshift);
    printf("Index file: %s\n",bfr);
    idx_out=fopen(bfr,"wb");  /* questo contiene gli indici per eventuale follow-up */
    strcpy(fltmp,"grp");
    for(i=1;i<Ngroups+1;i++) 
    {
        if(MassMin < hist[i] && hist[i] < MassMax) 
        {
            fwrite( &(G[i].npart), sizeof(G[i].npart), 1, idx_out); 
            fwrite( &i, sizeof(i), 1, idx_out); 
            snprintf(flnm, 256, "z%4.2f.%s%d.DM.dat",header1.redshift,fltmp,i); 
            dm_out = fopen(flnm,"w");
            snprintf(flnm, 256, "z%4.2f.%s%d.GAS.dat",header1.redshift,fltmp,i); 
            gas_out = fopen(flnm,"w");
            for(j=1;j<NumPart+1;j++) {
                if(idx[j-1]==i) {
                    /* converto in fisico: kpc/h, km/s, Msol */
#if 1
                    x=P[j].Pos[0] -G[i].Pos[0]; 
                    y=P[j].Pos[1] -G[i].Pos[1]; 
                    z=P[j].Pos[2] -G[i].Pos[2];
                    if(x>header1.BoxSize/2.0) x -= header1.BoxSize;
                    if(y>header1.BoxSize/2.0) y -= header1.BoxSize;
                    if(z>header1.BoxSize/2.0) z -= header1.BoxSize;
                    if(x< -header1.BoxSize/2.0) x += header1.BoxSize;
                    if(y< -header1.BoxSize/2.0) y += header1.BoxSize;
                    if(z< -header1.BoxSize/2.0) z += header1.BoxSize;
                    x *= header1.time;
                    y *= header1.time;
                    z *= header1.time;
                    vx=P[j].Vel[0] - G[i].Vel[0]; 
                    vy=P[j].Vel[1] - G[i].Vel[1]; 
                    vz=P[j].Vel[2] - G[i].Vel[2];
                    vx *= sqrt(header1.time);
                    vy *= sqrt(header1.time);
                    vz *= sqrt(header1.time);
                    m=P[j].Mass*UdM;

                    fprintf(P[j].Type == GAS ? gas_out : dm_out,
                        "%f %f %f %f %f %f %f %d %f %f %f\n",
                        x,y,z,vx,vy,vz,m,P[j].Type,P[j].Pos[0],P[j].Pos[1],P[j].Pos[2]);
#endif
                    fwrite( &(Id[j]), sizeof(Id[j]), 1, idx_out);  
                }
            }
            fprintf(dm_out,"# Group N.: %d Mass: %e N. Part: %d\n",i,hist[i]*UdM,G[i].npart);
            fprintf(dm_out,"# CM coords (COMOV): %f %f %f\n",G[i].Pos[0],G[i].Pos[1],G[i].Pos[2]);
            fprintf(dm_out,"# Here particle pos,vels are PHYSICAL and the group is CM centered\n");
            fprintf(dm_out,"# After part types you find COMOVING coordinates.\n"); 
            fclose(dm_out);
            fclose(gas_out);
        }
    }
    fclose(idx_out); 
}


int EvalGroups1(int *idx)
{
    float TotMass;
    float *hist,PMassMin;
    int i,j,Ngroups,NunderMin;
    char bfr[256];
    float *hhist;
    FILE *grps_out;

/* ---
    printf("\n\nTrovo gruppi fof con estremi di massa dati.\n");
    printf("Introdurre gli estremi di massa (in Msol):\n");
    scanf("%f %f",&MassMin,&MassMax);
    printf("Estremi in massa: %e %e\n",MassMin,MassMax);
    MassMin /= UdM; MassMax /= UdM;
   --- */
/* LE UNITA' NATURALI GADGET LAVORANO IN 10^10 MSOL. SE SONO STATE 
CAMBIATE, OCCORRE MODIFICARE QUI */

    hist=calloc((NumPart+1), sizeof(float));
    hhist=calloc((NumPart+1), sizeof(float));
    if(hist==NULL) { printf("Impossibile allocare memoria per l'istogramma delle masse.\n"); exit(-1); }
    if(hhist==NULL) { printf("Impossibile allocare memoria per l'istogramma delle masse.\n"); exit(-1); }

    PMassMin=1.e33;
    for(i=1;i<NumPart+1;i++) 
        if(PMassMin>P[i].Mass&&P[i].Type!=GAS) PMassMin=P[i].Mass;
    printf("HiRes Part Mass: %e Msol\n",PMassMin*UdM);
/* ATTENIONE! i numeri dei gruppi partono da 1 non da 0, idem le
   particelle; idx parte da zero quindi la particella i ha indice idx[i-1] */
/* idx[] e' il vettore che associa ogni particella i al gruppo idx[i-1];
   se idx[i]=0 la particella non appartiene ai gruppi */
/* siccome e' indicizzato da idx[], anche hist[] parte da 1 */
    TotMass=0.0;
    Ngroups=0;
    for(i=1;i<NumPart+1;i++) 
    {
// if(idx[i-1]>0&&P[i].Mass<=PMassMin+0.001) {
        if(idx[i-1]>0) 
        {
            hist[idx[i-1]] += P[i].Mass;
            if(P[i].Type==DM)  hhist[idx[i-1]] += P[i].Mass;
            if(P[i].Type==GAS) hhist[idx[i-1]] += P[i].Mass;
            TotMass += P[i].Mass;
            if(idx[i-1]>Ngroups) Ngroups=idx[i-1];
        }
    }

    printf("Total fof groups mass: %e Msol.\n",TotMass*UdM);
    printf("Number of fof groups: %d\n",Ngroups);
    NunderMin=0; hist[0]=0.0;
    for(i=1;i<Ngroups+1;i++) {
        if(hist[i]/PMassMin<nMin) {
            hist[i]=0.0;
            NunderMin++;
        }
    }
    printf("Groups under min. num. of particles (%d, mass %e): %d\n", nMin,nMin*PMassMin*UdM,NunderMin);
    printf("Groups over threshold: %d\n",Ngroups-NunderMin);
/* rimetto le particelle dei gruppi tolti nel campo */
    for (j=1; j<NumPart+1; j++) if(hist[idx[j-1]]==0.0) idx[j-1]=0;

    G=malloc( (Ngroups+1)*sizeof(struct group_data));
    if(G==NULL) { printf("Impossibile allocare la memoria per i dati dei gruppi.\n"); exit(-1); }

    for(i=1;i<Ngroups+1;i++) { G[i].Mass=hist[i]; G[i].HiResMass=hhist[i]; }
 
    CenterGroups(Ngroups,idx);

    sprintf(bfr,"Groups.z%6.4f.MERGTR.dat",header1.redshift);
    grps_out=fopen(bfr,"w");
    for(i=1;i<Ngroups+1;i++) 
    {
        if(G[i].Mass>0.0) 
        {
            fprintf(grps_out,
                "%d %f %f %f %f %f %f %f %e %d\n",
                i,
                G[i].Pos[0], G[i].Pos[1], G[i].Pos[2],
                G[i].Vel[0], G[i].Vel[1], G[i].Vel[2],
                G[i].rad, G[i].Mass*UdM, G[i].npart);
        }
    }
    fclose(grps_out);

    return Ngroups;
}

void EvalProgenitors(char *fname, int *idx, int *sonidx, int npson, int Ngroups, gmap_t *gmap, int gmapLen)
{
    int Ncho, nprog;
    FILE *idx_out, *dm_out, *gas_out;
    FILE *f1, *f2;
    int i, j, k;
    char flnm[256], bfr[256];
    char *fltmp = "grp";

    nprog=find_progenitors(idx,sonidx,npson,Ngroups);

/*
    mbig=0.0; ibig=0;
    for(i=1; i<Ngroups+1; i++)
    if(G[i].HiResMass>mbig) { mbig=G[i].HiResMass; ibig=i;}
    printf("Bigger HiRes group (%d) mass %e\n",ibig,mbig);
    printf("Coordinates %f %f %f\n",G[ibig].Pos[0],G[ibig].Pos[1],G[ibig].Pos[2]);
*/

    sprintf(bfr,"idxfofMERGTR.z%6.4f.dat",header1.redshift);
    idx_out=fopen(bfr,"wb"); /*  questo contiene gli indici per eventuale follow-up */

/* counting the total number of part written in idx&ascii files */
    Ncho=0;
    for(i=0;i<Ngroups+1;i++) 
        if (G[i].isProgenitor) Ncho += G[i].npart;

    printf("Index file: %s  # PROGENITORS: %d  # of selected particles: %d\n", 
           bfr, nprog, Ncho);

    fwrite( &Ncho,  sizeof(Ncho),  1, idx_out); 
    fwrite( &nprog, sizeof(nprog), 1, idx_out); 

    snprintf(flnm, 256,"z%6.4f.%s%d.MERGTR.DM.dat",header1.redshift,fltmp, Ngroups+1); 
    dm_out=fopen(flnm,"w");
    snprintf(flnm, 256,"z%6.4f.%s%d.MERGTR.GAS.dat",header1.redshift,fltmp, Ngroups+1); 
    gas_out=fopen(flnm,"w");
#if 1
    for (i=0; i < Ngroups+1; i++) 
    {
        if (G[i].isProgenitor) 
        {
            for (j=1; j < NumPart+1; j++) 
            {
                if ( (i!=0&&idx[j-1]==i)||(i==0&&idx[j-1] == -10) ) 
                {
                    fprintf(P[j].Type == GAS ? gas_out : dm_out,
                        "%f %f %f %f %f %f %d %f %d\n",
                        P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],
                        P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],P[j].Type,P[j].Mass,i);

                    fwrite( &(Id[j]), sizeof(Id[j]), 1, idx_out);   
                }
            }
        }
    }
#else

    for (i=0; i < Ngroups+1; i++) 
    {
        if (G[i].isProgenitor) 
        {
            if (i==0)
            {
                for (j=1; j < NumPart+1; j++) 
                {
                    if (idx[j-1] < 0) 
                    {
                        fprintf(P[j].Type == GAS ? gas_out : dm_out,
                            "%f %f %f %f %f %f %d %f %d\n",
                            P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],
                            P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],P[j].Type,P[j].Mass,i);

                        fwrite( &(Id[j]), sizeof(Id[j]), 1, idx_out);   
                    }
                }
            }
            else
            {
                fprintf(stderr, "npart: %i\n", gmap[i].npart);
                for (k=0; k < gmap[i].npart; k++)
                {
                    j = gmap[i].ps[k];

                    fprintf(P[j].Type == GAS ? gas_out : dm_out,
                        "%f %f %f %f %f %f %d %f %d\n",
                        P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],
                        P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],P[j].Type,P[j].Mass,i);

                    fwrite( &(Id[j]), sizeof(Id[j]), 1, idx_out);   
                }
            }
        }
    }
#endif

    fprintf(dm_out,"# Here particle pos,vels are GADGET (COMOVING) \n");
    fclose(dm_out);
    fclose(gas_out);
    fclose(idx_out); 

    f2=fopen("j.inp","w"); /* questo file serve a DensitySphere */
    fprintf(f2,
        "%s\n"
        "%s\n" 
        "jzgas.z%6.4f.vir\n"        
        "gas.z%6.4f.vir\n"          
        "jzdm.z%6.4f.vir\n"       
        "dm.z%6.4f.vir\n"           
        "JL.z%6.4f.vir\n", 
        fname,
        bfr,
        header1.redshift,
        header1.redshift,
        header1.redshift,
        header1.redshift,
        header1.redshift);
    fclose(f2); 



    /* Undo marking (see find_progenitors()) */
    for(j=1;j<NumPart+1;j++) 
    {
        if (idx[j-1] == -10) idx[j-1] = 0;
        //if (idx[j-1] < 0) idx[j-1] = -(idx[j-1] + 1);
    }

    //snprintf(flnm, 256, "cp %s idxtmp.dat",bfr);
    //system(flnm);
}

void CenterGroups(int Ngroups,int* idx)
{
    int i,j,k;
    float* histr,dx,offs,d[3],r;

    histr=calloc(Ngroups+1, sizeof(float));
    if(histr==NULL) 
    { 
        printf("Impossibile allocare memoria per istogramma running delle masse (nel calcolo del baricentro)\n"); 
        exit(-1); 
    }
    
    printf("Calcolo i baricentri\n");
    for(i=0;i<Ngroups+1; i++) {
        G[i].rad= -1.e33;
        G[i].npart=0;
        for(j=0;j<3;j++) {
            G[i].Pos[j]=0.0;
            G[i].Vel[j]=0.0;
        }

        assert(G[i].Mass >= 0);
    }

/* Calcolo il baricentro, con condizioni periodiche al contorno */
    for(i=1;i<NumPart+1; i++) {
        k=idx[i-1];
        if (k==0) continue;           /* particella isolata, non considerare */
        if (G[k].Mass == 0) continue; /* gruppo troppo piccolo, non considerare */
        G[k].npart++;
        if(G[k].npart==1) 
        { 
            for(j=0;j<3;j++) 
            { 
                G[k].Pos[j] = P[i].Pos[j]*P[i].Mass; 
                G[k].Vel[j] = P[i].Vel[j]*P[i].Mass; 
            }
        } 
        else 
        {
            for(j=0;j<3;j++) 
            {
                G[k].Vel[j] += P[i].Vel[j]*P[i].Mass;
                dx = (G[k].Pos[j]/histr[k])-P[i].Pos[j];
                offs = (dx < 0) ? -header1.BoxSize : header1.BoxSize;
                if(dx>header1.BoxSize/2.0||dx< -header1.BoxSize/2.0) 
                    G[k].Pos[j] += (P[i].Pos[j]+offs)*P[i].Mass;
                else 
                    G[k].Pos[j] += P[i].Pos[j]*P[i].Mass;
            }
        }

        histr[k] += P[i].Mass;
    }

    for(i=1;i<Ngroups+1;i++) {
        if(histr[i]<=0.0) continue;
        for(j=0;j<3;j++) {
            G[i].Pos[j] /= histr[i];
            G[i].Vel[j] /= histr[i];
            if(G[i].Pos[j]<0.0) G[i].Pos[j] += header1.BoxSize;
            if(G[i].Pos[j]>header1.BoxSize) G[i].Pos[j] -= header1.BoxSize;
        }
    }

    for(i=1;i<NumPart+1;i++) {
        k=idx[i-1];
        for(j=0;j<3;j++) {
            d[j]= P[i].Pos[j]-G[k].Pos[j];
            if (d[j] > header1.BoxSize/2.0) 
                d[j] = P[i].Pos[j] - G[k].Pos[j] - header1.BoxSize;
            else if (d[j] < -header1.BoxSize/2.0) 
                d[j] = P[i].Pos[j] - G[k].Pos[j] + header1.BoxSize;
        }
        r=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
        if(r > G[k].rad) G[k].rad=r;
    } 

    free(histr);
}

/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int files)
{
    FILE *fd;
    char   buf[200];
    int    i,k,dummy,ntot_withmasses;
    int    n,pc,pc_new,pc_sph;
    
#define SKIP do { fread(&dummy, sizeof(dummy), 1, fd); } while (0)
//#define SKIP do { fseek(fd, sizeof(dummy), SEEK_CUR); } while (0)

    for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
        if(files>1)
	        sprintf(buf,"%s.%d",fname,i);
        else
	        sprintf(buf,"%s",fname);

        if(!(fd=fopen(buf,"r")))
	    {
	        printf("can`t open file %s\n",buf);
	        exit(0);
	    }

        printf("reading %s ...\n",buf); fflush(stdout);
    
        SKIP;
        fread(&header1, sizeof(header1), 1, fd);
        SKIP;

        if(files==1)
	    {
	        for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	            NumPart+= header1.npart[k];
	        Ngas= header1.npart[0];
	    }
        else
	    {
	        for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	            NumPart+= header1.npartTotal[k];
	        Ngas= header1.npartTotal[0];
	    }

        for(k=0, ntot_withmasses=0; k<5; k++)
	    {
	        if(header1.mass[k]==0)
	            ntot_withmasses+= header1.npart[k];
	    }
	    printf("N: %d\n",ntot_withmasses);

        if (i==0)
	        allocate_memory();

        SKIP;
        for(k=0,pc_new=pc;k<6;k++)
	    {
	        for(n=0;n<header1.npart[k];n++)
	        {
	            fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	            pc_new++;
	        }
	    }
        SKIP;
    
        SKIP;
        for(k=0,pc_new=pc;k<6;k++)
	    {
	        for(n=0;n<header1.npart[k];n++)
	        {
	            fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	            pc_new++;
	        }
	    }
        SKIP;
    

        SKIP;
        for(k=0,pc_new=pc;k<6;k++)
	    {
	        for(n=0;n<header1.npart[k];n++)
	        {
	            fread(&Id[pc_new], sizeof(int), 1, fd);
	            pc_new++;
	        }
	    }
        SKIP;


        if(ntot_withmasses>0)
	        SKIP;

        for(k=0, pc_new=pc; k<6; k++)
	    {
	        for(n=0;n<header1.npart[k];n++)
	        {
	            P[pc_new].Type=k;

	            if(header1.mass[k]==0)
		            fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	            else
		            P[pc_new].Mass= header1.mass[k];
	            pc_new++;
	        }
	    }
        if(ntot_withmasses>0)
	        SKIP;
      

        if(header1.npart[0]>0)
	    {
	        SKIP;
	        for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	        {
	            fread(&P[pc_sph].U, sizeof(float), 1, fd);
	            pc_sph++;
	        }
	        SKIP;

	        SKIP;
	        for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	        {
	            fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	            pc_sph++;
	        }
	        SKIP;

	        if(header1.flag_cooling)
	        {
	            SKIP;
	            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		        {
		            fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		            pc_sph++;
		        }
	            SKIP;
	        }
	        else
            {
	            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	            {
		            P[pc_sph].Ne= 1.0;
		            pc_sph++;
	            }
            }
	    }

        fclose(fd);
    }


    Time= header1.time;
    Redshift= header1.time;

    return 0;
}




/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
    printf("allocating memory...\n");

    if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
        fprintf(stderr,"failed to allocate memory.\n");
        exit(0);
    }
  
    P--;   /* start with offset 1 */

  
    if(!(Id=malloc(NumPart*sizeof(int))))
    {
        fprintf(stderr,"failed to allocate memory.\n");
        exit(0);
    }
  
    Id--;   /* start with offset 1 */

    printf("allocating memory...done\n");

    return 0;
}




/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
int reordering(void)
{
    int i;
    int idsource, idsave, dest;
    struct particle_data psave, psource;


    printf("reordering....\n");

    for(i=1; i<=NumPart; i++)
    {
        if(Id[i] != i)
        {
            psource= P[i];
            idsource=Id[i];
            dest=Id[i];

            do
            {
                psave= P[dest];
                idsave=Id[dest];

                P[dest]= psource;
                Id[dest]= idsource;
      
                if(dest == i) 
                    break;

                psource= psave;
                idsource=idsave;

                dest=idsource;
            }
            while(1);
        }
    }

    printf("done.\n");

    Id++;   
    free(Id);

    printf("space for particle ID freed\n");

    return 0;
}

int find_progenitors(int *idx, int *sonidx, int npson, int Ngroups)
{
    int i,j,nprog,nfield;
    float *hrpmass,frac;

    printf("Beginning search...\n");

    hrpmass= (float *)calloc(Ngroups+1, sizeof(float));
    if(hrpmass==NULL) {
        printf("Not enough mem for vector hrpmass! Exiting. \n");
        exit(-4);
    }
    nfield=0;

#if 0
    int ds=npson/10;
    for(j=0; j<npson; j++) 
    {
        if ((j % ds) == 0) printf("...%d\n",j); 
        for(i=1; i< NumPart+1; i++) 
        {
            //printf("Id[%i] = %i\n", i, Id[i]);
            assert(Id[i] == i);
            if (Id[i] == 0x1648)
            {
                fprintf(stderr, "-- %i %i %i %i\n", i, Id[i], sonidx[j], idx[i-1]);
            }

            if(Id[i]==sonidx[j])   /* this part will be in our obj */
            {
                if (idx[i-1]>0) 
                    hrpmass[idx[i-1]] += P[i].Mass;
                else 
                {
                    nfield++;
                    idx[i-1] = -10;  
                }
                break;
            }
        }
    }
#else
    for(j=0; j<npson; j++) 
    {
        i = sonidx[j];
        assert(Id[i] == i);

        if (idx[i-1]>0) 
            hrpmass[idx[i-1]] += P[i].Mass;
        else 
        {
            nfield++;
            //idx[i-1] = -(idx[i-1] + 1);  /* Mark for later (see EvalProgenitors()) */
            idx[i-1] = -10;  /* Mark for later (see EvalProgenitors()) */
        }
    }
#endif

    printf("Field particles: %d\n",nfield);
    G[0].npart        = nfield;
    G[0].isProgenitor = (nfield != 0);

    nprog=0;
    for(i=1; i<Ngroups+1;i++) {
        G[i].isProgenitor=0;
        frac= hrpmass[i]/G[i].HiResMass;
        if(frac>PROGPERC) { 
            G[i].isProgenitor=1; 
            nprog++; 
        }
    }

    return(nprog);

}  

