#include <stdio.h>

#define GRID_RESOLUTION 4
#define EPSILON 0.005
#define EVAPORATE 4
#define NUM_TRIALS 1
#define MIN_LATTICE 0.2
#define USE_LIST 1
#define USE_CELL_LIST 0
#define NUM_KEYS 5
#define USE_INCLUSIONS 0
#define NUM_CELLS 5
#define USE_SWAP 0
#define MAX_DOMAIN_PARTS 1000

/* symmetries  */
#define VECTOR 1
#define DIRECTOR 0

/*boundary types*/
#define NONE 0
#define MEAN_FIELD 1
#define PERIODIC 2
#define PLANE 3 /*periodic except in Z direction*/
#define BOUNDARY PERIODIC

#define PARTICLE_MOVE 0
#define VOLUME_MOVE 1
#define EDGE 1
#define BULK 0
#define ALL  -1
#define DIMENSION 3
#define NUM_BINS 16
#define VAR_TEST 500
#define PRECISION 1E-30
#define PI M_PI
#define MAX_HOURS 23
#define RESOLUTION 3
#define PRINT_RESOLUTION 3
#define MAX_CONST 5
#define MAX_BONDS 20
#define X 0
#define Y 1
#define Z 2
#define ROTATION 3
#define FLIP 4
#define VOLUME 5
#define CLUSTER 7
#define PATCH 6
#define INCLUSION 2
#define LIPID  1
#define NUM_DOF 8
#define VARIATION 2
#define XPRINT 2
#define YPRINT 2
#define ZPRINT 2
#define INFNTY 10000000
#define	ERROR	1
#define CONST 3
#define DIM     2               /* Dimension of points */
#define PMAX    1000            /* Max # of pts in polygon */
#define OVERLAP 1
#define NO_OVERLAP 0
#define OUT_OF_BOUNDS 1
#define IN_BOUNDS 0
#define CENTER 0
#define ARROW 1
#define HINTERFACE 0
#define HTAIL 1
#define HREPULSE 3
#define HSOFTCORE 2
#define HEAD 0
#define INTERFACE 1
#define TAIL 2
#define TAIL2 3
#define SLOW 0
#define FAST 1
#define YES 1
#define NO 0
#define NOT_FOUND  -1
#define FOUND  1
#define UPDATE 1
#define NO_UPDATE 0
#define HARDCORE 0
#define SOFTCORE 1
#define NO_ERROR 10
#define HARMONIC 0
#define LJ 1
#define FIXED_POINT 0
#define FIXED_ANGLE 1
#define PERTURB_WHOLE_CHAIN 0
#define PERTURB_BEAD 1 
#define PERTURB_BOND 2

/*************************type definitions***************************/


typedef double vector[DIMENSION];


struct movesizestruct
{
   int attempts; 
   int accepts; 
   double size; 
};

typedef struct movesizestruct movesize; 

struct intnodestruct
{
    int contents;
    struct intnodestruct *previous;
    struct intnodestruct *next;
};

typedef struct intnodestruct intnode;

struct Btypepotential
{
	double range; 
	double coefficient; 
	double exponent; 
	double shift; 
};

typedef struct Btypepotential btypepotential; 

struct Mtablesite
{
	btypepotential attract, repulse; 
};

typedef struct Mtablesite mtablesite; 

struct nnstruct
{
	int index;
 	struct nn *next;
};

typedef struct nnstruct nn;

typedef struct Chainindex
{ int m;
	int b;
} chainindex;


struct cintnodestruct
{
    chainindex contents;
    struct cintnodestruct *previous;
    struct cintnodestruct *next;
};

typedef struct cintnodestruct cintnode;

typedef struct Site
{
	int *checklist; 
	cintnode *bead_list; 
}site; 

struct cellstruct
{
	cintnode *head;
};

typedef struct cellstruct cell;

typedef struct Lrnode
{
	chainindex item;
	int tag;
	struct Lrnode *left;
	struct Lrnode *right;
	int height;
} lrnode;



typedef struct Tree
{ lrnode *root;
	int size; } tree;

typedef struct Pair{
	lrnode *parent;
	lrnode *child; } pair;
	
	
struct mbeadstruct{
/*bead for a model structure*/
    int chainpos;  
	int type; 
    double radius;
	int user;  
};

typedef struct mbeadstruct mbead;

	
struct beadstruct{
/*bead for an actual lipid*/
    vector position;
    chainindex ci;
	int site_index; 
	int num_ghosts; 
	struct beadstruct *ghosts;  
	mbead *typeptr; 
    tree *ntree;
    vector init;
	int user;  
};

typedef struct beadstruct bead;



struct Domain
{
/*inclusion or fixed region*/
	int size; 
	vector center; 
	vector arrow; 
	double radius; 
	double spacing; 
	int partlist[MAX_DOMAIN_PARTS]; 
	//int chain_length; 
	int shells; 
}; 
typedef struct Domain domain; 

struct bondstruct{
int b1; 
int b2; 
double cbond; 
double equil; 
int constrained; };

typedef struct bondstruct bond;

struct anglestruct{
    int b1; 
    int b2; 
    int b3; 
    double cbend; 
    double equil; 
};

typedef struct anglestruct bondangle;

struct  particlestruct
{
/*actual lipid particle*/
	bead chain[MAX_BONDS];
	double cosangle[MAX_BONDS];
	vector virials[MAX_BONDS + 2];
	double energy;
	int chain_length;
	vector crossings;
	int cluster;
	int domain_index;
	domain *domainptr; 
	int model_index; 
};

typedef struct particlestruct particle;

struct  mparticlestruct
{
/*model particle for an actual lipid*/
	mbead chain[MAX_BONDS];
	bondangle bondangles[MAX_BONDS];
	bond bonds[MAX_BONDS]; 
	int num_bonds; 
	int num_angles;  
	double cbend;
	double length;
	int ibindex; 
	int chain1start, chain2start; /*index where each chain starts- chain2start = -1 if only one chain*/
	int chain_length;
	int num_chains; 
};

typedef struct mparticlestruct mparticle; 


typedef particle *ensemble;


struct potentialstruct
{
	double range;
	double falloff;
	int num_terms;
	double coefficients[5][3];
};

typedef struct potentialstruct potential;

struct Upotential
{	
	int type; 
	double r0;
	double coefficient; 
	double range; 
	chainindex ci1, ci2; 
	double exponent;  
}; 

typedef struct Upotential upotential; 

struct Uconstraint
{	
	int type; 
	chainindex ci; 
	vector point; 
	vector arrow; 	
}; 

typedef struct Uconstraint uconstraint;

struct hamiltonianstruct
{
	double head;
 	double center;
  	double tail;
   	double bend;
    double total;
    double samples;
} ;

typedef struct hamiltonianstruct hamiltonian;

struct averagestruct
{
	double indep;
	double depend;
	int num_samples;
};

typedef struct averagestruct average;

struct datasetstruct
{
	average *points;
	int num_points;
};

typedef struct datasetstruct dataset;

typedef enum { FALSE, TRUE }	bool;


/***************************function headers****************************/

double ran3(long *idum);

inline static void vsub_periodic(vector vr, vector v1, vector v2, vector box_length)
{
	vr[X] = v1[X] - v2[X];
	vr[Y] = v1[Y] - v2[Y];
	vr[Z] = v1[Z] - v2[Z];
	
	if (BOUNDARY == PERIODIC){
		vr[X] = (vr[X] > box_length[X]/2) ? -(box_length[X] - vr[X]):vr[X];
		vr[X] = (vr[X] < -box_length[X]/2) ? box_length[X] + vr[X] : vr[X];
		vr[Y] = (vr[Y] > box_length[Y]/2) ? -(box_length[Y] - vr[Y]):vr[Y];
		vr[Y] = (vr[Y] < -box_length[Y]/2) ? box_length[Y] + vr[Y] : vr[Y];
		vr[Z] = (vr[Z] > box_length[Z]/2) ? -(box_length[Z] - vr[Z]):vr[Z];
		vr[Z] = (vr[Z] < -box_length[Z]/2) ? box_length[Z] + vr[Z] : vr[Z];  }
}
