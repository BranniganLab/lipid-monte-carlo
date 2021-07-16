/***************************************************************************
 *            fourier.h
 *
 *  Mon Oct 17 17:35:23 2005
 *  Copyright  2005  User
 *  Email
 ****************************************************************************/

double fourier_step(  particle *the_membrane,   site *the_lattice, mparticle *the_types, int *accept,
                       movesize *fdelta, double previous_energy, vector box_length);
double generate_fourier_move(int *p, int *q, double *phase, movesize *fdelta);
void fourier_shift(particle *the_membrane, int p, int q, double phase, double shift, vector box_length);
void write_fdelta(movesize *fdelta);
void generate_fdelta(movesize *fdelta);
int read_in_fdelta(movesize *fdelta);
void adjust_fdelta(movesize *fdelta);
void set_fourier_constants(vector box_length);