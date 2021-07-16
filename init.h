void finit(particle *the_membrane, cell *cells, double *total_energy, int total_accept[], double move_size[], long step);
int frecenter_membrane(particle *the_membrane, vector box_length);
void init_potential(potential U);
void free_membrane(particle *the_membrane);
void init(particle *the_membrane, site *the_lattice, mparticle *the_types,  double *total_energy, int total_accept[], double move_size[], long step);
int set_chainpos(mparticle *the_types);