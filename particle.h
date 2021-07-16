double particle_move(particle *the_membrane,  site *the_lattice, mparticle *the_types, int attempts[],int accepts[])
;
double step_dof(int pert_index, particle *the_membrane, site *the_lattice, mparticle *the_types, int *accept, double *delta, vector box_length);
double step_bead(int pert_index, particle *the_membrane, site *the_lattice, mparticle *the_types, int *accept, double *delta,  vector box_length);
double step_bond(int pert_index, particle *the_membrane, site *the_lattice, mparticle *the_types,  int *accept, double *delta,  vector box_length);