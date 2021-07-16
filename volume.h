#define DEFINITELY_CONSTRAINED 1
#define POSSIBLY_CONSTRAINED 0

double volume_step(  particle *the_membrane,   site *the_lattice, mparticle *the_types,  int *accept,
					 double delta, double previous_energy, vector box_length);
double generate_box_shape_change(vector new_box, vector old_box, double delta);
void reshape_membrane(particle *the_membrane, site *the_lattice, mparticle *the_types,  vector vold, vector vnew);
int stretch_molecule(particle *part, particle *oldpart, mparticle *type, vector vnew, vector vold);