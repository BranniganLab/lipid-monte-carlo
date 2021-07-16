int perturb_chain(particle *pert_part, mparticle *pert_part_type, int bead_index, double delta, vector box_length);
void uniform_perturb_position(particle *pert_part, mparticle *pert_type, vector translation, double delta,
							  vector box_length );
int pivot_end(particle *pert_part, mparticle *pert_part_type, double delta, int end);
int rotate_bead(particle *pert_part, mparticle *pert_part_type, double delta, int bead);
int GetPertIndex();
int perturb_domain(domain *the_domain, particle *the_membrane, mparticle *the_types, vector translation, vector rotation, double delta, vector box_length);
void uniform_perturb_radius(particle *pert_part,mparticle *pert_type, vector translation, double delta,
							  vector box_length, uconstraint the_constraint);
int perturb_particle(particle *pert_part, mparticle *pert_part_type, chainindex ci, vector translation, double delta, vector box_length, int code);
void uniform_perturb_z(particle *pert_part, mparticle *pert_part_type,vector translation, double delta,
							  vector box_length );
int pivot_interface(particle *pert_part, mparticle *pert_part_type, double delta);
int rotate_interface(particle *pert_part, mparticle *pert_part_type, double delta);