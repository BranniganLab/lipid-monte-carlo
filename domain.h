#define MAX_DOMAINS 10

int fill_domain(domain *the_domain, particle *the_membrane, site *the_lattice, vector box_length);
double get_domain_energy(domain *the_domain, particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length);
double step_domain(domain *the_domain, particle *the_membrane, site *the_lattice, mparticle *the_types, double delta, vector box_length);
int initialize_domains(particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length);
int rotate_domain(domain *the_domain, particle *the_membrane, vector  rotation, vector box_length);
int translate_domain(domain *the_domain, particle *the_membrane, mparticle *the_types, vector  translation, vector box_length);
double get_intra_domain_energy(domain *the_domain, particle *the_membrane, site *the_lattice, vector box_length);
int get_domain_center_and_tilt(domain *the_domain, particle *the_membrane, mparticle *the_types,  vector box_length);
particle *store_domain_positions(domain *the_domain, particle *the_membrane);
int retrieve_domain_positions(domain *the_domain, particle *the_membrane, particle *domain_molecules);
double inter_domain_energy(particle *the_membrane, vector box_length);