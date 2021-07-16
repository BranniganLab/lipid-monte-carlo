double energy_tests(particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length);
double slow_total_energy(particle *the_membrane, mparticle *the_types, vector box_length);
double slow_total_energy_with_ghosts(particle *the_membrane, mparticle *the_types, vector box_length);
int check_all_ghosts(particle *the_membrane, site **the_lattice, vector box_length);
double check_molecule_energy(particle *the_membrane, site *the_lattice, mparticle *the_types, int index, vector box_length);
int check_bond_lengths(particle *part, mparticle *type);
double get_hypothetical_da_with_bond(particle *part, mparticle *type, vector box_length, double dA);
double get_hypothetical_bond_stretch(particle *part, mparticle *type, vector box_length, double dA);