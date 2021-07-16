double get_bead_energy(bead the_bead, particle *the_membrane, site *the_lattice, vector box_length);
double get_bead_pair_energy(bead bead1, bead bead2,vector box_length);
double get_molecule_energy(int index, particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length);
double get_total_energy(particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length);
double get_molecule_pair_energy(particle part1, particle part2, site *the_lattice, vector box_length);
double calculate_user_energy(double radius, upotential the_potential);
double get_bead_user_energy(int index,particle *the_membrane, vector box_length);
double get_molecule_self_energy(particle part, vector box_length);