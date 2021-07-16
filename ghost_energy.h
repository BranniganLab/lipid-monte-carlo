double get_bead_energy_with_ghosts(bead the_bead, particle *the_membrane, site *the_lattice,  vector box_length);
double get_molecule_energy_with_ghosts(int index, particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length);
double get_total_energy_with_ghosts(particle *the_membrane, site *the_lattice,mparticle *the_types, vector box_length);
double get_molecule_pair_energy_with_ghosts(particle part1, particle part2, site *the_lattice, vector box_length);
double get_total_energy_with_ghosts_slow(particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length);
double get_bead_pair_energy_with_ghosts(bead bead1, bead bead2, vector box_length);