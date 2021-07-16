struct stepstruct{
 	vector vector; 
	double scalar, scalar2, acceptance_ran, energy_diff, boltzmann_factor; 
	int domain_move_type, regular_move_type, bond_index, m, b, coord, accept, particle_or_volume; 
};

typedef struct stepstruct step;
