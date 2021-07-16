double get_bent(particle *part, mparticle *type);
double poly_energy(double radius, potential model, double core);
double mypow(double x, double y);
double get_angle(particle part, int i, int j, int k);
double poly_energy_noshift(double radius, potential model, double core);
double poly_attract_and_repulse(double radius, mtablesite model, double core);
int adjacent(bead bead1, bead bead2);
double get_bond_energy(particle *part, mparticle *type, vector box_length);