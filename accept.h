int test_acceptance(double energy_diff);
void CopyParticleContents(  particle *part_to_copy,   particle *part_to_hold);
void copy_bead(bead *part_to_copy, bead *part_to_hold);
double CheckAcceptanceRate(int accept, double differential, int interval);
double GetAcceptanceRate(int accept, int interval);
double GetTotalAcceptanceRate(int accept, int i);