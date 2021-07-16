int wrap_around_vector(vector wrap, vector position, vector box_length);
int translate_periodic(particle *the_particle, mparticle *the_type, vector translate, vector box_length);
int wrap_around(vector wrap, particle *the_particle, mparticle *the_type, vector box_length);
void print_out_of_bounds(particle *the_membrane, mparticle *the_types, vector box_length);
int check_bounds(particle part, mparticle type, vector box_length);