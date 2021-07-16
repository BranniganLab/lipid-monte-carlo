/***************************************************************************
 *            cluster.h
 *
 *  Wed Aug 31 14:16:36 2005
 *  Copyright  2005  User
 *  Email
 ****************************************************************************/
#define MAX_CLUSTER_SIZE 200
#define MAX_SHELL_SIZE 300

double cluster_move(domain *the_domain, particle *the_membrane, site *the_lattice, mparticle *the_types, int domain_index, vector box_length);
int cluster_swap(particle *the_membrane, mparticle *the_types, domain *the_domain, int the_cluster[], int cluster_size, int the_shell[], int shell_size,  vector box_length);
int find_cluster_slow(particle *the_membrane, int center_index, double cradius, int the_cluster[MAX_CLUSTER_SIZE], vector box_length);
int test_same_cluster(particle part1, particle part2, double cradius, vector box_length);
int find_coordination_shell_slow(particle *the_membrane, domain *the_domain, double cradius, int the_shell[MAX_SHELL_SIZE], vector box_length);
double hybrid_cluster_volume_move( domain *the_domain, particle *the_membrane,   site *the_lattice, mparticle *the_types, int *accept,
                       double delta,  vector box_length);
