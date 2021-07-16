int check_checklists(site *the_lattice, vector old_box_length, vector new_box_length);
int get_checklists(site *the_lattice, vector box_length);
int initialize_lattice(particle *the_membrane, site *the_lattice, vector box_length);
int add_bead_to_lattice(chainindex ci,site *new_site);
int update_beads_lattice_site(bead *the_bead, site *the_lattice, vector box_length);
int move_bead_in_lattice(chainindex ci,site *old_site, site *new_site);
void free_node_list(cintnode *head);