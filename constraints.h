/***************************************************************************
 *            constraints.h
 *
 *  Tue Apr  5 16:50:32 2005
 *  Copyright  2005  User
 *  Email
 ****************************************************************************/

int move_molecule_to_point(particle *the_membrane,chainindex ci, vector new_point, vector box_length);
int enforce_constraints(particle *the_membrane, vector box_length);
int print_constraints(particle *the_membrane, mparticle *the_types, vector box_length);
