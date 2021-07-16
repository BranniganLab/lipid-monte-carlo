void WriteLog(char *outfile,   particle *the_membrane, double total_energy,
              int total_accept[], double move_size[], long step);
void write_heights(char *outfile1,char *outfile2,particle *the_membrane);
void WriteXYZ(char *outfile,   particle *the_membrane, mparticle *the_types);
void write_crd(char *outfile,   particle *the_membrane, char *code);
void write_frame(char *outfile,   particle *the_membrane, mparticle *the_types, char *code);
void write_simple_psf(char *outfile,   particle *the_membrane, mparticle *the_types);
void print_translation(int move);
void WriteDiffusion(char *outfile,   vector diffusion, long step);
void write_double_array(char *outfile, double data[], int num_rows, int num_columns, char *code);
void write_double_array2(char *outfile, double data[][2], int num_rows, char *code);
void write_time_double(char *outfile, double data);
void write_xydata(char *outfile, dataset xydata);
void print_potential(potential U);
void write_interval(particle *the_membrane, mparticle *the_types, double total_energy, int total_accept[],
                    double move_size[], long i, vector box_length);
void write_trajectory(char *outfile,particle *the_membrane,mparticle *the_types, 
                      double total_energy, double move_size[], long step);
void write_version();
void write_domains(char *outfile, particle *the_membrane);