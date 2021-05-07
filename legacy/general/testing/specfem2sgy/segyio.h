/*---SEG-Y writing routines---*/

#ifndef SEGYIO_H
#define SEGYIO_H

#define SEGY_FILE_TEXT_HEADER_SIZE 3200
#define SEGY_FILE_HEADER_SIZE 3600
#define SEGY_TRACE_HEADER_SIZE 240

typedef struct
	{
	  char text_header[SEGY_FILE_TEXT_HEADER_SIZE];
	  int job_id;
	  int line;
	  int reel_n;
	  short Ntr_pa;
	  short Nat_pa;
	  short sample_interval;
	  short original_sample;
	  short ns;
	  short ns_original;
	  short code;
	  short ensemble_fold;
	  short sorting_code;
	  short vertical_sum_code;
	  short sw__freq_0;
	  short sw_freq_1;
	  short sw_l;
	  short sw_type;
	  short sw_ch;
	  short sw_taper0;
	  short sw_taper_l;
	  short taper_type;
	  short correlated;
	  short gain_recov;
	  short ampl_recov;
	  short mes_sys;
	  short imp_polarity;
	  short polarity_code;
	  char unassigned1[240];
	  short rev;
	  short fixed_length;
	  short num_ext_headers;
	  char unassigned2[94];
	}
	SEGY_file_header;
	
typedef struct
	{
	  int n_in_line;
	  int n_in_file;
	  int field_record_n;
	  int trace_field_record_n;
	  int esp_n;
	  int ensemble_n;
	  int trace_in_ensemple_n;
	  short trace_code;
	  short num_vert_sum;
	  short num_horiz_sum;
	  short data_use;
	  int dist_from_center;
	  int rec_elevation;
	  int src_surface_elevation;
	  int src_depth;
	  int datum_rec_elevation;
	  int datum_src_elevation;
	  int src_water_depth;
	  int group_water_depth;
	  short el_scalar;
	  short coord_scalar;
	  int src_X;
	  int src_Y;
	  int grp_X;
	  int grp_Y;
	  short coord_units;
	  short weathering_vel;
	  short subweathering_vel;
	  short uph_t_src;
	  short uph_t_grp;
	  short src_static;
	  short grp_static;
	  short tot_static;
	  short lag_A;
	  short lag_B;
	  short recording_delay;
	  short mute_start;
	  short mute_end;
	  short num_samples;
	  short sample_interval;
	  short gain_type;
	  short gain_constant;
	  short early_gain;
	  short correlated;
	  short sw__freq_0;
	  short sw_freq_1;
	  short sw_l;
	  short sw_type;
	  short sw_taper0;
	  short sw_taper_l;
	  short taper_type;
	  short alias_freq;
	  short alias_slope;
	  short notch_freq;
	  short notch_slope;
	  short low_cut_freq;
	  short high_cut_freq;
	  short low_cut_slope;
	  short high_cut_slope;
	  short year;
	  short day;
	  short hour;
	  short min;
	  short sec;
	  short time_basis_code;
	  short wf;
	  char chroup1[10];
	  int CDP_X;
	  int CDP_Y;
	  int in_line;
	  int cross_line;
	  int sp_n;
	  short scalar_spn;
	  short unit;
	  char transduction[6];
	  short transduction_unit;
	  short device_id;
	  short sc215;
	  short src_type;
	  char src_direction[6];
	  char scr_measure[6];
	  short src_measure_unit;
	  char unassigned[8];
	}
	SEGY_trace_header;
	
typedef struct
	{
	  int nx, ny, nz;
	  float x0,y0,z0;
	  float dx,dy,dz;
	  float ***data;
	}
	cube_type;
	
typedef struct
	{
	  int n_cr; //Number of crosslines, corrsponds to the third data dimension
	  int n_in; //Number of inlines, corresponds to the second data dimension 
	  int nz; //Number of trace elements, corresponds to the first data dimension
	  int cr_0; //Minimum crossline
	  int in_0; //Minimum inline
	  float z0; //Minimum Z (or time)
	  int d_cr,d_in; //Steps of crossline,inline
	  float dz; //Step of depth (time)
	  float *transform; //6 transform matrix coefficients, namely a11 a12 a13 a21 a22 a23
	  float ***data; //cube data[nz][n_in][n_cr]
	}
	cube_incr_type;

void swap2(short *x);
void swap4(void *x);

/*Simply rotate and shift coordinate system from inleine/crossline*/
/*to x (to the East) and y (to the North)*/
/*a must contain 6 elements treated as transform matrix coefficients*/
/*namely a11 a12 a13 a21 a22 a23*/
void incr2XY(float in, float cr, float *a, float *x, float *y);
	
void write3D_to_SEGY(cube_type &cube, const char *filename);
void write3D_to_SEGY_swapped(cube_type &cube, const char *filename);

/*Treating cube_incr info as inline/crossline and compute X,Y from transform*/
/*coefficient passed in transform member of cube_incr*/
/*Inlines must corrspond to the third data index*/
/*Crosslines must corrspond to the second data index*/
void write3D_to_SEGY_incr_swapped
(cube_incr_type &cube_incr, const char *filename, bool report_progress=0);

char *read_SEGY_text_header(FILE *f);

SEGY_file_header *read_SEGY_file_header(FILE *f, bool swap=0);

SEGY_trace_header *read_SEGY_trace_header(FILE *f, SEGY_trace_header &th, bool swap=0);

float *read_SEGY_trace_data(FILE *f, SEGY_trace_header &th, float *buf, bool swap_data=0);

float *read_SEGY_trace(FILE *f, SEGY_trace_header &th, float *buf, bool swap=0, bool swap_data=0);
	  
#endif
