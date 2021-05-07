#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "segyio.h"

#include "subio.h"

void swapwrite4(void *x, FILE *f)
{
  char *p = (char*)x;
  char buf;
  
  buf=p[3]; p[3]=p[0]; p[0]=buf;
  buf=p[2]; p[2]=p[1]; p[1]=buf;
  fwrite (p,4,1,f);
}

void swapwrite2(void *x, FILE *f)
{
  char *p = (char*)x;
  char buf;
  
  buf=p[1]; p[1]=p[0]; p[0]=buf;
  fwrite (p,2,1,f);
}

void swap2(short *x)
{
  char *p = (char*)x;
  char buf;
  
  buf=p[1]; p[1]=p[0]; p[0]=buf;
}

void swap4(void *x)
{
  char *p = (char*)x;
  char buf;
  
  buf=p[3]; p[3]=p[0]; p[0]=buf;
  buf=p[2]; p[2]=p[1]; p[1]=buf;
}

void SEGY_swap_file_header(SEGY_file_header &fh)
{
  swap4(&(fh.job_id));
  swap4(&(fh.line));
  swap4(&(fh.reel_n));
  swap2(&(fh.Ntr_pa));
  swap2(&(fh.Nat_pa));
  swap2(&(fh.sample_interval));
  swap2(&(fh.original_sample));
  swap2(&(fh.ns));
  swap2(&(fh.ns_original));
  swap2(&(fh.code));
  swap2(&(fh.ensemble_fold));
  swap2(&(fh.sorting_code));
  swap2(&(fh.vertical_sum_code));
  swap2(&(fh.sw__freq_0));
  swap2(&(fh.sw_freq_1));
  swap2(&(fh.sw_l));
  swap2(&(fh.sw_type));
  swap2(&(fh.sw_ch));
  swap2(&(fh.sw_taper0));
  swap2(&(fh.sw_taper_l));
  swap2(&(fh.taper_type));
  swap2(&(fh.correlated));
  swap2(&(fh.gain_recov));
  swap2(&(fh.ampl_recov));
  swap2(&(fh.mes_sys));
  swap2(&(fh.imp_polarity));
  swap2(&(fh.polarity_code));
  swap2(&(fh.rev));
  swap2(&(fh.fixed_length));
  swap2(&(fh.num_ext_headers));
}

void SEGY_swap_trace_header(SEGY_trace_header &th)
{
  swap4(&(th.n_in_line));					
  swap4(&(th.n_in_file));
  swap4(&(th.field_record_n));
  swap4(&(th.trace_field_record_n));
  swap4(&(th.esp_n));
  swap4(&(th.ensemble_n));
  swap4(&(th.trace_in_ensemple_n));
  swap2(&(th.trace_code));
  swap2(&(th.num_vert_sum));
  swap2(&(th.num_horiz_sum));
  swap2(&(th.data_use));
  swap4(&(th.dist_from_center));
  swap4(&(th.rec_elevation));
  swap4(&(th.src_surface_elevation));
  swap4(&(th.src_depth));
  swap4(&(th.datum_rec_elevation));
  swap4(&(th.datum_src_elevation));
  swap4(&(th.src_water_depth));
  swap4(&(th.group_water_depth));
  swap2(&(th.el_scalar));
  swap2(&(th.coord_scalar));
  swap4(&(th.src_X));
  swap4(&(th.src_Y));
  swap4(&(th.grp_X));
  swap4(&(th.grp_Y));
  swap2(&(th.coord_units));
  swap2(&(th.weathering_vel));
  swap2(&(th.subweathering_vel));
  swap2(&(th.uph_t_src));
  swap2(&(th.uph_t_grp));
  swap2(&(th.src_static));
  swap2(&(th.grp_static));
  swap2(&(th.tot_static));
  swap2(&(th.lag_A));
  swap2(&(th.lag_B));
  swap2(&(th.recording_delay));
  swap2(&(th.mute_start));
  swap2(&(th.mute_end));
  swap2(&(th.num_samples));
  swap2(&(th.sample_interval));
  swap2(&(th.gain_type));
  swap2(&(th.gain_constant));
  swap2(&(th.early_gain));
  swap2(&(th.correlated));
  swap2(&(th.sw__freq_0));
  swap2(&(th.sw_freq_1));
  swap2(&(th.sw_l));
  swap2(&(th.sw_type));
  swap2(&(th.sw_taper0));
  swap2(&(th.sw_taper_l));
  swap2(&(th.taper_type));
  swap2(&(th.alias_freq));
  swap2(&(th.alias_slope));
  swap2(&(th.notch_freq));
  swap2(&(th.notch_slope));
  swap2(&(th.low_cut_freq));
  swap2(&(th.high_cut_freq));
  swap2(&(th.low_cut_slope));
  swap2(&(th.high_cut_slope));
  swap2(&(th.year));
  swap2(&(th.day));
  swap2(&(th.hour));
  swap2(&(th.min));
  swap2(&(th.sec));
  swap2(&(th.time_basis_code));
  swap2(&(th.wf));
  swap4(&(th.CDP_X));
  swap4(&(th.CDP_Y));
  swap4(&(th.in_line));
  swap4(&(th.cross_line));
  swap4(&(th.sp_n));
  swap2(&(th.scalar_spn));
  swap2(&(th.unit));
  swap2(&(th.transduction_unit));
  swap2(&(th.device_id));
  swap2(&(th.sc215));
  swap2(&(th.src_type));
  swap2(&(th.src_measure_unit));
}

void SEGY_file_header_swapwrite(SEGY_file_header &fh, FILE *fo)
{
  fwrite(fh.text_header,1,3200,fo);
  swapwrite4(&(fh.job_id),fo);
  swapwrite4(&(fh.line),fo);
  swapwrite4(&(fh.reel_n),fo);
  swapwrite2(&(fh.Ntr_pa),fo);
  swapwrite2(&(fh.Nat_pa),fo);
  swapwrite2(&(fh.sample_interval),fo);
  swapwrite2(&(fh.original_sample),fo);
  swapwrite2(&(fh.ns),fo);
  swapwrite2(&(fh.ns_original),fo);
  swapwrite2(&(fh.code),fo);
  swapwrite2(&(fh.ensemble_fold),fo);
  swapwrite2(&(fh.sorting_code),fo);
  swapwrite2(&(fh.vertical_sum_code),fo);
  swapwrite2(&(fh.sw__freq_0),fo);
  swapwrite2(&(fh.sw_freq_1),fo);
  swapwrite2(&(fh.sw_l),fo);
  swapwrite2(&(fh.sw_type),fo);
  swapwrite2(&(fh.sw_ch),fo);
  swapwrite2(&(fh.sw_taper0),fo);
  swapwrite2(&(fh.sw_taper_l),fo);
  swapwrite2(&(fh.taper_type),fo);
  swapwrite2(&(fh.correlated),fo);
  swapwrite2(&(fh.gain_recov),fo);
  swapwrite2(&(fh.ampl_recov),fo);
  swapwrite2(&(fh.mes_sys),fo);
  swapwrite2(&(fh.imp_polarity),fo);
  swapwrite2(&(fh.polarity_code),fo);
  fwrite(fh.unassigned1,1,240,fo);
  swapwrite2(&(fh.rev),fo);
  swapwrite2(&(fh.fixed_length),fo);
  swapwrite2(&(fh.num_ext_headers),fo);
  fwrite(fh.unassigned2,1,94,fo);
}

void SEGY_trace_header_swapwrite(SEGY_trace_header &th, FILE *fo)
{
  //char buf[256]="1234567890";
  swapwrite4(&(th.n_in_line),fo);					
  swapwrite4(&(th.n_in_file),fo);
  swapwrite4(&(th.field_record_n),fo);
  swapwrite4(&(th.trace_field_record_n),fo);
  swapwrite4(&(th.esp_n),fo);
  swapwrite4(&(th.ensemble_n),fo);
  swapwrite4(&(th.trace_in_ensemple_n),fo);
  swapwrite2(&(th.trace_code),fo);
  swapwrite2(&(th.num_vert_sum),fo);
  swapwrite2(&(th.num_horiz_sum),fo);
  swapwrite2(&(th.data_use),fo);
  swapwrite4(&(th.dist_from_center),fo);
  swapwrite4(&(th.rec_elevation),fo);
  swapwrite4(&(th.src_surface_elevation),fo);
  swapwrite4(&(th.src_depth),fo);
  swapwrite4(&(th.datum_rec_elevation),fo);
  swapwrite4(&(th.datum_src_elevation),fo);
  swapwrite4(&(th.src_water_depth),fo);
  swapwrite4(&(th.group_water_depth),fo);
  swapwrite2(&(th.el_scalar),fo);
  swapwrite2(&(th.coord_scalar),fo);
  swapwrite4(&(th.src_X),fo);
  swapwrite4(&(th.src_Y),fo);
  swapwrite4(&(th.grp_X),fo);
  swapwrite4(&(th.grp_Y),fo);
  swapwrite2(&(th.coord_units),fo);
  swapwrite2(&(th.weathering_vel),fo);
  swapwrite2(&(th.subweathering_vel),fo);
  swapwrite2(&(th.uph_t_src),fo);
  swapwrite2(&(th.uph_t_grp),fo);
  swapwrite2(&(th.src_static),fo);
  swapwrite2(&(th.grp_static),fo);
  swapwrite2(&(th.tot_static),fo);
  swapwrite2(&(th.lag_A),fo);
  swapwrite2(&(th.lag_B),fo);
  swapwrite2(&(th.recording_delay),fo);
  swapwrite2(&(th.mute_start),fo);
  swapwrite2(&(th.mute_end),fo);
  swapwrite2(&(th.num_samples),fo);
  swapwrite2(&(th.sample_interval),fo);
  swapwrite2(&(th.gain_type),fo);
  swapwrite2(&(th.gain_constant),fo);
  swapwrite2(&(th.early_gain),fo);
  swapwrite2(&(th.correlated),fo);
  swapwrite2(&(th.sw__freq_0),fo);
  swapwrite2(&(th.sw_freq_1),fo);
  swapwrite2(&(th.sw_l),fo);
  swapwrite2(&(th.sw_type),fo);
  swapwrite2(&(th.sw_taper0),fo);
  swapwrite2(&(th.sw_taper_l),fo);
  swapwrite2(&(th.taper_type),fo);
  swapwrite2(&(th.alias_freq),fo);
  swapwrite2(&(th.alias_slope),fo);
  swapwrite2(&(th.notch_freq),fo);
  swapwrite2(&(th.notch_slope),fo);
  swapwrite2(&(th.low_cut_freq),fo);
  swapwrite2(&(th.high_cut_freq),fo);
  swapwrite2(&(th.low_cut_slope),fo);
  swapwrite2(&(th.high_cut_slope),fo);
  swapwrite2(&(th.year),fo);
  swapwrite2(&(th.day),fo);
  swapwrite2(&(th.hour),fo);
  swapwrite2(&(th.min),fo);
  swapwrite2(&(th.sec),fo);
  swapwrite2(&(th.time_basis_code),fo);
  swapwrite2(&(th.wf),fo);
  fwrite(th.chroup1,1,10,fo);
  swapwrite4(&(th.CDP_X),fo);
  swapwrite4(&(th.CDP_Y),fo);
  swapwrite4(&(th.in_line),fo);
  swapwrite4(&(th.cross_line),fo);
  swapwrite4(&(th.sp_n),fo);
  swapwrite2(&(th.scalar_spn),fo);
  swapwrite2(&(th.unit),fo);
  fwrite(th.transduction,1,6,fo);
  swapwrite2(&(th.transduction_unit),fo);
  swapwrite2(&(th.device_id),fo);
  swapwrite2(&(th.sc215),fo);
  swapwrite2(&(th.src_type),fo);
  fwrite(th.src_direction,1,6,fo);
  fwrite(th.scr_measure,1,6,fo);
  swapwrite2(&(th.src_measure_unit),fo);
  fwrite(th.unassigned,1,8,fo);
  //fwrite(buf,1,8,fo);
}

/*Simply rotate and shift coordinate system from inleine/crossline*/
/*to x (to the East) and y (to the North)*/
/*a must contain 6 elements treated as transform matrix coefficients*/
/*namely a11 a12 a13 a21 a22 a23*/
void incr2XY(float in, float cr, float *a, float *x, float *y)
{
  *x=a[0]*in+a[1]*cr+a[2];
  *y=a[3]*in+a[4]*cr+a[5];
}

void write3D_to_SEGY(cube_type &cube, const char *filename)
{
  int i,j,k;
  SEGY_file_header fh;
  SEGY_trace_header th;
  char *buf;
  
  FILE *ftest;
  
  FILE *fo=file_open(filename,"wb");
  
  
  if (sizeof(fh)!=SEGY_FILE_HEADER_SIZE)
  {
    fprintf(stderr,"SEG-Y file header size %lu != %d\n",sizeof(fh),SEGY_FILE_HEADER_SIZE);
    exit(EXIT_FAILURE);
  }
  if (sizeof(th)!=SEGY_TRACE_HEADER_SIZE)
  {
    fprintf(stderr,"SEG-Y file header size %lu != %d\n",sizeof(th),SEGY_TRACE_HEADER_SIZE);
    exit(EXIT_FAILURE);
  }
    
  buf=(char*)(&fh);
  for (i=0;i<SEGY_FILE_HEADER_SIZE;i++)
    buf[i]=0;
  buf=(char*)(&th);
  for (i=0;i<SEGY_TRACE_HEADER_SIZE;i++)
    buf[i]=0;
  
  strcpy(fh.text_header,"3D cube");
  //fh.text_header[3200]='q';
  fh.job_id=fh.line=fh.reel_n=1;
  fh.sample_interval = cube.dz;
  fh.ns = cube.nz;
  fh.code=5; //IEEE floating-point format
  fh.ensemble_fold=cube.nx*cube.ny;
  fh.sorting_code=1; //Horizontaly stacked
  fh.mes_sys=1; //Meters
  //fh.rev=0x100; // rev 1.0
  fh.rev=0; // rev 0
  fh.fixed_length=1;
  fh.num_ext_headers=0;
  
  fwrite(&fh,sizeof(fh),1,fo);
  
  th.field_record_n=1;
  th.coord_units=1;
  th.num_samples=cube.nz;
  th.sample_interval=cube.dz;
  for (i=0;i<cube.ny;i++)
  {
    th.src_X=th.grp_X=cube.y0+i*cube.dy;
    for (j=0;j<cube.nx;j++)
    {
      th.trace_field_record_n=th.n_in_file=th.n_in_line=i*cube.nx+j;
      th.trace_code=1;
      th.src_Y=th.grp_Y=cube.x0+j*cube.dx;
      fwrite(&th,sizeof(th),1,fo);
      if (i==cube.ny/2 && j==cube.nx/2)
      {
        //-------------
        sprintf(buf,"%0d-%0d.dat",i,j);
        ftest=file_open(buf,"wt");
        //-------------
      }
      for (k=0;k<cube.nz;k++)
      {
	//-------------------
	//cube.data[k][j][i]=10.0*pow(sin((float)k*cube.dz*(5.0*M_PI/((cube.nz-1)*cube.dz))),2.0);
	fwrite(&(cube.data[k][j][i]),sizeof(float),1,fo);
	//-------------------
	if (i==cube.ny/2 && j==cube.nx/2)
        {
	  //-------------
	  fprintf(ftest,"%d\t%g\n",k,cube.data[k][j][i]);
	  //-------------
	}
      }
      if (i==cube.ny/2 && j==cube.nx/2)
      {
        fclose(ftest);
      }
    }
  }
  fclose(fo);
}

void write3D_to_SEGY_swapped(cube_type &cube, const char *filename)
{
  int i,j,k;
  SEGY_file_header fh;
  SEGY_trace_header th;
  char *buf;
  
  FILE *ftest;
  
  FILE *fo=file_open(filename,"wb");
  
  
  if (sizeof(fh)!=SEGY_FILE_HEADER_SIZE)
  {
    fprintf(stderr,"SEG-Y file header size %lu != %d\n",sizeof(fh),SEGY_FILE_HEADER_SIZE);
    exit(EXIT_FAILURE);
  }
  if (sizeof(th)!=SEGY_TRACE_HEADER_SIZE)
  {
    fprintf(stderr,"SEG-Y file header size %lu != %d\n",sizeof(th),SEGY_TRACE_HEADER_SIZE);
    exit(EXIT_FAILURE);
  }
    
  buf=(char*)(&fh);
  for (i=0;i<SEGY_FILE_HEADER_SIZE;i++)
    buf[i]=0;
  buf=(char*)(&th);
  for (i=0;i<SEGY_TRACE_HEADER_SIZE;i++)
    buf[i]=0;
  
  strcpy(fh.text_header,"3D cube");
  //fh.text_header[3200]='q';
  fh.job_id=fh.line=fh.reel_n=1;
  fh.sample_interval = cube.dz;
  fh.ns = cube.nz;
  fh.code=5; //IEEE floating-point format
  fh.ensemble_fold=cube.nx*cube.ny;
  fh.sorting_code=1; //Horizontaly stacked
  fh.mes_sys=1; //Meters
  //fh.rev=0x100; // rev 1.0
  fh.rev=0; // rev 0
  fh.fixed_length=1;
  fh.num_ext_headers=0;
  
  //SEGY_file_header_swapwrite(fh,fo);
  fwrite(&fh,sizeof(fh),1,fo);
  
  th.field_record_n=1;
  th.coord_units=1;
  th.num_samples=cube.nz;
  th.sample_interval=cube.dz;
  th.src_measure_unit=1;
  for (i=0;i<cube.ny;i++)
  {
    th.src_X=th.grp_X=cube.y0+i*cube.dy;
    for (j=0;j<cube.nx;j++)
    {
      th.trace_field_record_n=th.n_in_file=th.n_in_line=i*cube.nx+j;
      th.trace_code=1;
      th.src_Y=th.grp_Y=cube.x0+j*cube.dx;
      //SEGY_trace_header_swapwrite(th,fo);
      fwrite(&th,sizeof(th),1,fo);
      
      if (i==cube.ny/2 && j==cube.nx/2)
      {
        //-------------
        sprintf(buf,"%0d-%0d.dat",i,j);
        ftest=file_open(buf,"wt");
        //-------------
      }
      for (k=0;k<cube.nz;k++)
      {
	//-------------------
	//cube.data[k][j][i]=10.0*pow(sin((float)k*cube.dz*(5.0*M_PI/((cube.nz-1)*cube.dz))),2.0);
	//fwrite(&(cube.data[k][j][i]),sizeof(float),1,fo);
	swapwrite4(&(cube.data[k][j][i]),fo);
	//-------------------
	if (i==cube.ny/2 && j==cube.nx/2)
        {
	  //-------------
	  fprintf(ftest,"%d\t%g\n",k,cube.data[k][j][i]);
	  //-------------
	}
      }
      if (i==cube.ny/2 && j==cube.nx/2)
      {
        fclose(ftest);
      }
    }
  }
  fclose(fo);
}

/*Treating cube_incr info as inline/crossline and compute X,Y from transform*/
/*coefficient passed in transform member of cube_incr*/
/*Inlines must corrspond to the third data index*/
/*Crosslines must corrspond to the second data index*/
void write3D_to_SEGY_incr_swapped
(cube_incr_type &cube_incr, const char *filename, bool report_progress)
{
  int i,j,k;
  SEGY_file_header fh;
  SEGY_trace_header th;
  float x,y;
  char *buf;
  
  FILE *fo=file_open(filename,"wb");
  
  
  if (sizeof(fh)!=SEGY_FILE_HEADER_SIZE)
  {
    fprintf(stderr,"SEG-Y file header size %lu != %d\n",sizeof(fh),SEGY_FILE_HEADER_SIZE);
    exit(EXIT_FAILURE);
  }
  if (sizeof(th)!=SEGY_TRACE_HEADER_SIZE)
  {
    fprintf(stderr,"SEG-Y file header size %lu != %d\n",sizeof(th),SEGY_TRACE_HEADER_SIZE);
    exit(EXIT_FAILURE);
  }
    
  buf=(char*)(&fh);
  for (i=0;i<SEGY_FILE_HEADER_SIZE;i++)
    buf[i]=0;
  buf=(char*)(&th);
  for (i=0;i<SEGY_TRACE_HEADER_SIZE;i++)
    buf[i]=0;
  
  strcpy(fh.text_header,"3D cube");
  fh.job_id=fh.line=fh.reel_n=1;
  fh.sample_interval = cube_incr.dz;
  fh.ns = cube_incr.nz;
  fh.code=5; //IEEE floating-point format
  fh.ensemble_fold=cube_incr.n_cr*cube_incr.n_in;
  fh.sorting_code=1; //Horizontaly stacked
  fh.mes_sys=1; //Meters
  //fh.rev=0x100; // rev 1.0
  fh.rev=0; // rev 0
  fh.fixed_length=1;
  fh.num_ext_headers=0;
  
  fwrite(&fh,sizeof(fh),1,fo);
  
  th.field_record_n=1;
  th.coord_units=1;
  th.num_samples=cube_incr.nz;
  th.sample_interval=cube_incr.dz;
  th.src_measure_unit=1;
  for (i=0;i<cube_incr.n_in;i++)
  {
    th.in_line=cube_incr.in_0+i*cube_incr.d_in;
    if (report_progress)
      fprintf(stderr,"Writing inline %d (%d max)...\n"
		    ,th.in_line,(int)(cube_incr.in_0+(cube_incr.n_in-1)*cube_incr.d_in));
    for (j=0;j<cube_incr.n_cr;j++)
    {
      th.trace_field_record_n=th.n_in_file=th.n_in_line=i*cube_incr.n_cr+j;
      th.trace_code=1;
      th.cross_line=cube_incr.cr_0+j*cube_incr.d_cr;
      incr2XY((float)th.in_line,(float)th.cross_line,cube_incr.transform,&x,&y);
      th.src_X=(int)x;
      th.src_Y=(int)y;
      fwrite(&th,sizeof(th),1,fo);
      for (k=0;k<cube_incr.nz;k++)
	swapwrite4(&(cube_incr.data[k][j][i]),fo);
    }
  }
  fclose(fo);
}

char *read_SEGY_text_header(FILE *f)
{
  char *buf = new char[SEGY_FILE_TEXT_HEADER_SIZE];
  if (!buf)
  {
    fprintf(stderr,"Error allocating memory in read_SEGY_text_header\n");
    exit(EXIT_FAILURE);
  }
  if (fread(buf,SEGY_FILE_TEXT_HEADER_SIZE,1,f) != 1) return NULL;
  else return buf;
}

SEGY_file_header *read_SEGY_file_header(FILE *f, bool swap)
{
  SEGY_file_header *h = new SEGY_file_header;
  if (!h)
  {
    fprintf(stderr,"Error allocating memory in read_SEGY_text_header\n");
    exit(EXIT_FAILURE);
  }
  if (fread(h,sizeof(SEGY_file_header),1,f) != 1) return NULL;
  if (swap) SEGY_swap_file_header(*h);
  return h;
}

SEGY_trace_header *read_SEGY_trace_header(FILE *f, SEGY_trace_header &th, bool swap)
{
  if (fread(&th,sizeof(SEGY_trace_header),1,f) != 1) return NULL;
  
  if (swap) SEGY_swap_trace_header(th);
  return &th;
}

float *read_SEGY_trace_data(FILE *f, SEGY_trace_header &th, float *buf, bool swap_data)
{
  if (fread(buf,sizeof(float),th.num_samples,f) != th.num_samples)
    return NULL;
  if (swap_data) for (int i=0;i<th.num_samples;i++) swap4(buf+i);
  return buf;
}

float *read_SEGY_trace(FILE *f, SEGY_trace_header &th, float *buf, bool swap, bool swap_data)
{
  if (!read_SEGY_trace_header(f,th,swap)) return NULL;
  return read_SEGY_trace_data(f,th,buf,swap_data);
}

/*
void print_SEGY_header(FILE *fo, SEGY_file_header &h)
{
  fprintf(fo,"SEGY file binary header:\n");
  fprintf(fo,"JOB_ID = %d\n",h.job_id);
  fprintf(fo,"line = %d\n",h.line);
  fprintf(fo,"reel_n=%d",h.reel_n);
  fprintf(fo,"Ntr_pa=%sd",h.Ntr_pa);
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
*/
