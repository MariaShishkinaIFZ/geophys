/*Reads and prints headers and the 2 traces from SEG-Y file*/
#include <cstdio>
#include <cmath>

#include "subio.h"
#include "submem.h"
#include "segyio.h"

#define NL_MAX 20
#define MAX_FNL 256
#define Z0_REF 500

/*
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
*/

int main(int narg, char **argv)
{
  FILE *f = file_open(argv[1],"rb");
  SEGY_file_header *fh = read_SEGY_file_header(f,1);
  SEGY_trace_header *th = new SEGY_trace_header;
  read_SEGY_trace_header(f,*th,1);
  fprintf(stdout,"\n\nnum_samples=%hd\n",th->num_samples);
  fprintf(stdout,"sample_interval=%hd\n",th->sample_interval);
  fprintf(stdout,"n_in_file=%d\n",th->n_in_file);
  fprintf(stdout,"n_in_line=%d\n",th->n_in_line);
  fprintf(stdout,"field_record_n=%d\n",th->field_record_n);
  fprintf(stdout,"trace_field_record_n=%d\n",th->trace_field_record_n);
  fprintf(stdout,"esp_n=%d\n",th->esp_n);
  fprintf(stdout,"ensemble_n=%d\n",th->ensemble_n);
  fprintf(stdout,"trace_in_ensemple_n=%d\n",th->trace_in_ensemple_n);
  fprintf(stdout,"Src_X=%d\n",th->src_X);
  fprintf(stdout,"Src_Y=%d\n",th->src_Y);
  fprintf(stdout,"Grp_X=%d\n",th->grp_X);
  fprintf(stdout,"Grp_Y=%d\n",th->grp_Y);
  fprintf(stdout,"CDP_X=%d\n",th->CDP_X);
  fprintf(stdout,"CDP_Y=%d\n",th->CDP_Y);
  fprintf(stdout,"inline=%d\n",th->in_line);
  fprintf(stdout,"crossline=%d\n",th->cross_line);
  float *buf = new float[th->num_samples];
  read_SEGY_trace_data(f,*th,buf,1);
  for (int i=0;i<th->num_samples;i++)
  {
    //swap4(buf+i);
    fprintf(stdout,"%d %g\n",i,buf[i]);
  }
  read_SEGY_trace(f,*th,buf,1,1);
  fprintf(stdout,"\n\nnum_samples=%hd\n",th->num_samples);
  fprintf(stdout,"sample_interval=%hd\n",th->sample_interval);
  fprintf(stdout,"n_in_file=%d\n",th->n_in_file);
  fprintf(stdout,"n_in_line=%d\n",th->n_in_line);
  fprintf(stdout,"field_record_n=%d\n",th->field_record_n);
  fprintf(stdout,"trace_field_record_n=%d\n",th->trace_field_record_n);
  fprintf(stdout,"esp_n=%d\n",th->esp_n);
  fprintf(stdout,"ensemble_n=%d\n",th->ensemble_n);
  fprintf(stdout,"trace_in_ensemple_n=%d\n",th->trace_in_ensemple_n);
  fprintf(stdout,"Src_X=%d\n",th->src_X);
  fprintf(stdout,"Src_Y=%d\n",th->src_Y);
  fprintf(stdout,"Grp_X=%d\n",th->grp_X);
  fprintf(stdout,"Grp_Y=%d\n",th->grp_Y);
  fprintf(stdout,"CDP_X=%d\n",th->CDP_X);
  fprintf(stdout,"CDP_Y=%d\n",th->CDP_Y);
  fprintf(stdout,"inline=%d\n",th->in_line);
  fprintf(stdout,"crossline=%d\n",th->cross_line);
  for (int i=0;i<th->num_samples;i++)
  {
    //swap4(buf+i);
    fprintf(stdout,"%d %g\n",i,buf[i]);
  }
  delete[] buf;
  delete th;
  delete fh;
  fclose(f);
}


