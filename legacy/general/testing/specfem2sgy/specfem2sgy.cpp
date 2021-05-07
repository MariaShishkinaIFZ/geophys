
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "nrutil.h"
#include "subio.h"
#include "submem.h"
#include "segyio.h"

#define FLNLN 256 	/*Maximum length of the file name*/
#define STRLN 256 	/*Maximum length of comments string*/

int main(int narg, char **argv)

{
  SEGY_file_header fh;			// The SEG-Y file header structure
  SEGY_trace_header th;			// The SEG-Y trace header structure
  FILE *fp, *fo, *fl, *ft;		// Files
  char *trace_name, *reg_type;		// The name of the string file
  float time, amp;			// Trace amplitude
  char buf[STRLN];
  int i, n_in_line, src_X, src_Y, grp_X;
  
  
 // fl=file_open(argv[0],"rt");	
  
  // The main function has 4 arguments. argv[0]=[name_of_programm] argv[1]=[SEG-Y_header_file_name] argv[2]=[trace_header_file_name] argv[3]=[Output_SEG-Y_file_name]
  
  if (narg<4)
  {
    fprintf(stderr,"Insufficient parameters.\n Use: specfem2sgy <SEG-Y_header_file_name> <Trace_header_file_name>, <Output_SEG-Y_file_name>\n");
    exit(EXIT_FAILURE);
  }

  
  fp=fopen(argv[1],"r"); 

// reading FILE HEADER parametrs
   
   
  fgets(buf,STRLN,fp);										/*1st line - comment */
  fscanf(fp,"%d",&(fh.job_id)); 			fgets(buf,STRLN,fp); 			/*2nd line - job_id */
  fscanf(fp,"%d",&(fh.line)); 				fgets(buf,STRLN,fp); 			/*3rd line - line_id */
  fscanf(fp,"%d",&(fh.reel_n)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.Ntr_pa)); 			fgets(buf,STRLN,fp); /* */
  fscanf(fp,"%d",&(fh.reel_n)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.Ntr_pa)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.Nat_pa)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.Ntr_pa)); 			fgets(buf,STRLN,fp); /* */
  fscanf(fp,"%hu",&(fh.sample_interval)); 		fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.original_sample)); 		fgets(buf,STRLN,fp); /* */
  fscanf(fp,"%hu",&(fh.ns)); 				fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.ns_original)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.code)); 				fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.ensemble_fold)); 		fgets(buf,STRLN,fp); /* */
  fscanf(fp,"%hu",&(fh.sorting_code)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.vertical_sum_code)); 		fgets(buf,STRLN,fp); /* */
  fscanf(fp,"%hu",&(fh.sw__freq_0)); 			fgets(buf,STRLN,fp); /* */
  fscanf(fp,"%hu",&(fh.sw_freq_1)); 			fgets(buf,STRLN,fp); /* */
  fscanf(fp,"%hu",&(fh.sw_l)); 				fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.sw_type)); 			fgets(buf,STRLN,fp); /* */
  fscanf(fp,"%hu",&(fh.sw_ch)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.sw_taper0)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.sw_taper_l)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.taper_type)); 			fgets(buf,STRLN,fp); /* */
  fscanf(fp,"%hu",&(fh.correlated)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.gain_recov)); 			fgets(buf,STRLN,fp); /* */
  fscanf(fp,"%hu",&(fh.ampl_recov)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.mes_sys)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.imp_polarity)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%c",&(fh.unassigned1[240])); 		fgets(buf,STRLN,fp); /* */
  fscanf(fp,"%hu",&(fh.rev)); 				fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.fixed_length)); 			fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%hu",&(fh.num_ext_headers)); 		fgets(buf,STRLN,fp); /* */ 
  fscanf(fp,"%c",&(fh.unassigned2[94])); 		fgets(buf,STRLN,fp); /* */
  
  fclose(fp);							//FILE HEADER closed. 
 
  fo=fopen(argv[3],"wb"); 					// We open SEG-Y file as an output file
  
  fwrite(&fh,sizeof(fh),1,fo);			// arguments: 1st - &fh - all the file header parrametres assighed above; 2nd - their size; 3rd - 1 ; 4th - The output file
   
  fl=fopen(argv[2],"rt"); 			// We open Trace header file and need to write all the necessary headers' data in SEG-Y
  
  fgets(buf,STRLN,fp);	
  
    while (!feof(fl))
      
  {
  
  fscanf(fp,"%s%s%d%d%d%d",reg_type, trace_name,&n_in_line,&src_X,&src_Y,&grp_X); 	fgets(buf,STRLN,fp); 	
  
  
  th.n_in_line=n_in_line;
  th.field_record_n=1;
  th.coord_units=1;
  th.num_samples=fh.ns;
  th.sample_interval=fh.sample_interval;
  th.src_measure_unit=1;
  
  
  fwrite(&th,sizeof(th),1,fo);								// We write all the TRACE HEADER data in output file
  
  ft=fopen(trace_name,"rt");								// We open each file with traces
  
  for (i=0; i< fh.ns; i++)
    
  {
   fscanf (ft,"%f%f", &time, &amp); 
   fwrite(&amp,sizeof(amp),1,fo);							// We write all the trace data in SEG-Y
   
  }
  
  fclose(ft);					//  We close each file with traces
  
  fclose(fo);					// We close SEG-Y file . TRACE data is written.
 
  }
  

     

     
   
  
  
  
  
}

  