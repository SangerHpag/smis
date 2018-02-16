/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2013, 2014  Genome Research Ltd.                          *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of smis_pipeline.                                     *
 *                                                                          *
 *  ssaha_pileup is free software: you can redistribute it and/or modify it *
 *  under the terms of the GNU General Public License as published by the   *
 *  Free Software Foundation, either version 3 of the License, or (at your  *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/
/****************************************************************************/

#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <sys/wait.h>
#include <sys/signal.h>
#include <errno.h>
#include "smissv.h"


static B64_long sBase;


/* SSAS default parameters   */


static int num_reads=0;
static int n_contigs = 0;
static int read_len = 2000;
static int step_len = 150;
static int min_len = 3000;
static float num_cover = 12.0;
static char ctgname[250];

void Align_Process_sh(char **argv,int args,int nRead);

fasta *seq;


int main(int argc, char **argv)
{
  int i,args;
  char line[2000]={0};
  fasta *seq; 

  seq=NULL;

  strcpy(ctgname,"fakemates");

  if(argc < 2)
    {
      printf(" %s <-rlength 1000> <-step 150> <-minlen 3000> <input fasta> <output fastq>\n",argv[0]);
      exit(1);
    }
     
  args=1;
  for(i=1;i<argc;i++)
    {
      if(!strcmp(argv[i],"-rlength")){
	sscanf(argv[++i],"%d",&read_len);
	args=args+2;
      }else if(!strcmp(argv[i],"-step")){
	sscanf(argv[++i],"%d",&step_len);
	args=args+2;
      }else if(!strcmp(argv[i],"-cover")){
	sscanf(argv[++i],"%f",&num_cover);
	args=args+2;
      }else if(!strcmp(argv[i],"-minlen")){
	sscanf(argv[++i],"%d",&min_len);
	args=args+2;
      }
    }
 
  num_reads = 0;
  Align_Process_sh(argv,args,num_reads);
  return EXIT_SUCCESS;
  
}
/* end of the main */

/*   Subroutine to process alignment information */
/* ====================================================  */
void Align_Process_sh(char **argv,int args,int nRead)
/* ====================================================  */
{
     int i,j,k,rc,nSeq = nRead;
     char **DBname,*ptr,RC,nametag1[250],nametag2[250];
     char *line,*st,*ed;
     FILE *fp,*namef,*namef2,*fpOutfast,*fpOutfast2;
     B64_long *read_offsets,big_num,sum_bases;
     int n_patch,idd,stopflag,num_gcs,*contig_index;
     float rate;
     fasta *segg,*seqp;
     B64_long Size_q_pdata;
     int num_seqque;
     char *pdata;

     if((fp=fopen(argv[args],"rb"))==NULL) printf("Cannot open file to shred\n");
       fseek(fp, 0, SEEK_END);
     Size_q_pdata = ftell(fp) + 1;
     fclose(fp);
     if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
       printf("calloc pdata\n");
     num_seqque = extractFastq(argv[args],pdata,Size_q_pdata);
     if((segg=(fasta*)calloc((num_seqque),sizeof(fasta)))==NULL)
       printf("calloc segg\n");
     if((seq=decodeFastq(argv[args],&num_seqque,&sBase,pdata,Size_q_pdata,segg))==NULL)
       printf("no query data found.\n");
     n_contigs = num_seqque;
     fastaUC(seq,n_contigs);

     nSeq = n_contigs;
     nRead = n_contigs; 
     pri=0;

     if((namef = fopen(argv[args+1],"w")) == NULL){
       printf("ERROR main:: reads group file \n");
       exit(1);
     }
     if((namef2 = fopen(argv[args+2],"w")) == NULL){
         printf("ERROR main:: reads group file \n");
         exit(1);
     }
     //printf("name1: %s\n",argv[args]);


     int survived=0;
     int pp=0;
     /*   output all the fastq files    */

     int n_contig = 0;
     for(i=0;i<nSeq;i++)  //nSeq
     {

       int orig=0;
       pri=0;
      
       char *dpp;
       int gen_len,num_reads,rdlen,rdidt,insert_size;

       seqp=seq+i;
       gen_len = seqp->length;	
       num_reads=num_cover*gen_len/read_len;
       seqp->finished = 1;
       st = seqp->name;
       ed = strchr(seqp->name,'_');
       memset(nametag1,'\0',250);
       memset(nametag2,'\0',250);

       if(orig){
	 strncpy(nametag1,seqp->name,ed-st);
	 ed = strrchr(seqp->name,'_');
	 rdlen = atoi(ed+1);
	 st = strchr(seqp->name,'_');
	 strncpy(nametag2,st+1,ed-st-1);
	 rdidt = atoi(nametag2);

	 
       }else{
	 strncpy(nametag1,ctgname,10);

	 //if(i<593 && i>590)pri=1;
	 //else if(i<400 && i>308)pri=1;
	 //else pri=0;


	 if(pri) std::cout << i << " "<< rdlen << " name1 " << nametag1 << " name2  "
                        << nametag2 << " ed-st " << ed-st << std::endl;
	 
	 sprintf(ed,"_%07d",gen_len);
	 
	 rdlen = gen_len; 
	 sprintf(st,"_%08d_%07d",n_contig,gen_len);
	 strncpy(nametag2,st+1,ed-st-1);
	 rdidt = n_contig;
	
 	//if(nametag2=="00000573_0006781")
	//std::cout << i << " "<< rdlen << " name1 " << nametag1 << " name2  "
          //              << nametag2 << std::endl;
 
      }
       
       insert_size = rdlen; 
       
       if(no)std::cout << insert_size << " "<< rdlen << " name1 " << nametag1 << " name2  "
			<< rdidt << " st " << st << std::endl;
       
      
       if(rdlen < min_len) continue;  

       if(gen_len > (2*read_len))  // len >2*set
	 {	  
	    if(seqp->finished)
	      {
		k = 0;
		j = gen_len;
			   
		for(insert_size=rdlen;insert_size>(2*read_len);)
		  {
		    //if(pp<7915 && pp> 7908) pri=1;
		    //else pri=0;
 
		    if(pri)std::cout << i << " " << nametag1
			//insert_size << " "<< rdlen << " ed " << ed << " name  "
			//      << seqp->length 
			<< std::endl;
		    fprintf(namef,">%s:%06d_%08d_%07d/1\n",nametag1,j,rdidt,insert_size);
		    if(pri)printf("file 1 >%s:%06d_%08d_%07d/1\n",nametag1,j,rdidt,insert_size);
		   
		    for(rc=0;rc<read_len;rc++){
		      fprintf(namef,"%c",seqp->data[rc+k]);
		      if(no)if(rc<10)printf("%c ",seqp->data[rc+k]);
		    }
		    fprintf(namef,"\n");
		    if(no)printf("\n");


		    k = k+step_len;
		    fprintf(namef2,">%s:%06d_%08d_%07d/2\n",nametag1,j,rdidt,insert_size);
		    if(pri)printf("file 2 >%s:%06d_%08d_%07d/2\n",nametag1,j,rdidt,insert_size);

		    dpp = seqp->data + j-1;
		    for(rc=0;rc<read_len;rc++)
		      {
			if(*dpp == 'A') fprintf(namef2,"%c",'T');
			else if(*dpp == 'C') fprintf(namef2,"%c",'G');
			else if(*dpp == 'G') fprintf(namef2,"%c",'C');
			else if(*dpp == 'T') fprintf(namef2,"%c",'A');
			else                 fprintf(namef2,"%c",*dpp);
			dpp--;
		      }
		    fprintf(namef2,"\n");
		    j = j-step_len;
		    insert_size = insert_size - step_len - step_len;

		    pp++;
		  }
		insert_size = rdlen;
		
	      }else{ // finished==0
           
	      fprintf(namef,"@%s/1\n",seqp->name);
	      for(rc=0;rc<read_len;rc++)
                fprintf(namef,"%c",seqp->data[rc]);
	      fprintf(namef,"\n");
	      
	      fprintf(namef,"+\n");
	      putc(0+041,namef);
	      for(rc=1;rc<read_len;rc++)
                putc(seqp->qual[rc]+041,namef);
	      
	      fprintf(namef,"\n");
	      fprintf(namef2,"@%s/2\n",seqp->name);
	      dpp = seqp->data + gen_len-1;
	      for(rc=0;rc<read_len;rc++)
		{
		  if(*dpp == 'A') fprintf(namef2,"%c",'T');
		  else if(*dpp == 'C') fprintf(namef2,"%c",'G');
		  else if(*dpp == 'G') fprintf(namef2,"%c",'C');
		  else if(*dpp == 'T') fprintf(namef2,"%c",'A');
		  else                 fprintf(namef2,"%c",*dpp);
		  dpp--;
		}
	      fprintf(namef2,"\n");
	      fprintf(namef2,"+\n");
	      putc(0+041,namef2);
	      for(rc=1;rc<read_len;rc++)
                putc(seqp->qual[read_len-rc]+041,namef2);
	      fprintf(namef2,"\n");
	    }
        }else if((gen_len > read_len)&&(gen_len > 1500)) {  // len >set && len>1500
          int read_len2 = 1500;
          if(seqp->finished)
          {
             j = 0; 
             if(pri)printf(">%s:%06d_%08d_%07d/1\n",nametag1,j,rdidt,insert_size);
             fprintf(namef,">%s:%06d_%08d_%07d/1\n",nametag1,j,rdidt,insert_size);
             for(rc=0;rc<read_len2;rc++)
                fprintf(namef,"%c",seqp->data[rc]);
             fprintf(namef,"\n");
             fprintf(namef2,">%s:%06d_%08d_%07d/2\n",nametag1,j,rdidt,insert_size);
             dpp = seqp->data + gen_len-1;
             for(rc=0;rc<read_len2;rc++)
             {
                if(*dpp == 'A') fprintf(namef2,"%c",'T');
                else if(*dpp == 'C') fprintf(namef2,"%c",'G');
                else if(*dpp == 'G') fprintf(namef2,"%c",'C');
                else if(*dpp == 'T') fprintf(namef2,"%c",'A');
                else                 fprintf(namef2,"%c",*dpp);
                dpp--;
             }
             fprintf(namef2,"\n");
          }else { // len < set or len < 1500
             j = 0; 
             fprintf(namef,"@%s:%06d%s/1\n",nametag1,j,ed);
             for(rc=0;rc<read_len2;rc++)
                fprintf(namef,"%c",seqp->data[rc]);
             fprintf(namef,"\n");

             fprintf(namef,"+\n");
             putc(0+041,namef);
             for(rc=1;rc<read_len2;rc++)
                putc(seqp->qual[rc]+041,namef);
        
             fprintf(namef,"\n");
             fprintf(namef2,"@%s:%06d%s/2\n",nametag1,j,ed);
             dpp = seqp->data + gen_len-1;
             for(rc=0;rc<read_len2;rc++)
             {
                if(*dpp == 'A') fprintf(namef2,"%c",'T');
                else if(*dpp == 'C') fprintf(namef2,"%c",'G');
                else if(*dpp == 'G') fprintf(namef2,"%c",'C');
                else if(*dpp == 'T') fprintf(namef2,"%c",'A');
                else                 fprintf(namef2,"%c",*dpp);
                dpp--;
             }
             fprintf(namef2,"\n");
             fprintf(namef2,"+\n");
             putc(0+041,namef2);
             for(rc=1;rc<read_len2;rc++)
                putc(seqp->qual[gen_len-rc]+041,namef2);
             fprintf(namef2,"\n");
          }
        }
     

       n_contig++;
     }// end loop
     fclose(namef);
     fclose(namef2);

}

