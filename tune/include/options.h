/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenko@ucdavis.edu
 **********************************************************************
 * Description:
 *  Automatic options processing in benchmarking executables.
 */

#include <stdio.h>

FILE * pFile;
#define MESSAGE(...)				\
  fprintf(pFile, __VA_ARGS__)

const int ntest = 12;

int process(int n_first, int n_last, int n_step, int L);

int main(int argc, char** argv)
{
  int i;
  int bUseFile = 0;
  int ret_val;
  int n_first = 32; 
  int n_last = 1024;
  int n_step = 128;
  int L = 2;
  char* strNlast=0;

  pFile = stdout; 
  if (argc != 1) {
    for(i = 1; i<argc; i++){
      if (strcmp("-n", argv[i])==0){
	strNlast = argv[++i];
	while(*strNlast != ':' && *strNlast != '\0') strNlast++;
	if(*strNlast == ':') {
	  *strNlast = '\0'; n_last = atoi(++strNlast);
	}
	n_first = atoi(argv[i]);
      }else if(strcmp("-f", argv[i])==0) {
	pFile = fopen (argv[++i],"w");
	bUseFile = 1;
      } else if(strcmp("-s", argv[i])==0)
	n_step  = atoi(argv[++i]);
      else if(strcmp("-L", argv[i])==0)
	L = atoi(argv[++i]);
    }
    /* { */
    /*   printf("\nUsage: \n"); */
    /*   printf("  <prog> -n %d -L %d\n\n", 1024, 10); */
    /*   exit(1); */
    /* } */
  }

  /* MESSAGE("\n!!%4d:%4d:%3d \n", n_first, n_last, n_step); */
  ret_val = process(n_first, n_last, n_step, L);
  if(bUseFile) fclose(pFile);
  return ret_val;
}
