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

const int ntest = 12;

int process(int n_first, int n_last, int n_step, int L);

int main(int argc, char** argv)
{
  int i;
  int n_first = 32; 
  int n_last = 1024;
  int n_step = 128;
  int L = 2;
 
  if (argc != 1) {
    for(i = 1; i<argc; i++){
      if (strcmp("-l", argv[i])==0)
	n_last  = atoi(argv[++i]);
      else if(strcmp("-f", argv[i])==0)
	n_first = atoi(argv[++i]);
      else if(strcmp("-s", argv[i])==0)
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
  
  return process(n_first, n_last, n_step, L);
}
