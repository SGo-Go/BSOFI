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

#define MSGFLUSH(...)				\
  fflush(pFile)

int process(int threads, int tests, int n, int L);

int main(int argc, char** argv)
{
  int bUseFile = 0, ret_val = 0;
  int n = 0, L = 10, i;
  int nthreads = 0, tests = 0;

  pFile = stdout; 
 
  if (argc != 1) {
    for(i = 1; i<argc; i++){
      if (strcmp("-n", argv[i])==0)
	n = atoi(argv[++i]);
      else if(strcmp("-L", argv[i])==0)
	L = atoi(argv[++i]);
      else if(strcmp("-t", argv[i])==0)
	tests = atoi(argv[++i]);
      else if(strcmp("-p", argv[i])==0)
	nthreads = atoi(argv[++i]);
      else if(strcmp("-f", argv[i])==0) {
	pFile = fopen (argv[++i],"w");
	bUseFile = 1;
      }
    } 

    if      (n>0 && !tests ) {tests = 1; }
    else if (n>0 && tests>0) {;}
    else if (!n  && !tests ) {tests = 1;}
    else if (!n  && tests>0) {;}
    else {
      printf("\nUsage: \n");
      printf("  <prog> -n %d -L %d\n\n", 1024, 10);
      exit(1);
    }
  }
  else {
    L = 10; tests = 1; n = 2;
  }

  if(n == 0)  n = 256;
  if(L == 0)  L = 4;
  if(nthreads == 0) nthreads = 4; 
  if(tests    == 0) tests = 2;

  ret_val = process(nthreads, tests, n, L);
  if(bUseFile) fclose(pFile);
  return ret_val;
}
