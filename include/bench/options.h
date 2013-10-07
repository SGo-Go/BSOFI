#include <stdio.h>

#define MAX_TEST 20
int sizen[MAX_TEST];
const int ntest = 12;

int process(int threads, int tests, int n, int L);

int main( int argc, char** argv)
{
  int n = 0, L = 10, i;
  int nthreads = 0, tests = 0;

  int k = 8;
  for(i = 0; i < MAX_TEST; i++, k+=2)
    sizen[i] = k*k;
 
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
    }

    if      (n>0 && !tests ) {sizen[0] = n; tests = 1; }
    else if (n>0 && tests>0) {
      k = n; 
      for(i = tests; i > 0; i--, k-=2){
	sizen[i-1] = k*k;
	//DBG(k <<i );
      }
    }
    else if (!n  && !tests ) {tests = ntest;}
    else if (!n  && tests>0) {;}
    else {
      printf("\nUsage: \n");
      printf("  <prog> -n %d -L %d\n\n", 1024, 10);
      exit(1);
    }
  }
  else {
    L = 10;
    tests = ntest;
  }

  //DBG(nthreads << tests << n << L);
  process(nthreads, tests, n, L);
}
