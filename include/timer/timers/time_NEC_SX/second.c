#include <time.h>
   void second_(double *time) {
   clock_t zeit;
   zeit = clock();
   *time = (double)zeit/CLOCKS_PER_SEC;
    }
