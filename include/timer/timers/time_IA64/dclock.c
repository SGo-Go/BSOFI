#include <sys/time.h>
#include <sys/resource.h>
double dclock_() /* change to dclock() for ifc */
{
/* This structure should be portable */
    struct rusage ru;
    double t;
    getrusage(RUSAGE_SELF, &ru);
    t = (double)ru.ru_utime.tv_sec +
        (double)ru.ru_stime.tv_sec;
    t += ( (double)ru.ru_utime.tv_usec +
           (double)ru.ru_stime.tv_usec ) *
          1.0e-6;
    return t;
}
