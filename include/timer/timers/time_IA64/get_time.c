#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#define BODY \
        struct timeval tv;\
        struct tms buf;\
\
        if(gettimeofday(&tv,NULL)!=0) {\
                perror("get_time: ");\
                exit(-1);\
        }\
        *real = (double)tv.tv_sec + ((double)tv.tv_usec/1000000.0);\
\
        times(&buf);\
        *user = (double)buf.tms_utime/(double)CLK_TCK ;\
        *sys  = (double)buf.tms_stime/(double)CLK_TCK;\


void get_time_(double *user, double *sys, double *real)
{
        BODY
}

void GET_TIME_(double *user, double *sys, double *real)
{
        BODY
}

void get_time__(double *user, double *sys, double *real)
{
        BODY
}

void get_time(double *user, double *sys, double *real)
{
        BODY
}

