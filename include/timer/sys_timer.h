/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   System timers (only for *nix)
 */

#ifndef _SYS_TIMER_H_
#define _SYS_TIMER_H_
#pragma once 

#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
//#include <stdio.h>
#include <errno.h>

#include <config.h>

typedef struct timeval timeval_t;

static FORCE_INLINE
timeval_t getwalltime() {
  struct timeval tv;
  /* struct tms buf; */
  if(gettimeofday(&tv,NULL)!=0) {
    //perror("get_time: ");
    exit(-1);
  }
  return tv;
  /* times(&buf); */
  /* user = (double)buf.tms_utime/(double)CLK_TCK; */
  /* sys  = (double)buf.tms_stime/(double)CLK_TCK; */
}

static FORCE_INLINE
double elapsed(timeval_t t1, timeval_t t0)
{
  return ((double)(t1.tv_sec - t0.tv_sec) + 
	  ((double)(t1.tv_usec - t0.tv_usec)/1000000.0));
}

#endif
