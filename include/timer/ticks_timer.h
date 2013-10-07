/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   CPU cycles timers
 */

#ifndef _TICKS_TIMER_H_
#define _TICKS_TIMER_H_
#pragma once 

#define HAVE_SYS_TIME_H 1
#undef TIME_WITH_SYS_TIME
#include <timer/timers/cycle.h>

typedef ticks timeval_t;

static FORCE_INLINE
timeval_t getwalltime() {
  return getticks();
}

#endif
