/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   Unified interface for timers
 */

#ifndef _TIMER_H_
#define _TIMER_H_
#pragma once 

#define USE_SYS_TIMER

#if   defined(USE_SYS_TIMER)
#  include <timer/sys_timer.h>
#elif defined(USE_TICKS_TIMER)
#  include <timer/ticks_timer.h>
#endif

#endif
