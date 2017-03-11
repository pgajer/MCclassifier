
/*****************************************************************
 *
 * sre_malloc() and sre_realloc() and the corresponding macros
 * were lifted from
 *
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

#if !defined(CUTILITIES_H)
#define CUTILITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>

void * sre_malloc(char *file, int line, size_t size);
void * sre_realloc(char *file, int line, void *p, size_t size);

#define MallocOrDie(x)     sre_malloc(__FILE__, __LINE__, (x))
#define ReallocOrDie(x,y)  sre_realloc(__FILE__, __LINE__, (x), (y))

#endif
