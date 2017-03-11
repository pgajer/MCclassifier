/*

 CUtilities.c

 These are memory allocation routines lifted from

 SQUID - a library of functions for biological sequence analysis
 Copyright (C) 1992-2002 Washington University School of Medicine

      This source code is freely distributed under the terms of the
      GNU General Public License. See the files COPYRIGHT and LICENSE
      for details.
*/

#include <stdlib.h>
#include <stdio.h>

void * sre_malloc(char *file, int line, size_t size)
{
  void *ptr;

  if ((ptr = malloc (size)) == NULL)
    printf("malloc of %ld bytes failed: file %s line %d", size, file, line);
  return ptr;
}

void * sre_realloc(char *file, int line, void *p, size_t size)
{
  void *ptr;

  if ((ptr = realloc(p, size)) == NULL)
    printf("realloc of %ld bytes failed: file %s line %d", size, file, line);
  return ptr;
}
