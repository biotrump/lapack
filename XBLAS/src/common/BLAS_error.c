#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "blas_extended.h"

void BLAS_error(const char *rname, int iflag, int ival, char *form, ...)
/*
 * Argument
 * ========
 * rname     (input) routine name
 *
 * iflag     (input) a negative value indicates that parameter number -IFLAG
 *                   caused the error; a nonnegative value is an
 *                   implementation-specific error code.
 *
 * ival      (input) the value of parameter number -IFLAG.
 */
{
  va_list argptr;

  va_start(argptr, form);
  fprintf(stderr, "Error #%d from routine %s:\n", iflag, rname);
  if (form)
    vfprintf(stderr, form, argptr);
  else if (iflag < 0)
    fprintf(stderr,
            "  Parameter number %d to routine %s had the illegal value %d\n",
            -iflag, rname, ival);
  else
    fprintf(stderr, "  Unknown error code %d from routine %s\n",
            iflag, rname);
  exit(iflag);
}
