#ifndef BLAS_EXTENDED_H
#define BLAS_EXTENDED_H

#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

#include "blas_enum.h"
#include "blas_malloc.h"

/* constants */

#define BITS_S  24
#define BITS_D  53
#define BITS_E  106

/* Split a double into 2 parts with at most 26 bits each. (2^27 + 1) */
#define split 	(134217729.0)


/* macros */

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))


/* prototypes */

#include "blas_extended_proto.h"
#include "blas_dense_proto.h"

#endif /* BLAS_EXTENDED_H */
