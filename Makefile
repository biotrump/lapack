#
#  Top Level Makefile for LAPACK
#  Version 3.0
#  June 30, 1999
#

include make.inc

all: lapack_install lib lapack_testing blas_testing

lib: lapacklib tmglib
#lib: blaslib lapacklib tmglib

clean: cleanlib cleantesting cleanblas_testing cleantiming

lapack_install:
	( cd INSTALL; $(MAKE); ./testlsame; ./testslamch; \
	  ./testdlamch; ./testsecond; ./testdsecnd; ./testversion \
	  cp lsame.f ../BLAS/SRC/; cp lsame.f ../SRC; \
	  cp slamch.f ../SRC/; cp dlamch.f ../SRC/; \
	  cp second.f ../SRC/; cp dsecnd.f ../SRC/; \
	  cp ilaver.f ../SRC )

blaslib:
	( cd BLAS/SRC; $(MAKE) )

xblaslib:
	( cd XBLAS;  $(MAKE) lib )

lapacklib:	lapack_install
	( cd SRC; $(MAKE) )

tmglib:
	( cd TESTING/MATGEN; $(MAKE) )

lapack_testing:	lib
	( cd TESTING ; $(MAKE) )

blas_testing:
	( cd BLAS/TESTING; $(MAKE) -f Makeblat1 )
	( cd BLAS; ./xblat1s > sblat1.out    ; \
	           ./xblat1d > dblat1.out    ; \
	           ./xblat1c > cblat1.out    ; \
	           ./xblat1z > zblat1.out    ) 
	( cd BLAS/TESTING; $(MAKE) -f Makeblat2 )
	( cd BLAS; ./xblat2s < sblat2.in     ; \
	           mv SBLAT2.SUMM sblat2.out ; \
	           ./xblat2d < dblat2.in     ; \
	           mv DBLAT2.SUMM dblat2.out ; \
	           ./xblat2c < cblat2.in     ; \
	           mv CBLAT2.SUMM cblat2.out ; \
	           ./xblat2z < zblat2.in     ; \
	           mv ZBLAT2.SUMM zblat2.out )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat3 )
	( cd BLAS; ./xblat3s < sblat3.in     ; \
	           mv SBLAT3.SUMM sblat3.out ; \
	           ./xblat3d < dblat3.in     ; \
	           mv DBLAT3.SUMM dblat3.out ; \
	           ./xblat3c < cblat3.in     ; \
	           mv CBLAT3.SUMM cblat3.out ; \
	           ./xblat3z < zblat3.in     ; \
	           mv ZBLAT3.SUMM zblat3.out )

xblas_testing:
	( cd XBLAS;  $(MAKE) tests )

lapack_timing:	lib lapack_testing blas_testing
	( cd TIMING; $(MAKE) )

blas_timing:	lapack_timing
	( cd TIMING/LIN; $(MAKE) )
	( cd TIMING; ./xlintims < sblasa_small.in > sblasa_small.out ; \
	             ./xlintims < sblasb_small.in > sblasb_small.out ; \
	             ./xlintims < sblasc_small.in > sblasc_small.out )
	( cd TIMING; ./xlintimd < dblasa_small.in > dblasa_small.out ; \
	             ./xlintimd < dblasb_small.in > dblasb_small.out ; \
	             ./xlintimd < dblasc_small.in > dblasc_small.out )
	( cd TIMING; ./xlintimc < cblasa_small.in > cblasa_small.out ; \
	             ./xlintimc < cblasb_small.in > cblasb_small.out ; \
	             ./xlintimc < cblasc_small.in > cblasc_small.out )
	( cd TIMING; ./xlintimz < zblasa_small.in > zblasa_small.out ; \
	             ./xlintimz < zblasb_small.in > zblasb_small.out ; \
	             ./xlintimz < zblasc_small.in > zblasc_small.out )

cleanlib:
	( cd INSTALL; $(MAKE) clean )
	( cd BLAS/SRC; $(MAKE) clean )
	( cd SRC; $(MAKE) clean )
	( cd TESTING/MATGEN; $(MAKE) clean )
	( cd XBLAS; $(MAKE) clean )

cleanblas_testing:	
	( cd BLAS/TESTING; $(MAKE) -f Makeblat1 clean )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat2 clean )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat3 clean )
	( cd BLAS; rm -f *.SUMM xblat* )

cleanxblas:
	( cd XBLAS; $(MAKE) clean ) 

cleantesting:
	( cd TESTING/LIN; $(MAKE) clean )
	( cd TESTING/EIG; $(MAKE) clean )
	( cd TESTING; rm -f xlin* xeig* )

cleantiming:
	( cd TIMING/LIN; $(MAKE) clean )
	( cd TIMING/LIN/LINSRC; $(MAKE) clean )
	( cd TIMING/EIG; $(MAKE) clean )
	( cd TIMING/EIG/EIGSRC; $(MAKE) clean )
	( cd TIMING; rm -f xlin* xeig* )

cleanall: cleanlib cleanblas_testing cleantesting cleantiming
	(cd XBLAS; $(MAKE) dist-clean )
	rm -f *.a TESTING/*.out TIMING/*.out INSTALL/test* \
                 BLAS/*.out TIMING/LIN/*.a TIMING/EIG/*.a

