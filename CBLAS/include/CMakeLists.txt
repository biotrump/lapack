SET (CBLAS_INCLUDE cblas.h cblas_f77.h cblas_test.h cblas_mangling.h cblas_mangling_with_flags.h)

#if(ANDROID_NDK)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cblas_mangling_with_flags.h ${CMAKE_CURRENT_SOURCE_DIR}/cblas_mangling.h COPYONLY)
#endif(ANDROID_NDK)

file(COPY ${CBLAS_INCLUDE} DESTINATION ${LAPACK_BINARY_DIR}/include)