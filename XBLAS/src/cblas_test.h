#ifndef __CBLAS_TEST_H_
#define __CBLAS_TEST_H_

#define TEST_PROBABILITY_DEFAULT 0.001
#define MAX_BAD_TESTS 100
#define TOTAL_FAILURE_THRESHOLD 1000

extern double xrand(int *);
extern double power(int, int);
extern int    FixedBits(double, double);

extern void ddmuld(double dda_l, double dda_t, double db,
		   double *ddc_l, double *ddc_t);
extern void ddadd(double dda_l, double dda_t, double ddb_l, double ddb_t,
		  double *ddc_l, double *ddc_t);
extern void dddiv(double dda_l, double dda_t, double ddb_l, double ddb_t,
		  double *ddc_l, double *ddc_t);
extern void z_ddmuld(double *dda_l, double *dda_t, double *db,
		     double *ddc_l, double *ddc_t);
extern void z_dddivd(double *dda_l, double *dda_t, double *db,
		     double *ddc_l, double *ddc_t);

/****************************
 *   Level 1 routines       *
 ****************************/

extern void testgen_BLAS_sdot(int n, int n_fix2, int n_mix, int norm, 
			      enum blas_conj_type conj,
			      float *alpha, int alpha_flag,
			      float *beta, int beta_flag,
			      float* x, float* y, int *seed,
			      float *r, double *r_true_l,
			      double *r_true_t);

extern void testgen_BLAS_ddot(int n, int n_fix2, int n_mix, int norm,
			      enum blas_conj_type conj,
			      double *alpha, int alpha_flag,
			      double *beta, int beta_flag,
			      double *x, double *y, int *seed, 
			      double *r, double *r_true_l, 
			      double *r_true_t);

extern void testgen_BLAS_cdot(int n, int n_fix2, int n_mix, int norm,
			      enum blas_conj_type conj,
			      void *alpha, int alpha_flag,
			      void *beta, int beta_flag,
			      void *x, void *y, int *seed,
			      void *r, double r_true_l[],
			      double r_true_t[]);

extern void testgen_BLAS_zdot(int n, int n_fix2, int n_mix, int norm,
			      enum blas_conj_type conj,
			      void *alpha, int alpha_flag,
			      void *beta, int beta_flag,
			      void *x, void *y, int *seed,
			      void *r, double r_true_l[],
			      double r_true_t[]);

extern void testgen_BLAS_sdot_x(int n, int n_fix2, int n_mix, 
				int norm, enum blas_conj_type conj,
				float *alpha, int alpha_flag, 
				float *beta, int beta_flag,
				double* x_l, double *x_t, float* y,
				int *seed, float *r,
				double *r_true_l, double *r_true_t);
  

extern void BLAS_sdot_testgen(int n, int n_fix2, int n_mix, int norm,
		      enum blas_conj_type conj,
		      float* alpha, int alpha_flag,
		      float* beta,  int beta_flag,
		      float* x, float* y, int *seed, float* r,
		      double *r_true_l, double *r_true_t);
extern void BLAS_ddot_testgen(int n, int n_fix2, int n_mix, int norm,
		      enum blas_conj_type conj,
		      double* alpha, int alpha_flag,
		      double* beta,  int beta_flag,
		      double* x, double* y, int *seed, double* r,
		      double *r_true_l, double *r_true_t);
extern void BLAS_ddot_s_s_testgen(int n, int n_fix2, int n_mix,
			  int norm, enum blas_conj_type conj,
			  double* alpha, int alpha_flag,
			  double* beta,  int beta_flag,
			  float* x, float* y, int *seed, double* r,
			  double *r_true_l, double *r_true_t);
extern void BLAS_ddot_s_d_testgen(int n, int n_fix2, int n_mix,
			  int norm, enum blas_conj_type conj,
			  double* alpha, int alpha_flag,
			  double* beta,  int beta_flag,
			  float* x, double* y, int *seed, double* r,
			  double *r_true_l, double *r_true_t);
extern void BLAS_ddot_d_s_testgen(int n, int n_fix2, int n_mix,
			int norm, enum blas_conj_type conj,
			double* alpha, int alpha_flag,
			double* beta,  int beta_flag,
			double* x, float* y, int *seed, double* r,
			double *r_true_l, double *r_true_t);
extern void BLAS_cdot_testgen(int n, int n_fix2, int n_mix, int norm, 
			enum blas_conj_type conj,
			void* alpha, int alpha_flag,
			void* beta,  int beta_flag,
			void* x, void* y, int *seed, void* r,
			double *r_true_l, double *r_true_t);
extern void BLAS_zdot_testgen(int n, int n_fix2, int n_mix, int norm, 
			enum blas_conj_type conj,
			void* alpha, int alpha_flag,
			void* beta,  int beta_flag,
			void* x, void* y, int *seed, void* r,
			double *r_true_l, double *r_true_t);
extern void BLAS_zdot_c_c_testgen(int n, int n_fix2, int n_mix, int norm, 
			enum blas_conj_type conj,
			void* alpha, int alpha_flag,
			void* beta,  int beta_flag,
			void* x, void* y, int *seed, void* r,
			double *r_true_l, double *r_true_t);
extern void BLAS_zdot_c_z_testgen(int n, int n_fix2, int n_mix, int norm, 
			enum blas_conj_type conj,
			void* alpha, int alpha_flag,
			void* beta,  int beta_flag,
			void* x, void* y, int *seed, void* r,
			double *r_true_l, double *r_true_t);
extern void BLAS_zdot_z_c_testgen(int n, int n_fix2, int n_mix, int norm, 
			enum blas_conj_type conj,
			void* alpha, int alpha_flag,
			void* beta,  int beta_flag,
			void* x, void* y, int *seed, void* r,
			double *r_true_l, double *r_true_t);
extern void BLAS_cdot_s_s_testgen(int n, int n_fix2, int n_mix,
                        int norm, enum blas_conj_type conj,
			void* alpha, int alpha_flag,
			void* beta,  int beta_flag,
			float* x, float* y, int *seed, void* r,
			double *r_true_l, double *r_true_t);
extern void BLAS_cdot_s_c_testgen(int n, int n_fix2, int n_mix, int norm, 
			enum blas_conj_type conj,
			void* alpha, int alpha_flag,
			void* beta,  int beta_flag,
			float* x, void* y, int *seed, void* r,
			double *r_true_l, double *r_true_t);
extern void BLAS_cdot_c_s_testgen(int n, int n_fix2, int n_mix,
			int norm, enum blas_conj_type conj,
			void* alpha, int alpha_flag,
			void* beta,  int beta_flag,
			void* x, float* y, int *seed, void* r,
			double *r_true_l, double *r_true_t);
extern void BLAS_zdot_d_d_testgen(int n, int n_fix2, int n_mix,
			int norm, enum blas_conj_type conj,
			void* alpha, int alpha_flag,
			void* beta,  int beta_flag,
			double* x, double* y, int *seed, void* r,
			double *r_true_l, double *r_true_t);
extern void BLAS_zdot_z_d_testgen(int n, int n_fix2, int n_mix,
  		int norm, enum blas_conj_type conj,
		void* alpha, int alpha_flag,
		void* beta,  int beta_flag,
		void* x, double* y, int *seed, void* r,
		double *r_true_l, double *r_true_t);
extern void BLAS_zdot_d_z_testgen(int n, int n_fix2, int n_mix,
		  int norm, enum blas_conj_type conj,
		  void* alpha, int alpha_flag,
		  void* beta,  int beta_flag,
		  double* x, void* y, int *seed, void* r,
		  double *r_true_l, double *r_true_t);




extern void test_BLAS_sdot(int n, enum blas_conj_type conj,
              float alpha, float beta, float rin, float rout,
              double r_true_l, double r_true_t,
              float* x, int incx, float* y, int incy,
              double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_ddot(int n, enum blas_conj_type conj,
              double alpha, double beta, double rin, double rout,
              double r_true_l, double r_true_t,
              double* x, int incx, double* y, int incy,
              double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_ddot_s_s(int n, enum blas_conj_type conj,
                double alpha, double beta, double rin, double rout,
                double r_true_l, double r_true_t,
                float* x, int incx, float* y, int incy,
		double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_ddot_s_d(int n, enum blas_conj_type conj,
                double alpha, double beta, double rin, double rout,
                double r_true_l, double r_true_t,
                float* x, int incx, double* y, int incy,
		double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_ddot_d_s(int n, enum blas_conj_type conj,
                double alpha, double beta, double rin, double rout,
                double r_true_l, double r_true_t,
                double* x, int incx, float* y, int incy,
		double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_cdot(int n, enum blas_conj_type conj,
              const void* alpha, const void* beta, const void* rin, const void* rout,
              double* r_true_l, double* r_true_t,
              void* x, int incx, void* y, int incy,
              double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_zdot(int n, enum blas_conj_type conj,
              const void* alpha, const void* beta, const void* rin, const void* rout,
              double* r_true_l, double* r_true_t,
              void* x, int incx, void* y, int incy,
              double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_zdot_c_c(int n, enum blas_conj_type conj,
                const void* alpha, const void* beta, const void* rin, const void* rout,
                double* r_true_l, double* r_true_t,
                void* x, int incx, void* y, int incy,
		double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_zdot_c_z(int n, enum blas_conj_type conj,
                const void* alpha, const void* beta, const void* rin, const void* rout,
                double* r_true_l, double* r_true_t,
                void* x, int incx, void* y, int incy,
		double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_zdot_z_c(int n, enum blas_conj_type conj,
                const void* alpha, const void* beta, const void* rin, const void* rout,
                double* r_true_l, double* r_true_t,
                void* x, int incx, void* y, int incy,
		double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_cdot_s_s(int n, enum blas_conj_type conj,
                const void* alpha, const void* beta, const void* rin, const void* rout,
                double* r_true_l, double* r_true_t,
                float* x, int incx, float* y, int incy,
		double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_cdot_s_c(int n, enum blas_conj_type conj,
                const void* alpha, const void* beta, const void* rin, const void* rout,
                double* r_true_l, double* r_true_t,
                float* x, int incx, void* y, int incy,
		double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_cdot_c_s(int n, enum blas_conj_type conj,
                const void* alpha, const void* beta, const void* rin, const void* rout,
                double* r_true_l, double* r_true_t,
                void* x, int incx, float* y, int incy,
		double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_zdot_d_d(int n, enum blas_conj_type conj,
                const void* alpha, const void* beta, const void* rin, const void* rout,
                double* r_true_l, double* r_true_t,
                double* x, int incx, double* y, int incy,
		double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_zdot_z_d(int n, enum blas_conj_type conj,
                const void* alpha, const void* beta, const void* rin, const void* rout,
                double* r_true_l, double* r_true_t,
                void* x, int incx, double* y, int incy,
		double eps_int, double un_int, double *test_ratio);

extern void test_BLAS_zdot_d_z(int n, enum blas_conj_type conj,
                const void* alpha, const void* beta, const void* rin, const void* rout,
                double* r_true_l, double* r_true_t,
                double* x, int incx, void* y, int incy,
		double eps_int, double un_int, double *test_ratio);


extern void BLAS_sdot_x_testgen(int n, int n_fix2, int n_mix, int norm,
			     enum blas_conj_type conj,
			     float *alpha, int alpha_flag,
			     float *beta, int beta_flag,
			     double *x_l, double *x_t,
			     float *y, int *seed, float *r,
			     double *r_true_l, double *r_true_t);

extern void BLAS_ddot_x_testgen(int n, int n_fix2, int n_mix, int norm,
			     enum blas_conj_type conj,
			     double *alpha, int alpha_flag,
			     double *beta, int beta_flag,
			     double *x_l, double *x_t,
			     double *y, int *seed, double *r,
			     double *r_true_l, double *r_true_t);

extern void BLAS_ddot_s_x_testgen(int n, int n_fix2, int n_mix, int norm,
			       enum blas_conj_type conj,
			       double *alpha, int alpha_flag,
			       double *beta, int beta_flag,
			       double *x_l, double *x_t,
			       float *y, int *seed, double *r,
			       double *r_true_l, double *r_true_t);


extern void test_BLAS_ssum(int n, float sum_comp, double sum_true_l, 
   double sum_true_t, float* x, int incx, 
   double eps_int, double un_int, double *test_ratio);
extern void test_BLAS_dsum(int n, double sum_comp, double sum_true_l, 
   double sum_true_t, double* x, int incx, 
   double eps_int, double un_int, double *test_ratio);
extern void test_BLAS_csum(int n, const void* sum_comp, double* sum_true_l, 
   double* sum_true_t, void* x, int incx, 
   double eps_int, double un_int, double *test_ratio);
extern void test_BLAS_zsum(int n, const void* sum_comp, double* sum_true_l, 
   double* sum_true_t, void* x, int incx, 
   double eps_int, double un_int, double *test_ratio);

extern void BLAS_ssum_testgen(int n, int norm, float* x, int *seed, 
   double *sum_true_l, double *sum_true_t);
extern void BLAS_dsum_testgen(int n, int norm, double* x, int *seed, 
   double *sum_true_l, double *sum_true_t);
extern void BLAS_csum_testgen(int n, int norm, void* x, int *seed, 
   double *sum_true_l, double *sum_true_t);
extern void BLAS_zsum_testgen(int n, int norm, void* x, int *seed, 
   double *sum_true_l, double *sum_true_t);


/* for DOT2 / GEMV2 testing */

extern void BLAS_sdot2_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      float* alpha, int alpha_flag,
      float* beta,  int beta_flag,
      float* head_x, float* tail_x, float* y, int *seed, float* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_ddot2_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      double* alpha, int alpha_flag,
      double* beta,  int beta_flag,
      double* head_x, double* tail_x, double* y, int *seed, double* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_ddot2_s_s_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      double* alpha, int alpha_flag,
      double* beta,  int beta_flag,
      float* head_x, float* tail_x, float* y, int *seed, double* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_ddot2_s_d_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      double* alpha, int alpha_flag,
      double* beta,  int beta_flag,
      float* head_x, float* tail_x, double* y, int *seed, double* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_ddot2_d_s_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      double* alpha, int alpha_flag,
      double* beta,  int beta_flag,
      double* head_x, double* tail_x, float* y, int *seed, double* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_cdot2_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      void* alpha, int alpha_flag,
      void* beta,  int beta_flag,
      void* head_x, void* tail_x, void* y, int *seed, void* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_zdot2_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      void* alpha, int alpha_flag,
      void* beta,  int beta_flag,
      void* head_x, void* tail_x, void* y, int *seed, void* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_zdot2_c_c_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      void* alpha, int alpha_flag,
      void* beta,  int beta_flag,
      void* head_x, void* tail_x, void* y, int *seed, void* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_zdot2_c_z_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      void* alpha, int alpha_flag,
      void* beta,  int beta_flag,
      void* head_x, void* tail_x, void* y, int *seed, void* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_zdot2_z_c_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      void* alpha, int alpha_flag,
      void* beta,  int beta_flag,
      void* head_x, void* tail_x, void* y, int *seed, void* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_cdot2_s_s_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      void* alpha, int alpha_flag,
      void* beta,  int beta_flag,
      float* head_x, float* tail_x, float* y, int *seed, void* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_cdot2_s_c_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      void* alpha, int alpha_flag,
      void* beta,  int beta_flag,
      float* head_x, float* tail_x, void* y, int *seed, void* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_cdot2_c_s_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      void* alpha, int alpha_flag,
      void* beta,  int beta_flag,
      void* head_x, void* tail_x, float* y, int *seed, void* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_zdot2_d_d_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      void* alpha, int alpha_flag,
      void* beta,  int beta_flag,
      double* head_x, double* tail_x, double* y, int *seed, void* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_zdot2_z_d_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      void* alpha, int alpha_flag,
      void* beta,  int beta_flag,
      void* head_x, void* tail_x, double* y, int *seed, void* r,
      double *r_true_l, double *r_true_t);
extern void BLAS_zdot2_d_z_testgen(int n, int n_fix2, int n_mix,
      int norm, enum blas_conj_type conj,
      void* alpha, int alpha_flag,
      void* beta,  int beta_flag,
      double* head_x, double* tail_x, void* y, int *seed, void* r,
      double *r_true_l, double *r_true_t);
extern void
testgen_BLAS_sdot2(int n, int n_fix2, int n_mix, int norm,
		   enum blas_conj_type conj,
		   float *alpha, int alpha_flag, float *beta, int beta_flag,
		   float *head_x, float *tail_x, float *y, int *seed,
		   float *r, double *r_true_l, double *r_true_t);
extern void
testgen_BLAS_ddot2(int n, int n_fix2, int n_mix, int norm,
		   enum blas_conj_type conj,
		   double *alpha, int alpha_flag, double *beta, int beta_flag,
		   double *head_x, double *tail_x, double *y, int *seed,
		   double *r, double *r_true_l, double *r_true_t);
extern void
testgen_BLAS_cdot2(int n, int n_fix2, int n_mix, int norm,
		   enum blas_conj_type conj,
		   void *alpha, int alpha_flag, void *beta, int beta_flag,
		   void *head_x, void *tail_x, void *y, int *seed,
		   void *r, double r_true_l[], double r_true_t[]);
extern void
testgen_BLAS_zdot2(int n, int n_fix2, int n_mix, int norm,
		   enum blas_conj_type conj,
		   void *alpha, int alpha_flag, void *beta, int beta_flag,
		   void *head_x, void *tail_x, void *y, int *seed,
		   void *r, double r_true_l[], double r_true_t[]);

extern void BLAS_sdot2_testgen(int n, int n_fix2, int n_mix, int norm,
			       enum blas_conj_type conj,
			       float *alpha, int alpha_flag,
			       float *beta, int beta_flag,
			       float *head_x, float *tail_x, float *y, int *seed,
			       float *r, double *r_true_l, double *r_true_t);
extern void BLAS_ddot2_testgen(int n, int n_fix2, int n_mix, int norm,
			       enum blas_conj_type conj,
			       double *alpha, int alpha_flag,
			       double *beta, int beta_flag,
			       double *head_x, double *tail_x, double *y, int *seed,
			       double *r, double *r_true_l, double *r_true_t);
extern void BLAS_cdot2_testgen(int n, int n_fix2, int n_mix, int norm,
			       enum blas_conj_type conj,
			       void *alpha, int alpha_flag,
			       void *beta, int beta_flag,
			       void *head_x, void *tail_x, void *y, int *seed,
			       void *r, double *r_true_l, double *r_true_t);
extern void BLAS_zdot2_testgen(int n, int n_fix2, int n_mix, int norm,
			       enum blas_conj_type conj,
			       void *alpha, int alpha_flag,
			       void *beta, int beta_flag,
			       void *head_x, void *tail_x, void *y, int *seed,
			       void *r, double *r_true_l, double *r_true_t);
extern void
test_BLAS_sdot2(int n, enum blas_conj_type conj,
		float alpha, float beta, float rin, float rout,
		double r_true_l, double r_true_t,
		float *x, int incx, float *head_y, float *tail_y, int incy,
		double eps_int, double un_int, double *test_ratio);
extern void
test_BLAS_ddot2(int n, enum blas_conj_type conj,
		double alpha, double beta, double rin, double rout,
		double r_true_l, double r_true_t,
		double *x, int incx, double *head_y, double *tail_y, int incy,
		double eps_int, double un_int, double *test_ratio);
extern void
test_BLAS_ddot2_s_s(int n, enum blas_conj_type conj,
		    double alpha, double beta, double rin, double rout,
		    double r_true_l, double r_true_t,
		    float *x, int incx, float *head_y, float *tail_y,
		    int incy, double eps_int, double un_int,
		    double *test_ratio);
extern void
test_BLAS_ddot2_s_d(int n, enum blas_conj_type conj,
		    double alpha, double beta, double rin, double rout,
		    double r_true_l, double r_true_t,
		    float *x, int incx, double *head_y, double *tail_y,
		    int incy, double eps_int, double un_int,
		    double *test_ratio);
extern void
test_BLAS_ddot2_d_s(int n, enum blas_conj_type conj,
		    double alpha, double beta, double rin, double rout,
		    double r_true_l, double r_true_t,
		    double *x, int incx, float *head_y, float *tail_y,
		    int incy, double eps_int, double un_int,
		    double *test_ratio);
extern void
test_BLAS_cdot2(int n, enum blas_conj_type conj,
		const void *alpha, const void *beta, const void *rin,
		const void *rout, double *r_true_l, double *r_true_t, void *x,
		int incx, void *head_y, void *tail_y, int incy,
		double eps_int, double un_int, double *test_ratio);
extern void
test_BLAS_zdot2(int n, enum blas_conj_type conj,
		const void *alpha, const void *beta, const void *rin,
		const void *rout, double *r_true_l, double *r_true_t, void *x,
		int incx, void *head_y, void *tail_y, int incy,
		double eps_int, double un_int, double *test_ratio);
extern void
test_BLAS_zdot2_c_c(int n, enum blas_conj_type conj,
		    const void *alpha, const void *beta, const void *rin,
		    const void *rout, double *r_true_l, double *r_true_t,
		    void *x, int incx, void *head_y, void *tail_y, int incy,
		    double eps_int, double un_int, double *test_ratio);
extern void
test_BLAS_zdot2_c_z(int n, enum blas_conj_type conj,
		    const void *alpha, const void *beta, const void *rin,
		    const void *rout, double *r_true_l, double *r_true_t,
		    void *x, int incx, void *head_y, void *tail_y, int incy,
		    double eps_int, double un_int, double *test_ratio);
extern void
test_BLAS_zdot2_z_c(int n, enum blas_conj_type conj,
		    const void *alpha, const void *beta, const void *rin,
		    const void *rout, double *r_true_l, double *r_true_t,
		    void *x, int incx, void *head_y, void *tail_y, int incy,
		    double eps_int, double un_int, double *test_ratio);
extern void
test_BLAS_cdot2_s_s(int n, enum blas_conj_type conj,
		    const void *alpha, const void *beta, const void *rin,
		    const void *rout, double *r_true_l, double *r_true_t,
		    float *x, int incx, float *head_y, float *tail_y,
		    int incy, double eps_int, double un_int,
		    double *test_ratio);
extern void
test_BLAS_cdot2_s_c(int n, enum blas_conj_type conj,
		    const void *alpha, const void *beta, const void *rin,
		    const void *rout, double *r_true_l, double *r_true_t,
		    float *x, int incx, void *head_y, void *tail_y, int incy,
		    double eps_int, double un_int, double *test_ratio);
extern void
test_BLAS_cdot2_c_s(int n, enum blas_conj_type conj,
		    const void *alpha, const void *beta, const void *rin,
		    const void *rout, double *r_true_l, double *r_true_t,
		    void *x, int incx, float *head_y, float *tail_y, int incy,
		    double eps_int, double un_int, double *test_ratio);
extern void
test_BLAS_zdot2_d_d(int n, enum blas_conj_type conj,
		    const void *alpha, const void *beta, const void *rin,
		    const void *rout, double *r_true_l, double *r_true_t,
		    double *x, int incx, double *head_y, double *tail_y,
		    int incy, double eps_int, double un_int,
		    double *test_ratio);
extern void
test_BLAS_zdot2_z_d(int n, enum blas_conj_type conj,
		    const void *alpha, const void *beta, const void *rin,
		    const void *rout, double *r_true_l, double *r_true_t,
		    void *x, int incx, double *head_y, double *tail_y,
		    int incy, double eps_int, double un_int,
		    double *test_ratio);
extern void
test_BLAS_zdot2_d_z(int n, enum blas_conj_type conj,
		    const void *alpha, const void *beta, const void *rin,
		    const void *rout, double *r_true_l, double *r_true_t,
		    double *x, int incx, void *head_y, void *tail_y, int incy,
		    double eps_int, double un_int, double *test_ratio);

void BLAS_sdot2_x(enum blas_conj_type conj, int n, float alpha,
		  const float *x, int incx, float beta,
		  const float *head_y, const float *tail_y, int incy,
		  float *r, enum blas_prec_type prec);
void BLAS_ddot2_x(enum blas_conj_type conj, int n, double alpha,
		  const double *x, int incx, double beta,
		  const double *head_y, const double *tail_y, int incy,
		  double *r, enum blas_prec_type prec);
void BLAS_cdot2_x(enum blas_conj_type conj, int n, const void *alpha,
		  const void *x, int incx, const void *beta,
		  const void *head_y, const void *tail_y, int incy,
		  void *r, enum blas_prec_type prec);
void BLAS_zdot2_x(enum blas_conj_type conj, int n, const void *alpha,
		  const void *x, int incx, const void *beta,
		  const void *head_y, const void *tail_y, int incy,
		  void *r, enum blas_prec_type prec);

void s_r_truth2(enum blas_conj_type conj, int n, float alpha,
		const float *x, int incx, float beta,
		const float *head_y, const float *tail_y, int incy, float *r,
		double *head_r_true, double *tail_r_true);
void d_r_truth2(enum blas_conj_type conj, int n, double alpha,
		const double *x, int incx, double beta,
		const double *head_y, const double *tail_y, int incy,
		double *r, double *head_r_true, double *tail_r_true);
void c_r_truth2(enum blas_conj_type conj, int n, const void *alpha,
		const void *x, int incx, const void *beta,
		const void *head_y, const void *tail_y, int incy,
		const void *r, double *head_r_true, double *tail_r_true);
void z_r_truth2(enum blas_conj_type conj, int n, const void *alpha,
		const void *x, int incx, const void *beta,
		const void *head_y, const void *tail_y, int incy,
		const void *r, double *head_r_true, double *tail_r_true);


/****************************
 *   Level 2 routines       *
 ****************************/

extern void BLAS_sspmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   float* alpha, int alpha_flag, float* beta, int beta_flag, 
   float* a, float* x, int incx, float* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dspmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   double* a, double* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dspmv_d_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   double* a, float* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dspmv_s_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   float* a, double* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dspmv_s_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   float* a, float* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_cspmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zspmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zspmv_c_z_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zspmv_z_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zspmv_c_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_zspmv_z_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, double* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zspmv_d_z_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   double* a, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zspmv_d_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   double* a, double* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_cspmv_c_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, float* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cspmv_s_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   float* a, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cspmv_s_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   float* a, float* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

extern void sspmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, const float* a, 
   float* a_vec, int row);
extern void dspmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, const double* a, 
   double* a_vec, int row);
extern void cspmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, const void* a, 
   void* a_vec, int row);
extern void zspmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, const void* a, 
   void* a_vec, int row);

extern void sspmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, float* a_packed, 
   float* a_full, int lda);
extern void dspmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, double* a_packed, 
   double* a_full, int lda);
extern void cspmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a_packed, 
   void* a_full, int lda);
extern void zspmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a_packed, 
   void* a_full, int lda);

extern void sprint_spmv_matrix (float* a, int n, 
     enum blas_order_type order, enum blas_uplo_type uplo);
extern void dprint_spmv_matrix (double* a, int n, 
     enum blas_order_type order, enum blas_uplo_type uplo);
extern void cprint_spmv_matrix (void* a, int n, 
     enum blas_order_type order, enum blas_uplo_type uplo);
extern void zprint_spmv_matrix (void* a, int n, 
     enum blas_order_type order, enum blas_uplo_type uplo);


extern void sgbmv_prepare(enum blas_order_type order, 
		     enum blas_trans_type trans, int m, int n, int kl, int ku,
                     float * AB, int lda, float* y, int row,
                     int *nfix2, int *nmix, int *ysize);

extern void sgbmv_commit(enum blas_order_type order, 
		     enum blas_trans_type trans, int m, int n, int kl, int ku,
                     float* AB, int lda, float* y, int row);

extern void sgbmv_copy(enum blas_order_type order, 
		  enum blas_trans_type trans, int m, int n, int kl, int ku,
                  const float* AB, int lda, float* y, int row);

extern void dgbmv_prepare(enum blas_order_type order, 
		     enum blas_trans_type trans, int m, int n, int kl, int ku,
                     double* AB, int lda, double* y, int row,
                     int *nfix2, int *nmix, int *ysize);

extern void dgbmv_commit(enum blas_order_type order, 
		     enum blas_trans_type trans, int m, int n, int kl, int ku,
                     double* AB, int lda, double* y, int row);

extern void dgbmv_copy(enum blas_order_type order, 
		  enum blas_trans_type trans, int m, int n, int kl, int ku,
                  const double* AB, int lda, double* y, int row);

extern void cgbmv_commit (enum blas_order_type order,
        enum blas_trans_type trans, int m, int n, int kl, int ku,
                void *AB, int lda, void *y, int row);

extern void cgbmv_prepare (enum blas_order_type order,
         enum blas_trans_type trans, int m, int n, int kl, int ku,
                  void *AB, int lda, void *y, int row,
                           int *nfix2, int *nmix, int *ysize);
extern void cgbmv_copy (enum blas_order_type order,
      enum blas_trans_type trans, int m, int n, int kl, int ku,
            const void *AB, int lda, void *y, int row);

extern void zgbmv_prepare (enum blas_order_type order,
         enum blas_trans_type trans, int m, int n, int kl, int ku,
                  void *AB, int lda, void *y, int row,
                           int *nfix2, int *nmix, int *ysize);

extern void zgbmv_commit (enum blas_order_type order,
        enum blas_trans_type trans, int m, int n, int kl, int ku,
                void *AB, int lda, void *y, int row);

extern void zgbmv_copy (enum blas_order_type order,
      enum blas_trans_type trans, int m, int n, int kl, int ku,
            const void *AB, int lda, void *y, int row);



extern void BLAS_sgbmv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      float* alpha, int alpha_flag, float* AB, int lda,
                      float* x, float* beta, int beta_flag, float* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_dgbmv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      double* alpha, int alpha_flag, double* AB, int lda,
                      double* x, double* beta, int beta_flag, double* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_dgbmv_d_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int ku, int kl,
                      double* alpha, int alpha_flag, double* AB, int lda,
                      float* x, double* beta, int beta_flag, double* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_dgbmv_s_d_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int ku, int kl,
                      double* alpha, int alpha_flag, float* AB, int lda,
                      double* x, double* beta, int beta_flag, double* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_dgbmv_s_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int ku, int kl,
                      double* alpha, int alpha_flag, float* AB, int lda,
                      float* x, double* beta, int beta_flag, double* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_cgbmv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      void* alpha, int alpha_flag, void* AB, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_cgbmv_s_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      void* alpha, int alpha_flag, float* AB, int lda,
                      float* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_cgbmv_c_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      void* alpha, int alpha_flag, void* AB, int lda,
                      float* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_cgbmv_s_c_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      void* alpha, int alpha_flag, float* AB, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);


extern void BLAS_zgbmv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      void* alpha, int alpha_flag, void* AB, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_zgbmv_d_d_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      void* alpha, int alpha_flag, double* AB, int lda,
                      double* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_zgbmv_z_d_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      void* alpha, int alpha_flag, void* AB, int lda,
                      double* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_zgbmv_d_z_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      void* alpha, int alpha_flag, double* AB, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_zgbmv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      void* alpha, int alpha_flag, void* AB, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_zgbmv_c_c_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      void* alpha, int alpha_flag, void* AB, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_zgbmv_z_c_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      void* alpha, int alpha_flag, void* AB, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_zgbmv_c_z_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n, int kl, int ku,
                      void* alpha, int alpha_flag, void* AB, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern
void BLAS_chbmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern
void BLAS_zhbmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern
void BLAS_zhbmv_c_z_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern
void BLAS_zhbmv_z_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern
void BLAS_zhbmv_c_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern
void BLAS_chbmv_c_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, float* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern
void BLAS_zhbmv_z_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, double* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

extern void sskew_commit_row_hbmv(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, float* a, int k, int lda, 
   float* a_vec, int row);
extern void dskew_commit_row_hbmv(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, double* a, int k, int lda, 
   double* a_vec, int row);
extern void chbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int k, int lda, 
   void* a_vec, int row);
extern void zhbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int k, int lda, 
   void* a_vec, int row);
extern void sskew_copy_row_hbmv(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, float* a, int k, int lda, 
   float* a_vec, int row);
extern void dskew_copy_row_hbmv(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, double* a, int k, int lda, 
   double* a_vec, int row);
extern void chbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int k, int lda, 
   void* a_vec, int row);
extern void zhbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int k, int lda, 
   void* a_vec, int row);
extern void cprint_hbmv_matrix (void* a, int n, int k, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);
extern void zprint_hbmv_matrix (void* a, int n, int k, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);

extern void BLAS_sskew_testgen_hbmv(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   float* alpha, float* beta, 
   float* a, int k, int lda, float* x, int incx, float* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dskew_testgen_hbmv(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, double* beta, 
   double* a, int k, int lda, double* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dskew_testgen_hbmv_d_s(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, double* beta, 
   double* a, int k, int lda, float* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dskew_testgen_hbmv_s_s(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, double* beta, 
   float* a, int k, int lda, float* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dskew_testgen_hbmv_s_s(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, double* beta, 
   float* a, int k, int lda, float* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

void BLAS_stpmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, 
   float* alpha, int alpha_flag,
   float* tp, float* x,   
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dtpmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, 
   double* alpha, int alpha_flag,
   double* tp, double* x,   
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dtpmv_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, 
   double* alpha, int alpha_flag,
   float* tp, double* x,   
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_ctpmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, 
   void* alpha, int alpha_flag,
   void* tp, void* x,   
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_ztpmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, 
   void* alpha, int alpha_flag,
   void* tp, void* x,   
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_ztpmv_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, 
   void* alpha, int alpha_flag,
   void* tp, void* x,   
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_ctpmv_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, 
   void* alpha, int alpha_flag,
   float* tp, void* x,   
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_ztpmv_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, 
   void* alpha, int alpha_flag,
   double* tp, void* x,   
   int *seed, double *r_true_l, double *r_true_t);
void stpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_trans_type trans, int n, float* a, 
   float* a_vec, int row);
void dtpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_trans_type trans, int n, double* a, 
   double* a_vec, int row);
void ctpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_trans_type trans, int n, void* a, 
   void* a_vec, int row);
void ztpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_trans_type trans, int n, void* a, 
   void* a_vec, int row);
void stpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_trans_type trans, int n, float* a, 
   float* a_vec, int row);
void dtpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_trans_type trans, int n, double* a, 
   double* a_vec, int row);
void ctpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_trans_type trans, int n, void* a, 
   void* a_vec, int row);
void ztpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_trans_type trans, int n, void* a, 
   void* a_vec, int row);
void stpmv_copy_vector(int n, float* x_to, int incx_to, 
   float* x_from, int incx_from);
void dtpmv_copy_vector(int n, double* x_to, int incx_to, 
   double* x_from, int incx_from);
void ctpmv_copy_vector(int n, void* x_to, int incx_to, 
   void* x_from, int incx_from);
void ztpmv_copy_vector(int n, void* x_to, int incx_to, 
   void* x_from, int incx_from);


extern void strsv_copy(enum blas_order_type order, 
	      	        enum blas_uplo_type uplo,
		        enum blas_trans_type trans, 
                        int n, const float* T, 
		        int lda, float* y, int row);

extern void strsv_commit(enum blas_order_type order, 
	         	  enum blas_uplo_type uplo, 
		          enum blas_trans_type trans,
                          int n, float* T, int lda,
		          const float* y, int row);

extern void dtrsv_copy(enum blas_order_type order, 
	      	        enum blas_uplo_type uplo,
		        enum blas_trans_type trans, 
                        int n, const double* T, 
		        int lda, double* y, int row);

extern void dtrsv_commit(enum blas_order_type order, 
	         	  enum blas_uplo_type uplo, 
		          enum blas_trans_type trans,
                          int n, double* T, int lda,
		          const double* y, int row);

extern void ctrsv_copy(enum blas_order_type order, 
		  enum blas_uplo_type uplo,
		  enum blas_trans_type trans, 
                  int n, const void* T, 
		  int lda, void* y, int row);

extern void ztrsv_copy(enum blas_order_type order, 
		  enum blas_uplo_type uplo,
		  enum blas_trans_type trans, 
                  int n, const void* T, 
		  int lda, void* y, int row);

extern void BLAS_strsv_testgen(int norm, enum blas_order_type order, 
			    enum blas_uplo_type uplo, enum blas_trans_type,
			    enum blas_diag_type diag, int n, float* alpha, 
			    int alpha_flag, float* T, int lda, float* x,  
			    int *seed, double *r_true_l, double *r_true_t, 
			    int row, enum blas_prec_type prec);

extern void BLAS_dtrsv_testgen(int norm, enum blas_order_type order, 
			    enum blas_uplo_type uplo, enum blas_trans_type,
			    enum blas_diag_type diag, int n, double* alpha, 
			    int alpha_flag, double* T, int lda, double* x,  
			    int *seed, double *r_true_l, double *r_true_t,
			    int row, enum blas_prec_type prec);

extern void BLAS_dtrsv_s_testgen(int norm, enum blas_order_type order, 
			      enum blas_uplo_type uplo, enum blas_trans_type,
			      enum blas_diag_type diag, int n, double* alpha, 
			      int alpha_flag, float* T, int lda, double* x,  
			      int *seed, double *r_true_l, double *r_true_t,
			      int row, enum blas_prec_type prec);

extern void BLAS_ctrsv_testgen(int norm, enum blas_order_type order, 
		      enum blas_uplo_type uplo, enum blas_trans_type trans,
		      enum blas_diag_type diag, int n, void* alpha, 
                      int alpha_flag, void* T, int lda, void* x,  
                      int *seed, double *r_true_l, double *r_true_t, 
		      int row, enum blas_prec_type prec);

extern void BLAS_ztrsv_c_testgen(int norm, enum blas_order_type order, 
		      enum blas_uplo_type uplo, enum blas_trans_type trans,
		      enum blas_diag_type diag, int n, void* alpha, 
                      int alpha_flag, void* T, int lda, void* x,  
                      int *seed, double *r_true_l, double *r_true_t, 
		      int row, enum blas_prec_type prec);

extern void BLAS_ztrsv_testgen(int norm, enum blas_order_type order, 
		      enum blas_uplo_type uplo, enum blas_trans_type trans,
		      enum blas_diag_type diag, int n, void* alpha, 
                      int alpha_flag, void* T, int lda, void* x,  
                      int *seed, double *r_true_l, double *r_true_t, 
		      int row, enum blas_prec_type prec);

extern void BLAS_ctrsv_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_uplo_type uplo, enum blas_trans_type trans,
		      enum blas_diag_type diag, int n, void* alpha, 
                      int alpha_flag, float* T, int lda, void* x,  
                      int *seed, double *r_true_l, double *r_true_t, 
		      int row, enum blas_prec_type prec);

extern void BLAS_ztrsv_d_testgen(int norm, enum blas_order_type order, 
		      enum blas_uplo_type uplo, enum blas_trans_type trans,
		      enum blas_diag_type diag, int n, void* alpha, 
                      int alpha_flag, double* T, int lda, void* x,  
                      int *seed, double *r_true_l, double *r_true_t, 
		      int row, enum blas_prec_type prec);


extern void BLAS_stbsv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, int k, int randomize, 
   float* alpha, int alpha_flag, float* T, int ldt, 
   float* x,int *seed, double *r_true_l, double *r_true_t, int row, enum blas_prec_type prec);
extern void BLAS_dtbsv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, int k, int randomize, 
   double* alpha, int alpha_flag, double* T, int ldt, 
   double* x,int *seed, double *r_true_l, double *r_true_t, int row, enum blas_prec_type prec);
extern void BLAS_dtbsv_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, int k, int randomize, 
   double* alpha, int alpha_flag, float* T, int ldt, 
   double* x,int *seed, double *r_true_l, double *r_true_t, int row, enum blas_prec_type prec);
extern void BLAS_ctbsv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, int k, int randomize, 
   void* alpha, int alpha_flag, void* T, int ldt, 
   void* x,int *seed, double *r_true_l, double *r_true_t, int row, enum blas_prec_type prec);
extern void BLAS_ztbsv_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, int k, int randomize, 
   void* alpha, int alpha_flag, void* T, int ldt, 
   void* x,int *seed, double *r_true_l, double *r_true_t, int row, enum blas_prec_type prec);
extern void BLAS_ztbsv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, int k, int randomize, 
   void* alpha, int alpha_flag, void* T, int ldt, 
   void* x,int *seed, double *r_true_l, double *r_true_t, int row, enum blas_prec_type prec);
extern void BLAS_ctbsv_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, int k, int randomize, 
   void* alpha, int alpha_flag, float* T, int ldt, 
   void* x,int *seed, double *r_true_l, double *r_true_t, int row, enum blas_prec_type prec);
extern void BLAS_ztbsv_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   enum blas_trans_type trans,
   enum blas_diag_type diag,
   int n, int k, int randomize, 
   void* alpha, int alpha_flag, double* T, int ldt, 
   void* x,int *seed, double *r_true_l, double *r_true_t, int row, enum blas_prec_type prec);

extern void stbsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
   enum blas_trans_type trans, 
   int n, int k, float* T, int ldt, 
   float* y, int row);
extern void dtbsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
   enum blas_trans_type trans, 
   int n, int k, double* T, int ldt, 
   double* y, int row);
extern void ctbsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
   enum blas_trans_type trans, 
   int n, int k, void* T, int ldt, 
   void* y, int row);
extern void ztbsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
   enum blas_trans_type trans, 
   int n, int k, void* T, int ldt, 
   void* y, int row);

extern void stbsv_copy(enum blas_order_type order, 
   enum blas_uplo_type uplo, 
   enum blas_trans_type trans,
   int n, int k, const float* T, int ldt, 
   float* y, int row);
extern void dtbsv_copy(enum blas_order_type order, 
   enum blas_uplo_type uplo, 
   enum blas_trans_type trans,
   int n, int k, const double* T, int ldt, 
   double* y, int row);
extern void ctbsv_copy(enum blas_order_type order, 
   enum blas_uplo_type uplo, 
   enum blas_trans_type trans,
   int n, int k, const void* T, int ldt, 
   void* y, int row);
extern void ztbsv_copy(enum blas_order_type order, 
   enum blas_uplo_type uplo, 
   enum blas_trans_type trans,
   int n, int k, const void* T, int ldt, 
   void* y, int row);

extern void sprint_tbsv_matrix (float* a, int n, int k, int ldt, 
     enum blas_order_type order, enum blas_uplo_type uplo, enum blas_trans_type);
extern void dprint_tbsv_matrix (double* a, int n, int k, int ldt, 
     enum blas_order_type order, enum blas_uplo_type uplo, enum blas_trans_type);
extern void cprint_tbsv_matrix (void* a, int n, int k, int ldt, 
     enum blas_order_type order, enum blas_uplo_type uplo, enum blas_trans_type);
extern void zprint_tbsv_matrix (void* a, int n, int k, int ldt, 
     enum blas_order_type order, enum blas_uplo_type uplo, enum blas_trans_type);


extern void sgemv_commit(enum blas_order_type order, 
		     enum blas_trans_type trans, int m, int n,
                     float* A, int lda, float* y, int row);

extern void sgemv_copy(enum blas_order_type order, 
		  enum blas_trans_type trans, int m, int n,
                  float* A, int lda, float* y, int row);
extern void dgemv_commit(enum blas_order_type order, 
		     enum blas_trans_type trans, int m, int n,
                     double* A, int lda, double* y, int row);

extern void dgemv_copy(enum blas_order_type order, 
		  enum blas_trans_type trans, int m, int n,
                  double* A, int lda, double* y, int row);
extern void cgemv_commit(enum blas_order_type order, 
		     enum blas_trans_type trans, int m, int n,
                     void* A, int lda, void* y, int row);

extern void cgemv_copy(enum blas_order_type order, 
		  enum blas_trans_type trans, int m, int n,
                  void* A, int lda, void* y, int row);
extern void zgemv_commit(enum blas_order_type order, 
		     enum blas_trans_type trans, int m, int n,
                     void* A, int lda, void* y, int row);

extern void zgemv_copy(enum blas_order_type order, 
		  enum blas_trans_type trans, int m, int n,
                  void* A, int lda, void* y, int row);


extern void BLAS_sgemv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      float* alpha, int alpha_flag, float* A, int lda,
                      float* x, float* beta, int beta_flag, float* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_dgemv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      double* alpha, int alpha_flag, double* A, int lda,
                      double* x, double* beta, int beta_flag, double* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_cgemv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_zgemv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dgemv_s_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      double* alpha, int alpha_flag, float* A, int lda,
                      float* x, double* beta, int beta_flag, double* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dgemv_s_d_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      double* alpha, int alpha_flag, float* A, int lda,
                      double* x, double* beta, int beta_flag, double* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dgemv_d_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      double* alpha, int alpha_flag, double* A, int lda,
                      float* x, double* beta, int beta_flag, double* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemv_c_c_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemv_c_z_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);         
extern void BLAS_zgemv_z_c_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cgemv_s_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, float* A, int lda,
                      float* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cgemv_s_c_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, float* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cgemv_c_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      float* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemv_d_d_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, double* A, int lda,
                      double* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemv_z_d_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      double* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemv_d_z_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, double* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_sgemv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      float* alpha, int alpha_flag, float* A, int lda,
                      float* x, float* beta, int beta_flag, float* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_dgemv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      double* alpha, int alpha_flag, double* A, int lda,
                      double* x, double* beta, int beta_flag, double* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_cgemv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_zgemv_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dgemv_s_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      double* alpha, int alpha_flag, float* A, int lda,
                      float* x, double* beta, int beta_flag, double* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dgemv_s_d_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      double* alpha, int alpha_flag, float* A, int lda,
                      double* x, double* beta, int beta_flag, double* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dgemv_d_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      double* alpha, int alpha_flag, double* A, int lda,
                      float* x, double* beta, int beta_flag, double* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemv_c_c_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemv_c_z_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemv_z_c_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cgemv_s_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, float* A, int lda,
                      float* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cgemv_s_c_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, float* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cgemv_c_s_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      float* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemv_d_d_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, double* A, int lda,
                      double* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemv_d_z_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, double* A, int lda,
                      void* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemv_z_d_testgen(int norm, enum blas_order_type order, 
		      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      double* x, void* beta, int beta_flag, void* y, 
                      int *seed, double *r_true_l, double *r_true_t);


extern void BLAS_sge_sum_mv_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   float* alpha, int alpha_flag, float* beta, int beta_flag, 
   float* a, int lda, float* b, int ldb, 
   float* x, int incx, 
   float* alpha_use_ptr, float* a_use, float* b_use, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dge_sum_mv_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   double* a, int lda, double* b, int ldb, 
   double* x, int incx, 
   double* alpha_use_ptr, double* a_use, double* b_use, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dge_sum_mv_d_s_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   double* a, int lda, double* b, int ldb, 
   float* x, int incx, 
   double* alpha_use_ptr, double* a_use, double* b_use, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dge_sum_mv_s_d_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   float* a, int lda, float* b, int ldb, 
   double* x, int incx, 
   double* alpha_use_ptr, float* a_use, float* b_use, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dge_sum_mv_s_s_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   float* a, int lda, float* b, int ldb, 
   float* x, int incx, 
   double* alpha_use_ptr, float* a_use, float* b_use, 
   int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_cge_sum_mv_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, 
   void* x, int incx, 
   void* alpha_use_ptr, void* a_use, void* b_use, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zge_sum_mv_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, 
   void* x, int incx, 
   void* alpha_use_ptr, void* a_use, void* b_use, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zge_sum_mv_c_z_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, 
   void* x, int incx, 
   void* alpha_use_ptr, void* a_use, void* b_use, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zge_sum_mv_z_c_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, 
   void* x, int incx, 
   void* alpha_use_ptr, void* a_use, void* b_use, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zge_sum_mv_c_c_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, 
   void* x, int incx, 
   void* alpha_use_ptr, void* a_use, void* b_use, 
   int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_cge_sum_mv_c_s_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, 
   float* x, int incx, 
   void* alpha_use_ptr, void* a_use, void* b_use, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cge_sum_mv_s_c_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   float* a, int lda, float* b, int ldb, 
   void* x, int incx, 
   void* alpha_use_ptr, float* a_use, float* b_use, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cge_sum_mv_s_s_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   float* a, int lda, float* b, int ldb, 
   float* x, int incx, 
   void* alpha_use_ptr, float* a_use, float* b_use, 
   int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_zge_sum_mv_z_d_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, 
   double* x, int incx, 
   void* alpha_use_ptr, void* a_use, void* b_use, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zge_sum_mv_d_z_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   double* a, int lda, double* b, int ldb, 
   void* x, int incx, 
   void* alpha_use_ptr, double* a_use, double* b_use, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zge_sum_mv_d_d_testgen(int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   double* a, int lda, double* b, int ldb, 
   double* x, int incx, 
   void* alpha_use_ptr, double* a_use, double* b_use, 
   int *seed, double *r_true_l, double *r_true_t);

extern void sge_sum_mv_commit(enum blas_order_type order, 
		     int m, int n,
                     float* A, int lda, float* y, int row);
extern void sge_sum_mv_copy(enum blas_order_type order, 
		  int m, int n,
                  float* A, int lda, float* y, int row);
extern void dge_sum_mv_commit(enum blas_order_type order, 
		     int m, int n,
                     double* A, int lda, double* y, int row);
extern void dge_sum_mv_copy(enum blas_order_type order, 
		  int m, int n,
                  double* A, int lda, double* y, int row);
extern void cge_sum_mv_commit(enum blas_order_type order, 
		     int m, int n,
                     void* A, int lda, void* y, int row);
extern void cge_sum_mv_copy(enum blas_order_type order, 
		  int m, int n,
                  void* A, int lda, void* y, int row);
extern void zge_sum_mv_commit(enum blas_order_type order, 
		     int m, int n,
                     void* A, int lda, void* y, int row);
extern void zge_sum_mv_copy(enum blas_order_type order, 
		  int m, int n,
                  void* A, int lda, void* y, int row);

/* for GEMV2 testing */

void BLAS_sgemv2_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      float* alpha, int alpha_flag, float* A, int lda,
                      float* head_x, float* tail_x, float* beta,
                      int beta_flag, float* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_dgemv2_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      double* alpha, int alpha_flag, double* A, int lda,
                      double* head_x, double* tail_x, double* beta,
                      int beta_flag, double* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_cgemv2_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* head_x, void* tail_x, void* beta,
                      int beta_flag, void* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_zgemv2_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* head_x, void* tail_x, void* beta,
                      int beta_flag, void* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_dgemv2_s_s_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      double* alpha, int alpha_flag, float* A, int lda,
                      float* head_x, float* tail_x, double* beta,
                      int beta_flag, double* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_dgemv2_s_d_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      double* alpha, int alpha_flag, float* A, int lda,
                      double* head_x, double* tail_x, double* beta,
                      int beta_flag, double* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_dgemv2_d_s_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      double* alpha, int alpha_flag, double* A, int lda,
                      float* head_x, float* tail_x, double* beta,
                      int beta_flag, double* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_zgemv2_c_c_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* head_x, void* tail_x, void* beta,
                      int beta_flag, void* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_zgemv2_c_z_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* head_x, void* tail_x, void* beta,
                      int beta_flag, void* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_zgemv2_z_c_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      void* head_x, void* tail_x, void* beta,
                      int beta_flag, void* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_cgemv2_s_s_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, float* A, int lda,
                      float* head_x, float* tail_x, void* beta,
                      int beta_flag, void* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_cgemv2_s_c_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, float* A, int lda,
                      void* head_x, void* tail_x, void* beta,
                      int beta_flag, void* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_cgemv2_c_s_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      float* head_x, float* tail_x, void* beta,
                      int beta_flag, void* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_zgemv2_d_d_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, double* A, int lda,
                      double* head_x, double* tail_x, void* beta,
                      int beta_flag, void* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_zgemv2_d_z_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, double* A, int lda,
                      void* head_x, void* tail_x, void* beta,
                      int beta_flag, void* y, int *seed,
                      double *r_true_l, double *r_true_t);
void BLAS_zgemv2_z_d_testgen(int norm, enum blas_order_type order, 
                      enum blas_trans_type trans, int m, int n,
                      void* alpha, int alpha_flag, void* A, int lda,
                      double* head_x, double* tail_x, void* beta,
                      int beta_flag, void* y, int *seed,
                      double *r_true_l, double *r_true_t);
void sgemv2_commit(enum blas_order_type order,
		   enum blas_trans_type trans, int m, int n,
		   float *A, int lda, float *y, int row);
void sgemv2_copy(enum blas_order_type order,
		 enum blas_trans_type trans, int m, int n,
		 float *A, int lda, float *y, int row);
void dgemv2_commit(enum blas_order_type order,
		   enum blas_trans_type trans, int m, int n,
		   double *A, int lda, double *y, int row);
void dgemv2_copy(enum blas_order_type order,
		 enum blas_trans_type trans, int m, int n,
		 double *A, int lda, double *y, int row);
void cgemv2_commit(enum blas_order_type order,
		   enum blas_trans_type trans, int m, int n,
		   void *A, int lda, void *y, int row);
void cgemv2_copy(enum blas_order_type order,
		 enum blas_trans_type trans, int m, int n,
		 void *A, int lda, void *y, int row);
void zgemv2_commit(enum blas_order_type order,
		   enum blas_trans_type trans, int m, int n,
		   void *A, int lda, void *y, int row);
void zgemv2_copy(enum blas_order_type order,
		 enum blas_trans_type trans, int m, int n,
		 void *A, int lda, void *y, int row);


/****************************
 *   Level 3 routines       *
 ****************************/

extern void BLAS_sgemm_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        float* alpha, int alpha_flag, float* a, int lda, 
        float* beta, int beta_flag, float* b, int ldb, 
        float* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dgemm_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        double* alpha, int alpha_flag, double* a, int lda, 
        double* beta, int beta_flag, double* b, int ldb, 
        double* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dgemm_d_s_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        double* alpha, int alpha_flag, double* a, int lda, 
        double* beta, int beta_flag, float* b, int ldb, 
        double* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dgemm_s_d_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        double* alpha, int alpha_flag, float* a, int lda, 
        double* beta, int beta_flag, double* b, int ldb, 
        double* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dgemm_s_s_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        double* alpha, int alpha_flag, float* a, int lda, 
        double* beta, int beta_flag, float* b, int ldb, 
        double* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cgemm_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        void* alpha, int alpha_flag, void* a, int lda, 
        void* beta, int beta_flag, void* b, int ldb, 
        void* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemm_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        void* alpha, int alpha_flag, void* a, int lda, 
        void* beta, int beta_flag, void* b, int ldb, 
        void* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemm_z_c_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        void* alpha, int alpha_flag, void* a, int lda, 
        void* beta, int beta_flag, void* b, int ldb, 
        void* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemm_c_z_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        void* alpha, int alpha_flag, void* a, int lda, 
        void* beta, int beta_flag, void* b, int ldb, 
        void* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemm_c_c_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        void* alpha, int alpha_flag, void* a, int lda, 
        void* beta, int beta_flag, void* b, int ldb, 
        void* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cgemm_c_s_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        void* alpha, int alpha_flag, void* a, int lda, 
        void* beta, int beta_flag, float* b, int ldb, 
        void* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cgemm_s_c_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        void* alpha, int alpha_flag, float* a, int lda, 
        void* beta, int beta_flag, void* b, int ldb, 
        void* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_cgemm_s_s_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        void* alpha, int alpha_flag, float* a, int lda, 
        void* beta, int beta_flag, float* b, int ldb, 
        void* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemm_z_d_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        void* alpha, int alpha_flag, void* a, int lda, 
        void* beta, int beta_flag, double* b, int ldb, 
        void* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemm_d_z_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        void* alpha, int alpha_flag, double* a, int lda, 
        void* beta, int beta_flag, void* b, int ldb, 
        void* c, int ldc, int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zgemm_d_d_testgen(int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        void* alpha, int alpha_flag, double* a, int lda, 
        void* beta, int beta_flag, double* b, int ldb, 
        void* c, int ldc, int *seed, double *r_true_l, double *r_true_t);

extern void sgemm_commit_row (enum blas_order_type order,
   enum blas_trans_type trans, 
   int m, int n, float* a, int lda, float* a_vec, int row);
extern void dgemm_commit_row (enum blas_order_type order,
   enum blas_trans_type trans, 
   int m, int n, double* a, int lda, double* a_vec, int row);
extern void cgemm_commit_row (enum blas_order_type order,
   enum blas_trans_type trans, 
   int m, int n, void* a, int lda, void* a_vec, int row);
extern void zgemm_commit_row (enum blas_order_type order,
   enum blas_trans_type trans, 
   int m, int n, void* a, int lda, void* a_vec, int row);

extern void sgemm_commit_col (enum blas_order_type order,
   enum blas_trans_type trans, 
   int m, int n, float* a, int lda, float* a_vec, int col);
extern void dgemm_commit_col (enum blas_order_type order,
   enum blas_trans_type trans, 
   int m, int n, double* a, int lda, double* a_vec, int col);
extern void cgemm_commit_col (enum blas_order_type order,
   enum blas_trans_type trans, 
   int m, int n, void* a, int lda, void* a_vec, int col);
extern void zgemm_commit_col (enum blas_order_type order,
   enum blas_trans_type trans, 
   int m, int n, void* a, int lda, void* a_vec, int col);

extern void sgemm_copy_row 
  (enum blas_order_type order, enum blas_trans_type trans, 
   int m, int n, float* a, int lda, float* a_vec, int row);
extern void dgemm_copy_row 
  (enum blas_order_type order, enum blas_trans_type trans, 
   int m, int n, double* a, int lda, double* a_vec, int row);
extern void cgemm_copy_row 
  (enum blas_order_type order, enum blas_trans_type trans, 
   int m, int n, void* a, int lda, void* a_vec, int row);
extern void zgemm_copy_row 
  (enum blas_order_type order, enum blas_trans_type trans, 
   int m, int n, void* a, int lda, void* a_vec, int row);

extern void sgemm_copy_col
  (enum blas_order_type order, enum blas_trans_type trans, 
   int m, int n, float* a, int lda, float* a_vec, int col);
extern void dgemm_copy_col
  (enum blas_order_type order, enum blas_trans_type trans, 
   int m, int n, double* a, int lda, double* a_vec, int col);
extern void cgemm_copy_col
  (enum blas_order_type order, enum blas_trans_type trans, 
   int m, int n, void* a, int lda, void* a_vec, int col);
extern void zgemm_copy_col
  (enum blas_order_type order, enum blas_trans_type trans, 
   int m, int n, void* a, int lda, void* a_vec, int col);

extern void sgemm_copy(enum blas_order_type order, int m, int n, 
    float* a, int lda, float* b, int ldb);
extern void dgemm_copy(enum blas_order_type order, int m, int n, 
    double* a, int lda, double* b, int ldb);
extern void cgemm_copy(enum blas_order_type order, int m, int n, 
    void* a, int lda, void* b, int ldb);
extern void zgemm_copy(enum blas_order_type order, int m, int n, 
    void* a, int lda, void* b, int ldb);

extern void sgemm_zero(enum blas_order_type order, int m, int n, float* a, int lda);
extern void dgemm_zero(enum blas_order_type order, int m, int n, double* a, int lda);
extern void cgemm_zero(enum blas_order_type order, int m, int n, void* a, int lda);
extern void zgemm_zero(enum blas_order_type order, int m, int n, void* a, int lda);

extern void sprint_matrix (float* a, int m, int n, int lda, 
     enum blas_order_type order);
extern void dprint_matrix (double* a, int m, int n, int lda, 
     enum blas_order_type order);
extern void cprint_matrix (void* a, int m, int n, int lda, 
     enum blas_order_type order);
extern void zprint_matrix (void* a, int m, int n, int lda, 
     enum blas_order_type order);


void BLAS_ssymm_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   float* alpha, int alpha_flag, float* beta, int beta_flag, 
   float* a, int lda, float* b, int ldb, float* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dsymm_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   double* a, int lda, double* b, int ldb, double* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dsymm_d_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   double* a, int lda, float* b, int ldb, double* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dsymm_s_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   float* a, int lda, double* b, int ldb, double* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dsymm_s_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   float* a, int lda, float* b, int ldb, double* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);

void BLAS_csymm_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymm_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymm_c_z_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymm_z_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymm_c_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);

void BLAS_csymm_c_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, float* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_csymm_s_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   float* a, int lda, void* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_csymm_s_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   float* a, int lda, float* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);

void BLAS_zsymm_z_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, double* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymm_d_z_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   double* a, int lda, void* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymm_d_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   double* a, int lda, double* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);

void ssymm_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, float* a, int lda, 
   float* a_vec, int row);
void dsymm_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, double* a, int lda, 
   double* a_vec, int row);
void csymm_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, void* a, int lda, 
   void* a_vec, int row);
void zsymm_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, void* a, int lda, 
   void* a_vec, int row);

void ssymm_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, float* a, int lda, 
   float* a_vec, int row);
void dsymm_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, double* a, int lda, 
   double* a_vec, int row);
void csymm_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, void* a, int lda, 
   void* a_vec, int row);
void zsymm_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, void* a, int lda, 
   void* a_vec, int row);

void sprint_symm_matrix (float* a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);
void dprint_symm_matrix (double* a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);
void cprint_symm_matrix (void* a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);
void zprint_symm_matrix (void* a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);



void BLAS_chemm_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhemm_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhemm_c_z_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhemm_z_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhemm_c_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);

void BLAS_chemm_c_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, float* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhemm_z_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, double* b, int ldb, void* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);

void chemm_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, void* a, int lda, 
   void* a_vec, int row);
void zhemm_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, void* a, int lda, 
   void* a_vec, int row);

void chemm_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, void* a, int lda, 
   void* a_vec, int row);
void zhemm_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, void* a, int lda, 
   void* a_vec, int row);

void sskew_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, float* a, int lda, 
   float* a_vec, int row);
void dskew_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, double* a, int lda, 
   double* a_vec, int row);

void sskew_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, float* a, int lda, 
   float* a_vec, int row);
void dskew_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, double* a, int lda, 
   double* a_vec, int row);

void BLAS_sskew_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n,  
   float* alpha, int alpha_flag, float* beta, int beta_flag, 
   float* a, int lda, float* b, int ldb, float* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dskew_s_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n,  
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   float* a, int lda, float* b, int ldb, double* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dskew_d_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n,  
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   double* a, int lda, float* b, int ldb, double* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dskew_s_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n,  
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   float* a, int lda, double* b, int ldb, double* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dskew_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo, enum blas_side_type side, 
   int m, int n,  
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   double* a, int lda, double* b, int ldb, double* c, int ldc, 
   int *seed, double *r_true_l, double *r_true_t);

void cprint_hemm_matrix (void* a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);
void zprint_hemm_matrix (void* a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);

   
void BLAS_ssymv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   float* alpha, int alpha_flag, float* beta, int beta_flag, 
   float* a, int lda, float* x, int incx, float* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dsymv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   double* a, int lda, double* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dsymv_d_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   double* a, int lda, float* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dsymv_s_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   float* a, int lda, double* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dsymv_s_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   float* a, int lda, float* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

void BLAS_csymv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymv_c_z_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymv_z_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymv_c_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

void BLAS_csymv_c_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, float* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_csymv_s_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   float* a, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_csymv_s_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   float* a, int lda, float* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

void BLAS_zsymv_z_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, double* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymv_d_z_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   double* a, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymv_d_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   double* a, int lda, double* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

void ssymv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, float* a, int lda, 
   float* a_vec, int row);
void dsymv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, double* a, int lda, 
   double* a_vec, int row);
void csymv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int lda, 
   void* a_vec, int row);
void zsymv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int lda, 
   void* a_vec, int row);

void ssymv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, float* a, int lda, 
   float* a_vec, int row);
void dsymv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, double* a, int lda, 
   double* a_vec, int row);
void csymv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int lda, 
   void* a_vec, int row);
void zsymv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int lda, 
   void* a_vec, int row);

void ssymv_copy_vector(int n, float* x_to, int incx_to, 
   float* x_from, int incx_from);
void dsymv_copy_vector(int n, double* x_to, int incx_to, 
   double* x_from, int incx_from);
void csymv_copy_vector(int n, void* x_to, int incx_to, 
   void* x_from, int incx_from);
void zsymv_copy_vector(int n, void* x_to, int incx_to, 
   void* x_from, int incx_from);

void sprint_symv_matrix (float* a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);
void dprint_symv_matrix (double* a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);
void cprint_symv_matrix (void* a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);
void zprint_symv_matrix (void* a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);

void sprint_vector (float* x, int n, int incx);
void dprint_vector (double* x, int n, int incx);
void cprint_vector (void* x, int n, int incx);
void zprint_vector (void* x, int n, int incx);


extern void BLAS_ssbmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   float* alpha, int alpha_flag, float* beta, int beta_flag, 
   float* a, int k, int lda, float* x, int incx, float* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dsbmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   double* a, int k, int lda, double* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dsbmv_d_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   double* a, int k, int lda, float* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dsbmv_s_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   float* a, int k, int lda, double* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_dsbmv_s_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, int alpha_flag, double* beta, int beta_flag, 
   float* a, int k, int lda, float* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_csbmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zsbmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zsbmv_c_z_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zsbmv_z_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zsbmv_c_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_zsbmv_z_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, double* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zsbmv_d_z_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   double* a, int k, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_zsbmv_d_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   double* a, int k, int lda, double* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

extern void BLAS_csbmv_c_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int k, int lda, float* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_csbmv_s_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   float* a, int k, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
extern void BLAS_csbmv_s_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   float* a, int k, int lda, float* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

extern void ssbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, float* a, int k, int lda,
   float* a_vec, int row);
extern void dsbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, double* a, int k, int lda,
   double* a_vec, int row);
extern void csbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int k, int lda,
   void* a_vec, int row);
extern void zsbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int k, int lda,
   void* a_vec, int row);

extern void ssbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, float* a, int k, int lda,
   float* a_vec, int row);
extern void dsbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, double* a, int k, int lda,
   double* a_vec, int row);
extern void csbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int k, int lda,
   void* a_vec, int row);
extern void zsbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int k, int lda,
   void* a_vec, int row);

extern void sprint_sbmv_matrix (float* a, int n,
	int k, int lda, enum blas_order_type order, 
	enum blas_uplo_type uplo);
extern void dprint_sbmv_matrix (double* a, int n,
	int k, int lda, enum blas_order_type order, 
	enum blas_uplo_type uplo);
extern void cprint_sbmv_matrix (void* a, int n,
	int k, int lda, enum blas_order_type order, 
	enum blas_uplo_type uplo);
extern void zprint_sbmv_matrix (void* a, int n,
	int k, int lda, enum blas_order_type order, 
	enum blas_uplo_type uplo);


void BLAS_chemv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhemv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhemv_c_z_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhemv_z_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhemv_c_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

void BLAS_chemv_c_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, float* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhemv_z_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, int lda, double* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

void sskew_commit_row_hemv(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, float* a, int lda, 
   float* a_vec, int row);
void dskew_commit_row_hemv(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, double* a, int lda, 
   double* a_vec, int row);

void sskew_copy_row_hemv(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, float* a, int lda, 
   float* a_vec, int row);
void dskew_copy_row_hemv(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, double* a, int lda, 
   double* a_vec, int row);
void chemv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int lda, 
   void* a_vec, int row);
void zhemv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, int lda, 
   void* a_vec, int row);

void cprint_hemv_matrix (void* a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);
void zprint_hemv_matrix (void* a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo);

void BLAS_sskew_testgen_hemv(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   float* alpha, float* beta, 
   float* a, int lda, float* x, int incx, float* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dskew_testgen_hemv(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, double* beta, 
   double* a, int lda, double* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dskew_testgen_hemv_d_s(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, double* beta, 
   double* a, int lda, float* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dskew_testgen_hemv_s_s(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, double* beta, 
   float* a, int lda, float* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_dskew_testgen_hemv_s_s(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   double* alpha, double* beta, 
   float* a, int lda, float* x, int incx, double* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);


void BLAS_chpmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhpmv_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhpmv_c_z_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhpmv_z_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhpmv_c_c_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, void* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

void BLAS_chpmv_c_s_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, float* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);
void BLAS_zhpmv_z_d_testgen(int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   void* alpha, int alpha_flag, void* beta, int beta_flag, 
   void* a, double* x, int incx, void* y, int incy, 
   int *seed, double *r_true_l, double *r_true_t);

void chpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, 
   void* a_vec, int row);
void zhpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a, 
   void* a_vec, int row);

void chpmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a_packed, 
   void* a_full, int lda);
void zhpmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, void* a_packed, 
   void* a_full, int lda);

void cprint_hpmv_matrix (void* a, int n, 
     enum blas_order_type order, enum blas_uplo_type uplo);
void zprint_hpmv_matrix (void* a, int n, 
     enum blas_order_type order, enum blas_uplo_type uplo);


extern void BLAS_strmv_testgen (int norm, enum blas_order_type order,
                                enum blas_uplo_type uplo,
                                enum blas_trans_type trans,
                                enum blas_diag_type diag, int n, float *alpha,
                                int alpha_flag, float *T, int lda, float *x,
                                int *seed, double *r_true_l,
                                double *r_true_t);
extern void BLAS_dtrmv_testgen (int norm, enum blas_order_type order,
                                enum blas_uplo_type uplo,
                                enum blas_trans_type trans,
                                enum blas_diag_type diag, int n,
                                double *alpha, int alpha_flag, double *T,
                                int lda, double *x, int *seed,
                                double *r_true_l, double *r_true_t);
extern void BLAS_dtrmv_s_testgen (int norm, enum blas_order_type order,
                                  enum blas_uplo_type uplo,
                                  enum blas_trans_type trans,
                                  enum blas_diag_type diag, int n,
                                  double *alpha, int alpha_flag, float *T,
                                  int lda, double *x, int *seed,
                                  double *r_true_l, double *r_true_t);

extern void BLAS_ctrmv_testgen (int norm, enum blas_order_type order,
                                enum blas_uplo_type uplo,
                                enum blas_trans_type trans,
                                enum blas_diag_type diag, int n, void *alpha,
                                int alpha_flag, void *T, int lda, void *x,
                                int *seed, double *r_true_l,
                                double *r_true_t);
extern void BLAS_ztrmv_testgen (int norm, enum blas_order_type order,
                                enum blas_uplo_type uplo,
                                enum blas_trans_type trans,
                                enum blas_diag_type diag, int n, void *alpha,
                                int alpha_flag, void *T, int lda, void *x,
                                int *seed, double *r_true_l,
                                double *r_true_t);
extern void BLAS_ztrmv_c_testgen (int norm, enum blas_order_type order,
                                  enum blas_uplo_type uplo,
                                  enum blas_trans_type trans,
                                  enum blas_diag_type diag, int n,
                                  void *alpha, int alpha_flag, void *T,
                                  int lda, void *x, int *seed,
                                  double *r_true_l, double *r_true_t);

extern void BLAS_ztrmv_d_testgen (int norm, enum blas_order_type order,
                                  enum blas_uplo_type uplo,
                                  enum blas_trans_type trans,
                                  enum blas_diag_type diag, int n,
                                  void *alpha, int alpha_flag, double *T,
                                  int lda, void *x, int *seed,
                                  double *r_true_l, double *r_true_t);
extern void BLAS_ctrmv_s_testgen (int norm, enum blas_order_type order,
                                  enum blas_uplo_type uplo,
                                  enum blas_trans_type trans,
                                  enum blas_diag_type diag, int n,
                                  void *alpha, int alpha_flag, float *T,
                                  int lda, void *x, int *seed,
                                  double *r_true_l, double *r_true_t);

extern void strmv_copy (enum blas_order_type order,
                        enum blas_uplo_type uplo,
                        enum blas_trans_type trans,
                        int n, const float *T, int lda, float *y, int row);
extern void dtrmv_copy (enum blas_order_type order,
                        enum blas_uplo_type uplo,
                        enum blas_trans_type trans,
                        int n, const double *T, int lda, double *y, int row);
extern void ctrmv_copy (enum blas_order_type order,
                        enum blas_uplo_type uplo,
                        enum blas_trans_type trans,
                        int n, const void *T, int lda, void *y, int row);
extern void ztrmv_copy (enum blas_order_type order,
                        enum blas_uplo_type uplo,
                        enum blas_trans_type trans,
                        int n, const void *T, int lda, void *y, int row);

extern void strmv_commit (enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans,
                          int n, float *T, int lda, const float *y, int row);
extern void dtrmv_commit (enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans,
                          int n, double *T, int lda,
                          const double *y, int row);
extern void ctrmv_commit (enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans,
                          int n, void *T, int lda, const void *y, int row);
extern void ztrmv_commit (enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans,
                          int n, void *T, int lda, const void *y, int row);

#endif
/* __CBLAS_TEST_H_ */

