#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>
#include <decomp_2d_iface.h>


int main(int argc, char **argv){
  int rank  = 0;
  int size  = 1;
  int niter = 20;
  double *titer = NULL;
  double *tfwd  = NULL;
  double *tbwd  = NULL;
  double *max_tfwd = NULL;
  double *max_tbwd = NULL;
  double *min_tfwd = NULL;
  double *min_tbwd = NULL;
  int nx = 8;
  int ny = 8;
  int nz = 8;
  int p_row, p_col;
  int i, j, k;
  double scaling = 0.0;
  MPI_Comm comm;

  if ( argc < 4 ) {
    fprintf ( stderr, "Error, the problem size must be given in the commandline\n" );
    return -1;
  }

  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  nz = atoi(argv[3]);
  
  nx = 32;

  p_row = 0; // The library will handle it
  p_col = 0; // The library will handle it

  scaling = 1 / (  (double)nx * ny * nz );

  MPI_Init(&argc, &argv);
  MPI_Comm_dup ( MPI_COMM_WORLD, &comm );

  decomp_2d_init( nx, ny, nz, p_row, p_col );
  MPI_Barrier( comm );
  decomp_2d_fft_init(PHYSICAL_IN_X);

  MPI_Comm_rank ( comm, &rank );
  MPI_Comm_size ( comm, &size );

  titer     = (double*) malloc ( niter * sizeof(double));
  tfwd      = (double*) malloc ( niter * sizeof(double));
  tbwd      = (double*) malloc ( niter * sizeof(double));
  max_tfwd  = (double*) malloc ( niter * sizeof(double));
  max_tbwd  = (double*) malloc ( niter * sizeof(double));
  min_tfwd  = (double*) malloc ( niter * sizeof(double));
  min_tbwd  = (double*) malloc ( niter * sizeof(double));

  /*get local data size and allocate*/
  int lxsize = 0;
  int lysize = 0;
  int lzsize = 0;
  double complex *data = NULL;
  double complex *ref = NULL;

  decomp_2d_get_local_sizes( &lxsize, &lysize, &lzsize);

  size_t localsize = lxsize * lysize * lzsize * sizeof(double complex);
  printf( "Local size = %zu\n", localsize );

  MPI_Barrier( comm );

#if 1
  data = (double complex*) malloc ( localsize );
  if ( ! data ) {
    fprintf( stderr, "Error, cannot allocate data of size %zu\n",
      localsize);
    goto finalize;
  }
  ref  = (double complex*) malloc ( localsize );
  if ( ! ref) {
    fprintf( stderr, "Error, cannot allocate ref of size %zu\n",
      localsize);
    goto finalize;
  }

  /*initialize rin to some functionmy_func(x,y,z) */
  for (i = 0; i < lxsize; ++i) {
    double tmp_x = i / nx;
    for (j = 0; j < lysize; ++j) {
      double tmp_y = j / ny;
      for (k = 0; k < lzsize; ++k) {
        double tmp_z = (k/nz); 
      //data[(i*M + j) * N + k] = rand() / ((double) RAND_MAX ) + 0 * I;
      //ref [(i*M + j) * N + k] = data[(i*M + j) * N + k];
        data[(i*lysize + j) * lzsize + k] = tmp_x * tmp_y * tmp_z;
        ref [(i*lysize + j) * lzsize + k] = data[(i*lysize + j) * lzsize + k];
      }
    }
  }

  MPI_Barrier( comm );

  /*compute transforms as many times as desired*/
  double tmp    = 0.0;
  double ti     = 0.0;
  double tn     = 0.0;
  double max_tn = 0.0;
  
  // Warmup
  decomp_2d_fft_3d_c2c( lxsize, lysize, lzsize, ref, 
      lxsize, lysize, lzsize, data, FORWARD );
  decomp_2d_fft_3d_c2c( lxsize, lysize, lzsize, data, 
      lxsize, lysize, lzsize, data, BACKWARD );

  // Normalization
  for (i = 0; i < lxsize; ++i)
    for (j = 0; j < lysize; ++j)
      for (k = 0; k < lzsize; ++k) {
        data[(i*lysize + j) * lzsize + k] *= scaling;
      }

  MPI_Barrier( comm );
  for ( int iter = 0; iter < niter; ++iter ) {
    ti = MPI_Wtime();
    decomp_2d_fft_3d_c2c( lxsize, lysize, lzsize, ref, 
        lxsize, lysize, lzsize, data, FORWARD );
    tfwd[iter] = MPI_Wtime() - ti;

    tmp = MPI_Wtime();
    decomp_2d_fft_3d_c2c( lxsize, lysize, lzsize, data, 
        lxsize, lysize, lzsize, data, BACKWARD );
    tbwd[iter] = MPI_Wtime() - tmp;
    titer[iter] = MPI_Wtime() - ti;
    tn += titer[iter];

    // Normalization
    for (i = 0; i < lxsize; ++i)
      for (j = 0; j < lysize; ++j)
        for (k = 0; k < lzsize; ++k) {
          data[(i*lysize + j) * lzsize + k] *= scaling;
        }
    
    MPI_Barrier ( comm );
  }

  // Reduce operation to compute min, max and avg
  MPI_Reduce ( &tn, &max_tn, 1, MPI_DOUBLE, MPI_MAX, 0, comm );

  for ( int i = 0; i < niter; ++i ) {
    MPI_Reduce ( tfwd + i, max_tfwd + i, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
    MPI_Reduce ( tbwd + i, max_tbwd + i, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
    MPI_Reduce ( tfwd + i, min_tfwd + i, 1, MPI_DOUBLE, MPI_MIN, 0, comm );
    MPI_Reduce ( tbwd + i, min_tbwd + i, 1, MPI_DOUBLE, MPI_MIN, 0, comm );
  }

  // Compute local error
  double err = 0.0;
  double ierr = 0.0;
  double max_ierr = 0.0;
  for (i = 0; i < lxsize; ++i)
    for (j = 0; j < lysize; ++j)
      for (k = 0; k < lzsize; ++k) {
        double complex cref_tmp   = ref[(i*lysize + j) * lzsize + k];
        double complex cdata_tmp  = data[(i*lysize + j) * lzsize + k];
        double tmp_r = creal(cref_tmp) - creal(cdata_tmp);
        double tmp_i = cimag(cref_tmp) - cimag(cdata_tmp);
        err += sqrt(tmp_r*tmp_r + tmp_i * tmp_i);
        double tmp_real = sqrt(tmp_r * tmp_r);
        if ( tmp_real > ierr )
          ierr = tmp_real;
      }

  MPI_Reduce ( &ierr, &max_ierr, 1, MPI_DOUBLE, MPI_MAX, 0, comm );

  // Computation of the flops
  double fwdAdd, fwdMul, fwdFma;
  double bwdAdd, bwdMul, bwdFma;
  double n1 = ((double)nx * ny * nz) / 3.;
  double flops = 5. * n1 * log(n1) / log(2);
  flops *= 3. * n1 * n1;

  if ( ! rank ) {
    printf ( "\n========================\nInput:" );
    for ( int i = 1; i < argc; ++i)
      printf ( " %s", argv[i] );
    printf ( "\n" );
    printf ( "2decomp&FFT test results\n---------\n"
        "Problem size:          %dx%dx%d\n"
        "MPI ranks:             %d\n"
        "Execution(niter=%3d):  %.8fs\n"
        "GFlops:                %.3e\n"
        "Error:                 %.2e\n",
        nx, ny, nz,
        size,
        niter, max_tn,
        ( 2 * flops ) * niter / max_tn,
        max_ierr );

    // Print some details
    printf ( "Fwd details: max_i=[" );
    for ( int i = 0; i < niter; ++i )
      printf ( "%.6e ", max_tfwd[i] );
    printf ( "]\n" );
    printf ( "Fwd details: min_i=[" );
    for ( int i = 0; i < niter; ++i )
      printf ( "%.6e ", min_tfwd[i] );
    printf ( "]\n" );
    printf ( "Bwd details: max_i=[" );
    for ( int i = 0; i < niter; ++i )
      printf ( "%.6e ", max_tbwd[i] );
    printf ( "]\n" );
    printf ( "Bwd details: min_i=[" );
    for ( int i = 0; i < niter; ++i )
      printf ( "%.6e ", min_tbwd[i] );
    printf ( "]\n" );
  }

finalize:

//decomp_2d_fft_finalize();
  decomp_2d_finalize();
#endif

  if ( data ) free ( data );
  if ( ref ) free ( ref );

  if ( titer ) free ( titer );
  if ( tfwd ) free ( tfwd );
  if ( tbwd ) free ( tbwd );
  if ( max_tfwd ) free ( max_tfwd );
  if ( max_tbwd ) free ( max_tbwd );
  if ( min_tfwd ) free ( min_tfwd );
  if ( min_tbwd ) free ( min_tbwd );

  MPI_Comm_free ( &comm );
  MPI_Finalize();

  return 0;
}
