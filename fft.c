/*
Copyright 2012 Jorge Velazquez
*/
#include <fftw3.h>
#include <stdio.h>
#include "libPP_5.0.h"


void CFFT(estado *es, Float2D_MP *correlacion)
{
	printf("Entra a fft\n");
fftw_complex *out;
fftw_plan p, plan2;
int NDX = es->NDX;
int NDY = es->NDY;
int **s = es->s;
double *in;
int i,j;

 in = fftw_alloc_real( sizeof ( double ) * NDX * NDY );
//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NDX * NDY);
int nyh = ( NDY / 2 ) + 1;

out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NDX * nyh);

p = fftw_plan_dft_r2c_2d ( NDX, NDY, in, out, FFTW_ESTIMATE );

	printf("Termina primer plan\n");
 
 for ( i = 0; i < NDX; i++ )
  {
    for ( j = 0; j < NDY; j++ )
    {
     // in[i*NDY+j][0] = s[i+1][j+1];
      //in[i*NDY+j][1] = 0;
      in[i*NDY + j] = (double)s[i+1][j+1];
    }
  }

	printf("Comienza FFT\n");
fftw_execute(p); /* repeat as needed */
	printf("Termina FFT\n");
 for ( i = 0; i < NDX; i++ )
  {
    for ( j = 0; j < nyh; j++ )
    {
      out[i*nyh + j][0] = out[i*nyh + j][0]*out[i*nyh + j][0] + out[i*nyh + j][1]*out[i*nyh + j][1];
      out[i*nyh + j][1] = 0;
    }
  }
		printf("2plan\n");
plan2 = fftw_plan_dft_c2r_2d ( NDX, NDY, out, in, FFTW_ESTIMATE );
	printf("2plan listo\n");
printf("Comienza FFT 2\n");
 fftw_execute ( plan2 );
 printf("Termina FFT 2\n");
 for ( i = 0; i < NDX; i++ )
  {
    for ( j = 0; j < NDY; j++ )
    {
      correlacion->array[i+1][j+1]=in[i*NDY + j];
    }
  }
 
fftw_destroy_plan(p);
fftw_destroy_plan(plan2);
fftw_free(in);
fftw_free(out);
}

