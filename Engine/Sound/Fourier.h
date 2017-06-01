// Fourier.h: interface for the Fourier class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

/*
 * fft.h
 *
 * loic fonteneau 15-feb-2001
 * Perform discrete FFT
 *
 * Original code : Don Cross <dcross@intersrv.com>
 * http://www.intersrv.com/~dcross/fft.html
 *
 */

typedef struct COMPLEX{
	float real;
	float imag;
}COMPLEX;

void FFT(unsigned int p_nSamples,bool p_bInverseTransform,float* p_lpRealIn,float* p_lpImagIn,float* p_lpRealOut,float* p_lpImagOut);

void FFT(unsigned int p_nSamples,bool p_bInverseTransform,double* p_lpRealIn,double* p_lpImagIn,double* p_lpRealOut,double* p_lpImagOut);

int FFT2D(COMPLEX *c,int nx,int ny,int dir);

bool IsPowerOfTwo(unsigned int p_nX);

unsigned int NumberOfBitsNeeded(unsigned int p_nSamples);

unsigned int ReverseBits(unsigned int p_nIndex, unsigned int p_nBits);

float Index_to_frequency(unsigned int p_nBaseFreq, unsigned int p_nSamples, unsigned int p_nIndex);
/////////////////////////////
//  constantes-definition  //
/////////////////////////////

#define  PI  (3.14159265358979323846)

/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems or
      the dimensions are not powers of 2
*/

