// Fourier.cpp: implementation of the Fourier class.
//
//////////////////////////////////////////////////////////////////////

#include "Fourier.h"
#include <math.h>
#include <stdlib.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

/*
 * fft.cpp
 *
 * loic fonteneau 15-feb-2001
 * Perform discrete FFT
 *
 * Original code : Don Cross <dcross@intersrv.com>
 * http://www.intersrv.com/~dcross/fft.html
 *
 */

#ifndef NULL
#define NULL '\0'
#endif

// dir=1 forward, -1 inverse
//#################################################
int FFT(int dir,int m,float*x,float *y)
{
   long nn,i,i1,j,k,i2,l,l1,l2;
   float c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   nn = 1;
   for (i=0;i<m;i++)
      nn *= 2;

   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0f;
   c2 = 0.0f;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<nn;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrtf((1.0f - c1) / 2.0f);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrtf((1.0f + c1) / 2.0f);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<nn;i++) {
         x[i] /= (float)nn;
         y[i] /= (float)nn;
      }
   }

   return(true);
}
int Powerof2( int n, int *m, int *twopm ) {
   if ( n <= 1 )
   {
      *m     = 0;
      *twopm = 1;
      return false;
   }

   *m     = 1;
   *twopm = 2;

   do
   {
      (*m)++;
      (*twopm) *= 2;
   } while ( 2 * (*twopm) <= n );

   if ( *twopm != n )
   {
      return false;
   }
   else
   {
      return true;
   }
}

int FFT2D(COMPLEX *c,int nx,int ny,int dir)
{
   int i,j;
   int m,twopm;
   float *real,*imag;

   /* Transform the rows */
   real=(float *)malloc(nx * sizeof(float));
   imag= (float *)malloc(nx * sizeof(float));
   if (real == NULL || imag == NULL)
      return(false);
   if (!Powerof2(nx,&m,&twopm) || twopm != nx)
      return(false);
   for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
		int index=i+nx*j;
         real[i] = c[index].real;
         imag[i] = c[index].imag;
      }
      FFT(dir,m,real,imag);
      for (i=0;i<nx;i++) {
		  int index=i+nx*j;
         c[index].real = real[i];
         c[index].imag = imag[i];
      }
   }
   free(real);
   free(imag);

   /* Transform the columns */
   real=(float *)malloc(ny*sizeof(float));
   imag=(float *)malloc(ny*sizeof(float));
   if (real == NULL || imag == NULL)
      return(false);
   if (!Powerof2(ny,&m,&twopm) || twopm != ny)
      return(false);
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
		  int index=i+nx*j;
         real[j] = c[index].real;
         imag[j] = c[index].imag;
      }
      FFT(dir,m,real,imag);
      for (j=0;j<ny;j++) {
		  int index=i+nx*j;
         c[i+nx*j].real = real[j];
         c[i+nx*j].imag = imag[j];
      }
   }
   free(real);
   free(imag);

   return(true);
}

/*-------------------------------------------------------------------------
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform

     Formula: forward
                  N-1
                  ---
              1   \          - j k 2 pi n / N
      X(n) = ---   >   x(k) e                    = forward transform
              N   /                                n=0..N-1
                  ---
                  k=0

      Formula: reverse
                  N-1
                  ---
                  \          j k 2 pi n / N
      X(n) =       >   x(k) e                    = forward transform
                  /                                n=0..N-1
                  ---
                  k=0
*/

//######################################33

/////////////////////////////////////////////////////////////////////////////////////
// do the fft for double numbers
//////////////////////////////////////////////////////////////////////////////////////
void FFT(unsigned int p_nSamples,bool p_bInverseTransform,double* p_lpRealIn,double*p_lpImagIn,double* p_lpRealOut,double* p_lpImagOut){
	if(!p_lpRealIn || !p_lpRealOut || !p_lpImagOut) return;
	unsigned int NumBits;
	unsigned int i,j,k,n;
	unsigned int BlockSize,BlockEnd;
	double angle_numerator=2.0*PI;
	double tr,ti;
	if(!IsPowerOfTwo(p_nSamples))return;
	if(p_bInverseTransform)angle_numerator=-angle_numerator;
	NumBits=NumberOfBitsNeeded(p_nSamples);
	for(i=0;i<p_nSamples;i++){
		j=ReverseBits(i,NumBits);
		p_lpRealOut[j]=p_lpRealIn[i];
		p_lpImagOut[j]=(p_lpImagIn==NULL)?0.0:p_lpImagIn[i];}

	BlockEnd=1;
	for(BlockSize=2;BlockSize<=p_nSamples;BlockSize<<=1){
		double delta_angle=angle_numerator/(double)BlockSize;
		double sm2=sin(-2*delta_angle);
		double sm1=sin(-delta_angle);
		double cm2=cos(-2*delta_angle);
		double cm1=cos(-delta_angle);
		double w=2*cm1;
		double ar[3],ai[3];

		for(i=0;i<p_nSamples;i+=BlockSize){
			ar[2]=cm2;
			ar[1]=cm1;
			ai[2]=sm2;
			ai[1]=sm1;
			for(j=i,n=0;n<BlockEnd;j++,n++){
				ar[0]=w*ar[1]-ar[2];
				ar[2]=ar[1];
				ar[1]=ar[0];
				ai[0]=w*ai[1]-ai[2];
				ai[2]=ai[1];
				ai[1]=ai[0];
				k=j+BlockEnd;
				tr=ar[0]*p_lpRealOut[k]-ai[0]*p_lpImagOut[k];
				ti=ar[0]*p_lpImagOut[k]+ai[0]*p_lpRealOut[k];
				p_lpRealOut[k]=p_lpRealOut[j]-tr;
				p_lpImagOut[k]=p_lpImagOut[j]-ti;
				p_lpRealOut[j]+=tr;
				p_lpImagOut[j]+=ti;}}
		BlockEnd = BlockSize;}

	if(p_bInverseTransform){
		double denom=(double)p_nSamples;
		for(i=0;i<p_nSamples;i++){
			p_lpRealOut[i]/=denom;
			p_lpImagOut[i]/=denom;}}
}

void FFT(unsigned int p_nSamples,bool p_bInverseTransform,float* p_lpRealIn,float*p_lpImagIn,float* p_lpRealOut,float* p_lpImagOut){
	if(!p_lpRealIn || !p_lpRealOut || !p_lpImagOut) return;
	unsigned int NumBits;
	unsigned int i,j,k,n;
	unsigned int BlockSize,BlockEnd;
	float angle_numerator=2.0f*(float)PI;
	float tr,ti;
	if(!IsPowerOfTwo(p_nSamples))return;
	if(p_bInverseTransform)angle_numerator=-angle_numerator;
	NumBits=NumberOfBitsNeeded(p_nSamples);
	for(i=0;i<p_nSamples;i++){
		j=ReverseBits(i,NumBits);
		p_lpRealOut[j]=p_lpRealIn[i];
		p_lpImagOut[j]=(p_lpImagIn==NULL)?0.0f:p_lpImagIn[i];}

	BlockEnd=1;
	for(BlockSize=2;BlockSize<=p_nSamples;BlockSize<<=1){
		float delta_angle=angle_numerator/(float)BlockSize;
		float sm2=sinf(-2*delta_angle);
		float sm1=sinf(-delta_angle);
		float cm2=cosf(-2*delta_angle);
		float cm1=cosf(-delta_angle);
		float w=2*cm1;
		float ar[3],ai[3];

		for(i=0;i<p_nSamples;i+=BlockSize){
			ar[2]=cm2;
			ar[1]=cm1;
			ai[2]=sm2;
			ai[1]=sm1;
			for(j=i,n=0;n<BlockEnd;j++,n++){
				ar[0]=w*ar[1]-ar[2];
				ar[2]=ar[1];
				ar[1]=ar[0];
				ai[0]=w*ai[1]-ai[2];
				ai[2]=ai[1];
				ai[1]=ai[0];
				k=j+BlockEnd;
				tr=ar[0]*p_lpRealOut[k]-ai[0]*p_lpImagOut[k];
				ti=ar[0]*p_lpImagOut[k]+ai[0]*p_lpRealOut[k];
				p_lpRealOut[k]=p_lpRealOut[j]-tr;
				p_lpImagOut[k]=p_lpImagOut[j]-ti;
				p_lpRealOut[j]+=tr;
				p_lpImagOut[j]+=ti;}}
		BlockEnd = BlockSize;}

	if(p_bInverseTransform){
		float denom=(float)p_nSamples;
		for(i=0;i<p_nSamples;i++){
			p_lpRealOut[i]/=denom;
			p_lpImagOut[i]/=denom;}}
}

//////////////////////////////////////////////////////////////////////////////////////
// check is a number is a power of 2
//////////////////////////////////////////////////////////////////////////////////////

bool IsPowerOfTwo(unsigned int p_nX){
	if(p_nX<2)return false;
	if(p_nX&(p_nX-1))return false;
    return true;
}
//////////////////////////////////////////////////////////////////////////////////////
// return needed bits for fft
//////////////////////////////////////////////////////////////////////////////////////
unsigned int NumberOfBitsNeeded(unsigned int p_nSamples){
	int i;
	if(p_nSamples<2)return 0;
	for (i=0;;i++){if(p_nSamples&(1<<i))return i;}
}

//////////////////////////////////////////////////////////////////////////////////////
// ?
//////////////////////////////////////////////////////////////////////////////////////
unsigned int ReverseBits(unsigned int p_nIndex,unsigned int p_nBits){
	unsigned int i, rev;
	for(i=rev=0;i<p_nBits;i++){
		rev=(rev<<1)|(p_nIndex&1);
		p_nIndex>>=1;}
	return rev;
}
//////////////////////////////////////////////////////////////////////////////////////
// return a frequency from the basefreq and num of samples
//////////////////////////////////////////////////////////////////////////////////////
float Index_to_frequency(unsigned int p_nBaseFreq,unsigned int p_nSamples,unsigned int p_nIndex){
	if(p_nIndex >= p_nSamples){return 0.0;}
	else if(p_nIndex <= p_nSamples/2){return ((float)p_nIndex/(float)p_nSamples*p_nBaseFreq );}
	else{return ( -(float)(p_nSamples-p_nIndex) / (float)p_nSamples * p_nBaseFreq );}
}
