/*
 * A program that reads a PGM file and compresses it with the given
 * JPEG parameters.
 *
 * Usage:
 *
 *   subjectivity input.pgm start end
 *
 * Ben Bunett
 *
 */

#include "image.h"
#include "measure.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;


long double
spatial_frequency(const RAW_IMAGE &source)
{
  unsigned int  i, j, offset, width, height;
  int           lhs;
  long double   value, row_frequency, column_frequency;
  unsigned char *buffer;

  /* shorten some names for convenience */
  width  = source.width();
  height = source.height();
  buffer = source.buffer;

  /* calculate the row frequency */
  row_frequency = 0.0;
  for (i = 0; i < height; ++i) {
    offset = i*width;
    for (j = 1; j < width; ++j) {
      lhs            = static_cast<int>(buffer[j+offset]);
      value          = lhs-buffer[j+offset-1];
      row_frequency += value*value;
    }
  }
  row_frequency *= 1.0/(width*height);
  row_frequency  = sqrt(row_frequency);
  row_frequency *= row_frequency;

  /* calculate the column frequency */
  column_frequency = 0.0;
  for (i = 1; i < height; ++i) {
    offset = i*width;
    for (j = 0; j < width; ++j) {
      lhs               = static_cast<int>(buffer[j+offset]);
      value             = lhs-buffer[j+offset-width];
      column_frequency += value*value;
    }
  }
  column_frequency *= 1.0/(width*height);
  column_frequency  = sqrt(column_frequency);
  column_frequency *= column_frequency;

  /* calculate the spatial frequency */
  value = sqrt(row_frequency+column_frequency);

  return value;
}

long double
average_distance(const RAW_IMAGE &original,
		 const RAW_IMAGE &compressed)
{
  unsigned int i, size;
  long double  value;

  /* shorten some names for convenience */
  size = original.width()*original.height();

  /* calculate the average distance */
  value = 0;
  for (i = 0; i < size; ++i) {
    value += abs(original.buffer[i]-compressed.buffer[i]);
  }

  return value/size;
}

long double
structural_content(const RAW_IMAGE &original, 
		   const RAW_IMAGE &compressed)
{
  unsigned int i, size;
  long double  u, v;

  /* shorten some names for convenience */
  size = original.width()*original.height();

  /* calculate the structural content */
  u = v = 0.0;
  for (i = 0; i < size; ++i) {
    u += original.buffer[i]*original.buffer[i];
    v += compressed.buffer[i]*compressed.buffer[i];
  }

  return u/v;
}

long double
cross_correlation(const RAW_IMAGE &original, 
		  const RAW_IMAGE &compressed)
{
  unsigned int i, size;
  long double  u, v;

  /* shorten some names for convenience */
  size = original.width()*original.height();

  /* calculate the cross correlation */
  u = v = 0.0;
  for (i = 0; i < size; ++i) {
    u += original.buffer[i]*compressed.buffer[i];
    v += original.buffer[i]*original.buffer[i];
  }

  return u/v;
}

long double
correlation_quality(const RAW_IMAGE &original, 
		    const RAW_IMAGE &compressed)
{
  unsigned int i, size;
  long double  u, v;

  /* shorten some names for convenience */
  size = original.width()*original.height();

  /* calculate the correlation quality */
  u = v = 0.0;
  for (i = 0; i < size; ++i) {
    u += original.buffer[i]*compressed.buffer[i];
    v += original.buffer[i];
  }

  return u/v;
}

long double
max_difference(const RAW_IMAGE &original, 
	       const RAW_IMAGE &compressed)
{
  unsigned int i, size;
  long double  diff, value;

  /* shorten some names for convenience */
  size = original.width()*original.height();

  /* calculate the maximum difference */
  diff = 0.0;
  for (i = 0; i < size; ++i) {
    diff  = abs(original.buffer[i]-compressed.buffer[i]);
    value = max(value, diff);
  }

  return value;
}

long double
image_fidelity(const RAW_IMAGE &original,
	       const RAW_IMAGE &compressed)
{
  unsigned int i, size;
  long double  u, v, value;

  /* shorten some names for convenience */
  size = original.width()*original.height();

  /* calculate the image fidelity */
  u = v = 0.0;
  for (i = 0; i < size; ++i) {
    value = original.buffer[i]-compressed.buffer[i];
    u    += value*value;
    v    += compressed.buffer[i]*compressed.buffer[i];
  }
  
  return 1-(u/v);
}

long double
laplacian(const RAW_IMAGE &original, int i)
{
  int          size, width;
  long double  u;

  /* shorten some names for convenience */
  width  = original.width();
  size   = width*original.height();

  /* calculate the laplacian */
  u = 0.0;
  if (i+1 < size) {
    u += original.buffer[i+1];
  } else if (i-1 >= 0) {
    u += original.buffer[i-1];
  } else if (i-width >= 0) {
    u += original.buffer[i-width];
  } else if (i+width < size) {
    u += original.buffer[i+width];
  }
  u -= 4*original.buffer[i];

  return u;
}

long double
laplacian_sum(const RAW_IMAGE &original)
{
  int          i, size;
  long double  u, tmp;

  /* shorten some names for convenience */
  size = original.width()*original.height();

  /* calculate the laplacian */
  u = 0.0;
  for (i = 0; i < size; ++i) {
    tmp    = laplacian(original, i);
    tmp   *= tmp;
    u += tmp;
  }

  return u;
}

long double
laplacian_diff(const RAW_IMAGE &original,
	       const RAW_IMAGE &compressed)
{
  int         i, size;
  long double u, diff;

  /* shorten some names for convenience */
  size = original.width()*original.height();

  /* calculate the laplacian */
  u = 0.0;
  for (i = 0; i < size; ++i) {
    diff  = laplacian(original, i)-laplacian(compressed, i);
    diff *= diff;
    u    += diff;
  }

  return u;
}

long double
laplacian_mean_square_error(const RAW_IMAGE &original, 
			    const RAW_IMAGE &compressed)
{
  long double u, v;

  /* calculate the leplatian mean square error */
  u = laplacian_diff(original, compressed);
  v = laplacian_sum(original);

  return u/v;
}

long double
peak_mean_square_error(const RAW_IMAGE &original, 
		       const RAW_IMAGE &compressed)
{
  unsigned char maximum;
  unsigned int  i, size;
  long double   diff, u;

  /* shorten some names for convenience */
  size = original.width()*original.height();

  /* find the maximum pixel value in the original image */
  maximum = 0;
  for (i = 0; i < size; ++i) {
    maximum = max(maximum, original.buffer[i]);
  }
  maximum *= maximum;

  /* calculate the peak mean square error */
  u = 0;
  for (i = 0; i < size; ++i) {
    diff   = original.buffer[i]-compressed.buffer[i];
    u += (diff*diff)/maximum;
  }

  return u/size;
}

/* http://en.wikipedia.org/wiki/Peak_signal-to-noise_ratio */
long double
peak_signal_to_noise_ratio(const RAW_IMAGE &original, 
			   const RAW_IMAGE &compressed)
{
  unsigned char maximum;
  unsigned int  i, size;

  /* shorten some names for convenience */
  size = original.width()*original.height();

  /* find the maximum pixel value in the original image */
  maximum = 0;
  for (i = 0; i < size; ++i) {
    maximum = max(maximum, original.buffer[i]);
  }
  maximum *= maximum;

  return 10*log(maximum/peak_mean_square_error(original, compressed));
}

long double
identity(unsigned char c)
{
  return static_cast<long double>(c);
}

long double
cube_root(unsigned char c)
{
  long double d = c;

  return cbrt(d);
}

/* adpated from "Hacker's Delight", Addison-Wesley, 2003 */
template <class T>
T next_power_of_two(T k) 
{
  unsigned int i;

  if (0 == k) {
    return 1;
  }

  k--;
  for (i = 1; i < sizeof(T)*CHAR_BIT; i <<= 1) {
    k = k | k >> i;
  }

  return ++k;
}

long double
hvs_measure(const RAW_IMAGE &original, 
	    const RAW_IMAGE &compressed)
{
#if 0
  unsigned int i, j, 
               width, height, offset, 
               side, size;
  double       *input, *output;
  //fftw_complex  *output;
  fftw_plan    plan;


  /* shorten some names for convenience */
  width  = original.width();
  height = original.height();
  side   = 2*next_power_of_two(max(width, height));
  size   = side*side;

  /* allocate the buffers */
  input  = static_cast<double*>(fftw_malloc(size*sizeof(double)));
  output = static_cast<double*>(fftw_malloc(size*sizeof(double)));

  /* create an FFTW plan for building the DCT (FFTW_REDFT10). To
     compute IDCT, use FFTW_REDFT01 instead. */
  plan = fftw_plan_r2r_2d(side, side, input, output, 
			  FFTW_REDFT10, FFTW_REDFT10,
			  FFTW_MEASURE|FFTW_DESTROY_INPUT);
  
  /* populate input buffer */
  for (i = 0; i < size; ++i) {
    input[i] = 0.0;
  }
  for (i = 0; i < width; ++i) {
    offset = i*side;
    for (j = 0; j < height; ++j) {
      input[j+offset] = compressed.buffer[j+offset];
    }
  }

  /* build the DCT */
  fftw_execute(plan);

  
  for (i = 0; i < size; ++i) {
    cout << output[i];
    if (0 == i%width) {
      cout << endl;
      i += side-width;
    } else {
      cout << ", ";
    }
  }

  /* clean up */
  fftw_destroy_plan(plan);
  fftw_free(input);
  fftw_free(output);



  fftw_plan plan;

  size_t allocSize = windowSize_ * sizeof(float);
  float *fltReturnArr = (float *)fftwf_malloc(allocSize);
  
  //FFT Optimization Options: FFTW_ESTIMATE,FFTW_PATIENT, FFTW_MEASURE
  ffpPlan = fftwf_plan_r2r_1d(windowSize_, pcmVals, fltReturnArr, FFTW_R2H
			      C, FFTW_ESTIMATE);     
  fftwf_execute_r2r(ffpPlan, pcmVals, fltReturnArr);
  
  vector<double> fftVec(windowSize_);
  for(int i = 0; i < windowSize_; i++){ //because fftw requires arrays not
    vectors... 
      fftVec[i] = static_cast<double>(fltReturnArr[i]);
  }
  
  fftwf_destroy_plan(ffpPlan);
  fftwf_free(fltReturnArr);
  return fftVec;
#endif
  return 0.0;
}

template<class Func>
long double
absolute_error(const RAW_IMAGE &original, 
	       const RAW_IMAGE &compressed,
	       Func            preprocess)
{
  unsigned int i, size;
  long double  u, v;

  /* shorten some names for convenience */
  size = original.width()*original.height();
  
  /* calculate the absolute error */
  u = v = 0.0;
  for (i = 0; i < size; ++i) {
    u += abs(preprocess(original.buffer[i])
	     -preprocess(compressed.buffer[i]));
    v += abs(preprocess(original.buffer[i]));
  }

  return u/v;
}

long double
absolute_error_1(const RAW_IMAGE &original, 
		 const RAW_IMAGE &compressed)
{
  return absolute_error(original, compressed, identity);
}

long double
absolute_error_2(const RAW_IMAGE &original, 
		 const RAW_IMAGE &compressed)
{
  return absolute_error(original, compressed, cube_root);
}

long double
absolute_error_3(const RAW_IMAGE &original, 
		 const RAW_IMAGE &compressed)
{
  return absolute_error(original, compressed, identity);
}


template<class Func>
long double
mean_square_error(const RAW_IMAGE &original, 
		  const RAW_IMAGE &compressed,
		  Func            preprocess)
{
  unsigned int i, size;
  long double  u, v, value;

  /* shorten some names for convenience */
  size = original.width()*original.height();

  /* calculate the image mean square error */
  u = v = 0.0;
  for (i = 0; i < size; ++i) {
    value  = preprocess(original.buffer[i]);
    v     += value*value;
    value -= preprocess(compressed.buffer[i]);
    u     += value*value;
  }
  
  return u/v;
}

long double
mean_square_error_1(const RAW_IMAGE &original, 
		    const RAW_IMAGE &compressed)
{
  return mean_square_error(original, compressed, identity);
}

long double
mean_square_error_2(const RAW_IMAGE &original, 
		    const RAW_IMAGE &compressed)
{
  return mean_square_error(original, compressed, cube_root);
}

long double
mean_square_error_3(const RAW_IMAGE &original, 
		    const RAW_IMAGE &compressed)
{
  return mean_square_error(original, compressed, identity);
}

template<class Func>
long double
lp_norm(const RAW_IMAGE    &original, 
	const RAW_IMAGE    &compressed,
	const unsigned int p,
	Func               preprocess)
{
  unsigned int i, j, size;
  long double  u, value;

  /* shorten some names for convenience */
  size = original.width()*original.height();

  /* calculate the image mean square error */
  u = 0.0;
  for (i = 0; i < size; ++i) {
    value = abs(preprocess(original.buffer[i])
		-preprocess(compressed.buffer[i]));
    for (j = 1; j < p; ++j) {
      value *= value;
    }
    u += value;
  }
  u /= size;
  switch (p) {
  case 2:
    u = sqrt(u);
    break;
  case 3:
    u = cbrt(u);
    break;
  }
  u /= 255.0; /* scale */
  
  return u;
}

long double
lp_1_norm(const RAW_IMAGE &original, 
	  const RAW_IMAGE &compressed)
{
  return lp_norm(original, compressed, 1, identity);
}

long double
lp_2_norm_1(const RAW_IMAGE &original, 
	    const RAW_IMAGE &compressed)
{
  return lp_norm(original, compressed, 2, identity);
}

long double
lp_2_norm_2(const RAW_IMAGE &original, 
	    const RAW_IMAGE &compressed)
{
  return lp_norm(original, compressed, 2, cube_root);
}

long double
lp_2_norm_3(const RAW_IMAGE &original, 
	    const RAW_IMAGE &compressed)
{
  return lp_norm(original, compressed, 2, identity);
}

long double
lp_3_norm(const RAW_IMAGE &original, 
	  const RAW_IMAGE &compressed)
{
  return lp_norm(original, compressed, 3, identity);
}

