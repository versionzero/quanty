
#ifndef _MEASURE_H_
#define _MEASURE_H_

#include "image.h"

long double ssim(char const *fn1, char const *fn2);

long double
spatial_frequency(const RAW_IMAGE &source);

typedef long double (*quality_measure_function)(const RAW_IMAGE &original, 
						const RAW_IMAGE &compressed);

long double
average_distance(const RAW_IMAGE &original,
		 const RAW_IMAGE &compressed);

long double
structural_content(const RAW_IMAGE &original, 
		   const RAW_IMAGE &compressed);

long double
cross_correlation(const RAW_IMAGE &original, 
		  const RAW_IMAGE &compressed);

long double
correlation_quality(const RAW_IMAGE &original, 
		    const RAW_IMAGE &compressed);

long double
max_difference(const RAW_IMAGE &original, 
	       const RAW_IMAGE &compressed);

long double
image_fidelity(const RAW_IMAGE &original, 
	       const RAW_IMAGE &compressed);

long double
laplacian_mean_square_error(const RAW_IMAGE &original, 
			    const RAW_IMAGE &compressed);

long double
peak_mean_square_error(const RAW_IMAGE &original, 
		       const RAW_IMAGE &compressed);

long double
absolute_error_1(const RAW_IMAGE &original, 
		 const RAW_IMAGE &compressed);

long double
absolute_error_2(const RAW_IMAGE &original, 
		 const RAW_IMAGE &compressed);

long double
absolute_error_3(const RAW_IMAGE &original, 
		 const RAW_IMAGE &compressed);

long double
mean_square_error_1(const RAW_IMAGE &original, 
		    const RAW_IMAGE &compressed);

long double
mean_square_error_2(const RAW_IMAGE &original, 
		    const RAW_IMAGE &compressed);

long double
mean_square_error_3(const RAW_IMAGE &original, 
		    const RAW_IMAGE &compressed);

long double
lp_1_norm(const RAW_IMAGE &original, 
	  const RAW_IMAGE &compressed);

long double
lp_2_norm_1(const RAW_IMAGE &original, 
	    const RAW_IMAGE &compressed);

long double
lp_2_norm_2(const RAW_IMAGE &original, 
	    const RAW_IMAGE &compressed);

long double
lp_2_norm_3(const RAW_IMAGE &original, 
	    const RAW_IMAGE &compressed);

long double
lp_3_norm(const RAW_IMAGE &original, 
	  const RAW_IMAGE &compressed);

#endif

