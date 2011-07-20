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
#include "macro.h"
#include "measure.h"
#include "random.h"

#include <fftw3.h>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

/* functions */
void read_image(char *filename, PGM_IMAGE &source);
int  random_between(int min, int max);
void generate_table(unsigned int table[8][8], unsigned int seed);
void compress_image(JPEG_IMAGE &destination, const RAW_IMAGE &source,
		    const unsigned int table[8][8], int quality, 
		    bool standard = false);
void decompress_image(RAW_IMAGE &destination, const JPEG_IMAGE &source);
void write_image(char *filename, const RAW_IMAGE &source, 
		 const unsigned int table[8][8], int quality, 
		 bool standard = false, bool compress = true);
long double spatial_frequency(const RAW_IMAGE &source);
void run_metrics(const RAW_IMAGE &original, const RAW_IMAGE &compressed);
void usage();

/* the program's name */
const char *program_name = NULL;

int
main(int argc, char *argv[])
{
  unsigned int i, j;
  unsigned int start, end;
  int          desired_quality, actual_quality, best_quality;
  unsigned int raw_compressed_size;
  int          previous_direction;
  unsigned int table[8][8];
  double       compression_ratio, best_compression_ratio;
  char         file_name[512], *base_name, *dot;
  bool         store;
  PGM_IMAGE    pgm_image;
  RAW_IMAGE    raw_image, decompressed;
  JPEG_IMAGE   compressed;  

  /* save the program name */
  program_name = argv[0];

  /* check arguments */
  if (argc < 4) {
    usage();
  }

  /* determine the range of quantization tables we will be using */
  start           = atoi(argv[2]);
  end             = atoi(argv[3]);

  /* did the user supply a quality scalar? */
  desired_quality = (argc >= 5) ? atoi(argv[4]) : 80;

  /* are we storing the image represenations we are creating? */
  store           = (argc >= 6) ? atoi(argv[5]) != 0 : false;

  /* initializes the (netpbm) library */
  pm_init(argv[0], 0);

  /* read the image */
  read_image(argv[1], pgm_image);

  /* convert the pgm image to a "raw" version */
  pgm2raw(raw_image, pgm_image);  

  /* determine the appropriate base_name for the name of the output
     files */
  base_name = strrchr(argv[1], '/');
  if (NULL != base_name && '/' == *base_name) {
    base_name++;
  } else {
    base_name = argv[1];
  }
  if (NULL != (dot = strrchr(argv[1], '.'))) {
    *dot = '\0';
  }

  /* save a copy of the raw image */
  sprintf(file_name, "%s.raw.jpeg", base_name);
  if (store) {
    write_image(file_name, raw_image, table, desired_quality, true, false);
  }

  /* create the destination buffer */
  compressed.create(pgm_image.width(), pgm_image.height());
  decompressed.create(pgm_image.width(), pgm_image.height());

  /* print some random numbers to convice us that are number generator
     truly is hardware independant */
  cout << "# ";
  seed_random(6+28+496+8128); /* perfect :) */
  for ( i = 0; i < 10; ++i ) {
    cout << random_between(0, 100) << " ";
  }
  cout << endl;

  /* print a header describing the file format */
  cout << "# seed, quality, spatial_frequency, original_size, "
    "compressed_size, compression_ratio, average_distance, "
    "structural_content, cross_correlation, image_fidelity, "
    "laplacian_mean_square_error, peak_mean_square_error, "
    "absolute_error_1, absolute_error_2, mean_square_error_1, "
    "mean_square_error_2, lp_1_norm, lp_2_norm_1, lp_2_norm_2, "
    "lp_3_norm" << endl;
  
  /* print out the spatial frequency for the image */
  cout << fixed << setprecision(10);
  cout << "# spatial frequency for " << base_name << " = " 
       << spatial_frequency(raw_image) << endl;

  /* save a copy of the image using the standard JPEG quantization
     table */
  if (store) {
    sprintf(file_name, "%s.%d.raw.jpeg", base_name, desired_quality);
    write_image(file_name, raw_image, table, desired_quality, true, true);    
  }

  /* we would like to know how the standard quantization matrix
     performs with respect to bit rate so that we can later pick a
     scalr value to use with the non-standard matricies such that the
     bit rate of the images compressed with the non-standard matrices
     closely matched that of the ones compressed with the standard
     matrix. */    
  compress_image(compressed, raw_image, table, desired_quality, true);
  raw_compressed_size = compressed.size();
  
  /* now decompress the image */
  decompress_image(decompressed, compressed);  
  
  /* start printing the output file: seed, spatial frequency,
     original image size, compressed image size */
  cout << "0 " 
       << desired_quality << " "
       << spatial_frequency(decompressed) << " "
       << raw_compressed_size << " " 
       << compressed.size() << " "
       << 1.0 << " ";
  
  /* run all of our quality metrics */
  run_metrics(raw_image, decompressed);

  /* adjust range numbers such that we don't collide with the 0 entry
     printed above */

  /* generate all quantization tables and try compressing the input
     image with them */  
  for (i = start; i < end; ++i) {

    /* create the quantization table */
    generate_table(table, i);
    
    /* compress the image - note that we want to pick a JPEG quality
       parameter here that create an image using Q_i such that it is
       close in bit rate to the image compressed with the standard
       quantization matrix (we stop after looping 10 times, or if we
       adjust the "actualy quality" in both directions, since it mean
       that we were not able to find a satisfactory integer value) */
    actual_quality         = desired_quality;
    best_quality           = desired_quality;
    best_compression_ratio = 100.0; /* unlikely that this will ever be
				       to considered be the best */
    previous_direction     = 0;     /* track the direction in which we
				       are taking the qulity measure,
				       if we ever switch direction,
				       then stop looking */
    for ( j = 1; actual_quality > 0 && actual_quality <= 100; ++j ) {
      
      /* compress the image to see if we come close to the same bit
	 rate attained by the JPEG standard quantization matrix */
      compress_image(compressed, raw_image, table, actual_quality);
      compression_ratio  = compressed.size();
      compression_ratio /= raw_compressed_size;
      
      cout << "# " << setw(6) << setfill(' ') 
	   << j << " : " << setw(2) << actual_quality
	   << " --> " << compressed.size() << " (" 
	   << compression_ratio << ") ";
      if (0 == previous_direction) {
	cout << "first";
      } else if (1 == previous_direction) {
	cout << "increased";
      } else {
	cout << "decreased";
      }
      cout << " quality." << endl;

      /* if we do find a "better" compression ratio, then remember it,
	 just in case we don't find a "good" one that meets the
	 criteria bellow. */
      if (abs(1.0-compression_ratio) >
	  abs(1.0-best_compression_ratio)) {
	best_compression_ratio = compression_ratio;
	best_quality           = actual_quality;
      }

      /* determine if the compression ratio we have attained is as
	 "good" as we can make it (~5% larger or smaller); otherwise,
	 increase or decrease the quality measure and try comressing
	 again. */
      if (compression_ratio > 1.05) {
	if (1 == previous_direction) {
	  break;
	}
	previous_direction = -1; /* subtracted */
	actual_quality--;
      } else if (compression_ratio < 0.95) {
	if (-1 == previous_direction) {
	  break;
	}
	previous_direction = 1; /* added */
	actual_quality++;
      } else {
	best_quality = actual_quality;
	break;
      }

    }

    /* now decompress the image */
    decompress_image(decompressed, compressed);

    /* save a copy of the image */
    if (store) {
      sprintf(file_name, "%s.%d.%d.jpeg", base_name, best_quality, i);
      write_image(file_name, raw_image, table, best_quality);
    }

    /* start printing the output file: seed, spatial frequency,
       original image size, compressed image size */
    cout << i+1 << " " 
	 << actual_quality << " "
	 << spatial_frequency(decompressed) << " "
	 << raw_compressed_size << " " 
	 << compressed.size() << " "
	 << compression_ratio << " ";
    
    /* run all of our quality metrics */
    run_metrics(raw_image, decompressed);
    
  }

  /* clean-up */
  fftw_cleanup();

  return 0;
}

/* prints the usage information for the application */
void
usage()
{
  cerr << "usage: " << program_name 
       << " input.pgm start end [quality] [save]" << endl;
  cerr << "  input.pgm    - input image." << endl;
  cerr << "  [start, end) - define the range of quantization tables "
       << "to use." << endl;
  cerr << "  quality      - (desired) scaling factor for the quantization table"
       << endl;
  cerr << "  save         - save a copy of the compressed image "
       << "(default: false)" << endl;
  exit(1);
}

/* reads the image into the netpbm structure */
void
read_image(char *file_name, PGM_IMAGE &source)
{  
  FILE  *f;
  tuple **A;
  pam   *I;

  if ((f = pm_openr(file_name)) == NULL) {
    cerr << "Cannot open file \"" << file_name << "\" for reading." << endl;
    exit(1);
  }

  I = &source.information;
  if ((A = pnm_readpam(f, I, PAM_STRUCT_SIZE(tuple_type))) == NULL) {
    cerr << "Cannot read image \"" << file_name << "\"." << endl;
    exit(1);
  }
  
  source.buffer = A;
  pm_close(f);
}

/* generate a random number in the range [min, max] */
int 
random_between(int min, int max)
{
  return marsaglia_random() % (max-min+1) + min;
}

/* create a quantization table from a given seed number */
void
generate_table(unsigned int table[8][8], unsigned int seed)
{
  int          i, r, c;
  unsigned int base[] = { 16, 11, 10, 16, 24, 40, 51, 61,
			  55, 56, 63, 77, 92, 101, 99 };
  int          value;

  /* seed the random number generator based on the seed number */
  seed_random(seed);

  /* use the seed number to generate modify our base numbers */
  for (i = 14; i >= 0; --i) {
    if (((seed>>(14-i)) & 1) == 1) {
      base[i] += random_between(-5,5);
    }
  }

  /* fill the table boarder with the base numbers. These were taken
     from the original JPEG standard. */
  for (i = 0; i < 8; ++i) {
    table[0][i] = base[i];   /* horizontal */
    table[i][7] = base[i+7]; /* vertical */
  }

  /* fill the interior of the table. We assume the table is symmetric
     in the same sense of a symmetric matrix. We fill the interior in
     by decending down the 15 right-to-left diagonals. */
  for (r = 1; r < 8; ++r) {
    for (c = 0; c < 7; ++c) {
      value = table[r-1][c+1];
      if (r != c+1) {
	value = table[r-1][c+1] + random_between(-5,5);
	if (value < 0) {
	  value = 0;
	}
	if (value > 255) {
	  value = 255;
	}
      }
      table[r][c] = value;
    }
  }

#if 0
  cerr << "quantization table #" << seed << ": " << endl;
  for (i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      cerr << setw(3) << setfill(' ') << table[i][j] << " & ";
    }
    cerr << " \\\\" << endl;
  }
  cerr << endl << endl;
#endif

}

/* compress an image */
void
compress_image(JPEG_IMAGE         &destination, 
	       const RAW_IMAGE    &source,
	       const unsigned int table[8][8],
	       int                quality,
	       bool               standard)
{
  jpeg_compress_struct *jpeg;
  jpeg_error_mgr       jerr;
  JSAMPROW             row[1];
  unsigned int         stride, height, width;
  const unsigned int   *p;

  /* shorten some names for convenience */
  stride = source.width();
  width  = source.width();
  height = source.height();

  /* use the default error handler */
  jpeg      = &destination.information;
  jpeg->err = jpeg_std_error(&jerr);
  jpeg_create_compress(jpeg);

  /* set the destination to be a memory buffer */
  jpeg_mem_dest(jpeg, &destination.buffer, &destination._size);

  /* set some basic JPEG parameters */  
  jpeg->image_width      = width;
  jpeg->image_height     = height;
  jpeg->input_components = 1;
  jpeg->in_color_space   = JCS_GRAYSCALE;

  /* set compression parameters (i.e. no compression) */
  jpeg_set_defaults(jpeg);  

  /* set the quantization table */
  if (standard) {
    /* use the IJG-recommended quantization table */
    jpeg_set_quality(jpeg, quality, TRUE);
  } else {      
    /* use our generated quantization table */
    quality = jpeg_quality_scaling(quality);
    p = reinterpret_cast<const unsigned int*>(table);  
    jpeg_add_quant_table(jpeg, 0, p, quality, TRUE);
  }


  /* compress the JPEG image */
  jpeg_start_compress(jpeg, TRUE);
  while (jpeg->next_scanline < height) {
    row[0] = &source.buffer[jpeg->next_scanline*stride];
    jpeg_write_scanlines(jpeg, row, 1);
  }

  /* clean up */
  jpeg_finish_compress(jpeg);
  jpeg_destroy_compress(jpeg);
}

void
decompress_image(RAW_IMAGE &destination, const JPEG_IMAGE &source)
{
  jpeg_decompress_struct dinfo;
  jpeg_decompress_struct *jpeg;
  jpeg_error_mgr         jerr;
  JSAMPROW               row[1];
  unsigned int           stride, height, width;

  /* shorten some names for convenience */
  stride = source.width();
  width  = source.width();
  height = source.height();

  /* use the default error handler */
  jpeg      = &dinfo;
  jpeg->err = jpeg_std_error(&jerr);
  jpeg_create_decompress(jpeg);

  /* set the source to be a memory buffer */
  jpeg_mem_src(jpeg, source.buffer, source.size());

  /* read  JPEG parameters */  
  jpeg_read_header(jpeg, TRUE);

  /* decompress the JPEG image */
  jpeg_start_decompress(jpeg);
  while (jpeg->output_scanline < height) {
    row[0] = &destination.buffer[jpeg->output_scanline*stride];
    jpeg_read_scanlines(jpeg, row, 1);
  }

  /* clean up */
  jpeg_finish_decompress(jpeg);
  jpeg_destroy_decompress(jpeg);
}

/* writes the image to the given file */
void 
write_image(char               *file_name, 
	    const RAW_IMAGE    &source, 
	    const unsigned int table[8][8], 
	    int                quality,
	    bool               standard,
	    bool               compress)
{
  FILE                 *file;
  jpeg_compress_struct information, *jpeg;
  jpeg_error_mgr       jerr;
  JSAMPROW             row[1];
  unsigned int         stride, height, width;
  const unsigned int   *p;

  /* shorten some names for convenience */
  stride = source.width();
  width  = source.width();
  height = source.height();
  
  /* open the output file */
  if ((file = fopen(file_name, "wb")) == NULL) {
    cerr << "Cannot open file \"" << file_name << "\" for writing." << endl;
    exit(1);
  }

  /* use the default error handler */
  jpeg      = &information;
  jpeg->err = jpeg_std_error(&jerr);
  jpeg_create_compress(jpeg);

  /* set some basic JPEG parameters */  
  jpeg->image_width      = width;
  jpeg->image_height     = height;
  jpeg->input_components = 1;
  jpeg->in_color_space   = JCS_GRAYSCALE;

  /* set compression parameters (i.e. no compression) */
  jpeg_set_defaults(jpeg);  

  /* set the percentage scaling factor for an quantization table
     (using IJG-recommended scaling curve). */
  if (compress) {    
    if (standard) {
      /* use the IJG-recommended quantization table */
      jpeg_set_quality(jpeg, quality, TRUE);
    } else {      
      /* use our generated quantization table */
      quality = jpeg_quality_scaling(quality);
      p = reinterpret_cast<const unsigned int*>(table);  
      jpeg_add_quant_table(jpeg, 0, p, quality, TRUE);
    }
  }

  /* set the destination to be a file */
  jpeg_stdio_dest(jpeg, file);

  /* compress the JPEG image */
  jpeg_start_compress(jpeg, TRUE);
  while (jpeg->next_scanline < height) {
    row[0] = &source.buffer[jpeg->next_scanline*stride];
    jpeg_write_scanlines(jpeg, row, 1);
  }

  /* clean up */
  jpeg_finish_compress(jpeg);
  jpeg_destroy_compress(jpeg);
  fclose(file);
}

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

void 
run_metrics(const RAW_IMAGE &original, 
	    const RAW_IMAGE &compressed)
{
  unsigned int             i;
  quality_measure_function measures[] = { average_distance, 
					  structural_content, 
					  cross_correlation,
					  image_fidelity,
					  laplacian_mean_square_error,
					  peak_mean_square_error,
					  absolute_error_1,
					  absolute_error_2,
					  //absolute_error_3,
					  mean_square_error_1,
					  mean_square_error_2,
					  //mean_square_error_3,
					  lp_1_norm,
					  lp_2_norm_1,
					  lp_2_norm_2,
					  //lp_2_norm_3,
					  lp_3_norm };
  
  cout << fixed << setprecision(10);
  for (i = 0; i < count_of(measures)-1; ++i) {
    cout << measures[i](original, compressed) << " ";
  }
  cout << measures[i](original, compressed) << endl;
}
