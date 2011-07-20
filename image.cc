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

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

void
read_image(char *file_name, RAW_IMAGE &destination)
{
  FILE                   *file;
  jpeg_decompress_struct information, *jpeg;
  jpeg_error_mgr         jerr;
  JSAMPROW               row[1];
  unsigned int           stride, height, width;
  
  /* open the output file */
  if ((file = fopen(file_name, "rb")) == NULL) {
    cerr << "Cannot open file \"" << file_name << "\" for reading." << endl;
    exit(1);
  }
  
  /* use the default error handler */
  jpeg      = &information;
  jpeg->err = jpeg_std_error(&jerr);
  jpeg_create_decompress(jpeg);
  
  /* set the source to be a the file */
  jpeg_stdio_src(jpeg, file);
  
  /* read  JPEG parameters */  
  jpeg_read_header(jpeg, TRUE);
  
  /* shorten some names for convenience */
  stride = jpeg->image_width;
  width  = jpeg->image_width;
  height = jpeg->image_height;
  
  /* create the actual image buffer */
  if (destination.size() <= 0) {
    destination.create(width, height);
  }
  
  /* decompress the JPEG image */
  jpeg_start_decompress(jpeg);
  while (jpeg->output_scanline < height) {
    row[0] = &destination.buffer[jpeg->output_scanline*stride];
    jpeg_read_scanlines(jpeg, row, 1);
  }
  
  /* clean up */
  jpeg_finish_decompress(jpeg);
  jpeg_destroy_decompress(jpeg);
  fclose(file);
}

void 
write_image(char *file_name, const RAW_IMAGE &source, const unsigned int table[8][8], int quality)
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
     (using IJG-recommended scaling curve) with the given quantization
     table */
  quality = jpeg_quality_scaling(quality);
  p = reinterpret_cast<const unsigned int*>(table);
  jpeg_add_quant_table(jpeg, 0, p, quality, TRUE);
  
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

