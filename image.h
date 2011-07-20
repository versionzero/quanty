
#ifndef _CONVERT_H_
#define _CONVERT_H_

#include "types.h"
#include <pam.h>
#include <jpeglib.h>
#include <iostream>

struct PGM_IMAGE {
  
  pam   information;
  tuple **buffer;
  
  PGM_IMAGE() : buffer(NULL) {
  }
  
  ~PGM_IMAGE() {
    
    if (buffer) {
      pnm_freepamarray(buffer, &information);
      buffer = NULL;
    }
    
  }
  
  uint width()  const { return information.width;  }
  uint height() const { return information.height; }
  uint size()   const { return width()*height(); }
  
};

struct RAW_IMAGE {
  
  uchar *buffer;
  uint  _width, _height;
  
  RAW_IMAGE() : buffer(NULL), _width(0), _height(0) {
  }
  
  ~RAW_IMAGE() {
    
    _width = _height = 0;  

    if (NULL != buffer) {
      delete [] buffer;      
      buffer = NULL;
    }

  }
  
  void create(uint width, uint height) {

    _width  = width;
    _height = height;

    if ((buffer = new uchar[width*height]) == NULL) {
      std::cerr << "Cannot allocate raw image buffer." << std::endl;
      exit(1);
    }

  }
  
  uint width()  const { return _width;  }
  uint height() const { return _height; }
  uint size()   const { return width()*height(); }

};

struct JPEG_IMAGE {

  jpeg_compress_struct   information;
  jpeg_decompress_struct dinfo;
  uchar          *buffer;
  uint      _size;
  bool                   _compressed;

  JPEG_IMAGE(bool compressed = true) : buffer(NULL), _size(0), _compressed(compressed) {
  }

  ~JPEG_IMAGE() {

    _size = 0;

    if (NULL != buffer) {
      delete [] buffer;
      buffer = NULL;
    }

  }

  void create(uint width, uint height) {
    
    /* set the buffer size */
    _size = width*height;
    
    /* allocate the image's buffer */
    if ((buffer = new uchar[_size]) == NULL) {
      std::cerr << "Cannot allocate JPEG image buffer." << std::endl;
      exit(1);
    }

  }

  uint width()  const { 
    return _compressed ? information.image_width  : dinfo.image_width;  
  }
  uint height() const { 
    return _compressed ? information.image_height : dinfo.image_height; 
  }
  uint size()   const { return _size; }

};

void read_image(char *file_name, RAW_IMAGE &destination);
void write_image(char *file_name, const RAW_IMAGE &source, const uint table[8][8], int quality);

#endif
