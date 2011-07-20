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

#include "error.h"
#include "image.h"
#include "measure.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

char              *tool_name;
bool              verbose;
verbosity::type_t noisiness;

int
main(int argc, char *argv[])
{
  long double m1, m2, ratio;
  char        c;
  char        *name;
  RAW_IMAGE   original, standard, generaged;
  struct stat status;
  
  /* save the program name */
  tool_name = argv[0];
  
  /* we will privide our own error messages */
  opterr = 0;
  
  /* extract any command-line options the user provided */
  while (-1 != (c = getopt(argc, argv, ":v"))) {
    switch (c) {
    case 'v': 
      verbose = !verbose;
      break;
    case ':':
      die("Option -%c requires an operand; that is, an integer or string value.\n", optopt);
      break;
    case '?':
      die("Unknown option: `-%c'\n", optopt);
      break;
    default:
      abort();
      break;
    }
  }
  
  name = argv[optind++];
  read_image(name, original);
    
  name = argv[optind++];
  read_image(name, standard);
  
  stat(name, &status);
  ratio = status.st_size;
  
  name = argv[optind++];
  read_image(name, generaged);
  
  stat(name, &status);
  ratio /= status.st_size;
  
  m1 = average_distance(original, standard);
  m2 = average_distance(original, generaged);
  
  printf("%Lf %Lf %Lf", m1, m2, ratio);
  
  return 0;
}
