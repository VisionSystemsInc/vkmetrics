#ifndef vkm_vrml_color_h
#define vkm_vrml_color_h

//:
// \file
// \brief A class with vrml utilities
// \author Isabel Restrepo mir@lems.brown.edu
// \date  Dec 8, 2009
//
// \verbatim
//  Modifications
//   <none yet>
// \endverbatim

#include <iostream>
#include <fstream>
#include <string>

class vkm_vrml_color
{
 public:
  //: store the color scheme 'classic' to generate a heatmap
  static unsigned heatmap_classic_size;
  static unsigned char heatmap_classic[256][3];

  static unsigned heatmap_custom_size;
  static unsigned char heatmap_custom[256][3];

  //apply a function to compute the custom map (one time only)
  //since the result is coded as an array initializer for future use
  static void compute_custom();
 private:
  static void custom_map(float t, float& r, float& g, float& b);// a function that defines the map
};

#endif
