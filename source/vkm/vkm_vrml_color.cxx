#include <iostream>
#include <cmath>
#include "vkm_vrml_color.h"

unsigned vkm_vrml_color::heatmap_classic_size = 256;
unsigned vkm_vrml_color::heatmap_custom_size = 256;
unsigned char vkm_vrml_color::heatmap_classic[256][3] = {
    {255, 237, 237},
    {255, 224, 224},
    {255, 209, 209},
    {255, 193, 193},
    {255, 176, 176},
    {255, 159, 159},
    {255, 142, 142},
    {255, 126, 126},
    {255, 110, 110},
    {255, 94, 94},
    {255, 81, 81},
    {255, 67, 67},
    {255, 56, 56},
    {255, 46, 46},
    {255, 37, 37},
    {255, 29, 29},
    {255, 23, 23},
    {255, 18, 18},
    {255, 14, 14},
    {255, 11, 11},
    {255, 8, 8},
    {255, 6, 6},
    {255, 5, 5},
    {255, 3, 3},
    {255, 2, 2},
    {255, 2, 2},
    {255, 1, 1},
    {255, 1, 1},
    {255, 0, 0},
    {255, 0, 0},
    {255, 0, 0},
    {255, 0, 0},
    {255, 0, 0},
    {255, 0, 0},
    {255, 0, 0},
    {255, 0, 0},
    {255, 1, 0},
    {255, 4, 0},
    {255, 6, 0},
    {255, 10, 0},
    {255, 14, 0},
    {255, 18, 0},
    {255, 22, 0},
    {255, 26, 0},
    {255, 31, 0},
    {255, 36, 0},
    {255, 41, 0},
    {255, 45, 0},
    {255, 51, 0},
    {255, 57, 0},
    {255, 62, 0},
    {255, 68, 0},
    {255, 74, 0},
    {255, 81, 0},
    {255, 86, 0},
    {255, 93, 0},
    {255, 99, 0},
    {255, 105, 0},
    {255, 111, 0},
    {255, 118, 0},
    {255, 124, 0},
    {255, 131, 0},
    {255, 137, 0},
    {255, 144, 0},
    {255, 150, 0},
    {255, 156, 0},
    {255, 163, 0},
    {255, 169, 0},
    {255, 175, 0},
    {255, 181, 0},
    {255, 187, 0},
    {255, 192, 0},
    {255, 198, 0},
    {255, 203, 0},
    {255, 208, 0},
    {255, 213, 0},
    {255, 218, 0},
    {255, 222, 0},
    {255, 227, 0},
    {255, 232, 0},
    {255, 235, 0},
    {255, 238, 0},
    {255, 242, 0},
    {255, 245, 0},
    {255, 247, 0},
    {255, 250, 0},
    {255, 251, 0},
    {253, 252, 0},
    {250, 252, 1},
    {248, 252, 2},
    {244, 252, 2},
    {241, 252, 3},
    {237, 252, 3},
    {233, 252, 3},
    {229, 252, 4},
    {225, 252, 4},
    {220, 252, 5},
    {216, 252, 5},
    {211, 252, 6},
    {206, 252, 7},
    {201, 252, 7},
    {197, 252, 8},
    {191, 251, 8},
    {185, 249, 9},
    {180, 247, 9},
    {174, 246, 10},
    {169, 244, 11},
    {164, 242, 11},
    {158, 240, 12},
    {151, 238, 13},
    {146, 236, 14},
    {140, 233, 14},
    {134, 231, 15},
    {128, 228, 16},
    {122, 226, 17},
    {116, 223, 18},
    {110, 221, 19},
    {105, 218, 20},
    {99, 216, 21},
    {93, 214, 22},
    {88, 211, 23},
    {82, 209, 24},
    {76, 207, 25},
    {71, 204, 26},
    {66, 202, 28},
    {60, 200, 30},
    {55, 198, 31},
    {50, 196, 33},
    {45, 194, 34},
    {40, 191, 35},
    {36, 190, 37},
    {31, 188, 39},
    {27, 187, 40},
    {23, 185, 43},
    {19, 184, 44},
    {15, 183, 46},
    {12, 182, 48},
    {9, 181, 51},
    {6, 181, 53},
    {3, 180, 55},
    {1, 180, 57},
    {0, 180, 60},
    {0, 180, 62},
    {0, 180, 65},
    {0, 181, 68},
    {0, 182, 70},
    {0, 182, 74},
    {0, 183, 77},
    {0, 184, 80},
    {0, 184, 84},
    {0, 186, 88},
    {0, 187, 92},
    {0, 188, 95},
    {0, 190, 99},
    {0, 191, 104},
    {0, 193, 108},
    {0, 194, 112},
    {0, 196, 116},
    {0, 198, 120},
    {0, 200, 125},
    {0, 201, 129},
    {0, 203, 134},
    {0, 205, 138},
    {0, 207, 143},
    {0, 209, 147},
    {0, 211, 151},
    {0, 213, 156},
    {0, 215, 160},
    {0, 216, 165},
    {0, 219, 171},
    {0, 222, 178},
    {0, 224, 184},
    {0, 227, 190},
    {0, 229, 197},
    {0, 231, 203},
    {0, 233, 209},
    {0, 234, 214},
    {0, 234, 220},
    {0, 234, 225},
    {0, 234, 230},
    {0, 234, 234},
    {0, 234, 238},
    {0, 234, 242},
    {0, 234, 246},
    {0, 234, 248},
    {0, 234, 251},
    {0, 234, 254},
    {0, 234, 255},
    {0, 232, 255},
    {0, 228, 255},
    {0, 224, 255},
    {0, 219, 255},
    {0, 214, 254},
    {0, 208, 252},
    {0, 202, 250},
    {0, 195, 247},
    {0, 188, 244},
    {0, 180, 240},
    {0, 173, 236},
    {0, 164, 232},
    {0, 156, 228},
    {0, 147, 222},
    {0, 139, 218},
    {0, 130, 213},
    {0, 122, 208},
    {0, 117, 205},
    {0, 112, 203},
    {0, 107, 199},
    {0, 99, 196},
    {0, 93, 193},
    {0, 86, 189},
    {0, 78, 184},
    {0, 71, 180},
    {0, 65, 175},
    {0, 58, 171},
    {0, 52, 167},
    {0, 46, 162},
    {0, 40, 157},
    {0, 35, 152},
    {0, 30, 147},
    {0, 26, 142},
    {0, 22, 136},
    {0, 18, 131},
    {0, 15, 126},
    {0, 12, 120},
    {0, 9, 115},
    {1, 8, 110},
    {1, 6, 106},
    {1, 5, 101},
    {2, 4, 97},
    {3, 4, 92},
    {4, 5, 89},
    {5, 5, 85},
    {6, 6, 82},
    {7, 7, 79},
    {8, 8, 77},
    {10, 10, 77},
    {12, 12, 77},
    {14, 14, 76},
    {16, 16, 74},
    {19, 19, 73},
    {21, 21, 72},
    {24, 24, 71},
    {26, 26, 69},
    {29, 29, 70},
    {32, 32, 69},
    {35, 35, 68},
    {37, 37, 67},
    {40, 40, 67},
    {42, 42, 65},
    {44, 44, 65},
    {46, 46, 64},
    {48, 48, 63},
    {49, 50, 62},
    {51, 51, 61},
    {53, 52, 61}
};
unsigned char vkm_vrml_color::heatmap_custom[256][3] = {
    {0,0,255},
    {0,2,252},
    {0,5,249},
    {0,8,246},
    {0,11,243},
    {0,14,240},
    {0,17,237},
    {0,20,234},
    {0,23,231},
    {0,26,228},
    {0,29,225},
    {0,32,222},
    {0,35,219},
    {0,38,216},
    {0,41,213},
    {0,44,210},
    {0,47,207},
    {0,50,204},
    {0,53,201},
    {0,56,198},
    {0,59,195},
    {0,62,192},
    {0,65,189},
    {0,68,186},
    {0,71,183},
    {0,74,180},
    {0,77,177},
    {0,80,174},
    {0,83,171},
    {0,86,168},
    {0,89,165},
    {0,92,162},
    {0,95,159},
    {0,98,156},
    {0,101,153},
    {0,104,150},
    {0,107,147},
    {0,110,144},
    {0,113,141},
    {0,116,138},
    {0,119,135},
    {0,122,132},
    {0,125,129},
    {0,128,126},
    {0,131,123},
    {0,134,120},
    {0,137,117},
    {0,140,114},
    {0,143,111},
    {0,146,108},
    {0,149,105},
    {0,152,102},
    {0,155,99},
    {0,158,96},
    {0,161,93},
    {0,164,90},
    {0,167,87},
    {0,170,84},
    {0,173,81},
    {0,176,78},
    {0,179,75},
    {0,182,72},
    {0,185,69},
    {0,188,66},
    {0,191,63},
    {0,194,60},
    {0,197,57},
    {0,200,54},
    {0,203,51},
    {0,206,48},
    {0,209,45},
    {0,212,42},
    {0,215,39},
    {0,218,36},
    {0,221,33},
    {0,224,30},
    {0,227,27},
    {0,230,24},
    {0,233,21},
    {0,236,18},
    {0,239,15},
    {0,242,12},
    {0,245,9},
    {0,248,6},
    {0,251,3},
    {0,254,0},
    {1,255,0},
    {4,255,0},
    {7,255,0},
    {10,255,0},
    {13,255,0},
    {16,255,0},
    {19,255,0},
    {22,255,0},
    {25,255,0},
    {28,255,0},
    {31,255,0},
    {34,255,0},
    {37,255,0},
    {40,255,0},
    {43,255,0},
    {46,255,0},
    {49,255,0},
    {52,255,0},
    {55,255,0},
    {58,255,0},
    {61,255,0},
    {64,255,0},
    {67,255,0},
    {70,255,0},
    {73,255,0},
    {76,255,0},
    {79,255,0},
    {82,255,0},
    {85,255,0},
    {88,255,0},
    {91,255,0},
    {94,255,0},
    {97,255,0},
    {100,255,0},
    {103,255,0},
    {106,255,0},
    {109,255,0},
    {112,255,0},
    {115,255,0},
    {118,255,0},
    {121,255,0},
    {124,255,0},
    {127,255,0},
    {130,255,0},
    {133,255,0},
    {136,255,0},
    {139,255,0},
    {142,255,0},
    {145,255,0},
    {148,255,0},
    {151,255,0},
    {154,255,0},
    {157,255,0},
    {160,255,0},
    {163,255,0},
    {166,255,0},
    {169,255,0},
    {172,255,0},
    {175,255,0},
    {178,255,0},
    {181,255,0},
    {184,255,0},
    {187,255,0},
    {190,255,0},
    {193,255,0},
    {196,255,0},
    {199,255,0},
    {202,255,0},
    {205,255,0},
    {208,255,0},
    {211,255,0},
    {214,255,0},
    {217,255,0},
    {220,255,0},
    {223,255,0},
    {226,255,0},
    {229,255,0},
    {232,255,0},
    {235,255,0},
    {238,255,0},
    {241,255,0},
    {244,255,0},
    {247,255,0},
    {250,255,0},
    {253,255,0},
    {255,254,0},
    {255,251,0},
    {255,248,0},
    {255,245,0},
    {255,242,0},
    {255,239,0},
    {255,236,0},
    {255,233,0},
    {255,230,0},
    {255,227,0},
    {255,224,0},
    {255,221,0},
    {255,218,0},
    {255,215,0},
    {255,212,0},
    {255,209,0},
    {255,206,0},
    {255,203,0},
    {255,200,0},
    {255,197,0},
    {255,194,0},
    {255,191,0},
    {255,188,0},
    {255,185,0},
    {255,182,0},
    {255,179,0},
    {255,176,0},
    {255,173,0},
    {255,170,0},
    {255,167,0},
    {255,164,0},
    {255,161,0},
    {255,158,0},
    {255,155,0},
    {255,152,0},
    {255,149,0},
    {255,146,0},
    {255,143,0},
    {255,140,0},
    {255,137,0},
    {255,134,0},
    {255,131,0},
    {255,128,0},
    {255,125,0},
    {255,122,0},
    {255,119,0},
    {255,116,0},
    {255,113,0},
    {255,110,0},
    {255,107,0},
    {255,104,0},
    {255,101,0},
    {255,98,0},
    {255,95,0},
    {255,92,0},
    {255,89,0},
    {255,86,0},
    {255,83,0},
    {255,80,0},
    {255,77,0},
    {255,74,0},
    {255,71,0},
    {255,68,0},
    {255,65,0},
    {255,62,0},
    {255,59,0},
    {255,56,0},
    {255,53,0},
    {255,50,0},
    {255,47,0},
    {255,44,0},
    {255,41,0},
    {255,38,0},
    {255,35,0},
    {255,32,0},
    {255,29,0},
    {255,26,0},
    {255,23,0},
    {255,20,0},
    {255,17,0},
    {255,14,0},
    {255,11,0},
    {255,8,0},
    {255,5,0},
    {255,2,0}
};

// not used except to generate the "custom" table
void vkm_vrml_color::compute_custom(){
  for(size_t i = 0; i<heatmap_custom_size; ++i){
    float t = static_cast<float>(i)/static_cast<float>(heatmap_custom_size);
    float r, g, b;
    custom_map(t, r, g, b);
    float s = 255.0f;
    unsigned rc = static_cast<unsigned char>(r*s);
    unsigned gc = static_cast<unsigned char>(g*s);
    unsigned bc = static_cast<unsigned char>(b*s);
    // print out lookup table of colors
    std::cout << '{' << rc << ',' << gc << ',' << bc <<"}," << std::endl;
  }
}

void vkm_vrml_color::custom_map(float value, float& r, float& g, float& b){
  const int NUM_COLORS = 4;
  static float color[NUM_COLORS][3] = { {0,0,1}, {0,1,0}, {1,1,0}, {1,0,0} };
    // A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.
   int idx1;        // |-- Our desired color will be between these two indexes in "color".
   int idx2;        // |
   float fractBetween = 0;  // Fraction between "idx1" and "idx2" where our value is.
  if(value <= 0)      {  idx1 = idx2 = 0;            }    // accounts for an input <=0
  else if(value >= 1)  {  idx1 = idx2 = NUM_COLORS-1; }    // accounts for an input >=0
  else
  {
    value = value * (NUM_COLORS-1);        // Will multiply value by 3.
    idx1  = floor(value);                  // Our desired color will be after this index.
    idx2  = idx1+1;                        // ... and before this index (inclusive).
    fractBetween = value - float(idx1);    // Distance between the two indexes (0-1).
  }
  r   = (color[idx2][0] - color[idx1][0])*fractBetween + color[idx1][0];
  g = (color[idx2][1] - color[idx1][1])*fractBetween + color[idx1][1];
  b  = (color[idx2][2] - color[idx1][2])*fractBetween + color[idx1][2];
}