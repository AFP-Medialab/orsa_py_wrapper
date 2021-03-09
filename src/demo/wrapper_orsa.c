#include "liborsa.h"

bool estimate_homography_py(float* image1, int w1, int h1, int c1, float* image2, int w2, int h2, int c2, double precision, float fSiftRatio, float* H, float* in, float* out, float* im1w, float* im2w, float* mosaic)
{
    return wrapper_estimate_homography(image1, w1, h1, c1, image2, w2, h2, c2, precision, fSiftRatio, H, in, out, im1w, im2w, mosaic);
}
