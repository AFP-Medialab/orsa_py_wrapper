/**
 * @file liborsa.cpp
 * @brief Homographic image registration library adaptation
 * @author Pascal Monasse, Pierre Moulon, Thibaud Ehret
 * 
 * Copyright (c) 2020 Pascal Monasse, Pierre Moulon, Thibaud Ehret
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <ctime>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "libImage/image_io.hpp"
#include "libImage/image_drawing.hpp"
#include "libOrsa/homography_model.hpp"

#include "extras/libNumerics/numerics.h"
#include "extras/sift/library.h"

#include "demo/siftMatch.hpp"
#include "demo/warping.hpp"

/// Number of random samples in ORSA
static const int ITER_ORSA=10000;

static void display_stats(const std::vector<Match>& vec_matchings,
                          const std::vector<int>& vec_inliers,
                          libNumerics::matrix<double>& H) {
  std::vector<int>::const_iterator it=vec_inliers.begin();
  double l2=0, linf=0;
  for(; it!=vec_inliers.end(); ++it) {
    const Match& m=vec_matchings[*it];
    double x1=m.x1, y1=m.y1;
    TransformH(H, x1, y1);
    double e = (m.x2-x1)*(m.x2-x1) + (m.y2-y1)*(m.y2-y1);
    l2 += e;
    if(linf < e)
      linf = e;
  }
  std::cout << "Average/max error: "
            << sqrt(l2/vec_inliers.size()) << "/"
            << sqrt(linf) <<std::endl;
}

/// Output inlier and oulier matches in image files.
void display_match(const std::vector<Match>& vec_matchings,
                   std::vector<int>& vec_inliers,
                   const libNumerics::matrix<double>* H,
                   const libNumerics::matrix<double>& H1,
                   const libNumerics::matrix<double>& H2,
                   Image<RGBColor>& in, Image<RGBColor>& out)
{
  std::sort(vec_inliers.begin(), vec_inliers.end());

  // For outliers, show vector (yellow) from prediction to observation
  const RGBColor col=YELLOW;
  std::vector<int>::const_iterator it = vec_inliers.begin();
  matchingslist::const_iterator m = vec_matchings.begin();
  if(H) // Otherwise, no prediction
    for(int i=0; m != vec_matchings.end(); ++m, ++i) {
      if(it != vec_inliers.end() && i==*it)
        ++it;
      else { //Outlier
          double x1=m->x1, y1=m->y1, x2=m->x2, y2=m->y2;
          TransformH(H2 * *H, x1, y1);
          TransformH(H2, x2, y2);
          libs::DrawLine((int)x1,(int)y1,(int)x2,(int)y2, col,&out);
      }
    }

  // Show link for inliers (green) and outliers (red)
  it = vec_inliers.begin();
  m = vec_matchings.begin();
  for(int i=0; m != vec_matchings.end(); ++m, ++i)
  {
    Image<RGBColor>* im=&out;
    RGBColor col=RED;
    if(it != vec_inliers.end() && i==*it) {
      ++it;
      im=&in;
      col=GREEN;
    }
    double x1=m->x1, y1=m->y1, x2=m->x2, y2=m->y2;
    TransformH(H1, x1, y1);
    TransformH(H2, x2, y2);
    libs::DrawLine((int)x1,(int)y1,(int)x2, (int)y2, col, im);
  }
}

/// ORSA homography estimation
bool ORSA(const std::vector<Match>& vec_matchings, int w1,int h1, int w2,int h2,
          double precision,
          libNumerics::matrix<double>& H, std::vector<int>& vec_inliers)
{
  const int n = static_cast<int>( vec_matchings.size() );
  if(n < 5)
  {
      std::cerr << "Error: ORSA needs 5 matches or more to proceed" <<std::endl;
      return false;
  }
  libNumerics::matrix<double> xA(2,n), xB(2,n);

  for (int i=0; i < n; ++i)
  {
    xA(0,i) = vec_matchings[i].x1;
    xA(1,i) = vec_matchings[i].y1;
    xB(0,i) = vec_matchings[i].x2;
    xB(1,i) = vec_matchings[i].y2;
  }

  orsa::HomographyModel model(xA, w1, h1, xB, w2, h2, true);
  //model.setConvergenceCheck(true);

  if(model.orsa(vec_inliers, ITER_ORSA, &precision, &H, true)>0.0)
    return false;
  //std::cout << "Before refinement: ";
  //display_stats(vec_matchings, vec_inliers, H);
  if( model.ComputeModel(vec_inliers,&H) ) // Re-estimate with all inliers
  {
    //std::cout << "After  refinement: ";
    //display_stats(vec_matchings, vec_inliers, H);
  } else
    std::cerr << "Warning: error in refinement, result is suspect" <<std::endl;
  return true;
}

/// Return 3x3 translation matrix
libNumerics::matrix<double> translation(double dx, double dy) {
    libNumerics::matrix<double> T = libNumerics::matrix<double>::eye(3);
    T(0,2) = dx;
    T(1,2) = dy;
    return T;
}

/// Return 3x3 zoom-translation matrix
libNumerics::matrix<double> zoomtrans(double z, double dx, double dy) {
    libNumerics::matrix<double> T = libNumerics::matrix<double>::eye(3);
    T(0,0)=T(1,1) = z;
    T(0,2) = dx;
    T(1,2) = dy;
    return T;
}


#ifdef __cplusplus
extern "C"
#endif
bool wrapper_estimate_homography(float* im1, int _w1, int _h1, int c1, float* im2, int _w2, int _h2, int c2, double precision, float fSiftRatio, float* _H, float* _in, float* _out, float* im1w, float* im2w, float* mosaic)
{
    // Init random seed
    time_t seed = time(0); // Replace by a fixed value to debug a reproducible run
    srand((unsigned int)seed);

    // Cast images into unsigned char
    std::vector<unsigned char> tim1(_w1 * _h1 * 3);
    for(int i=0; i<_w1*_h1; ++i)
    {
        tim1[i*3    ] = (c1 == 1) ? (unsigned char) im1[i] : (unsigned char) im1[i*3    ];
        tim1[i*3 + 1] = (c1 == 1) ? (unsigned char) im1[i] : (unsigned char) im1[i*3 + 1];
        tim1[i*3 + 2] = (c1 == 1) ? (unsigned char) im1[i] : (unsigned char) im1[i*3 + 2];
    }
    Image<RGBColor> image1(_w1, _h1, (const RGBColor*) &tim1[0]);

    std::vector<unsigned char> tim2(_w2 * _h2 * 3);
    for(int i=0; i<_w2*_h2; ++i)
    {
        tim2[i*3    ] = (c2 == 1) ? (unsigned char) im2[i] : (unsigned char) im2[i*3    ];
        tim2[i*3 + 1] = (c2 == 1) ? (unsigned char) im2[i] : (unsigned char) im2[i*3 + 1];
        tim2[i*3 + 2] = (c2 == 1) ? (unsigned char) im2[i] : (unsigned char) im2[i*3 + 2];
    }
    Image<RGBColor> image2(_w2, _h2, (const RGBColor*) &tim2[0]);

    Image<unsigned char> image1Gray, image2Gray;
    libs::convertImage(image1, &image1Gray);
    libs::convertImage(image2, &image2Gray);

    std::vector<Match> vec_matchings;
    SIFT(image1Gray, image2Gray, vec_matchings, fSiftRatio, 0);

    // Remove duplicates (frequent with SIFT)
    rm_duplicates(vec_matchings);

    const int w1=int(image1Gray.Width()), h1=int(image1Gray.Height());
    const int w2=int(image2Gray.Width()), h2=int(image2Gray.Height());

    // Estimation of homography with ORSA
    libNumerics::matrix<double> H(3,3);
    std::vector<int> vec_inliers;
    bool ok = ORSA(vec_matchings, w1, h1, w2, h2, precision, H, vec_inliers);
    if(ok)
        H /= H(2,2);

    _H[0] = H(0,0);
    _H[1] = H(0,1);
    _H[2] = H(0,2);
    _H[3] = H(1,0);
    _H[4] = H(1,1);
    _H[5] = H(1,2);
    _H[6] = H(2,0);
    _H[7] = H(2,1);
    _H[8] = H(2,2);

    if(_in != NULL && _out != NULL)
    {
        // Create a visual representation of the inliers and outliers
        int w = std::max(w1,w2);
        float z = w/(float)(w1+w2); //Set width as max of two images
        int h = int(z*std::max(h1,h2));
        Image<unsigned char> concat(w, h, 255);
        libNumerics::matrix<double> T1=zoomtrans(z, 0,   (concat.Height()-h1*z)/2);
        libNumerics::matrix<double> T2=zoomtrans(z, w1*z,(concat.Height()-h2*z)/2);
        Warp(image1Gray, T1, image2Gray, T2, concat);
        Image<RGBColor> in;
        libs::convertImage(concat, &in);
        libs::DrawLine(int(w1*z),0, int(w1*z),int(concat.Height()), BLUE, &in);
        Image<RGBColor> out(in);
        display_match(vec_matchings, vec_inliers, ok? &H:0, T1, T2, in, out);

        for(int y = 0; y < h; ++y)
        for(int x = 0; x < w; ++x)
        {
            RGBColor pixel = in(y, x);
            _in[      y*w+x] = pixel.r;
            _in[1*w*h+y*w+x] = pixel.g;
            _in[2*w*h+y*w+x] = pixel.b;
        }

        for(int y = 0; y < h; ++y)
        for(int x = 0; x < w; ++x)
        {
            RGBColor pixel = out(y, x);
            _out[      y*w+x] = pixel.r;
            _out[1*w*h+y*w+x] = pixel.g;
            _out[2*w*h+y*w+x] = pixel.b;
        }
    }

    if(im1w != NULL || im2w != NULL || mosaic != NULL)
    {
        Rect intersection;
        if(IntersectionBox(w1, h1, w2, h2, H, intersection) &&
                intersection.Width() > 0 && intersection.Height() > 0)
        {
            int xc=(intersection.left+intersection.right)/2;
            int yc=(intersection.top+intersection.bottom)/2;
            size_t wM = std::max(image1.Width(), image2.Width());
            size_t hM = std::max(image1.Height(), image2.Height());
            int xo=static_cast<int>(wM/2);
            int yo=static_cast<int>(hM/2);
            libNumerics::matrix<double> T = translation(xo-xc, yo-yc);

            if(mosaic != NULL)
            {
                Image<RGBColor> imageMosaic(wM, hM);
                imageMosaic.fill(WHITE);
                Warp(image1, T*H, image2, T, imageMosaic);

                for(int y = 0; y < hM; ++y)
                for(int x = 0; x < wM; ++x)
                {
                    RGBColor pixel = imageMosaic(y, x);
                    mosaic[        y*wM+x] = pixel.r;
                    mosaic[1*wM*hM+y*wM+x] = pixel.g;
                    mosaic[2*wM*hM+y*wM+x] = pixel.b;
                }
            }

            if(im1w != NULL)
            {
                Image<RGBColor> WarpingA(wM, hM, WHITE);
                Warp(image1, T*H, WarpingA);

                for(int y = 0; y < hM; ++y)
                for(int x = 0; x < wM; ++x)
                {
                    RGBColor pixel = WarpingA(y, x);
                    im1w[        y*wM+x] = pixel.r;
                    im1w[1*wM*hM+y*wM+x] = pixel.g;
                    im1w[2*wM*hM+y*wM+x] = pixel.b;
                }
            }

            if(im2w != NULL)
            {
                Image<RGBColor> WarpingB(wM, hM, WHITE);
                Warp(image2, T, WarpingB);

                for(int y = 0; y < hM; ++y)
                for(int x = 0; x < wM; ++x)
                {
                    RGBColor pixel = WarpingB(y, x);
                    im2w[        y*wM+x] = pixel.r;
                    im2w[1*wM*hM+y*wM+x] = pixel.g;
                    im2w[2*wM*hM+y*wM+x] = pixel.b;
                }
            }
        }
    }

    return ok;
}
