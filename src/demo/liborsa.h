/*
 * Copyright (c) 2019, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU Affero General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LIBORSA_H
#define LIBORSA_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C"
#endif
bool wrapper_estimate_homography(float* image1, int w1, int h1, int c1, float* image2, int w2, int h2, int c2, double precision, float fSiftRatio, float* H, float* in, float* out, float* im1w, float* im2w, float* mosaic);

#endif // LIBORSA_H
