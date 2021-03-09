import iio
from orsa import orsa, NULL, STR, wrap
import numpy as np

im1 = iio.read('data/carc1.jpg').astype(np.float32)
im2 = iio.read('data/carc2.jpg').astype(np.float32)
H = np.zeros(9, dtype=np.float32)
h1, w1, c1 = im1.shape
h2, w2, c2 = im2.shape

im1g = np.mean(im1, axis=2)
im2g = np.mean(im2, axis=2)
print(im1.shape)

w = max(w1,w2)
z = w / (w1 + w2)
h = int(z*max(h1,h2))

inl = np.zeros((3, h, w), dtype = np.float32)
outl = np.zeros((3, h, w), dtype = np.float32)

h = max(h1,h2)
im1w = np.zeros((3, h, w), dtype = np.float32)
im2w = np.zeros((3, h, w), dtype = np.float32)
mosaic = np.zeros((3, h, w), dtype = np.float32)

detection = orsa.estimate_homography_py(im1, w1, h1, c1, im2, w2, h2, c2, 0, 0.6, H, inl, outl, im1w, im2w, mosaic)

iio.write('in.png', inl.transpose(1,2,0))
iio.write('out.png', outl.transpose(1,2,0))

iio.write("im1_warped.png", im1w.transpose(1,2,0))
iio.write("im2_warped.png", im2w.transpose(1,2,0))
iio.write("mosaic.png", mosaic.transpose(1,2,0))
