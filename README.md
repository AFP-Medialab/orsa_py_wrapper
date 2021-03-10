PYTHON WRAPPER FOR ORSA
=======================

OVERVIEW
--------

This source code is a fork of an [IPOL](http://www.ipol.im/) publication. Please cite 
[http://www.ipol.im/pub/art/2012/mmm-oh/](http://www.ipol.im/pub/art/2012/mmm-oh/) if you use this code
as part of your research. It provides an additional python wrapper to use the code natively in python.

COMPILATION
-----------

The code is compilable on Unix/Linux and on Mac OS (requires a simple manual fix though!). 

**Compilation:** requires the cmake and make programs.

**Dependencies:**  Requires libpng, libtiff and libjpeg.
 
Configure and compile the source code using cmake and make.  It is recommended
that you create a folder for building:

UNIX/LINUX/MAC:
```
$ mkdir build; cd build
$ cmake ../src
$ make
```

/!\ If you're on mac, the libraries are compiled as *.dlib* instead of *.so*. You'll need to change the extension 
for the provided python scripts to work out of the box (this will be fixed later)
```
$ mv demo/liborsa.dlib demo/liborsa.so
$ mv demo/libworsa.dlib demo/libworsa.so
```

USAGE
-----

The wrapper can be accessed by importing orsa in python:
```
from orsa import orsa
```
This library provides a single function `orsa.estimate_homography_py(im1, w1, h1, c1, im2, w2, h2, c2, precision, siftRatio, H, inl, outl, im1w, im2w, mosaic)`:
* `im1`: Numpy array corresponding to the first image. It has a shape of (`c1`, `h1`, `w1`) and be of type `np.float32`.
* `im2`: Numpy array corresponding to the first image. It has a shape of (`c2`, `h2`, `w2`) and be of type `np.float32`.
* `precision`: Precision of the processing. Recommended value is `0`.
* `siftRatio`: Sift ratio used during the matching part. Recommended value is `0.6`.
* `H`: Numpy array that will contain the estimated homography (it needs to be of shape `9` and of type `np.float32`).
* `inl` and `outl`: Numpy arrays that will contain a visual represensation of the inliers for `inl` and the outliers for `outl`. They both need to be of shape (`3`, `max(h1,h2) * max(w1,w2)/(w1+w2)`, `max(w1,w2)`) and be of type `np.float32`. It is not be computed when set to NULL.
* `im1w`, `im2w`, `mosaic`: Numpy arrays that will contain respectively the warped `im1`, the warped `im2` and a mosaic made from the two. They all need to be of shape (`3`, `max(h1,h2)`, `max(w1,w2)`) and be of type `np.float32`. It is not be computed when set to NULL.

A test code `main.py` is also provided to show how this can be used. Dependencies for this test code can be installed using the provided requirements file `requirements.txt`.
