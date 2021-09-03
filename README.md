# LassoHyper
This is a collection of MATLAB codes to reproduce all the figures in the paper "Lasso hyperinterpolation over general regions" by Congpei An and Hao-Ning Wu.

This collection consists of four foulders, including Lasso hyperinterpolation on the interval, on the disc, on the sphere, and in the cube.
# Interval
To run the demo in the folder entitled "Interval", you may wish to install [Chebfun](http://www.chebfun.org/) in advance, which is a great function computing toolbox developed by Oxford's Numerical Analysis Group. See http://www.chebfun.org/ for Chebfun information. We also provide Chenfun in this folder! You are recommended to ADD IT TO PATH first.

# Disc
We are very grateful to Professor Kendall E. Atkinson of the University of Iowa for providing us with MATLAB codes of disc-related experiments, which were conducted in their paper entitled "On the norm of the hyperinterpolation operator on the unit disc and its use for the solution of the nonlinear Poisson equation" (by O. Hansen, K. Atkinson, and D. Chien). Our main function hyperdisc.m is written based on their M-file NonlinearElliptic.m.

# Sphere
Dr. Congpei An has developed and maintained a toolbox "sphere_approx_toolbox_v3.0" since 2017. This toolbox is also provided in this folder. Enjoy the spherical world!

# Cube
Codes in this folder is based on a package "Hyper3", which is available at https://www.math.unipd.it/~demarchi/software.html. This package corresponds to a paper entitled "New cubature formulae and hyperinterpolation in three variables" (by S. De Marchi, M. Vianello, and Y. Xu). We put this package in a subfolder called "utilities", and the credit of codes in this subfolder must go to S. De Marchi and M. Vianello. Based on their M-files, we write a function CX_Interpolation.m for hyperinterpolation with cosideration of noisy data. Our main function is CX_Interpolationlasso.m.






