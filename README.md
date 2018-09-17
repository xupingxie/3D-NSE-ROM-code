# 3D-NSE-ROM-code

This repo contains 3D ROM code for Navier-Stokes-equations.  It contains standard Galerkin projection based ROM and the new AD-ROM. All the detailed algorithms are in the paper "approximate deconvolution reduced order modeling".  The main driver script is ADROM.m 
The only thing missing in this repo is the offline data matrix (POD matrix and tensors). The data is available upon request.  Please let me know if you need.

## Navier-Stokes equations
![equation](https://latex.codecogs.com/gif.latex?%5Cmathbf%7Bu%7D_t%20&plus;%5Cmathbf%7Bu%7D%5Ccdot%5Cnabla%5Cmathbf%7Bu%7D-Re%5CDelta%20%5Cmathbf%7Bu%7D%20&plus;%5Cnabla%20p%20%3D%20f%20%5C%5C%20%5Cnabla%5Ccdot%5Cmathbf%7Bu%7D%20%3D%200)
## Citation:
    @article{xie2017approximate,
      title={Approximate deconvolution reduced order modeling},
      author={Xie, Xuping and Wells, David and Wang, Zhu and Iliescu, Traian},
      journal={Computer Methods in Applied Mechanics and Engineering},
      volume={313},
      pages={512--534},
      year={2017},
      publisher={Elsevier}
    }
