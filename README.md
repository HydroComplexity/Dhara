
<p align="center">
<img src="http://hydrocomplexity.net/images/dharaD.png" width="150px" hspace=5 /> <br/>
</p>



**DHARA** is an open source software for distributed hydrology and river analysis. We designed it to be optimal for hybrid CPU-GPU parallel computing. 

This project aims to provide a consistent framework for scientists to use and futher develop  modules for their modeling purposes. By participating, you are expected to uphold this code and contribute your model to the community.

Visit [http://hydrocomplexity.github.io/Dhara](http://hydrocomplexity.github.io/Dhara) to learn more information.

## Getting Started

### Prerequisites

The following packages are required to compile and run DHARA
```
mpich or openmpi
CUDA (5.0 or later)
netcdf4
```


### Download

Go to the folder you where you want to download Dhara, do:
```
$ git clone git@github.com:HydroComplexity/Dhara.git
```

### Installation
```
$ cd Dhara
$ make
```

## Documentation
Dhara model is developed based on a set of following research work: [Drewry et al. 2010a](http://onlinelibrary.wiley.com/doi/10.1029/2010JG001340/abstract), [Drewry et al. 2010b](http://onlinelibrary.wiley.com/doi/10.1029/2010JG001341/abstract), [Le et al. 2011](http://www.pnas.org/content/108/37/15085.abstract), [Quijano et al. 2012](http://onlinelibrary.wiley.com/doi/10.1029/2011WR011416/abstract), [Le et al. 2015](http://www.sciencedirect.com/science/article/pii/S1364815215300207).

See our tutorials [here](https://github.com/HydroComplexity/Dhara/blob/master/docs/notebooks/Dhara_model.ipynb) to learn about DHARA and run examples.

## Citation
If you use Dhara model, please use the following citations:
* Le, P. V. V., & Kumar, P. (2017). Interaction between ecohydrologic dynamics and microtopographic variability under climate change. *Water Resources Research*, 53, 8383–8403. https://doi.org/10.1002/2017WR020377.
* Woo, D. K. & Kumar, P. (2017). Role of micro-topographic variability on the distribution of inorganic soil-nitrogen age in intensively managed landscape. *Water Resources Research*, 53, 8404–8422. https://doi.org/10.1002/2017WR021053
* Le, P.V.V., Kumar, P., Valocchi, A.J., and Dang, H.-V. (2015): GPU-based high-performance computing for integrated surface–sub-surface flow modeling. *Environmental Modelling & Software*. doi: 10.1016/j.envsoft.2015.07.015
* Drewry, D.T., P. Kumar, S. Long, C. Bernacchi, X.-Z. Liang, and M. Sivapalan (2010), Ecohydrological responses of dense canopies to environmental variability: 1. Interplay between vertical structure and photosynthetic pathway, *J. Geophys. Res.*, 115, G04022, doi:10.1029/2010JG001340.

## License
This software is freeware and is released under restricted licences. See LICENSE.txt for more information. 

## Contributions
Patches, bug fixes, developing new modules, and pull request are welcome on the GitHub page for this project.

## Acknowledgements
Big thanks to the following people for contributing to this project in myriad ways:
* Darren Drewry: Multi-layer canopy module development
* Juan Quijano: Multi-species modeling
* Venkatraman Srinivasan: Multi-layer canopy modeling 
* Hoang-Vu Dang: Parallel computing integration and development

## Contact Authors
* Phong Le: <mailto:levuvietphong@gmail.com>
* Praveen Kumar: <mailto:kumar1@illinois.edu>