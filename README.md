
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
* mpich or openmpi
* CUDA (5.0 or later)
* netcdf4
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

See our tutorials [here](https://github.com/HydroComplexity/Dhara/blob/master/docs/notebooks/Dhara_model.ipynb) to learn about DHARA and run examples.


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