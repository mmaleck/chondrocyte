# chondrocyte

<p align="center">
  <img src="https://github.com/mmaleck/chondrocyte/blob/master/docs/body_nankotsu.png" width="200">
  <img src="https://github.com/mmaleck/chondrocyte/blob/master/docs/body_nankotsu_suriheri.png" width="200">
</p>

## Description
Model of the fundamental electrophysiological characteristics of human articular chondrocytes, focused on K+ transport mechanisms, and their ability to regulate to the resting membrane potential, Em.

## Installation

The followinfg python modules are needed to use Choondrocyte, 

* Numpy (developed with Version: 1.18.1)
* Scipy (developed with Version: 1.6.0)

To use chondrocyte, please clone the repository.

```
git clone https://github.com/mmaleck/chondrocyte.git
cd chondrocyte
```

## Usage

chondrocyte was first developed with MATLAB and thus we have a folder [matlab_code](matlab_code). Then, it was reimplemented in python with some improvements and model expansions. \
So, we enocurage users to use [python_code](python_code) and also we only write instructions how to use chondrocyte in python for that reason. \

Inside [python_code](python_code), there are four python files, `chondrocyte.py`, `functions.py`, `params.py`, and `postprocess.py` in addition to four jupter notebooks `demo_*.ipynb`. \
`chondrocyte.py` is the core part of the chdondrocyte. It imports necessary functions from `functions.py` and parameters from `params.py` to execute the program, and `postprocess.py` is merely to demonstrate how to postporocess the results. 

For detailed instruction aboout the code, we ask users to have a look on the code as it contains the explanation as comments.

## Authors
chondrocyte is developed by

  * Molly Maleckar 
  * Kei Yamamoto 
  * Sophie Fischer 
  * Maria Perona Fjeldstad 

## License

chondrocyte is licensed under the MIT License.

chondrocyte is Copyright (2021) by the authors.

## Citation

If you use chondrocyte and want to cite it, you can either cite this web site or the following paper.

```
@article{Maleckar2018,
  title = {The resting potential and K+ currents in Primary Human articular chondrocytes},
  author = {Maleckar, Mary M. and Clark, Robert B. and Votta, Bartholomew and Giles, Wayne R.},
  doi = {10.3389/fphys.2018.00974},
  journal = {Frontiers in Physiology},
  volume = {9},
  year = {2018}
}
```

## Contact

Please report bugs and other issues through the [issue tracker](https://github.com/mmaleck/chondrocyte/issues)

## How to contribute ?

See [guideline how to contribute](docs/CONTRIBUTING.md)



