# chondrocyte

<p align="center">
  <img src="https://github.com/mmaleck/chondrocyte/blob/master/docs/body_nankotsu.png" width="200">
  <img src="https://github.com/mmaleck/chondrocyte/blob/master/docs/body_nankotsu_suriheri.png" width="200">
</p>

## Description
Model of the fundamental electrophysiological characteristics of human articular chondrocytes, focused on K+ transport mechanisms, and their ability to regulate to the resting membrane potential, Em.

## Installation

The following python modules are needed to use chondrocyte, 

* Numpy (developed with Version: 1.18.1)
* Scipy (developed with Version: 1.6.0)
* matplotlib (developed with Version: 3.4.1)

To use chondrocyte, please clone the repository.

```
git clone https://github.com/mmaleck/chondrocyte.git
cd chondrocyte
```

## Usage

chondrocyte was first developed with MATLAB and thus we have the folder [matlab_code](matlab_code). Then, it was reimplemented in python with some improvements and model expansions.
Therefore, we encourage users to use [python_code](python_code), and also we only write instructions on how to use chondrocyte in python for that matter. 

Inside [python_code](python_code), there are four python files, `chondrocyte.py`, `functions.py`, `params.py`, and `postprocess.py` in addition to four jupyter notebooks `demo_*.ipynb`. 

These jupyter notebooks are great examples of how to use the chondrocyte interactively and the demos include replications of figures from the published paper. 

`chondrocyte.py` is the core part of the chondrocyte. It imports necessary functions from `functions.py` and parameters from `params.py` to execute the program, and `postprocess.py` is merely to demonstrate how to post-process the results. 

For detailed instructions about the code, we ask users to have a look at the code as it contains the explanation as comments.

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

If you use chondrocyte and want to cite it, you can either cite this website or the following paper.

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



