# chondrocyte

<p align="center">
  <img src="https://github.com/mmaleck/chondrocyte/blob/master/docs/body_nankotsu.png" width="200">
  <img src="https://github.com/mmaleck/chondrocyte/blob/master/docs/body_nankotsu_suriheri.png" width="200">
</p>

## Description
Model of the fundamental electrophysiological characteristics of human articular chondrocytes, focused on K+ selective ion channels and ion pumps, but also including the contributions of Na+, Ca2+, and Mg2+ on the chondrocyte membrane potential (Em). 

## Installation

The following python modules are needed to use chondrocyte, 

* NumPy (developed with Version: 1.18.1)
* SciPy (developed with Version: 1.6.0)
* Matplotlib (developed with Version: 3.4.1)

To use chondrocyte, please clone the repository:

```
git clone https://github.com/mmaleck/chondrocyte.git
cd chondrocyte
```

## Usage

chondrocyte was originally developed with MATLAB and can be found in the [matlab_code](matlab_code) folder. It has been reimplemented in python with some improvements and model expansions.
Therefore, we encourage users to use [python_code](python_code), and also we only write instructions on how to use chondrocyte in python for that matter. 

Inside [python_code](python_code), there are four python files, `chondrocyte.py`, `functions.py`, `params.py`, and `postprocess.py` in addition to four jupyter notebooks `demo_*.ipynb`. 

These jupyter notebooks are great examples of how to use the chondrocyte interactively, and the demos include replications of figures from the published papers. 

`chondrocyte.py` is the core part of the chondrocyte model. 
It imports necessary functions from `functions.py` and parameters from `params.py` to execute the program, and `postprocess.py` demonstrates one way to process the results. 

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

If you use chondrocyte and want to cite it, you can either cite this website or one of the following papers.

```
@article{Maleckar2018,
  title = {{The resting potential and K+ currents in Primary Human articular chondrocytes}},
  author = {Maleckar, Mary M. and Clark, Robert B. and Votta, Bartholomew and Giles, Wayne R.},
  doi = {10.3389/fphys.2018.00974},
  journal = {Frontiers in Physiology},
  volume = {9},
  year = {2018}
}
```

```
@article{Maleckar2020,
  title = {{ Physiological Effects of the Electrogenic Current Generated by the Na + /K + Pump in Mammalian Articular Chondrocytes }},
  author = {Maleckar, Mary M. and Mart{\'{i}}n-Vasallo, Pablo and Giles, Wayne R. and Mobasheri, Ali},
  doi = {10.1089/bioe.2020.0036},
  journal = {Bioelectricity},
  volume = {2},
  year = {2020}
}
```

## Contact

Please report bugs and other issues through the [issue tracker](https://github.com/mmaleck/chondrocyte/issues)

## How to contribute ?

See [guideline how to contribute](docs/CONTRIBUTING.md)



