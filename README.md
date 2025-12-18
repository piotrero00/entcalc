# entcalc
entcalcpy is a Matlab package for computing the geometric entanglement of a given quantum state. 
It works by computing lower and upper bounds for the geometric entanglement.
These bounds are often close to each other, allowing us to estimate the value of the geometric entanglement in a given quantum state.
All functions of the package are in entcalc directory.
## Table of Contents
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Getting started](#getting-started)
- [Issues](#issues)
- [Acknowledgment](#acknowledgment)
- [License](#license)
- [Citing](#citing)

## Dependencies
This package requires the following Python packages:
- `cvx` ≥ 2.2.2  
- `qetlab` ≥ 1.0 
## Getting started
To get started with entcalc, we recommend reading readme file. The documentation is written in docstrings.
For example, the following code computes the upper bound of the geometric entanglement of a random quantum state.
```matlab
% MATLAB usage example
rho = RandomDensityMatrix(8);
rho = rho*rho';
rho = rho/trace(rho);

disp(uppersame(rho, [2 2 2]))
disp(upperbip(rho, [4 2]))
```
Similarly, we can compute the lower bound
```matlab
rho = RandomDensityMatrix(8);
rho = rho*rho';
rho = rho/trace(rho);
disp(ge_mixed_gr(rho,[2,2,2])) %second argument specifies that we have a 3-qubit state
disp(ge_mixed_gr(rho,[4,2])) %here we have a 4x2 bipartite state
```
The lower bound can be computed using four different methods.
## ge_mixed_gr vs ge_mixed_ra_gr
ge_mixed_gr and ge_mixed_sm as a first step computes a purification of the input state. They do it by so-called canonical purification. 
ge_mixed_ra_gr takes input state in the form of orthogonal decomposition. Thanks to it, it can compute lower-dimensional purification and compute
 lower boudnd more efficiently.
 For example, if ge_mixed_gr takes n-qubit state as input, it must optimize over matrix with dimesnions $2^{2n}\times 2^{2n}$ matrices. If the state is rank-2, ge_mixed_ra_gr optimizes
 over matrix with dimesnions $2^{n+1}\times 2^{n+1}$. Thanks to this, it can compute entanglement for more qubits and can make it much faster.
## Solver
entcalc used SPDT3 to solve semi-definite programming problems. For a detailed information of solver setting visit CVX page.
## Issues
If you find any issues, we encourage you to report them via GitHub or by emailing maspiotr00@gmail.com.
## Acknowledgment
Special thanks to Krystyna Mordzińska for coming up with the package name after a creative brainstorming session.
## License
This project is licensed under the MIT License - see the LICENSE file for details.
## Citing
If you use entcalc in academic work, please cite this package. Below we give bibtex version of citation:
```latex
@misc{masajada2025entcalctoolkitcalculatinggeometric,
      title={ENTCALC: Toolkit for calculating geometric entanglement in multipartite quantum systems}, 
      author={Piotr Masajada and Aby Philip and Alexander Streltsov},
      year={2025},
      eprint={2512.10884},
      archivePrefix={arXiv},
      primaryClass={quant-ph},
      url={https://arxiv.org/abs/2512.10884}, 
}
```
