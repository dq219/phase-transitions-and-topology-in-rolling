# Phase transitions and topology in rolling

PNAS link: https://www.pnas.org/doi/10.1073/pnas.2417161122
ArXiv link: https://arxiv.org/abs/2407.19861

This code repository contains details of constructing irregular 2D and 3D bodies using Fourier modes, and simulation scripts for their rolling motion under an external force on an inclined ramp. Non-slip constraints are assumed, and rotational viscosity and moment of inertia are isotropic for simplicity. The Jupyter notebooks contain sample code showing how the libraries work.

- 2D rolling scripts
	- [2D_overdamped.ipynb](./2D/2D_overdamped.ipynb): 2D rolling without inertia. Second-order halting-cruising phase transition is observed.
	- [2D_inertial.ipynb](./2D/2D_inertial.ipynb): first-order halting-cruising transition in 2D rolling with inertia, and an associated first-order cruising-halting transition as a hysteresis effect.
	
- 3D rolling scripts
	- [3D_overdamped.ipynb](./3D/3D_overdamped.ipynb): simulation details for 3D overdamped rolling, generates various plots for data visualisation and numerical error quality control. The 'Q factor' is the quality factor as defined in the SI, not the quaternion.
	- [3D_inertial.ipynb](./3D/3D_inertial.ipynb): simulations for 3D rolling with non-zero mass, with a focus on comparing long-time trajectories at different masses. Period-doubling due to the doubly-connectedness of *SO(3)* is observed.
	- 3D_advanced series: analysis script that first performs many simulations for overdamped/inertial rolling, then computes:
		- [3D_advanced_Poincare.ipynb](./3D/3D_advanced_Poincare.ipynb): computation of the Poincare section in overdamped rolling;
		- [3D_advanced_double.ipynb](./3D/3D_advanced_double.ipynb): quaternion illustration for period-doubling;
		- [3D_advanced_tube.ipynb](./3D/3D_advanced_tube.ipynb): how a point surface perturbation divides a 3-ball into two parts.

In addition, curated experimental data and the analysis/plotting scripts are provided in the [Experiments folder](./Experiments).

Please cite the [publication](https://arxiv.org/abs/2407.19861](https://www.pnas.org/doi/10.1073/pnas.2417161122 ) if you find the code helpful for your work. Get in touch with dq219@cam.ac.uk or 000.gdz@gmail.com for any questions.

Thanks for your interest! :)

## Dependencies
Libraries used:
- `numpy`
- `scipy`
- `sympy`
- `matplotlib`
- `pandas`
- `numpy-quaternion` (if installing using `conda` instead of `pip`, the package is named just `quaternion` - I think)

## License
This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.

