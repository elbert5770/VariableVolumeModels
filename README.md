# VariableVolumeModels
Variable volume compartmental models in Julia - improving computational accuracy of exact solutions.

The model is inspired by the following work: http://www.physics.uoi.gr/assimakopoulos/files/scientific_papers/39_Health_Physics_61_(1991)_245-253.pdf

The paper is in need of some corrections.  In equation 1, the $x_i$ in the first term on the right hand side should be $x_j$.  In equation 19, $G^{-1}C$ is unnecessary due to the use of a definite integral and should be just $y_0$.  Equation 23 should be zero.

At a fundamental level, the volumetric flow of liquid into the udder will carry radioisotope from the body to the milk.  The mass of radioisotope will decrease as milk leaves the sheep.  These flows really need to be in the equations and are added here.

x = concentration of cesium (Bq/vol)

y = amount of cesium (Bq)

V = volumes of compartments, a diagonal matrix

P = production rate (rate of input of cesium)

R = various reactions and flows between compartments

Amount is related to the concentration through the volumes
$\mathbf{y} = \mathbf{Vx}$

The form of the exact solution affects the accuracy of the solution, which involves numerical integration, matrix exponentiation and matrix inversion.  The system of equations can be solved in three ways:



The equation written as:

$\frac{d\mathbf{y}}{dt} = \mathbf{P} + \mathbf{RV}^{-1}\mathbf{y}$

may be solved with an integrating factor, and

$\frac{d\mathbf{Vx}}{dt} = \mathbf{P} + \mathbf{R}\mathbf{x}$

may be solved after applying the product rule to the LHS as a separable equation or solved with an integrating factor.

The interesting result is that during milking where the volume of the udder is decreasing, the solution derived from solving $\frac{d\mathbf{Vx}}{dt}$ with an integrating factor is more accurate than the other solutions.  During refilling, where the volume is increasing, the solution of $\frac{d\mathbf{y}}{dt}$ with an integrating factor provides a higher accuracy calculation of the solution.


