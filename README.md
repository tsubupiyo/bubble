[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)
[![experimental](http://badges.github.io/stability-badges/dist/experimental.svg)](http://github.com/badges/stability-badges)
[![unstable](http://badges.github.io/stability-badges/dist/unstable.svg)](http://github.com/badges/stability-badges)


# bubble
3-Dimensional Voronoi Diagram Divided by Cones.

## Definition of Voronoi Diagram
In a metric space ![X](doc/fig/X.svg), the Voronoi diagram associated with a set of subset ![P](doc/fig/P.svg) is defined as the set of 
![Rk](doc/fig/Rk.svg).

In our calculation, the distance function is Euclidean.

## Starting Point
Although bubbling is the elegant algorithm, the data structure and algorithm are complex and heavy.
However, contacted bubbles hinted starting point of an algorithm described later.

![bubbling](doc/fig/Voronoi_growth_euclidean.gif)


## Algorithm
Our approach accepts a partial reconstruction.
Maybe a novel approach. 
Please contact me if you know a similar one.

### Discretization and Data Structure
Voronoi cells are represented as a set of cones in our discretization.
The center of the bottom of each cone is adjusted to match the boundary of the cell.
In the discretization, overlaps and gaps of cones occur.
The property that the actual cell and the discretized cell coincide with each other as the number of divisions increases is the same as when using a cubic voxel.
The compressibility of information is higher when using a cone.

![bubbling](doc/fig/discretization.jpeg)

Our algorithm requires a grid to determine the direction of the cone.
The grid points are distributed almost uniformly on the spherical surface.
In addition, all of the grid points are tethered to create a closed network, in order to speed up the calculation that determines the length of the vector extending from ![Pk](doc/fig/Pk.svg) to the bottom of each cone.
The grid point is managed as ![theta_phi](doc/fig/theta_phi.svg), which are same with angular coordinates in the spherical coordinate system.
The length between ![Pk](doc/fig/Pk.svg) and the center of bottom of cone is a function of ![theta_phi](doc/fig/theta_phi.svg), ![u_func](doc/fig/u_func.svg).

### Outline of Tesselation Procedure
#### 1. Find minimum ![uk](doc/fig/uk.svg)
- The minimum ![uk](doc/fig/uk.svg) is approximately equal to half of the distance to ![Pj](doc/fig/Pj.svg) closest to ![Pj](doc/fig/Pj.svg). Also, the direction ![theta_phi](doc/fig/theta_phi.svg) is direction to ![Pj](doc/fig/Pj.svg).

#### 2. Solve the others ![u_func](doc/fig/u_func.svg)

- Reference values for a time-optimization:

    - ![u_func_dash](doc/fig/u_func_dash.svg) : u for the neighbor of ![theta_phi](doc/fig/theta_phi.svg)
    - ![betaset_dash](doc/fig/betaset_dash.svg) : the set of parameter beta of distance function for each ![P](doc/fig/P.svg) and ![theta_phi_dash](doc/fig/theta_phi_dash.svg) direction
- Solve recursively
	1. Pop ![theta_phi](doc/fig/theta_phi.svg) from the stack
    2. If ![u_func](doc/fig/u_func.svg) is already solved, pop again
    	- When the stack size is zero, procedure is terminated
    3. Solve the distance functions (the initial parameter sets is set to ![betaset_dash](doc/fig/betaset_dash.svg) )
    4. Solve the first cross point of the distance functions (the inital cross point is set to ![u_func_dash](doc/fig/u_func_dash.svg) )
    5. Stack neighbor grid points of ![theta_phi](doc/fig/theta_phi.svg) using the network of grid points
    6. Stock the parameters of distance function as ![u_func_dash](doc/fig/u_func_dash.svg)
        - When ![theta_phi](doc/fig/theta_phi.svg) taken out from the stack, the stocked parameters are a parameters of distance function for neighbor of ![theta_phi](doc/fig/theta_phi.svg)
    
    Note: If ![u_func](doc/fig/u_func.svg) is negative or infinite, the ![Rk](doc/fig/Rk_simple.svg) is open.
    <img src="doc/fig/cross_point.png" width="600px">

## Contribution
Pull Request

1. Fork it ( https://github.com/toyaku-phys/bubble/fork )
2. Create your feature branch (git checkout -b my-new-feature)
3. Commit your changes (git commit -am 'Add some feature')
4. Push to the branch (git push origin my-new-feature)
5. Create a new Pull Request to the bubble/master branch

Issue

1. Write your new feature or bug report

## Submodules
- [misteltein/Levenberg-Marquardt](https://github.com/misteltein/Levenberg-Marquardt)
    - [eigenteam/eigen-git-mirror](https://github.com/eigenteam/eigen-git-mirror)
- [toyaku-phys/Chaperone](https://github.com/toyaku-phys/Chaperone)


## Versioning
We use [SemVer](http://semver.org/) for versioning. 
For the versions available, see the tags on this repository.

## Authors
* [**Hibiki Itoga**](https://github.com/misteltein) -Key programmer-
* [**yde**](https://github.com/master-yde) -Discussion partner-

## License
[MIT-license](LICENSE)