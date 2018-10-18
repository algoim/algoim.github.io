_Algoim_ is a collection of high-order accurate numerical methods and C++ algorithms for working with implicitly defined geometry and level set methods. Motivated by multi-phase multi-physics applications, particularly those with evolving dynamic interfaces, these algorithms target core, fundamental techniques in level set methods. They have been designed with a view to standard finite difference implementations as well as more advanced finite element and discontinuous Galerkin implementations, multi-threading and massively parallel MPI computation. The collection includes:
- High-order accurate quadrature algorithms for implicitly defined domains in hyperrectangles (such as for computing integrals in curved domains or over curved surfaces)
- High-order accurate closest point calculations to implicitly defined surfaces _(To be added)_
- k-d trees optimized for codimension-one point clouds _(To be added)_
- Accurate level set reinitialization and extension velocity schemes _(To be added)_
- Voronoi implicit interface methods for multi-phase interconnected interface dynamics _(To be added)_

## Prerequisites

Many of the numerical algorithms implemented in Algoim are templated on the spatial dimension `N`, allowing one to develop numerical schemes for any number of spatial dimensions. To assist with this functionality, Algoim makes use of the open source C++ library [blitz++](https://github.com/blitzpp/blitz) for high-performance fixed-length vector arithmetic, though `template<typename T, int N> blitz::TinyVector<T,N>`. Therefore, **[blitz++](https://github.com/blitzpp/blitz) should first be downloaded, configured for your compiler, and installed** such that it can be found in the appropriate include directories, i.e., so that `#include <blitz/array.h>` resolves correctly.

In addition, a C++ compiler supporting standard C++14 or higher is required.

## Installation

Algoim is a header-only C++ library. Except for small example/demonstration applications, all of the files are C++ `.hpp` header files. As such, it requires minimal installation effort: simply download and configure so that the appropriate header driver can be found by your compiler when you include it in your C++ program, e.g., `#include "algoim/src/algoim_quad.hpp"`.

## High-Order Quadrature Algorithms for Implicitly Defined Domains

An implicitly defined domain is either a volumetric region or codimension-one surface whose shape is characterised implicitly by an isosurface of a continuous scalar function. A variety of applications involving implicitly defined geometry require the evaluation of integrals over such domains, including level set methods for propagating interfaces in computational physics, embedded boundary methods for solving partial differential equations on curved domains, and in treating jump conditions and singular source terms in weak formulations. A specific example is that of [implicit mesh discontinuous Galerkin methods](https://doi.org/10.1016/j.jcp.2017.04.076), which has been developed to facilitate high-order accurate modelling of interfacial fluid dynamics.

In practice, to calculate integrals over implicitly defined domains, a quadrature scheme must be computed. In [R. I. Saye, _High-Order Quadrature Methods for Implicitly Defined Surfaces and Volumes in Hyperrectangles_, SIAM Journal on Scientific Computing, 37(2), A993-A1019 (2015)](http://dx.doi.org/10.1137/140966290), a general purpose, high-order accurate quadrature algorithm has been developed. Based on the idea of converting the implicitly defined geometry into the graph of an implicitly defined height function, the approach requires only one-dimensional root finding and simple one-dimensional Gaussian quadrature schemes. These algorithms produce quadrature schemes with strictly positive quadrature weights and inherits the high-order accuracy of Gaussian quadrature, e.g., with 4 quadrature points (per dimension), 8th order accuracy can be achieved. Examples of generated quadrature schemes are shown in the figure.

<div style="width:502px; margin: 0 auto; font-size: 70%"><img src="img-quad.png"/><br/>Some examples illustrating the quadrature schemes constructed by the algorithm (based on a one-dimensional Gaussian quadrature scheme of order 4) for surface integrals (left column) and volume integrals (middle/right columns). The weights are coloured according to a scale that is normalised for each particular case: pale indicates a low-valued weight and dark blue indicates a high-valued weight.</div>

### Citing

If you make use of these quadrature algorithms in your research, in any documentation or publication please cite (in addition to _Algoim_) the original paper describing these algorithms:
- [R. I. Saye, _High-Order Quadrature Methods for Implicitly Defined Surfaces and Volumes in Hyperrectangles_, SIAM Journal on Scientific Computing, 37(2), A993-A1019 (2015)](http://dx.doi.org/10.1137/140966290)

### Usage

The driver header file for the quadrature algorithms is located at `algoim/src/algoim_quad.hpp`. There is only one driver routine, `quadGen` which is templated on the level set function object `phi` and the dimension `N`:

```C++
template<typename F, int N>
Algoim::QuadratureRule<N> Algoim::quadGen(const F& phi, const Algoim::BoundingBox<Real,N>& xrange,
                                          int dim, int side, int qo);
```
The provided parameters are as follows:
- `phi` is a user-defined function object which evaluates the level set function and its gradient. It must implement both `template<typename T> operator() (const blitz::TinyVector<T,N>& x) const` and `template<typename T> grad(const blitz::TinyVector<T,N>& x) const`. In the simplest case, `T = double` and the role of `phi` is to simply evaluate the level set function (e.g., `return x(0)*x(0) + x(1)*x(1) - 1;` for a unit circle) and its gradient (e.g., `return TinyVector<double,2>(2.0*x(0), 2.0*x(1));`). However, it is crucial that these two member functions be templated on `T` in order to enable the interval arithmetic underlying the algorithms of the paper cited above. In essence, the interval arithmetic automatically computes a first-order Taylor series (with bounded remainder) of the given level set function, and uses that to make decisions concerning the existence of the interface and what direction to use when converting the implicitly defined geometry into the graph of an implicitly defined height function. This requirement on `phi` being able to correctly perform interval arithmetic places restrictions on the type of level set functions `quadGen` can be applied to; these restrictions are discussed later.
- `xrange` is a user-specified bounding box, indicating the extent of the hyperrectangle in `N` dimensions to which the quadrature algorithm is applied.
- `dim` is used to specify the type of quadrature:
    - If `dim < 0`, compute a **volumetric quadrature scheme**, whose domain is implicitly defined by `{phi < 0}` intersected with `xrange`.
    - If `dim == N`, compute a **curved surface quadrature scheme**, whose domain is implicitly defined by `{phi == 0}` intersected with `xrange`.
    - If `0 <= dim && dim < N`, compute a **flat surface quadrature schme** for one of the sides of the hyperrectangle, i.e., `{phi < 0}` intersected with `xrange` intersected with the side `{x(dim) == xrange(side)(dim)`.
- `side` is used only when `0 <= dim && dim < N` and specifies which side, either `side == 0` or `side == 1` for the "left" or "right", respectively.
- `qo` specifies the degree of the underlying one-dimensional Gaussian quadrature scheme, e.g., `qo = 4` in the figure above. `qo` must satisfy `1 <= qo && qo <= 10`.

The output of `quadGen` is an `Algoim::QuadratureRule<N>` object. This object is essentially a `std::vector` listing the set of quadrature points (as `blitz::TinyVector<Real,N>`) and their associated weights. A `QuadratureRule` object is also a function object -- its associated `operator()(const F& f)` can be applied to any user-specified function to compute the result of applying the quadrature rule to a functional.

**Requirements on `phi`.** As mentioned above, the `phi` function object must be able to be applied to arguments of type `blitz::TinyVector<Algoim::Interval<N>,N>`. Here, `Algoim::Interval<N>` is a special type whose purpose is to calculate a first-order Taylor series with bounded remainder, and shares concepts in common with _automatic differentiation_. `Interval<N>` implements `operator+`, `operator*`, `operator/`, etc., and can be used in a variety of ways, most commonly in evaluating polynomial expressions. Common unary operators are also implemented, e.g., `sin(Interval<N>)`. However, one cannot straightforwardly apply `Interval<N>` arithmetic to max or min statements, if conditions, and so forth. This relates to the fact that it is very difficult and subtle to compute quadrature schemes for shapes which have corners, holes, or cusps, etc. In using Algoim quadrature schemes for the first time, it is recommended they be applied to smooth level set functions, made out of polynomial expressions and common smooth functions like `sin`, `exp`, etc.

### Examples

The quadrature algorithms of Algoim are here demonstrated with a level set function describing an ellipse (in `N = 2` dimensions) or ellipsoid (in `N = 3` dimensions). First, we define a function object implementing the requirements above for the level set function:

```
template<int N>
struct Ellipsoid
{
    template<typename T>
    T operator() (const blitz::TinyVector<T,N>& x) const
    {
        if (N == 2)
            return x(0)*x(0) + 4.0*x(1)*x(1) - 1.0;
        else
            return x(0)*x(0) + 4.0*x(1)*x(1) + 9.0*x(2)*x(2) - 1.0;
    }

    template<typename T>
    blitz::TinyVector<T,N> grad(const blitz::TinyVector<T,N>& x) const
    {
        if (N == 2)
            return blitz::TinyVector<T,N>(2.0*x(0), 8.0*x(1));
        else
            return blitz::TinyVector<T,N>(2.0*x(0), 8.0*x(1), 18.0*x(2));
    }
};
```

## Advanced

QD
