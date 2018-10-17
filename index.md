_Algoim_ is a collection of high-order accurate numerical methods and C++ algorithms for working with implicitly defined geometry and level set methods. Motivated by multi-phase multi-physics applications, particularly those with evolving dynamic interfaces, these algorithms target core, fundamental techniques in level set methods. They have been designed with a view to standard finite difference implementations as well as more advanced finite element and discontinuous Galerkin implementations, multi-threading and massively parallel MPI computation. The collection includes:
- High-order accurate quadrature algorithms for implicitly defined domains in hyperrectangles (such as for computing integrals in curved domains or over curved surfaces)
- High-order accurate closest point calculations to implicitly defined surfaces _(To be added)_
- k-d trees optimized for codimension-one point clouds _(To be added)_
- Accurate level set reinitialization and extension velocity schemes _(To be added)_
- Voronoi implicit interface methods for multi-phase interconnected interface dynamics _(To be added)_

## Prerequisites

Many of the numerical algorithms implemented in Algoim are templated on the spatial dimension `N`, allowing one to develop numerical schemes for any number of spatial dimensions. To assist with this functionality, Algoim makes use of the open source C++ library [blitz++](https://github.com/blitzpp/blitz) for high-performance fixed-length vector arithmetic, though `template<typename T, int N> blitz::TinyVector<T,N>`. Therefore, **[blitz++](https://github.com/blitzpp/blitz) should first be downloaded, configured for your compiler, and installed** such that it can be found in the appropriate include directories, i.e., so that `#include <blitz/array.h>` resolves correctly.

## Installation

Algoim is a header-only C++ library. Except for small example/demonstration applications, all of the files are C++ `.hpp` header files. As such, it requires minimal installation effort: simply download and configure so that the appropriate header driver can be found by your compiler when you include it in your C++ program, e.g., `#include "algoim/src/algoim_quad.hpp"`.

## High-Order Quadrature Algorithms for Implicitly Defined Domains

An implicitly defined domain is either a volumetric region or codimension-one surface whose shape is characterised implicitly by an isosurface of a continuous scalar function. A variety of applications involving implicitly defined geometry require the evaluation of integrals over such domains, including level set methods for propagating interfaces in computational physics, embedded boundary methods for solving partial differential equations on curved domains, and in treating jump conditions and singular source terms in weak formulations. A specific example is that of [implicit mesh discontinuous Galerkin methods](https://doi.org/10.1016/j.jcp.2017.04.076), which has been developed to facilitate high-order accurate modelling of interfacial fluid dynamics.

In practice, to calculate integrals over implicitly defined domains, a quadrature scheme must be computed. In [R. I. Saye, _High-Order Quadrature Methods for Implicitly Defined Surfaces and Volumes in Hyperrectangles_, SIAM Journal on Scientific Computing, 37(2), A993-A1019 (2015)](http://dx.doi.org/10.1137/140966290), a general purpose, high-order accurate quadrature algorithm has been developed. Based on the idea of converting the implicitly defined geometry into the graph of an implicitly defined height function, the approach requires only one-dimensional root finding and simple one-dimensional Gaussian quadrature schemes. These algorithms produce quadrature schemes with strictly positive quadrature weights and inherits the high-order accuracy of Gaussian quadrature, e.g., with 4 quadrature points (per dimension), 8th order accuracy can be achieved. Examples of generated quadrature schemes are shown in the figure.

<div style="width:502px; margin: 0 auto; font-size: 70%"><img src="img-quad.png"/><br/>Some examples illustrating the quadrature schemes constructed by the algorithm (based on a one-dimensional Gaussian quadrature scheme of order 4) for surface integrals (left column) and volume integrals (middle/right columns). The weights are coloured according to a scale that is normalised for each particular case: pale indicates a low-valued weight and dark blue indicates a high-valued weight.</div>

### Citing

If you make use of these quadrature algorithms in your research, in any documentation or publication please cite (in addition to _Algoim_) the original paper describing these algorithms:
- [R. I. Saye, _High-Order Quadrature Methods for Implicitly Defined Surfaces and Volumes in Hyperrectangles_, SIAM Journal on Scientific Computing, 37(2), A993-A1019 (2015)](http://dx.doi.org/10.1137/140966290)

### Using

main driver
templated on function object that implements operator and grad, both templated on type to do interval arithmetic

### Examples


## Advanced

QD



## Welcome to GitHub Pages

You can use the [editor on GitHub](https://github.com/algoim/algoim.github.io/edit/master/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/algoim/algoim.github.io/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
