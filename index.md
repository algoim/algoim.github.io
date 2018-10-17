_Algoim_ is a collection of high-order accurate numerical methods and C++ algorithms for working with implicitly defined geometry and level set methods. Motivated by multi-phase multi-physics applications, particularly those with evolving dynamic interfaces, these algorithms target core, fundamental techniques in level set methods. They have been designed with a view to standard finite difference implementations as well as more advanced finite element and discontinuous Galerkin implementations, multi-threading and massively parallel MPI computation. The collection includes:
- High-order accurate quadrature algorithms for implicitly defined domains in hyperrectangles (such as for computing integrals in curved domains or over curved surfaces)
- High-order accurate closest point calculations to implicitly defined surfaces _(To be added)_
- k-d trees optimized for codimension-one point clouds _(To be added)_
- Accurate level set reinitialization and extension velocity schemes _(To be added)_
- Voronoi implicit interface methods for multi-phase interconnected interface dynamics _(To be added)_

## Prerequisites

Many of the numerical algorithms implemented in Algoim are templated on the spatial dimension `N`, allowing one to develop numerical schemes which operate in any number of spatial dimensions. To assist with this functionality, Algoim makes heavy use of the open source C++ library [blitz++](https://github.com/blitzpp/blitz) for high-performance fixed-length vector arithmetic, via `template<typename T, int N> blitz::TinyVector<T,N>`. In order to use Algoim, [blitz++](https://github.com/blitzpp/blitz) should first be downloaded and configured for your compiler, and installed such that it can be found in the appropriate include directories, i.e., so that `#include <blitz/array.h>` resolves correctly.

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
