# Nonbacktracking Spectral Clustering of Nonuniform Hypergraphs

*[Phil Chodrow](https://www.philchodrow.com), [Nicole Eikmeier](https://eikmeier.sites.grinnell.edu), and [Jamie Haddock](https://jamiehaddock.com).*

This repository houses code used in the manuscript "Nonbacktracking spectral clustering of nonuniform hypergraphs." Julia is used for primary computations, while R is used for data visualization. 

If you use this code in published work, please include the following bibtex entries: 

```
@article{chodrow2022nonbacktracking,
  title = {Nonbacktracking Spectral Clustering of Nonuniform Hypergraphs},
  author = {Chodrow, Philip S. and Eikmeier, Nicole and Haddock, Jamie L.},
  year = {2022},
  journal = {TBD},
  volume = {TBD},
  number = {TBD},
  pages = {TBD},
  publisher = {{TBD}}
}

@article{bezanson2017julia,
  title = {Julia: {{A}} Fresh Approach to Numerical Computing},
  author = {Bezanson, Jeff and Edelman, Alan and Karpinski, Stefan and Shah, Viral B},
  year = {2017},
  journal = {SIAM Review},
  volume = {59},
  number = {1},
  pages = {65--98},
  publisher = {{SIAM}}
}
```

## Repository Outline

This repository contains the following directories: 

- `src`: the primary source code of the project, organized as a Julia package called `HypergraphNB` (`NB` = "NonBacktracking"). The source code implements a `hypergraph` struct, functions for computing various nonbacktracking matrices, inference of parameters, and simple spectral algorithms.  
- `test`: Unit tests for the package. The primary purpose of these unit tests is to verify identities related to matrix reduction formulae. To run unit tests: 
```
>>> cd HypergraphNB
>>> julia

julia> ]

pkg> activate .
pkg> test
```
- `scripts`: these files perform the primary computations. The `makefile` provides an automated description of the sequence in which they should be run in order to produce the figures used in the manuscript. To produce all such figures, run `make figs` in the top-level directory after installing all necessary Julia and R packages. **Warning**: as written, several of the scripts require hours or days of computation time. Running as-is on personal equipment is not recommended. 
- `data`: contains data used in empirical experiments. These data sets were retrieved from [Austin Benson's data repository](https://www.cs.cornell.edu/~arb/data/). 
- `notebooks`: this directory primarily contains deprecated Jupyter notebooks related to early versions of code development and experimentation. The notebooks in this directory are not expected to run or be comprehensible, and may be broken to an arbitrary extent. 


## Requirements

The file `Manifest.toml` describes the Julia requirements of `HypergraphNB.jl` and the files in the `scripts` directory. In order to create the visualizations, an installation of R with the following R packages is also required: 

- `tidyverse`
- `patchwork`
- `ggrepel`
- `colorspace`
- `ggforce`