# Legendrian links

An interactive application for analyzing Legendrian links. Or it will be some day, maybe.

Much of this recreates [Sivek's lch.sage](https://www.ma.imperial.ac.uk/~ssivek/code/lch.sage) so that it will be applicable to generalizations of the LCH algebra. We also want to make analysis of knots accessible through a web application (with minimal dependencies) so that users don't need to know how to code.

# Installation

To install and run the server from your terminal from the root of the project:

```
$ python3 -m pip install pipenv
$ pipenv install
$ pipenv shell
$ cd main
$ python app.py
```

# Web interface

When you run `app.py` as above, a URL should appear which you can access from your web browser. Parameters `n_strands` and `crossings` (indicating front crossings) can be added. Here is a screenshot for `http://127.0.0.1:5000/?n_strands=6&crossings=3,1,2,2,1,3,3`, giving a [polyfillable link](https://arxiv.org/abs/1307.7998):

![image info](./main/static/screenshot.png)

# Python interface & batch processing

Only use the web interface to compute augmentations for links with small numbers (say, < 30) of crossings. Computing augmentations of DGAs with large numbers of generators can take hours or be impossible due to memory constraints. We use batch processing to deal with these cases.

Here is an example. Suppose we want to compute augmentations of the RSFT algebra of the 2-copy of `m(9_46)`. Run the following code from inside the `main/` folder to setup the `PlatDiagram` object:

```
$ python
>>> import legendrian_links as ll
>>> front = [1,2,3,4,2,0,1,4,1,2,0,3,1,3,2,4,1,3]
>>> pd = ll.PlatDiagram(n_strands=6, front_crossings=front, n_copy=2, lazy_disks=False, lazy_lch=True, lazy_rsft=True)
>>> pd.set_rsft(lazy_augs=True, lazy_bilin=True)
```
The `lazy...` options prevent the program from jumping into too many heavy calculations. A bunch of data will be logged to the command line:
```
2022-03-08 10:10:00,488|utils|INFO|DGA has 180 generators with 0 grading
2022-03-08 10:10:00,492|utils|INFO|Differentials required to compute augs: 304
...
2022-03-08 10:13:31,488|utils|INFO|Aug analysis table:
{'symbol': {22}, 'freq': 140, 'n_cum_polys': 0, 'poly_to_sym_ratio': 0.0}
...
{'symbol': {13}, 'freq': 40, 'n_cum_polys': 7, 'poly_to_sym_ratio': 0.5}
...
{'symbol': {4,23}, 'freq': 6, 'n_cum_polys': 100, 'poly_to_sym_ratio': 1.2345679012345678}
...
{'symbol': {43,4}, 'freq': 4, 'n_cum_polys': 251, 'poly_to_sym_ratio': 1.56875}
...
{'symbol': {71,32}, 'freq': 2, 'n_cum_polys': 304, 'poly_to_sym_ratio': 1.6888888888888889}
```
This tells us that in order to enumerate all augmentations of the DGA, we need to find the common zero set of 304 polynomials in 180 variables over Z/2Z. Batch processing allows us to solve a subset of the polynomials first, and then iteratively add more polynomials. Looking at 'poly_to_sym_ratio' in the logs, if we try work with the first 7 polynomials, we will have a very sparse problem (with more variables then polynomials). If we work with the first 251 polynomials, the problem will be less sparse but this might be too much data to process at once. We'll go with 100:
```
>>> pd.rsft_dga.set_augmentations(batch_size=100)
2022-03-08 10:15:20,450|utils|INFO|Starting execution set_augmentations
...
2022-03-08 11:25:07,513|utils|INFO|95th loop of zero_set search: 95 nodes, 95 spawned nodes, 0 solution nodes
2022-03-08 11:25:07,584|utils|INFO|Ending execution zero_set
2022-03-08 11:25:07,584|utils|INFO|Ending execution batch_zero_set
2022-03-08 11:25:07,584|utils|INFO|Found 0 augmentations of DGA
```
After about 70 minutes we see that the RSFT DGA has no augmentations! This computation would have been impossible without batch processing due to memory constraints of my laptop.

# To do list

## Features

- Some threads used for Groebner appear to never die.
- During set spawn, we take a cartesian product over all of the sunset variables. We should avoid this so we don't spawn billions of nodes. Currently implemented but needs to be undoable from `algebra.py` or `polynomials.py`.
- Application is maintaining some variables between requests.
- Is there any way to speed up the computations of poincare polynomials? This should boil down to speeding up `rref` computations.
- Grid -> plat algorithm. From grids could import the knot atlas or do algorithmic exploration. Difficult to enumerate links using plat presentations.
- Copy knot tables. Have to remember how to translate Sivek front crossing notation to mine.
- UI: Ordering of generators is annoyingly out of place. Should also count numbers of augs.
- Check if two augmentations are homotopic or not, by seeing if the bilinearized homology has non-zero hom to the base field.
- Should be able to use Groebner bases to tell if the commutative version of a DGA is trivial with Z/2 ceoffs. Implement in dga.py.
- Make tables nicer using some JS library. Big tables can be condensed.
- Introduce t coordinate in differentials. (This is not so important for augmentations where we can use t=1 for Z/2Z coeffs).
- Ability to flip orientations of link components.
- Ability to reverse orientations on link components.
- Orientations: Process disks into differentials for LCH with Z coefficients.
- Carry out orientation processing for RSFT differentials.
- Algorithmically determine a plat diagram from a grid diagram.
- Compute differentials for 2-copies and twisted 2-copies.
- In the we interface, users can do basic input links, flip components, etc with HTML / javascript input.

## Cleanup, performance, and testing

- Tests for different modules should have their own files.
- More test cases for matrices.
- dga.py test cases needed for more polynomials.
- dga.py needs tests for differentials and DGA class methods.
- Add pylint or something to ensure code cleanliness.
- Review relationships between data structures.
- Make `legendrian_links` importable as a python library.

## Data sets

- Port knots from the Legendrian knot atlas or other resource into `links.json`.
- Get two component links from the link atlas.
