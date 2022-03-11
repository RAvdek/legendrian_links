# Legendrian links

An interactive application for analyzing Legendrian links. Or it will be some day, maybe.

Much of this recreates [Sivek's lch.sage](https://www.ma.imperial.ac.uk/~ssivek/code/lch.sage) so that it will be applicable to generalizations of the LCH algebra. We also want to make analysis of knots accessible through a web application (with minimal dependencies) for easy visual inspection of link diagrams.

The program computes augmentations of LCH algebras as well as for a version of Legendrian RSFT (which is not on the arXiv yet). The program uses plat diagrams to represent Legendrians in R3. Plats make [holomorphic disks particularly nice](https://arxiv.org/abs/2104.00505) (helping with algorithmic computation), although putting a Legendrian in plat position will typically introduce extra crossings (hurting our ability to algorithmically compute).

# Installation

To install and run the server from your terminal from the root of the project:

```
$ python3 -m pip install pipenv
$ python3 -m pipenv install
$ python3 -m pipenv shell
$ cd main
$ python app.py
```

You may need to modify the version of python in Pipfile so that this is compatible with your system.

# Web interface

When you run `app.py` as above, a URL should appear which you can access from your web browser. Parameters `n_strands` and `crossings` (indicating front crossings) can be added. Here is a screenshot for `http://127.0.0.1:5000/?n_strands=6&crossings=3,1,2,2,1,3,3&auto_dgas=rsft`, giving a [polyfillable link](https://arxiv.org/abs/1307.7998):

![image info](./main/static/screenshot.png)

Only use the web interface to compute augmentations for links with small numbers (say, < 30) of crossings. To visualize a plat diagram without computing any holomorphic disks, use a `lazy_disks=True` flag in your URL. For example `http://127.0.0.1:5000/?n_strands=6&crossings=3,1,2,2,1,3,3&lazy_disks=True`. You can visualize links of any size without putting much strain on your computer.

# Python interface & batch processing

Computing augmentations of DGAs with large numbers of generators can take hours or be impossible due to memory constraints. Even computing gradings for DGAs can be time consuming.

Run the following code from inside the `main/` folder to setup the `PlatDiagram` object:
```
$ python
>>> import legendrian_links as ll
>>> front = [1 for _ in range(21)]
>>> pd = ll.PlatDiagram(n_strands=4, front_crossings=front, n_copy=2, lazy_disks=False, lazy_lch=True, lazy_rsft=True)
>>> pd.set_lch(lazy_augs=True, lazy_bilin=True)
>>> pd.lch_dga.set_augmentations()
>>> pd.lch_dga.set_all_bilin()
```
The `lazy...` options can prevent the kick-off of some potentially heavy computations. Some experimental options for `set_augmentations` are trying to help speed up computations (see the code). We can also fun the above with less commands:
```
$ python
>>> import legendrian_links as ll
>>> front = [1 for _ in range(21)]
>>> pd = ll.PlatDiagram(n_strands=4, front_crossings=front, n_copy=2, lazy_disks=False, lazy_lch=False, lazy_rsft=True)
```

# Technical notes on threading and Groebner bases

Groebner basis computations used to search for augmentations can be very heavy and they are skipped if they take too long. This timeout functionality is very difficult to implement when using the web app (for threading reasons). In general, Groebner computations will be skipped whenever `pd.rsft_dga.set_augmentations(...)` or `polynomials.zero_set(...)` are called outside of the main thread.

# To do list

## Features

- Get dual betti numbers of chain complexes by transposing matrices. Currently broken.
- Command line interface would make it easy to run scripts.
- Way behind on testing...
- Application is maintaining some variables between requests? I think this is due to mutable function args. Should be solved now.
- Is there any way to speed up the computations of poincare polynomials? This should boil down to speeding up `rref` computations.
- Grid -> plat algorithm. From grids could import the knot atlas or do algorithmic exploration. Difficult to enumerate links using plat presentations.
- Copy knot tables. Have to remember how to translate Sivek front crossing notation to mine.
- UI: Ordering of generators is annoyingly out of place. Should also count numbers of augs.
- Check if two augmentations are homotopic or not, by seeing if the bilinearized homology has non-zero hom to the base field.
- Make tables nicer using some JS library. Big tables can be condensed. It would also be nice to sort data.
- Introduce t coordinate in differentials. (This is not so important for augmentations where we can use t=1 for Z/2Z coeffs).
- Ability to flip orientations of link components. Currently upper-left corner always points right.
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
