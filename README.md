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

When you run `app.py` as above, a URL should appear which you can access from your web browser. Parameters `n_strands` and `crossings` (indicating front crossings) can be added. The `crossings` parameter gives a plat presentation of the link in the front projection, with crossing indices ranging from 0 to `n_strands` - 2, ordered from left to right.

Here is a screenshot for `http://127.0.0.1:5000/?n_strands=6&crossings=3,1,2,2,1,3,3&auto_dgas=rsft`, giving a [polyfillable link](https://arxiv.org/abs/1307.7998):

![image info](./main/static/screenshot.png)

Only use the web interface to compute augmentations for links with small numbers of crossings (say, < 30). You can monitor the terminal to see what computations are happening. For small numbers of crossings, computations of bilinearized Poincare polynomials may still take a while.

To visualize a plat diagram without computing any holomorphic disks, use a `lazy_disks=True` flag in your URL. For example `http://127.0.0.1:5000/?n_strands=6&crossings=3,1,2,2,1,3,3&lazy_disks=True`. You can visualize links of any size without putting much strain on your computer.

# Python interface

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

It can take a long time to compute augmentations and bilinearized homologies of DGAs. To store a DGA for later analysis, we use the `pickle` functions. Here is an example which assumes we have `pd` as above:
```
>>> lch = pd.lch_dga
>>> lch.pickle('my_favorite_dga.pk')
```
Now we can reload the dga with all of its computed data in a new session:
```
$ python
>>> import legendrian_links as ll
>>> lch = ll.DGA.from_pickle('my_favorite_dga.pk')
>>> lch.augmentations
... will show list of all augs ...
```
The pickle functionality only stores the data of a DGA, so that we can recover old DGAs even after our code has been updated.

It's often the case that the number of augmentations is too large to be stored in memory by any computer (eg. 2**100). If the number of augmentations is very large, they can still be computed and stored using ``compressed representations''. In the following example, we'll see that there are a small number of augmentations, meaning it is ok to ``decompress'' them to usual augmentations.
```
$ python
>>> import legendrian_links as ll
>>> front = [1 for _ in range(5)]
>>> pd = ll.PlatDiagram(n_strands=4, front_crossings=front, n_copy=1, lazy_disks=False, lazy_lch=True, lazy_rsft=True)
>>> pd.set_lch(lazy_augs=True, lazy_bilin=True)
>>> pd.lch_dga.set_augmentations(decompress=False)
2022-03-16 10:09:56,784|utils|INFO|Starting execution set_augmentations
...
2022-03-16 10:11:02,702|utils|INFO|Found 8 compressed augmentations of DGA
2022-03-16 10:11:02,702|utils|INFO|Found 21 (uncompressed) augmentations of DGA
2022-03-16 10:11:02,702|utils|INFO|Ending execution set_augmentations
>>> pd.lch_dga.decompress_augmentations()
```
Now we can proceed with computing bilinearized homologies, etc. We're working on more functionality to deal with large numbers of augmentations.

# Technical notes

## Threading and Groebner bases

Groebner basis computations used to search for augmentations can be very heavy and they are skipped if they take too long. This timeout functionality is very difficult to implement when using the web app (for threading reasons). In general, Groebner computations will be skipped whenever `pd.rsft_dga.set_augmentations(...)` or `polynomials.zero_set(...)` are called outside of the main thread.

## Recursion limits

We use [sympy](https://www.sympy.org/en/index.html) to encode and manipulate polynomials. The functionality for performing variable substitutions in `sympy` relies on a [recursive method](https://github.com/sympy/sympy/blob/master/sympy/polys/densebasic.py#L452) which we've seen exceed python's [default recursion depth limit](https://stackoverflow.com/questions/3323001/what-is-the-maximum-recursion-depth-in-python-and-how-to-increase-it) for polynomials in ~300 variables. The error message looks like
```
  File ".../lib/python3.6/site-packages/sympy/polys/densebasic.py", line 474, in dmp_to_tuple
    return tuple(dmp_to_tuple(c, v) for c in f)
  File ".../lib/python3.6/site-packages/sympy/polys/densebasic.py", line 474, in <genexpr>
    return tuple(dmp_to_tuple(c, v) for c in f)
  File ".../lib/python3.6/site-packages/sympy/polys/densebasic.py", line 474, in dmp_to_tuple
    return tuple(dmp_to_tuple(c, v) for c in f)
  File ".../lib/python3.6/site-packages/sympy/polys/densebasic.py", line 474, in <genexpr>
    return tuple(dmp_to_tuple(c, v) for c in f)
RecursionError: maximum recursion depth exceeded
```
To manually override this, add the following lines to your code:
```
import sys
sys.setrecursionlimit(10000)
```
I cannot guarantee that this is not a bad idea!

# To do list

## Features

- How do we more efficiently store augmentations? Maybe as `SubsNode`s? Currently it appears that compressed augs can expand enormously.
- Capping paths is storing too much info. We really only need this for rotation numbers.
- Get dual betti numbers of chain complexes by transposing matrices. Currently broken.
- Command line interface would make it easy to run scripts.
- Way behind on testing...
- Application is maintaining some variables between requests? I think this is due to mutable function args. Should be solved now.
- Is there any way to speed up the computations of poincare polynomials? This should boil down to speeding up `rref` computations.
- Grid -> plat algorithm. From grids could import the knot atlas or do algorithmic exploration. Difficult to enumerate links using plat presentations.
- Copy knot tables. To translate Sivek, do (y - 1).
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
