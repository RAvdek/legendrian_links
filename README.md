# Legendrian links

An interactive application for analyzing Legendrian links. Or it will be some day, maybe.

Much of this recreates [Sivek's lch.sage](https://www.ma.imperial.ac.uk/~ssivek/code/lch.sage) so that it will be applicable to generalizations of the LCH algebra. We also want to make analysis of knots accessible through a web application (with minimal dependencies) so that users don't need to know how to code.

# Installation

To install and run the server from your terminal from the root of the project:

```
python3 -m pip install pipenv
pipenv install
pipenv shell
cd main
python app.py
```

A URL should appear which you can access from your web browser. Parameters `n_strands` and `crossings` (indicating front crossings) can be added. Here is a screenshot for `http://127.0.0.1:5000/?n_strands=6&crossings=3,1,2,2,1,3,3`, giving a [polyfillable link](https://arxiv.org/abs/1307.7998):

![image info](./main/static/screenshot.png)

# To do list

## Features

- Process disks into differentials for RSFT with Z/2Z coefficients and grading. Then put Z/2Z augs into HTML table.
- Create ChainComplex class which we'll use to store bilinearized invariants. We only care to work with bounded complexes over QQ or finite fields so that homology groups will be determined by their ranks. In this case, we can store kernels of differentials and compute betti numbers as b_n = dim ker del_n - dim im del_{n-1} = dim ker del_n + rank C_{n-1} - dim ker del_{n-1}. Note that computing dimensions of kernels in sympy is not already built. It seems we'd need to develop an rref algorithm by hand for finite fields.
- Check if two augmentations are homotopic or not, by seeing if the bilinearized homology has non-zero hom to the base field.
- Put bilinearized homologies of Poincare polynomials into an HTML table.
- Introduce t coordinate in differentials. (This is not so important for augmentations where we can use t=1 for Z/2Z coeffs).
- Ability to flip orientations of link components.
- Ability to reverse orientations on link components.
- Orientations: Process disks into differentials for LCH with Z coefficients.
- Carry out orientation processing for RSFT differentials.
- Algorithmically determine a plat diagram from a grid diagram.
- Compute differentials for 2-copies and twisted 2-copies.
- In the we interface, users can do basic input links, flip components, etc with HTML / javascript input.

## Cleanup, performance, and testing

- dga.py test cases needed for more polynomials.
- dga.py needs tests for differentials and DGA class methods.
- Add pylint or something to ensure code cleanliness.
- Review relationships between data structures.
- Make `legendrian_links` importable as a python library.

## Data sets

- Port knots from the Legendrian knot atlas or other resource into `links.json`.
- Get two component links from the link atlas.
