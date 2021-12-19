# Legendrian links

An interactive application for analyzing Legendrian links. Or it will be some day, maybe.

# Installation

To install and run the server from your terminal from the root of the project:

```
python3 -m pip install pipenv
pipenv install
pipenv shell
cd main
python app.py
```

A URL should appear which you can access from your web browser. Parameters `n_strands` and `crossings` (indicating front crossings) can be added. Here is a screenshot for `http://127.0.0.1:5000/?n_strands=4&crossings=1,1`, giving a Hopf link:

![image info](./main/static/screenshot.png)

# To do list

## Features

- Introduce t coordinate in differentials.
- Process disks into differentials for RSFT with Z/2Z coefficients and grading.
- Compute numbers of Z/2Z augmentations. Should be in a new module.
- Process disks into differentials for LCH and RSFT with Z coefficients.
- Ability to reverse orientations on link components.
- Compute differentials for 2-copies and twisted 2-copies.

## Data sets

- Port knots from the Legendrian knot atlas or other resource into `links.json`.
- Get two component links from the link atlas.

## Code cleanup

- Add pylint or something to ensure code cleanliness.
- Review relationships between data structures.
- Make `legendrian_links` importable as a python library.
