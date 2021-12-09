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

- Compute rot and tb of knots.
- Compute Z indices of chords and their cyclic words.
- Process disks into differentials for LCH and RSFT.
- Compute numbers of Z/2Z augmentations.
