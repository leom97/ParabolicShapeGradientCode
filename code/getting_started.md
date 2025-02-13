# Code Usage

The main purpose of this code is to produce the results presented in Table 1 in the article.

The main script is `./shape_gradients_ooc/shape_gradients_ooc_verification.py`.
It relies on `./shape_gradients_ooc/configuration.py`, where the data necessary to set up each run is collected.

Running `shape_gradients_ooc_verification.py` with the default `configuration.py` yields the first row of Table 1,
where the computed orders of convergence are printed to the console.

In `configuration.py` further instructions are given on how to obtain the other rows.
