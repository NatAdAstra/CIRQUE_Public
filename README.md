# CIRQUE

## Completeness Information for high Redshift QUEsts

Currently in closed development, CIRQUE will be a python toolset to allow for the calculation of completeness corrections in astrophysical galaxy surveys. A public version should be available in late 2018.

CIRQUE has the same foundations as GLACiAR ( https://github.com/danielacarrasco/GLACiAR , https://arxiv.org/pdf/1805.08985.pdf ). However it takes a slightly different approach to some of the optimisations and assumptions that get made.

### CIRQUE when complete will do the following:

-Generate a library of galaxy luminosity profile templates, these can have a variety of tuned parameters from sersec index to angle of observation. These templates are generated on the pixel scale of the science images.

-Iteratively populate your science images with these templates and test extraction success rates. Priors will be available to be tuned to allow the population of galaxies tested for extraction to be fine tuned. e.g. a prior on how many galaxies are of particular sersec indicies, size and redshift.

-Extract inserted galaxies with SExtractor and record success rates depending on the input parameters. Allowing for completeness corrections to be calculated.

### Compared to GLACiAR, CIRQUE will do the following differently.

+generate non-centered galaxies. CIRQUE will produce galaxy templates where the center of the luminosity profile is not centered in the middle of a pixel. Allowing for greater realism and to greater test source extraction of faint objects when the central bright peak of the profile is spread over many pixels.

+'Color' galaxies with SED's chosen from a folder of SED templates. Instead of assuming galaxies are a lyman break and power law continuum. CIRQUE can use real or simulated SED's along with chosen filter transparancies to calculate brightnesses in any filtered image. Useful for work on Photometric Redshifts.

-The 2 above processes consequently make CIRQUE slower and more computationally intensive than GLACiAR. To compensate, an option to parallelise where possible will be available. Instead of generating galaxy templates on the fly like GLACiAR, CIRQUE saves them in order to reduce repeat calculations and allow one time to play more with priors and image incersion.

+CIRQUE will allow for full control over the priors for how galaxies are selected from the template library to be placed into the images. Only care about elliptical galaxies? select only high sersic indicies. Want to see how many galaxies would drop out if the luminosity function was slightly different? You can do that too.
