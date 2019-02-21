#! /usr/bin/env python

# Convert fits and VO catalogues to Andre's sky model format

import os, sys
from string import replace

import numpy as np
# tables and votables
import astropy.io.fits as fits
from astropy.io.votable import parse_single_table

from optparse import OptionParser

__author__ = ['Natasha Hurley-Walker', 'Paul Hancock']
__date__ = '2019-02-21'

usage = "Usage: %prog [options] <file>\n"
parser = OptionParser(usage=usage)
parser.add_option('--catalogue', type="string", dest="catalogue",
                  help="The filename of the catalogue you want to read in.", default=None)
parser.add_option('--output', type="string", dest="output",
                  help="The filename of the output (default=test.txt).", default='test.txt')
parser.add_option('--namecol', type="string", dest="namecol",
                  help="The name of the Name column (default=Name).", default="Name")
parser.add_option('--racol', type="string", dest="racol",
                  help="The name of the RA column (default=ra_str).", default="ra_str")
parser.add_option('--decol', type="string", dest="decol",
                  help="The name of the Dec column (default=dec_str).", default="dec_str")
parser.add_option('--acol', type="string", dest="acol",
                  help="The name of the major axis column (default=a_wide).", default="a")
parser.add_option('--bcol', type="string", dest="bcol",
                  help="The name of the minor axis column (default=b_wide).", default="b")
parser.add_option('--pacol', type="string", dest="pacol",
                  help="The name of the position angle column (default=pa_wide).", default="pa")
parser.add_option('--fluxcol', type="string", dest="fluxcol",
                  help="The name of the flux density column (default=peak_flux).", default="peak_flux")
parser.add_option('--fcol', type='string', dest="fcol",
                  help="The name of the frequecy column (default = freq).",
                  default='freq')
parser.add_option('--alphacol', type="string", dest="alphacol",
                  help="The name of the spectral index alpha-term column (default=alpha).", default="alpha")
# parser.add_option('--betacol',type="string", dest="betacol",
#                    help="The name of the spectral index beta-term column (default=beta).", default="beta")
(options, args) = parser.parse_args()

if options.catalogue is None:
    print "must specify input catalogue"
    sys.exit(1)
else:
    filename, file_extension = os.path.splitext(options.catalogue)
    if file_extension == ".fits":
        temp = fits.open(options.catalogue)
        data = temp[1].data
    elif file_extension == ".vot":
        temp = parse_single_table(options.catalogue)
        data = temp.array

# Generate an output in Andre's sky model format
gformatter = "source {{\n  name \"{Name:s}\"\n  component {{\n    type {shape:s}\n    position {RA:s} {Dec:s}\n    shape {a:2.1f} {b:2.1f} {pa:4.1f}\n    sed {{\n      frequency {freq:3.0f} MHz\n      fluxdensity Jy {flux:4.7f} 0 0 0\n      spectral-index {{ {alpha:2.2f} {beta:2.2f} }}\n    }}\n  }}\n}}\n"
pformatter = "source {{\n  name \"{Name:s}\"\n  component {{\n    type {shape:s}\n    position {RA:s} {Dec:s}\n    sed {{\n      frequency {freq:3.0f} MHz\n      fluxdensity Jy {flux:4.7f} 0 0 0\n      spectral-index {{ {alpha:2.2f} {beta:2.2f} }}\n    }}\n  }}\n}}\n"

shape = np.empty(shape=data[options.namecol].shape, dtype="S8")
shape.fill("gaussian")

indices = np.where(data[options.acol]==0.)
shape[indices] = "point"

bigzip = zip(data[options.namecol], data[options.racol], data[options.decol], data[options.acol], data[options.bcol],
             data[options.pacol], data[options.fluxcol], data[options.alphacol], shape, data[options.fcol])

with open(options.output, "w") as f:
    f.write("skymodel fileformat 1.1\n")
    for Name, RA, Dec, a, b, pa, flux, alpha, shape, freq in bigzip:
        RA = replace(RA, ":", "h", 1)
        RA = replace(RA, ":", "m", 1)
        RA = RA + "s"
        Dec = replace(Dec, ":", "d", 1)
        Dec = replace(Dec, ":", "m", 1)
        Dec = Dec + "s"
        if shape == "gaussian":
            f.write(gformatter.format(Name=Name, RA=RA, Dec=Dec, a=a, b=b, pa=pa, flux=flux, alpha=alpha, beta=0.0,
                                      freq=freq, shape=shape))
        elif shape == "point":
            f.write(pformatter.format(Name=Name, RA=RA, Dec=Dec, flux=flux, alpha=alpha, beta=0.0, freq=freq,
                                      shape=shape))
