# Example of making an additive table model with polarization
#
# Makes an additive table model with a power-law energy spectrum having an exponential cut-off and 
# 100% linear polarization. The result should be idential to the analytical 'cutoffpl' model of Xspec, 
# which you may verify yourself loading this table model to XSPEC:
# XSPEC> model atable{cutoffpl}
#
# This example uses an analytical function for spectra, but in a similar way you can build a model 
# using external data that you read from text or binary files (e.g. a result of a simulation).
#
# Each grid point, i.e. each distinct combination of model parameters `alpha`` and `beta`, has three spectra 
# that are fitted simultaneously in XSPEC - one for each of the Stokes parameters `I`, `Q`, `U` 
# (`V` is not included as only linear polarization is assumed)

import sys
import math
import numpy as np

# add the path to xspec_table_model.py and import the class from the module
sys.path.append('..')
from xspec_table_models import XspecTableModelAdditiveWithLinearPolarization


# Define a function that will serve spectra for given parameters of your model.
def spectrum(energies, params):
    """ 
    Return energy (not photon) spectrum in [erg/s/cm2/keV]. 
    energies: the energy earray [keV]
    params: array of parameters
    Note: in this function you implement a method of getting your spectra. You can provide a formula,
    you can read and process external numerical data, you can call external service or application, etc.
    """

    # In our table model, we have defined two parameters (alpha, beta; in this order), so params 
    # is an array with two items, first is current value of alpha parameter taken from the list of values, 
    # second is beta. 
    # We return a power-law spectrum with a cut-off, where alpha is the photon index and beta is 
    # the cutoff energy. Normalization is such that the model gives exactly 1 photon at 1keV.
    # 6.241507e+08 is erg-to-keV conversion factor (erg/(1e3*electronvolt))
    param_alpha = params[0]
    param_beta = params[1]
    # note: as the result has to be an energy spectrum (not photon spectrum), the photon index alpha 
    # shall be increased by 1 to have spectral index
    norm = 6.241507e+08 * 1.0**(-param_alpha+1) * math.exp(-1.0/param_beta)
    
    # let's define Stokes spectra for a 100% vertically polarized light
    Iv = 1./norm * np.power(energies, -param_alpha+1) * np.exp(-energies/param_beta)
    Qv = Iv
    Uv = 0.0 * energies
    
    return Iv, Qv, Uv
#end def





# define the name of the model (should only contain letters and _) and 
# the name of the FITS file the model fill be saved to
model_name = 'cutoffplpol'
model_fits_file = 'cutoffplpol.fits'

# define the energy grid for your specra in [keV]
en_min = 1e-2
en_max = 1e+2
en_bins = 500
energies = np.geomspace(en_min, en_max, en_bins)

# define the parameters of your model
# each parameter is a tuple with the following 4 items: 
#  - parameter name
#  - list of parameter values 
#  - interpolation method: 0/False for linear, 1/True for logarithmic
#  - flag saying if the parameter shall be initially frozen in Xspec (True/False)
# In this example, we have two parameters:
#  - alpha is the spectral index of the power law
#  - beta is the cut-off energy [keV]
param1 = ('alpha', np.linspace(-5.0, 5.0, 50), False, False)
param2 = ('beta',  np.geomspace(0.01, 500.0, 50), True, False)

# set up the fits file
fits = XspecTableModelAdditiveWithLinearPolarization(model_fits_file, model_name, energies, [param1, param2])

# fill the model with spectra
# the loop body is called once for each combination of model parameters; in our case it will be
# called 50x50=2500 times
for g in fits.generator():
    # g tuple contains 4 items:
    #  - index: the row index of the spectral table (starts at 0)
    #  - param_values: a tuple containing values of the current parameter combination
    #  - param_indexes: a tuple containing indices of the current parameter combination
    #  - energies: the energy grid [keV]
    # the paremeters are in the same order as they have been given to XspecTableModelAdditive(); 
    # the last parameter changes the fastest
    index, param_values, param_indexes, energies = g
    sys.stderr.write("\r")
    sys.stderr.write("> spectrum for row index %d, params:%s" % (index, str(param_values)))
    sys.stderr.write(" "*20)

    # get the I, Q, U spectrum for the given combination of parameters;
    # the returned spectrum units must be [erg/s/cm2/keV]
    Iv, Qv, Uv = spectrum(energies, param_values)
    
    # write the spectrum to the model table
    fits.write(index, Iv, Qv, Uv, False)
#end if
sys.stderr.write("\n")

# finally, save the fits file
sys.stderr.write("> Saving FITS file %s\n" % (model_fits_file))
fits.save();

sys.stderr.write("> Done\n")


