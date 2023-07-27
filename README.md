


# XSPEC Table Model Generator

Create XSPEC table models with the help of this Python class.

## Purpose

Should you want to create your own spectral model for [XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) (an X-ray Spectral Fitting Package), the easiest way is to make a FITS table in [OGIP](https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/summary/ogip_92_009_summary.html) standard as it is described in the XSPEC manual "Appendix C: Adding Models to XSPEC".  The basic concept of a table model is that the file contains an N-dimensional grid of model spectra with each point on the grid having a spectrum for a particular combination of values of the N parameters in the model.

This python module provides a class that helps to create the FITS file in the right format, so you only need to focus on getting your model spectra.

## Limitations

* `XspecTableModelAdditive` class does not support so called "additional parameters" at the moment
* only additive models can be created at the moment

## Dependencies

The module uses [AstroPy](https://www.astropy.org/) for writing the FITS files.

## Usage

The basic usage involves importing the helper class from the module, creating an instace of `XspecTableModel` and filling the model with spectra. The skeleton looks like this:

```python
from xspec_table_models import XspecTableModelAdditive

def spectrum(energies, param1, param2):
  ''' This function provides the spectrum for the given set of parameters. '''
	return [] # specific fluxes [erg/s/cm2/keV] for the grid of energies
#end def

# array of energies [keV]
energies = [...]
# define parameters:
# each parameter is a tuple of: parameter name, array of values, 
#    False/True for linear/logarithimic interpolation, False/True for free/frozen parameter
param1 = ('param1', [...], False, False)
param2 = ...

fits = XspecTableModelAdditive('mymodel.fits', 'mymodel', energies, [param1, param2, ...])

for g in fits.generator():
    # recover parameters for the current grid point
    index, param_values, param_indexes, energies = g
    param1 = param_values[0]
    param2 = param_values[1]
    # get the spectrum
    Iv = spectrum(energies, param1, param2)
    # write to FITS file
    fits.write(index, Iv, False)
#end if

fits.save();
```

Then within the XSPEC environment, you simply load your table model with
```
XSPEC> model atable{mymodel.fits}
```

#### Polarization

In case you work with spectro-polarimetric data (eg. X-ray polarization), you can use `XspecTableModelAdditiveWithLinearPolarization` class that is a descendant of `XspecTableModelAdditive` and produces FITS file with three spectra (for `I`,`Q`,`U` Stokes parameters) for each parameter grid point. Only linear polarization is assumed, so only spectra for `I`,`Q`,`U` are given and `V` is ignored.

## Examples

You can refer to the [examples](tree/main/examples) folder for complete and commented examples.

## Documentation

### XspecTableModelAdditive class

Use `XspecTableModelAdditive` class to create an additive table model.

```
class XspecTableModelAdditive(file_name, model_name, energies, params, redshift=False)
```

Creates the class instance and opens the FITS file from the filesystem if it exists already or creates a new one.

**file_name**  
The path to the FITS file that will be created in the filesystem.  
**model_name**  
The name of the model (letters only, 12 characters max).  
**energies**  
Array of energies to consider for spectra in [keV].  
**params**  
Array of parameters of the model. Each parameter is a tupe with 4 items:
* **name** - name of the parameter (letters only, 12 characters max)
* **grid** - array of parameter values (in increasing order).
* **logarithmic** - True/False flag saying if Xspec shall interpolate in the parameter values linearly (False) or logarithmicaly (True)
* **frozen** - True/False flag saying if Xspec shall initially set this parameter as frozen 
Example: `par1 = ('mass', [10,20,30,40,50], False, False)`
**redshift**
If `redshift` parameter shall be added  by Xspec to the model (boolean). The `redshift` parameter will shift the model in energy space and divide by (1+z) factor.

  Do not include a normalization parameter, it will be added by Xspec automaticaly.

Example:
```python
energies = np.geomspace(1e-2, 1e+2, 100)
param1 = ('alpha', np.linspace(-5,+5,20), False, False)
param1 = ('beta', np.geomspace(1e-1,1e+1, 50), True, False)
model = XspecTableModelAdditive('mymodel.fits', 'mymodel', energies, [param1,param2], False)
```
<br>

```
def XspecTableModelAdditive.generator()
```
Gives an iterator that loops over all combinations of parameter values and allows to provide a spectrum for each row of the spectral table. The iterator returns a tuple with 4 items:
* **index** - index of the row in the spectral table (shall be passed to `write()`)
* **param_values** - array of parameter values for the current spectral row (values are in the same order in which parameters have been passed to the class constructor)
* **param_indexes** - array of parameter indexes for the current spectral row (not really needed, but provided for completeness; you can use the index to get the parameter value from the parameter value grid)
* **energies** - array of energies in [keV] (a copy of the energy grid passed to the class constructor)

The iterator skips any rows that have the spectra filled already. In that way, if the FITS file have existed before, only the missing spectra are computed and so the script allows for a recovery from an interrupted run. On the other hand, if you want to start over, you need to remove the existing FITS file before starting the script.

Example:
```python
model = XspecTableModelAdditive(...)
for g in fits.generator():
    index, param_values, param_indexes, energies = g
    Iv = spectrum(energies, param_values)
    model.write(index, Iv, False)
#end if
```
<br>

```
def XspecTableModelAdditive.write(index, Iv, flush=False)
```
Write a single spectrum to the table. 
**index**  
Row index of the spectrum (given by the generator).  
**Iv**  
Energy spectrum (specific flux) in [erg/s/cm2/keV] given at each point of the energy grid.  
**flush**  
If True, the model table is saved to the file system after the spectrum is written.

**Note**: XSPEC requires the table model to contain spectra in units of photons/cm/s (photon spectrum integrated over the energy bin). The spectrum that is passed to `write()` method, however, must be an energy spectrum (specific flux) given at each energy point (not integrated). The integration and conversion to photon spectrum is done inside the function.
<br>

```
def XspecTableModelAdditive.save()
```
Save the content of the FITS to the filesystem.



### XspecTableModelAdditiveWithLinearPolarization class

```
def XspecTableModelAdditiveWithPolarization.write(index, Iv, Qv, Uv, flush=False)
```
Write a single spectrum to the table. 
**index**  
Row index of the spectrum (given by the generator).  
**Iv**  
Energy spectrum (specific flux) in [erg/s/cm2/keV] given at each point of the energy grid.  
**Qv**  
Energy spectrum (specific flux) in [erg/s/cm2/keV] of the Q Stokes parameter given at each point of the energy grid.  
**Uv**  
Energy spectrum (specific flux) in [erg/s/cm2/keV] of the U Stokes parameter given at each point of the energy grid.  
**flush**  
If True, the model table is saved to the file system after the spectrum is written.

**Note**: XSPEC requires the table model to contain spectra in units of photons/cm/s (photon spectrum integrated over the energy bin). The spectrum that is passed to `write()` method, however, must be an energy spectrum (specific flux) given at each energy point (not integrated). The integration and conversion to photon spectrum is done inside the function.
<br>


## History

* 2021/07/22 - intial release
* 2023/07/26 - added support for spectro-polarimetric data (linear polarization)