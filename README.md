
# Xspec Table Models Generator

Python class that helps to create additive table models for XSPEC.

## Purpose

Should you want to create your own spectral model for XSPEC, the easiest way is to make a FITS table in [OGIP](https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/summary/ogip_92_009_summary.html) standard as it is described in the XSPEC manual "Appendix C: Adding Models to XSPEC".  The basic concept of a table model is that the file contains an N-dimensional grid of model spectra with each point on the grid having a spectrum for a particular combination of values of the N parameters in the model.

This python module provides a class that helps to create the FITS file in the right format, so you only need to focus on getting your model spectra.

## Limitations

* `XspecTableModelAdditive` class does not support so called "additional parameters" at the moment
* only additive models can be created at the moment

## Usage

The basic usage involves importing the helper class from the module, creating an instace of `XspecTableModel` and filling the model with spectra. The skeleton looks like this:

```python
from xspec_table_models import XspecTableModel

def spectrum(energies, params):
	return [] # specific fluxes [erg/s/cm2/keV] for the grid of energies
#end def

fits = XspecTableModelAdditive('mymodel.fits', 'mymodel', energies, [param1, param2, ...])

for g in fits.generator():
    index, param_indexes, param_values, energies = g
    Iv = spectrum(energies, param_values)
    fits.write(index, Iv, False)
#end if

fits.save();
```

Then within the XSPEC environment, you simply load your table model with
```
XSPEC> model atable{mymodel.fits}
```

You can refer to the [examples](tree/main/examples) folder for complete and commented examples.

<!--
## Documentation

The  **atable**  model requires the table model to contain spectra in units of photons/cm/s (the standard for additive models models). XSPEC will include a normalization parameter. 
-->
