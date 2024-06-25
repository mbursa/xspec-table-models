# A helper class to make XSPEC table models
#
# It creates/opens a FITS file that contains an Xspec table model, 
# (a set of tabulated spectra that Xspec knows how to work with)
# and fills it with given spectra. 
#
# The table model FITS file is constructed according to the following spec:
# ftp://legacy.gsfc.nasa.gov/fits_info/fits_formats/docs/general/ogip_92_009
#
# If the FITS file has existed before, only missing spectra are computed. 
# In that way the script allows recovery from an interrupted run.
#
# See a provided description and examples that demonstrate the usage of the class.
#
# Dependency: astropy (https://www.astropy.org/)
# 
# (c) 2017 Michal Bursa
# Astronomical Institute of the Czech Academy of Sciences
# http://astro.cas.cz/bursa
#
# MIT Licence


import os
import sys
import math
import numpy as np

import astropy.io.fits as fits
import astropy.table as fits_table



# constants
kev2erg = 1.602177e-09  # keV-erg ((1e3*electronvolt)/erg)
erg2kev = 6.241507e+08  # erg->keV (erg/(1e3*electronvolt))


class XspecTableModelAdditive:

    # define the number of model spectra for each parameter grid
    # (in most situations there is exactly one spectrum per grid point)
    NXFLTEXP = 1

    def __init__(self, file_name, model_name, energies, params, redshift=False):
        """
        Create the target FITS file containing the tabulated model.
        
        Args:
            file_name: full path to the target file
            model_name: name of the model (should contain only letters and numbers)
            energies: grid of energies [keV]
            params: array of tuples (name, grid, logarithmic, frozen), where
                    `name` is the name of the parameter (max 11 characters),
s                    `grid` is an array of floats containing discrete values for 
                    the parameter, `logarithmic` is to tell Xspec whether to treat
                    the parameter linearly (0/False) or logarithmicly (1/True), 
                    and `frozen` tells if the parameter is frozen by default
            redshift: if the redshift parameter shall be added by Xspec
        """
        
        # remember some inputs
        self.model_name = model_name
        self.redshift   = redshift
        self.filename   = file_name
        self.energies   = energies
        self.params     = params
       
        # Xspec uses bin-integrated spectra, so the number of table bins is one less 
        # than the number of energy points
        en_bins = len(energies)-1

        # number of parameters
        n_params = len(params)
        
        # dimension of the parameter space (parameters x values)
        # (this gives the number of rows of the fits table)
        self.paramspace = 1
        for (param_name, param_grid, param_log, param_frozen) in params: self.paramspace *= len(param_grid)
        

        # create the primary HDU
        prihdr = fits.Header()
        self.define_hdu_primary(prihdr)
        self.primary_hdu = fits.PrimaryHDU(header=prihdr)


        # create PARAMETERS table
        params_hdr = fits.Header()
        params_hdr['NINTPARM'] = len(params)
        params_hdr['NADDPARM'] = 0
        params_hdr['HDUCLASS'] = "OGIP"
        params_hdr['HDUCLAS1'] = "XSPEC TABLE MODEL"
        params_hdr['HDUCLAS2'] = "PARAMETERS"
        params_hdr['HDUVERS1'] = "1.0.0"
        # determine the maximul length of param grid
        max_param_len = 0
        for (param_name, param_grid, param_log, param_frozen) in params: max_param_len = max(max_param_len, len(param_grid))
        # define table structure
        self.params_hdu = fits.BinTableHDU.from_columns([
            fits.Column(name='NAME', format='12A'),
            fits.Column(name='METHOD', format='J'),
            fits.Column(name='INITIAL', format='E'),
            fits.Column(name='DELTA', format='E'),
            fits.Column(name='MINIMUM', format='E'),
            fits.Column(name='BOTTOM', format='E'),
            fits.Column(name='TOP', format='E'),
            fits.Column(name='MAXIMUM', format='E'),
            fits.Column(name='NUMBVALS', format='J'),
            fits.Column(name='VALUE', format=str(max_param_len)+'E')
        ], params_hdr, nrows=len(params))
        self.params_hdu.name = 'PARAMETERS'
        # add rows with parameters
        param_index = 0
        for (param_name, param_grid, param_log, param_frozen) in params: 
            param_init  = 0.5*(min(param_grid)+max(param_grid))
            param_delta = (max(param_grid)-min(param_grid))/(5.*len(param_grid))
            self.params_hdu.data[param_index] = (
                param_name, 1 if param_log else 0, 
                param_init, 
                param_delta if (not param_frozen) else -1,
                min(param_grid), min(param_grid), 
                max(param_grid), max(param_grid), 
                len(param_grid), 
                np.pad(param_grid, (0,max_param_len-len(param_grid)), 'constant')
            )
            param_index += 1
        #end for


        # create ENERGIES table
        energs_hdr = fits.Header()
        energs_hdr['HDUCLASS'] = "OGIP"
        energs_hdr['HDUCLAS1'] = "XSPEC TABLE MODEL"
        energs_hdr['HDUCLAS2'] = "ENERGIES"
        energs_hdr['HDUVERS1'] = "1.0.0"
        self.energs_hdu = fits.BinTableHDU.from_columns([
            fits.Column(name='ENERG_LO', format='E'),
            fits.Column(name='ENERG_HI', format='E'),
        ], energs_hdr, nrows=en_bins)
        self.energs_hdu.name = 'ENERGIES'
        for i in range(en_bins):
            self.energs_hdu.data[i] = (energies[i], energies[i+1])
        #end of for


        # create SPECTRA table
        spectra_hdr = fits.Header()
        spectra_hdr['HDUCLASS'] = "OGIP"
        spectra_hdr['HDUCLAS1'] = "XSPEC TABLE MODEL"
        spectra_hdr['HDUCLAS2'] = "MODEL SPECTRA"
        spectra_hdr['HDUVERS1'] = "1.0.0"
        self.spectra_hdu = fits.BinTableHDU.from_columns([
            fits.Column(name='PARAMVAL', format=str(n_params)+'E'),
            fits.Column(name='INTPSPEC', format=str(en_bins)+'E'),
        ], spectra_hdr, nrows=self.NXFLTEXP*self.paramspace)
        self.spectra_hdu.name = 'SPECTRA'
    #end of def


    def define_hdu_primary(self, hdr):
        hdr['HDUCLASS'] = "OGIP"
        hdr['HDUCLAS1'] = "XSPEC TABLE MODEL"
        hdr['HDUVERS1'] = "1.0.0"
        hdr['MODLNAME'] = self.model_name
        hdr['MODLUNIT'] = "photons/cm2/s"
        hdr['REDSHIFT'] = bool(self.redshift)
        hdr['ADDMODEL'] = True
    #end def


    def generator(self):
        """
        Return a generator that iterates over the whole space of parameters.
        Rows that are already filed with data are skipped.

        Returns:
            Tuple of four items: the global row index, array of parameter indexes, array of parameter values, energy array
        """
        global_index = 0
        n_params = len(self.params)
        while (global_index < self.NXFLTEXP*self.paramspace):
            # skip row that have data already
            while (np.sum(self.spectra_hdu.data[global_index][1]) > 0.0): 
                global_index += self.NXFLTEXP
                if (global_index >= self.NXFLTEXP*self.paramspace): break
            if (global_index >= self.NXFLTEXP*self.paramspace): break

            # get indexes in each grid; the last grid changing the fastest
            param_indexes = np.zeros(n_params, dtype=int)
            param_values  = np.zeros(n_params)
            N0 = self.paramspace
            for i in range(n_params):
                (p_name, p_grid, p_log, p_frozen) = self.params[i]
                N = len(p_grid)
                N0 /= N
                p_index = int((global_index/self.NXFLTEXP)//N0 % N)
                #print('global_index',global_index)
                #print('p_index',p_index)
                #print('p_grid[p_index]',p_grid[p_index])
                #print('p_grid',p_grid)
                param_indexes[i] = p_index
                param_values[i]  = p_grid[p_index]

            # write parameter values (repeat the same parameters for each spectrum of the set)    
            for i in range(self.NXFLTEXP):
                self.spectra_hdu.data[global_index+i][0] = param_values
            #end for

            # return total index, array of grid indexes, and array of grid values
            #sys.stderr.write("> generator: passing spectrum index %d (%s %s)\n" % (global_index, str(param_indexes), str(param_values)))
            yield (global_index, param_values, param_indexes, self.energies)
            global_index += self.NXFLTEXP
        #end while
    #end def


    def write(self, index, Iv, flush=False):
        """
        Write a spectrum at a row position of `index`.
        
        Iv is a spectrum (specific flux) in units [erg/s/cm2/keV].
        If `flush` is set, the FITS file content will be saved to disk.
        """
        if (len(self.energies) != len(Iv)):
            sys.stderr.write("ERROR: invalid spectrum dimension (%d, expected %d)\n" % (len(Iv), len(self.energies)))
            return

        # convert spectrum from [erg/s/cm2/keV] to [photons/s/cm2]
        Iv_integrated = np.zeros(self.energies.size-1)
        for i in range(len(Iv_integrated)): 
            e1 = self.energies[i]*kev2erg
            e2 = self.energies[i+1]*kev2erg
            Iv_integrated[i] = 0.5*(Iv[i]/e1+Iv[i+1]/e2) * (e2-e1)*erg2kev
        #end for

        self.spectra_hdu.data[index][1] = Iv_integrated
        
        if flush: self.save()
    #end of def


    def save(self):
        """
        Write the FITS table to a file.
        """
        if os.path.isfile(self.filename): os.remove(self.filename)
        fits.HDUList([self.primary_hdu, self.energs_hdu, self.params_hdu, self.spectra_hdu]).writeto(self.filename)
    #end of def


#end of class (XspecTableModelAdditive)



class XspecTableModelAdditiveWithLinearPolarization(XspecTableModelAdditive):

    # define the number of model spectra for each parameter grid
    # for polarimetry we have 3 spectra for each grid point that are fitted simoultaneously (I, Q, U)
    NXFLTEXP = 3

    def define_hdu_primary(self, hdr):
        XspecTableModelAdditive.define_hdu_primary(self,hdr)

        # add polarization-specific headers
        hdr['NXFLTEXP'] = self.NXFLTEXP
        hdr['XFXP0001'] = 'Stokes:0'
        hdr['XFXP0002'] = 'Stokes:1'
        hdr['XFXP0003'] = 'Stokes:2'
    #end def


    def write(self, index, Iv, Qv, Uv, flush=False):
        """
        Write a spectrum at a row position of `index`.
        
        Spectrum is supposed to be in units [erg/s/cm2/keV].
        If `flush` is set, the FITS file content will be saved to disk.
        """
        if (len(self.energies) != len(Iv)) or (len(self.energies) != len(Qv)) or (len(self.energies) != len(Uv)):
            sys.stderr.write("ERROR: invalid spectrum dimension (%d, expected %d)\n" % (len(Iv), len(Qv), len(Uv), len(self.energies)))
            return

        # convert spectrum from [erg/s/cm2/keV] to [photons/s/cm2]
        Iv_integrated = np.zeros(self.energies.size-1)
        Qv_integrated = np.zeros(self.energies.size-1)
        Uv_integrated = np.zeros(self.energies.size-1)
        for i in range(len(Iv_integrated)): 
            e1 = self.energies[i]*kev2erg
            e2 = self.energies[i+1]*kev2erg
            Iv_integrated[i] = 0.5*(Iv[i]/e1+Iv[i+1]/e2) * (e2-e1)*erg2kev
            Qv_integrated[i] = 0.5*(Qv[i]/e1+Qv[i+1]/e2) * (e2-e1)*erg2kev
            Uv_integrated[i] = 0.5*(Uv[i]/e1+Uv[i+1]/e2) * (e2-e1)*erg2kev
        #end for

        self.spectra_hdu.data[index+0][1] = Iv_integrated
        self.spectra_hdu.data[index+1][1] = Qv_integrated
        self.spectra_hdu.data[index+2][1] = Uv_integrated
        
        if flush: self.save()
    #end of def


#end of class





# if executed as main file then
# run the main function and exit with returned error code
if __name__ == "__main__": 
    print("Provides a classes to help to create a table model for XSPEC. See the documentation for a usage guide.")
    sys.exit(0)
#end if

