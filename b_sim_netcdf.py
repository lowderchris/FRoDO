import struct
import numpy as np
from scipy import ndimage
import sphb
#import grid
from scipy.io import netcdf
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import matplotlib.cm as cm


class SphB_sim(sphb.SphB):
    """Sub-class for 3d magnetic field from MF simulation
    -----
    This version reads in the new netcdf output format (then maps it to uniform
    grid).
    -----
    modified 15/9/14 --ary
    """

    def __init__(self, fileName, nr, nth, nph):
        """Class constructor.
        The size of the uniform grid is given by nr, nth, nph.
        """

        # Read in magnetic field from netcdf file:
        self.bv, self.ri, self.thi, self.phi = self.readB(fileName)

        # Create uniform spherical grid:
        dph = 2*np.pi/nph
        dth = np.pi/nth
        dr = (self.ri[-1] - self.ri[0])/nr
        ph = np.linspace(0, 2*np.pi, nph)
        #th = np.linspace(self.thi[-1], self.thi[0], nth)
        th = np.linspace(0, np.pi, nth)
        r = np.linspace(self.ri[0], self.ri[-1], nr) 
        r3, th3, ph3 = np.meshgrid(r,th,ph)
        r3 = np.swapaxes(r3,0,1)
        th3 = np.swapaxes(th3,0,1)
        ph3 = np.swapaxes(ph3,0,1)

        self.r = r
        self.th = th
        self.ph = ph
        self.nr = nr
        self.nth = nth
        self.nph = nph
        self.r0 = r[0]
        self.r1 = r[-1]
        self.th0 = th[0]
        self.th1 = th[-1]
        self.ph0 = ph[0]
        self.ph1 = ph[-1]

        # Regrid to uniform spherical grid:
        self._variableToUnif_(self.bv)

        # Compute magnitude, needed for field line tracing:
        self.bbg = np.sqrt(self.brg**2 + self.bthg**2 + self.bphg**2)
        self.bbg[self.bbg < 1e-12] = 1e-12
        self.jjg = np.sqrt(self.jrg**2 + self.jthg**2 + self.jphg**2)     

    def readB(self, fileName):
        """Read magnetic field and grid from file "bv[snap].dat"
        ---
        modified 15/9/14 --ary
        """

        fh = netcdf.netcdf_file(fileName, 'r', mmap=False)
        r = fh.variables['r'][:]
        th = fh.variables['th'][:]
        ph = fh.variables['ph'][:]
        br = fh.variables['br'][:]
        bth = fh.variables['bth'][:]
        bph = fh.variables['bph'][:]
        jr = fh.variables['jr'][:]
        jth = fh.variables['jth'][:]
        jph = fh.variables['jph'][:]
        ar = fh.variables['ar'][:]
        ath = fh.variables['ath'][:]
        aph = fh.variables['aph'][:]
        fh.close()
        
        br = np.swapaxes(br,0,2)
        bth = np.swapaxes(bth,0,2)
        bph = np.swapaxes(bph,0,2)
        jr = np.swapaxes(jr,0,2)
        jth = np.swapaxes(jth,0,2)
        jph = np.swapaxes(jph,0,2)
        ar = np.swapaxes(ar,0,2)
        ath = np.swapaxes(ath,0,2)
        aph = np.swapaxes(aph,0,2)

        bv = SubB(br, bth, bph, jr, jth, jph, ar, ath, aph)
        return bv, r, th, ph


    def _variableToUnif_(self, bv):
        """Interpolate magnetic field from netcdf grid to uniform spherical grid.
        Uses scipy.ndimage.map_coordinates function for speed.
        -----
        modified 15/9/14 --ary
        """
        interpOrder = 3
        
        # Interpolation indices in each direction:
        fr = interp1d(self.ri,np.arange(np.size(self.ri)))
        fth = interp1d(self.thi[::-1],(np.arange(np.size(self.thi)))[::-1], \
                       bounds_error=False, fill_value=len(self.thi)+10)
        fph = interp1d(self.phi,np.arange(np.size(self.phi)))
       
        jr0 = fr(self.r)
        jth0 = fth(self.th)
        jph0 = fph(self.ph)
      
        jr, jth, jph = np.meshgrid(jr0,jth0,jph0)
        jr = np.swapaxes(jr,0,1)
        jth = np.swapaxes(jth,0,1)
        jph = np.swapaxes(jph,0,1)
        self.bphg = ndimage.map_coordinates(bv.bph,[jr, jth, jph], order=interpOrder)
        self.bthg = ndimage.map_coordinates(bv.bth,[jr, jth, jph], order=interpOrder)
        self.brg = ndimage.map_coordinates(bv.br,[jr, jth, jph], order=interpOrder)
        self.jphg = ndimage.map_coordinates(bv.jph,[jr, jth, jph], order=interpOrder)
        self.jthg = ndimage.map_coordinates(bv.jth,[jr, jth, jph], order=interpOrder)
        self.jrg = ndimage.map_coordinates(bv.jr,[jr, jth, jph], order=interpOrder)
        self.aphg = ndimage.map_coordinates(bv.aph,[jr, jth, jph], order=interpOrder)
        self.athg = ndimage.map_coordinates(bv.ath,[jr, jth, jph], order=interpOrder)
        self.arg = ndimage.map_coordinates(bv.ar,[jr, jth, jph], order=interpOrder)

class SubB:
     def __init__(self, br, bth, bph, jr, jth, jph, ar, ath, aph):
         """Structure to store magnetic field vector.
         ---
         modified 15/9/14 --ary
         """
         self.br = br
         self.bth = bth
         self.bph = bph
         self.jr = jr
         self.jth = jth
         self.jph = jph
         self.ar = ar
         self.ath = ath
         self.aph = aph
