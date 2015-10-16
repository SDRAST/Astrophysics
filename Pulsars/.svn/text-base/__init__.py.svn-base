# -*- coding: utf-8 -*-
"""
package Astrophysics.Pulsars - analyze pulsar data

Contents
========

1. Effects related to propagation of pulsed signals (main module).
2. Pulsar data from the Parkes Catalogue (module pulsar_data)
3. Fast folding algorithm for DM searching using numpy (module np_ffa).
"""

import Astronomy
import Physics
import math
import pickle
import pulsar_data
import np_ffa
from math import log10

def propagation_delay(electron_density, path_length, frequency):
    """
    EM radiation propagation time

    Time-of-flight is calculated from physical parameters

    @type electron_density : float
    @param electron_density : in cm^{-3}
    
    @type path_length : float
    @param path_length : in pc

    @type frequency : float
    @param frequency : in MHz

    @return: time in seconds
    """
    L = path_length*Astronomy.pc
    nu = 2*math.pi*frequency*1e6 # in radians/sec
    nu_p = Physics.plasma_frequency(electron_density) # in radians
    prop_delay = L*pow(nu_p/nu,2)/Physics.c/2
    # print "Propagation delay = %15.10e" % (prop_delay)
    return prop_delay

def drift_rate( f, electron_density, path_length ):
    """Given the frequency f in Hz, the mean electron density, and the
       path length, returns the drift rate in Hz/s"""
    f_p = Physics.plasma_frequency(electron_density)/(2*math.pi)
    return -Physics.c*pow(f,3)/f_p**2/path_length/Astronomy.pc

def DfDt(f_GHz,DM):
    """Given the frequency in GHz, and the dispersion measure in pc/cm^3,
    returns the slope of f vs t in MHz/microsec"""
    return drift_rate(f_GHz*1e9,1,DM)/1e6/1e6

def emission_measure (electron_density, path_length):
    """ electron density in electrons/cc """
    """ path length in pc """
    """ returns emission measure in cm^-2 """
    return electron_density*path_length*Astronomy.pc

def differential_time_delay(DM, f, BW):
    """ Given the dispersion measure in cm^-3/pc, f in GHz, BW in MHz
        returns time_delay in microseconds"""
    f   = f*1.e9   # Hz
    BW  = BW*1.e6  # Hz
    factor = 1./(4.*math.pow(math.pi,2))
    factor = factor * math.pow(Physics.e,2) \
                    *DM/(Physics.c*Physics.m_e*Physics.eps_0)
    factor = 1e6*factor*Astronomy.pc
    return 1e6*factor*BW/pow(f,3)

def time_delay(DM,f):
    """
    EM transit time from dispersion measure
    
    @type DM : float
    @param DM : dispersion measure in cm^-3 x pc

    @type f : float
    @param f : frequency in GHz

    @return: the time delay relative in microseconds
    """
    # The factor of 1e6 converts the electron density in the dispersion
    # measure from cm^-3 to m^-3
    factor = 1e6*DM*pow(Physics.e,2)*Astronomy.pc/\
             (8*pow(math.pi,2)*Physics.c*Physics.m_e*Physics.eps_0)
    # This converts seconds to microseconds
    f = f*1e9
    return factor*DM/math.pow(f,2)

def broadening(f,dm,model="B",dalpha=0):
  """
  Pulse broadeing due to interstellar scattering

  The model is taken from Bhat et al., Ap.J. 605, 759 (2004), eqn. 7.
  The coefﬁcients, a = -6.46, b = 0.154, and c = 1.07, are only slightly
  different from those of Cordes & Lazio (2003), a = -6.59, b = 0.129, and
  c = 1.02.  The frequency exponent alpha was determined by Bhat et al.
  to be 3.86 +/- 0.16.

  @param f : frequency (GHz)
  @type  f : float

  @param dm : dispersion measure (cm^3 / pc)
  @type  dm : float

  @return: the broadening timescale in millisec
  """
  if model[0].upper() == 'B':
    a = -6.46
    b = 0.154
    c = 1.07
    alpha = 3.86
  elif model[0].upper() == "C":
    a = -6.59
    b = 0.129
    c = 1.02
    alpha = 4.4
  else:
    raise "Unknown scattering model",model

  log_tau = a + b*log10(dm) + c*(log10(dm)**2) - (alpha+dalpha)*log10(f)
  return 10**(log_tau)
