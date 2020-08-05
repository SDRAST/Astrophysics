# -*- coding: utf-8 -*-
"""
Astrophysics.Pulsars.pulsar_data - data from the Parkes pulsar catalog

This module creates a dictionary 'data' of pulsar data keyed to the pulsar's
JNAME. Each pulsar's data is a dictionary keyed to the keys used in the
psrcat.db from the Parkes Pulsar Survey website,
        http://www.atnf.csiro.au/research/pulsar/.

An example of getting specific data
>>> import Astrophysics.Pulsars.pulsar_data as PD
>>> data = PD.data
>>> PD.names['B0531+21']
'J0534+2200'
>>> data['J0534+2200']['RAJ'], data['J0534+2200']['DECJ']
('05:34:31.973', '+22:00:52.06')
>>> data['J0538+2817']['RAJ'], data['J0538+2817']['DECJ'], data['J0538+2817']['S1400']
>>> ('05:38:25.0572', '+28:17:09.161', '1.9')

To get a sorted list of data keys, do
>>> K = list(help.keys())
>>> K.sort()
>>> print K
['A1', 'ASSOC', 'BINARY', 'CLK', 'DECJ', 'DIST_AMN', 'DIST_AMX', 'DIST_DM',
 'DIST_DM1', 'DM', 'DM1', 'ECC', 'EPHEM', 'EPS1', 'EPS2', 'F0', 'F1', 'F2',
 'FINISH', 'OM', 'OMDOT', 'P0', 'P1', 'PB', 'PEOCH', 'PMA', 'PMDEC',
 'POSEPOCH', 'PSRB', 'PSRJ', 'RAJ', 'S1400', 'S400', 'S600', 'SPINDX',
 'START', 'SURVEY', 'T0', 'TASC', 'TRES', 'TZRFRQ', 'TZRMJD', 'TZRSITE',
 'W10', 'W50']

To get a (very long!) list of pulsar names, do
>>> P = data.keys()
>>> P.sort()
>>> P[0]
'J0006+1834'
>>> P[1]
'J0014+4746'
>>> len(P)
1627
 """

import pickle
import Astronomy as A
import ephem as E
import sys
import math
import logging

module_logger = logging.getLogger(__name__)

pyver=sys.version[0:3]
dbfile=open(
  '/usr/local/RATools/Astrophysics/Pulsars/pulsars-pickle','rb')
data=pickle.load(dbfile)
dbfile.close()

namefile=open(
  '/usr/local/RATools/Astrophysics/Pulsars/BtoJnames','rb')
names=pickle.load(namefile)
namefile.close()

"""This is has incomplete list of the dictionary keys"""
help={'PSRB'  :'B1950 name',\
    'PSRJ'    :'J2000 name',\
    'RAJ'     :'J2000 right ascension',\
    'DECJ'    :'J2000 declination',\
    'ELONG'   :'Ecliptic longitude (degrees)',\
    'ELAT'    :'Ecliptic latitude (degrees)',\
    'POSEPOCH':'Epoch at which the position is measured (MJD)',\
    'PMA'     :'Proper motion in right ascension (mas/yr)',\
    'PMDEC'   :'Proper motion in declination (mas/yr)',\
    'F0':'Barycentric rotation frequency (Hz)',\
    'F1':'Time derivative of barycentric rotation frequency (s^-2)',\
    'F2':'Second derivative of barycentric rotation frequency (s^-3)',\
    'P0':'Barycentric period of the pulsar (s)',\
    'P1':'Time derivative of barcycentric period (dimensionless)',\
    'PEPOCH':'Epoch of period or frequency (MJD)',\
    'DM':'Dispersion measure (cm^-3 pc)',\
    'DM1':'First time derivative of dispersion measure (cm^-3 pc yr^-1)',\
    'BINARY':'Binary model (normally one of several recognised by the\n\
    pulsar timing program TEMPO ',\
    'TASC':'Epoch of ascending node(MJD)',\
    'T0':'Epoch of periastron (MJD)',\
    'PB':'Binary period of pulsar (days)',
    'A1':'Projected semi-major axis of orbit (lt s)',\
    'ECC':'Eccentricity',\
    'OM':'Longitude of periastron (degrees)',\
    'OMDOT':'',\
    'EPS1':'Ecc * sin(OM) - ELL1 binary model',\
    'EPS2':'Ecc * cos(OM) - ELL1 binary model',\
    'START':'',\
    'FINISH':'',\
    'TRES':'',\
    'CLK':'',\
    'EPHEM':'',\
    'TZRMJD':'',\
    'TZRFRQ':'',\
    'TZRSITE':'',\
    'S400':'Mean flux density at 400 MHz (mJy)',\
    'S600':'Mean flux density at 600 MHz (mJy)',\
    'S1400':'Mean flux density at 1400 MHz (mJy)',\
    'SPINDX':'Measured spectral index',\
    'W50':'Width of pulse at 50% of peak (ms).  Note, pulse widths are a\n\
    function of both observing frequency and observational time\n\
    resolution,so quoted widths are indicative only. Refer to the\n\
    original reference for details. ',\
    'W10':'Width of pulse at 10% (ms). Note the comments above for W50.',\
    'DIST_DM':'Distance based on the Taylor & Cordes (1993) electron\n\
    density model. In LONG or PUBLICATION QUALITY modes,\n\
    lower limits from the distance model are preceded by a\n\
    "+" sign.',\
    'DIST_DM1':'',\
    'DIST_AMN':'',\
    'DIST_AMX':'',\
    'ASSOC':'Names of other objects, for example, supernova remnants or\n\
    globular clusters, associated with the pulsar',\
    'SURVEY':'Surveys that detected the pulsar (discovery survey first).'}
    
def key_help(key):
    return help[key]

def get(data,key):
    """
    Returns the value corresponding to a key, if it exists.
    
    Otherwise it returns an empty string.
    """
    try:
        value = data[key]
    except KeyError:
        value = ''
    return value

def equatorial(data):
    """
    Returns the J2000 right ascension and declination
    
    Taken from keyed data or computed from from ecliptic or galactic 
    coordinates, whatever is given.
    
    @return: hours (float), degrees (float)
    """
    ra = get(data,'RAJ')
    if ra != '':
        decl = get(data,'DECJ')
        ra   = E.hours(data['RAJ']) # E.Angle(data['RAJ'], A.u.hourangle)
        decl = E.degrees (data['DECJ']) # E.Angle(data['DECJ'], A.u.deg)
        return ra, decl
    else:
        elong = get(data,'ELONG')
        module_logger.debug("equatorial: elong is %s", elong)
        if elong != '':
            elat = float(get(data,'ELAT'))
            try:
                mjd = float(get(data,'POSEPOCH'))
            except (IndexError, ValueError):
                try:
                    mjd = float(get(data,'PEPOCH'))
                except (IndexError, ValueError):
                    # default to J2000
                    mjd = 51544
            ra,decl = A.ecliptic_to_J2000(elong,elat,mjd)
            return ra*12/math.pi, decl*180/math.pi
        else:
            # try galactic coordinates
            return "","" # for now
    return None
        
def period(data):
    """
    Returns the pulsar period in millisec
    """
    # There are two ways to record the pulse period (in ms)
    period = get(data,'P0')
    if period == '':
        freq = get(data,'F0')
        if freq:
          period = 1000./float(freq)
        else:
          pass
    else:
        period = 1000*float(period)
    return period

def period_change_rate(data):
    """Returns the time rate of change of the period in s/s"""
    rate = get(data,'P1')
    if rate == '':
        f0 = get(data,'F0')
        f1 = get(data,'F1')
        if f1 == '':
            rate = 0
        else:
            rate = -math.pow(float(f0),-2)*float(f1)
    else:
        rate = float(rate)
    return rate
