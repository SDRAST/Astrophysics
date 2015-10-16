"""Astrophysics functions and modules"""

def Tb_galaxy_NCP(f):
    """Returns the sky brightness in K due to galactic synchrotron radiation
    at the north celestial pole given the frequency f in MHz"""
    return 280*pow(f/150.,-2.5)
