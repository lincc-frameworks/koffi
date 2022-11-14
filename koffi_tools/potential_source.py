from astropy.time import Time
from astropy.io import fits
from astropy.wcs import WCS

class PotentialSource:
    def __init__(self, form = "pixel", x=None, y=None, vx=None, vy=None,
                ra=None, dec=None, vra=None, vdec=None):

        if form != "pixel" and form != "radec":
            raise ValueError("value 'format' must be either 'pixel' or 'radec'")

        self.format = form
            
        # pixel format
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy

        # ra/dec format
        self.ra = ra
        self.dec = dec
        self.vra = vra
        self.vdec = vdec