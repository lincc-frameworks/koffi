import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS


def create_fake_fits_file(fname, x_dim, y_dim):
    # Create a primary HDU with just the date/time
    # and observatory info in the header.
    hdr0 = fits.Header()
    hdr0["DATE-AVG"] = "2022-08-15T06:00:00.000000000"
    hdr0["OBS-LAT"] = -30.166060
    hdr0["OBS-LONG"] = 70.814890
    hdr0["OBS-ELEV"] = 2215.000000
    hdr0["NAXIS"] = 2
    hdr0["NAXIS1"] = x_dim
    hdr0["NAXIS2"] = y_dim
    hdu0 = fits.PrimaryHDU(header=hdr0)

    # Create and image HDU with a header containing
    # minimal celestial information for the WCS.
    data1 = np.ones((y_dim, x_dim))
    hdr1 = fits.Header()
    hdr1["WCSAXES"] = 2

    # (0,0) corner is at RA=201.614 and Dec=-10.788
    # with 0.001 degrees per pixel.
    hdr1["CRPIX1"] = 1.0
    hdr1["CRVAL1"] = 201.614
    hdr1["CDELT1"] = 0.001
    hdr1["CTYPE1"] = "RA"

    hdr1["CRPIX2"] = 1.0
    hdr1["CRVAL2"] = -10.788
    hdr1["CDELT2"] = 0.001
    hdr1["CTYPE2"] = "DEC"
    hdu1 = fits.ImageHDU(data1, header=hdr1)

    # Write both HDUs to the given file.
    h = fits.HDUList([hdu0, hdu1])
    h.writeto(fname, overwrite=True)
