import tempfile
import unittest

import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

from koffi import ImageMetadata, ImageMetadataStack

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


class TestImageMetadata(unittest.TestCase):
    def setUp(self):
        self.metadata = ImageMetadata()

    def test_init_file_none(self):
        self.assertEqual(self.metadata.hdu_list, None)
        self.assertEqual(self.metadata.center, None)
        self.assertFalse(self.metadata.obs_loc_set)
        self.assertFalse(self.metadata.epoch_set)

    def test_item_access(self):
        with tempfile.TemporaryDirectory() as dir_name:
            # Create two fake files in the temporary directory.
            fname = "%s/tmp.fits" % dir_name
            create_fake_fits_file(fname, 20, 30)
            self.metadata.populate_from_fits_file(fname)
            self.assertEqual(self.metadata["NAXIS1"], 20)
            self.assertEqual(self.metadata["FAKE_COEFF"], None)



    def test_load_from_file(self):
        with tempfile.TemporaryDirectory() as dir_name:
            # Create two fake files in the temporary directory.
            fname = "%s/tmp.fits" % dir_name
            create_fake_fits_file(fname, 20, 30)
            self.metadata.populate_from_fits_file(fname)
            self.assertEqual(self.metadata["NAXIS1"], 20)
            self.assertEqual(self.metadata["NAXIS2"], 30)
            self.assertEqual(self.metadata["CTYPE1"], "RA")
            self.assertEqual(self.metadata["CTYPE2"], "DEC")

    def test_epoch(self):
        mjd_time = 57130.25769920395
        epoch = Time(mjd_time, format='mjd')

        self.metadata.set_epoch(epoch)

        self.assertTrue(self.metadata.epoch_set)
        self.assertEqual(self.metadata.get_epoch().mjd, mjd_time)

    def test_set_obs_position(self):
        lattitude = -30.16606
        longitude = 70.81489
        altitude = 2215.0

        self.metadata.set_obs_position(lattitude, longitude, altitude)

        self.assertTrue(self.metadata.obs_loc_set)
        self.assertEqual(self.metadata.obs_lat, lattitude)
        self.assertEqual(self.metadata.obs_long, longitude)
        self.assertEqual(self.metadata.obs_alt, altitude)

    def test_set_obs_code(self):
        obs_code = 'W84'

        self.metadata.set_obs_code(obs_code)

        self.assertTrue(self.metadata.obs_loc_set)
        self.assertEqual(self.metadata.obs_code, obs_code)

    def test_pixels_to_skycoords(self):
        with tempfile.TemporaryDirectory() as dir_name:
            # Create two fake files in the temporary directory.
            fname = "%s/tmp.fits" % dir_name
            create_fake_fits_file(fname, 20, 30)
            self.metadata.populate_from_fits_file(fname)
            sc = self.metadata.pixels_to_skycoords(0,0)
            self.assertEqual(sc, SkyCoord(201.614,-10.788, unit="deg"))

    def test_approximate_radius(self):
        with tempfile.TemporaryDirectory() as dir_name:
            # Create two fake files in the temporary directory.
            fname = "%s/tmp.fits" % dir_name
            create_fake_fits_file(fname, 20, 30)
            self.metadata.populate_from_fits_file(fname)
            radius = self.metadata.approximate_radius()
            self.assertAlmostEqual(radius.arcsecond, Angle(0.01793046, unit='deg').arcsecond, delta = 0.005)

    def test_ra_radius(self):
        with tempfile.TemporaryDirectory() as dir_name:
            # Create two fake files in the temporary directory.
            fname = "%s/tmp.fits" % dir_name
            create_fake_fits_file(fname, 20, 30)
            self.metadata.populate_from_fits_file(fname)
            radius = self.metadata.ra_radius()
            self.assertAlmostEqual(radius.arcsecond, Angle(0.00982375, unit='deg').arcsecond, delta = 0.005)
        
    def test_approximate_radius(self):
        with tempfile.TemporaryDirectory() as dir_name:
            # Create two fake files in the temporary directory.
            fname = "%s/tmp.fits" % dir_name
            create_fake_fits_file(fname, 20, 30)
            self.metadata.populate_from_fits_file(fname)
            radius = self.metadata.dec_radius()
            self.assertAlmostEqual(radius.arcsecond, Angle(0.015, unit='deg').arcsecond, delta = 0.005)

class TestImageMetadataStack(unittest.TestCase):
    def test_init_none(self):
        images = ImageMetadataStack()
        self.assertEqual(images.image_metadatas, [])

    def test_build_files(self):
        with tempfile.TemporaryDirectory() as dir_name:
            # Create two fake files in the temporary directory.
            fname1 = "%s/tmp1.fits" % dir_name
            fname2 = "%s/tmp2.fits" % dir_name
            create_fake_fits_file(fname1, 10, 20)
            create_fake_fits_file(fname2, 20, 30)
            images = ImageMetadataStack()
            images.build_from_filenames([fname1, fname2])
            
            self.assertEqual(images[0]["NAXIS1"], 10)
            self.assertEqual(images[1]["NAXIS2"], 30)
            self.assertEqual(len(images), 2)

if __name__ == '__main__':
    unittest.main()
