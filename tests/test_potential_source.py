import unittest
import tempfile

from koffi import ImageMetadataStack
from .koffi_test_helpers import *

from koffi import PotentialSource

class TestImageMetadata(unittest.TestCase):
    def test_init(self):
        ps = PotentialSource()

        self.assertEqual(ps.position_at, {})
        self.assertEqual(ps.times, None)

    def test_build_from_times_and_known_positions(self):
        pos = [
            [200.501433, -14.166194],
            [200.502433, -14.166194],
            [200.503433, -14.166194],
            [200.504433, -14.166194],
            [200.505433, -14.166194]
        ]

        mjds = [
            57131.25239706122,
            57131.25888067302,
            57133.24996839132,
            57133.25643622691,
            57133.26293793291
        ]

        ps = PotentialSource()
        ps.build_from_times_and_known_positions(pos, mjds)

        self.assertEqual(ps[57133.25643622691][0], 200.504433)

    def test_build_from_images_and_xy_position(self):
        pos = [
            [0,0]
        ]

        with tempfile.TemporaryDirectory() as dir_name:
            fname = "%s/tmp.fits" % dir_name
            create_fake_fits_file(fname, 20, 30)
            images = ImageMetadataStack([fname])

            ps = PotentialSource()
            ps.build_from_images_and_xy_positions(pos, images)
            self.assertEqual(ps[59806.25][0], 201.614)
            self.assertEqual(ps[59806.25][1], -10.788)    

    def test_build_from_images_and_xy_position_len_mismatch(self):
        pos = [
            [0,0]
        ]

        with tempfile.TemporaryDirectory() as dir_name:
            fname1 = "%s/tmp1.fits" % dir_name
            fname2 = "%s/tmp2.fits" % dir_name
            create_fake_fits_file(fname1, 20, 30)
            create_fake_fits_file(fname2, 20, 30)
            images = ImageMetadataStack([fname1, fname2])
            images[1].set_epoch(Time(59806.30, format='mjd'))

            ps = PotentialSource()
            with self.assertRaises(ValueError) as context:
                ps.build_from_images_and_xy_positions(pos, images)

                self.assertTrue(
                    'number of positions does not match number of images provided' in context.exception
                )


if __name__ == '__main__':
    unittest.main()