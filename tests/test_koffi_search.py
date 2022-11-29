import io
import json
import tempfile
import unittest
import urllib.request as libreq
from unittest import mock

from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.table import QTable
from astropy.time import Time
from astropy.wcs import WCS

from koffi import *

from .koffi_test_helpers import *


class TestKoffiSearch(unittest.TestCase):
    def setUp(self):
        self.mock_skybot_results = QTable()
        self.mock_skybot_results["Number"] = [240334]
        self.mock_skybot_results["Name"] = ["2003 OH28"]
        self.mock_skybot_results["RA"] = [240334]
        self.mock_skybot_results["DEC"] = [-10.78208]

        self.jpl_data = io.BytesIO(
            b'{"n_second_pass": 1,"data_second_pass": [["(2013 FD28)", "13:22:00.37", "-14 09\'58.3"]]}'
        )

        self.jpl_data_second = io.BytesIO(
            b'{"n_second_pass": 1,"data_second_pass": [["(2013 FD28)", "13:22:00.37", "-14 09\'58.3"]]}'
        )

        self.source_skybot = PotentialSource()
        self.source_skybot.build_from_times_and_known_positions([[240334, -10.78208]], [59806.25])

        self.source_jpl = PotentialSource()
        self.source_jpl.build_from_times_and_known_positions([[200.501562, -14.166247]], [59806.25])

        self.source_jpl_mult = PotentialSource()
        self.source_jpl_mult.build_from_times_and_known_positions(
            [[200.501562, -14.166247], [200.501562, -14.166247]], [59806.25, 59806.30]
        )

    @mock.patch("astroquery.imcce.Skybot.cone_search")
    def test_skybot_query_known_objects(self, mock_cone_search):
        mock_cone_search.return_value = self.mock_skybot_results

        with tempfile.TemporaryDirectory() as dir_name:
            fname = "%s/tmp.fits" % dir_name
            create_fake_fits_file(fname, 20, 30)
            metadata = ImageMetadata(fname)
            results = skybot_query_known_objects([self.source_skybot], metadata)

            # mock_cone_search.assert_called()
            self.assertEqual(results[0][0], 0)
            self.assertEqual(results[0][1][0], "2003 OH28")
            self.assertEqual(results[0][1][1].ra.arcsecond, 770400.0)
            self.assertEqual(results[0][1][1].dec.arcsecond, -38815.488000000005)

    @mock.patch("astroquery.imcce.Skybot.cone_search")
    def test_skybot_query_no_sources(self, mock_cone_search):
        mock_cone_search.return_value = self.mock_skybot_results

        with tempfile.TemporaryDirectory() as dir_name:
            fname = "%s/tmp.fits" % dir_name
            create_fake_fits_file(fname, 20, 30)
            metadata = ImageMetadata(fname)
            results = skybot_query_known_objects([], metadata)

            # mock_cone_search.assert_called()
            self.assertEqual(results, [])

    @mock.patch("astroquery.imcce.Skybot.cone_search")
    def test_skybot_query_known_objects_stack(self, mock_cone_search):
        mock_cone_search.return_value = self.mock_skybot_results

        with tempfile.TemporaryDirectory() as dir_name:
            fname1 = "%s/tmp1.fits" % dir_name
            fname2 = "%s/tmp2.fits" % dir_name
            create_fake_fits_file(fname1, 20, 30)
            create_fake_fits_file(fname2, 20, 30)

            images = ImageMetadataStack([fname1, fname2])
            images[1].set_epoch(Time(59806.30, format="mjd"))

            source = PotentialSource()
            source.build_from_times_and_known_positions(
                [
                    [240334, -10.78208],
                    [240334, -10.78208],
                ],
                [59806.25, 59806.30],
            )

            results = skybot_query_known_objects_stack([source], images)

            # mock_cone_search.assert_called()
            self.assertEqual(mock_cone_search.call_count, 2)
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0]["2003 OH28"], 2)

    @mock.patch("astroquery.imcce.Skybot.cone_search")
    def test_skybot_query_stack_no_source(self, mock_cone_search):
        mock_cone_search.return_value = self.mock_skybot_results

        with tempfile.TemporaryDirectory() as dir_name:
            fname1 = "%s/tmp1.fits" % dir_name
            fname2 = "%s/tmp2.fits" % dir_name
            create_fake_fits_file(fname1, 20, 30)
            create_fake_fits_file(fname2, 20, 30)

            images = ImageMetadataStack([fname1, fname2])
            images[1].set_epoch(Time(59806.30, format="mjd"))

            results = skybot_query_known_objects_stack([], images)

            # mock_cone_search.assert_called()
            self.assertEqual(mock_cone_search.call_count, 2)
            self.assertEqual(len(results), 0)

    def test_create_jpl_query_string(self):
        regex_loc = ".*lat=-30.166060&lon=70.814890&alt=2215.000000&obs-time=2459806.750000.*"
        regex_fov = ".*fov-ra-lim=13-26-27.40,13-26-32.12&fov-dec-lim=M10-47-16.80,M10-45-28.80"
        with tempfile.TemporaryDirectory() as dir_name:
            fname = "%s/tmp.fits" % dir_name
            create_fake_fits_file(fname, 20, 30)
            image = ImageMetadata(fname)

            jpl_query_string = create_jpl_query_string(image)
            self.assertRegex(jpl_query_string, regex_loc)
            self.assertRegex(jpl_query_string, regex_fov)

    def test_jpl_query_known_objects(self):
        with mock.patch.object(libreq, "urlopen", return_value=self.jpl_data):
            with tempfile.TemporaryDirectory() as dir_name:
                fname = "%s/tmp.fits" % dir_name
                create_fake_fits_file(fname, 20, 30)
                image = ImageMetadata(fname)

                results = jpl_query_known_objects([self.source_jpl], image)
                self.assertEqual(results[0][0], 0)
                self.assertEqual(results[0][1][0], "(2013 FD28)")
                self.assertEqual(results[0][1][1].ra.arcsecond, 721805.55)
                self.assertEqual(results[0][1][1].dec.arcsecond, -50998.3)

    def test_jpl_query_known_objects_stack(self):
        with mock.patch.object(libreq, "urlopen", side_effect=[self.jpl_data, self.jpl_data_second]):
            with tempfile.TemporaryDirectory() as dir_name:
                fname1 = "%s/tmp1.fits" % dir_name
                fname2 = "%s/tmp2.fits" % dir_name
                create_fake_fits_file(fname1, 20, 30)
                create_fake_fits_file(fname2, 20, 30)
                images = ImageMetadataStack([fname1, fname2])
                images[1].set_epoch(Time(59806.30, format="mjd"))

                results = jpl_query_known_objects_stack([self.source_jpl_mult], images)
                self.assertEqual(len(results), 1)
                self.assertEqual(results[0]["(2013 FD28)"], 2)

    @mock.patch("astroquery.imcce.Skybot.cone_search")
    def test_skybot_search_frame(self, mock_cone_search):
        mock_cone_search.return_value = self.mock_skybot_results

        with tempfile.TemporaryDirectory() as dir_name:
            fname = "%s/tmp.fits" % dir_name
            create_fake_fits_file(fname, 20, 30)
            metadata = ImageMetadata(fname)
            results = skybot_search_frame(metadata)

            mock_cone_search.assert_called()
            self.assertEqual(results[0][0], "2003 OH28")
            self.assertEqual(results[0][1].ra.arcsecond, 770400.0)
            self.assertEqual(results[0][1].dec.arcsecond, -38815.488000000005)

    def test_jpl_query_search_frame(self):
        with mock.patch.object(libreq, "urlopen", return_value=self.jpl_data):
            with tempfile.TemporaryDirectory() as dir_name:
                fname = "%s/tmp.fits" % dir_name
                create_fake_fits_file(fname, 20, 30)
                image = ImageMetadata(fname)

                results = jpl_search_frame(image)
                self.assertTrue(results is not None)
                self.assertEqual(results[0][0], "(2013 FD28)")
                self.assertEqual(results[0][1].ra.arcsecond, 721805.55)
                self.assertEqual(results[0][1].dec.arcsecond, -50998.3)
