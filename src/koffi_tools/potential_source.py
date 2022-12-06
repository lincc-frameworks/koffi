from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS


class PotentialSource:
    def __init__(self):
        self.position_at = {}
        self.times = None

    def __getitem__(self, time):
        """given an MJD value return the [ra, dec] position"""
        if time not in self.position_at.keys():
            raise ValueError("no location associated with provided time")
        return self.position_at[time]

    def build_from_times_and_known_positions(self, positions, times):
        """
        Build out the position_at dict using times and positions.

        Arguments:
            positions : a list containing lists of [ra, dec] positions.
            times : a list of doubles containing MJD times of given positions.
                NOTE: to avoid headaches, use the MJDs given from the
                ImageMetadata. You can conveniently get such a list from
                an ImageMetadataStack using .get_mjds().
        """
        if len(positions) != len(times):
            raise ValueError("number of positions does not match number of times provided")

        for i in range(len(times)):
            self.position_at[times[i]] = positions[i]

    def build_from_images_and_xy_positions(self, positions, images):
        """
        Build out the position_at dict using x&y coordinates and an
            ImageMetadataStack

        Arguments:
            positions : a list containing lists of [x, y] positions.
            images : an ImageMetadataStack.
        """
        if len(positions) != len(images):
            raise ValueError("number of positions does not match number of images provided")

        mjds = images.get_mjds()
        for i in range(len(mjds)):
            x, y = positions[i]
            position = images[i].pixels_to_skycoords(x, y)
            self.position_at[mjds[i]] = [position.ra.degree, position.dec.degree]
