from astropy.time import Time
from astropy.io import fits
from astropy.wcs import WCS

class PotentialSource:
    def __init__(self):
        self.position_at = {}
        self.times = None

    def __getitem__(self, time):
        """ given an MJD value return the [ra, dec] position """
        if time not in self.position_at.keys():
            raise ValueError('no location associated with provided time')
        return self.position_at[time]

    def build_from_times_and_known_positions(self, positions, times):
        """
        Build out the position_at dict using times and positions.

        Arguments:
            positions : a list containing lists of [ra, dec] positions.
            times : a list of doubles containing MJD times of given positions.
                NOTE: to avoid headaches, use the MJDs given from the
                ImageMetadata epoch.
        """
        if len(positions) != len(times):
            raise ValueError('number of positions does not match number of times provided')

        for i in range(len(times)):
            self.position_at[times[i]] = positions[i]