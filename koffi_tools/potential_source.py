from astropy.time import Time
from astropy.io import fits
from astropy.wcs import WCS

class PotentialSource:
    def __init__(self, times = None):
        self.position_at = {}
        self.times = times

    def __getitem__(self, time):
        if time not in self.position_at.keys():
            raise ValueError('no location associated with provided time')
        return self.position_at[time]

    def build_from_times_and_known_positions(self, positions, times):
        if len(positions) != len(times):
            raise ValueError('number of positions does not match number of times provided')

        for i in range(len(times)):
            self.position_at[times[i]] = positions[i]