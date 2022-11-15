from astropy.time import Time
from astropy.io import fits
from astropy.wcs import WCS

class ImageMetadata:
    def __init__(self, file = None, mjd_key = 'MJD_OBS', mjd_val = None):
        self.obs_loc_set = False
        self.epoch_set_ = False
        self.wcs = None
        self.center = None
        self.obs_code = ""
        self.hdu_list = None

        if file is not None:
            self.populate_from_fits_file(file, mjd_key, mjd_val)


    def __getitem__(self, element):
        return self.get_header_element(element)

    def populate_from_fits_file(self, filename, mjd_key = 'MJD_OBS', mjd_val = None):
        """
        Read the file stats information from a FITS file.

        Arguments:
            filename : string
                The path and name of the FITS file.
        """
        self.filename = filename
        with fits.open(filename) as hdu_list:
            self.hdu_list = hdu_list.copy()

            # a fix for DECam's weird format, where they don't
            # put the required elements in the primary header
            if "NAXIS1" not in hdu_list[0].header:
                self.wcs = WCS(hdu_list[1].header)
            else:
                self.wcs = WCS(hdu_list[0].header)
            self.width = self.get_header_element("NAXIS1")
            self.height = self.get_header_element("NAXIS2")

            # get the time and set the epoch
            if mjd_val is not None:
                self.set_epoch(Time(mjd_val, format="mjd"))
            else:
                # Patch for DECam data
                if "DATE-AVG" in hdu_list[0].header:
                    self.set_epoch(Time(hdu_list[0].header["DATE-AVG"], format="isot"))
                elif self.get_header_element(mjd_key) is not None:
                    self.set_epoch(Time(self.get_header_element(mjd_key), format="mjd"))

            # Extract information about the location of the observatory.
            # Since this doesn't seem to be standardized, we try some
            # documented versions.
            observat = self.get_header_element("OBSERVAT")
            obs_lat = self.get_header_element("OBS-LAT")
            lat_obs = self.get_header_element("LAT_OBS")
            if observat is not None:
                self.obs_code = observat
                self.obs_loc_set = True
            elif obs_lat is not None:
                self.obs_lat = float(obs_lat)
                self.obs_long = float(self.get_header_element("OBS-LONG"))
                self.obs_alt = float(self.get_header_element("OBS-ELEV"))
                self.obs_loc_set = True
            elif lat_obs is not None:
                self.obs_lat = float(lat_obs)
                self.obs_long = float(self.get_header_element("LONG_OBS"))
                self.obs_alt = float(self.get_header_element("ALT_OBS"))
                self.obs_loc_set = True
            else:
                self.obs_loc_set = False

            # Compute the center of the image in sky coordinates.
            self.center = self.wcs.pixel_to_world(self.width / 2, self.height / 2)

    def set_obs_code(self, obs_code):
        """
        Manually set the observatory code.

        Arguments:
            obs_code : string
               The Observatory code.
        """
        self.obs_code = obs_code
        self.obs_loc_set = True

    def set_obs_position(self, lat, long, alt):
        """
        Manually set the observatory location and clear
        the observatory code.

        Arguments:
            lat : float - Observatory latitude.
            long : float - Observatory longitude.
            alt : float - Observatory altitude.
        """
        self.obs_code = ""
        self.obs_lat = lat
        self.obs_long = long
        self.obs_alt = alt
        self.obs_loc_set = True

    def set_epoch(self, epoch):
        """
        Manually set the epoch for this image.

        Arguments:
            epoch : astropy Time object.
        """
        self.epoch_ = epoch
        self.epoch_set_ = True

    def get_epoch(self):
        """
        Get the epoch for this image.

        Returns:
            epoch : astropy Time object.
        """
        if not self.epoch_set_:
            raise ValueError('epoch (astropy.Time object) has not been set. do so with set_epoch()')
        return self.epoch_

    def pixels_to_skycoords(self, pos):
        """
        Transform the pixel position within an image
        to a SkyCoord.

        Arguments:
            pos : pixel_pos
                A pixel_pos object containing the x and y
                coordinates on the pixel.

        Returns:
            A SkyCoord with the transformed location.
        """
        return self.wcs.pixel_to_world(pos.x, pos.y)

    def approximate_radius(self):
        """
        Compute an approximate radius of the image.

        Arguments: NONE

        Returns:
            A radius in degrees.
        """
        corner = self.wcs.pixel_to_world(0.0, 0.0)
        radius = self.center.separation(corner)
        return radius

    def ra_radius(self):
        edge = self.wcs.pixel_to_world(0.0, self.height / 2)
        radius = self.center.separation(edge)
        return radius

    def dec_radius(self):
        edge = self.wcs.pixel_to_world(self.width / 2, 0.0)
        radius = self.center.separation(edge)
        return radius

    def get_header_element(self, element):
        # small wrapper element to grab a certain element from an
        # hdu_list, since in the case where a FITS file has multiple
        # hdus we can't be sure exactly which one it's in.
        for hdu in self.hdu_list:
            if element in hdu.header:
                return hdu.header[element]
        return None

class ImageMetadataStack:
    def __init__(self, files=None):
        self.image_metadatas = []
        if files is not None:
            self.build_from_filenames(files)

    def __getitem__(self, index):
        return self.image_metadatas[index]

    def __len__(self):
        return len(self.image_metadatas)

    def build_from_filenames(self, filenames):
        self.image_metadatas = []
        for f in filenames:
            img = ImageMetadata()
            img.populate_from_fits_file(f)
            self.image_metadatas.append(img)

