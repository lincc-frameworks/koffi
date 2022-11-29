"""
koffi - Known Objects From Fits Indices
A tool that provides a unified interface to query SkyBoT and JPL's Horizons
APIs using metadata pulled from FITS files and user provided sources.
Developed by Max West (https://github.com/maxwest-uw) as part of work for
the LINCC Frameworks team at the University of Washington.
Based on code written by Jeremy Kubica (https://github.com/jeremykubica) for
the KBMOD package (https://github.com/dirac-institute/kbmod).
Matching is approximate and results might not be comprehensive.
All scientific studies should conduct their own matching analysis.

Acknowledgements: The known object matching uses the IMCCE's SkyBoT VO tool
(Berthier et. al. 2006) and JPL’s SSD (Solar System Dynamics) API service.
"""

import json
import urllib.request as libreq

import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astroquery.imcce import Skybot
from astroquery.jplhorizons import Horizons
from tqdm import tqdm

# Import the data_tools packages directly in koffi so that they are
# accesible at the top level of the package.
from koffi_tools.image_metadata import *
from koffi_tools.potential_source import *


def skybot_search_frame(image):
    """
    Gets all known objects within the frame of a single FITS image from SkyBoT.

    Arguments:
        image : an ImageMetadata object.
    Returns:
        A list of objects, where each object row is in
            the form [Name, SkyCoord].
    """
    image_center = image.center
    image_radius = image.approximate_radius()
    image_epoch = image.get_epoch()

    skybot_results = Skybot.cone_search(image_center, image_radius, image_epoch)

    objects = []
    for row in range(len(skybot_results)):
        name = skybot_results["Name"][row]
        ra = skybot_results["RA"][row]
        dec = skybot_results["DEC"][row]
        row_coord = SkyCoord(ra, dec, unit="deg")
        objects.append([name, row_coord])

    return objects


def skybot_query_known_objects(potential_sources, image, tolerance=0.5):
    """
    Finds all known objects that should appear in an image
    given meta data from a FITS file in the form of an
    ImageMetadata and checks them against a list of PotentialSources
    using the IMCCE's SkyBoT VO tool.

    Arguments:
        potential_sources : List of PotentialSource objects
        image : ImageMetadata object
           The metadata for the current image.
        tolerance : the allowed separation between objects in arcseconds
    Returns:
        matches : a matrix of matches, where row[0] is the index of the
            potential source, row[1] a list with the object name and a SkyCoord
            denoting the location of the given object.
    """

    # Use SkyBoT to look up the known objects with a conesearch.
    # The function returns a list of objects attached to an index of
    # the potential_sources list.
    skybot_frame_objects = skybot_search_frame(image)
    matches = []

    for i in range(len(potential_sources)):
        ps = potential_sources[i]
        time = image.get_epoch().mjd
        ps_coord = SkyCoord(ps[time][0], ps[time][1], unit="deg")
        for row in skybot_frame_objects:
            name = row[0]
            row_coord = row[1]
            sep = ps_coord.separation(row_coord)
            if sep.arcsecond < tolerance:
                matches.append([i, [name, row_coord]])

    return matches


def skybot_query_known_objects_stack(potential_sources, images, tolerance=0.5, min_observations=1):
    """
    Finds all known objects that should appear in a series of
    images given the meta data from the corresponding FITS files
    using the IMCCE's SkyBoT VO tool.

    Arguments:
        potential_sources : List of PotentialSource objects
        images : an ImageMetadataSet holding the metadata for a stack of images.
        tolerance : the allowed separation between objects in arcseconds
        min_obeservations : minimum number of times a source has to be found throughout
            the frames to be returned in the results.
    Returns:
        matches : dictionary using the potential source id as a key that
            contains another dict of found objects with a corresponding count
            of observations (object name is the key).
    """
    matches = {}
    for i in range(len(potential_sources)):
        matches[i] = {}

    for image in tqdm(images):
        frame_sources = skybot_query_known_objects(potential_sources, image, tolerance)
        for res in frame_sources:
            ps_id = res[0]
            obj_name = res[1][0]
            if obj_name in matches[ps_id].keys():
                matches[ps_id][obj_name] += 1
            else:
                matches[ps_id][obj_name] = 1

    for ps_id in matches.keys():
        bad_obs = []
        for obj in matches[ps_id].keys():
            if matches[ps_id][obj] < min_observations:
                bad_obs.append(obj)
        for rem in bad_obs:
            matches[ps_id].pop(rem)

    return matches


def create_jpl_query_string(image):
    """
    Create JPL query string out of the component
    information.

    Argument:
        image : An ImageMetadata object holding the
                metadata for the current image.

    Returns:
        The query string for JPL conesearch queries or None
        if the ImageInfo object does not have sufficient
        information.
    """
    if not image.obs_loc_set or image.center is None:
        return None

    base_url = "https://ssd-api.jpl.nasa.gov/sb_ident.api?sb-kind=a&mag-required=true&req-elem=false"

    # Format the time query and MPC string.
    t_str = "obs-time=%f" % image.get_epoch().jd

    # Create a string of data for the observatory.
    if image.obs_code:
        obs_str = "mpc-code=%s" % self.obs_code
    else:
        obs_str = "lat=%f&lon=%f&alt=%f" % (image.obs_lat, image.obs_long, image.obs_alt)

    # Format the RA query including half width.
    if image.center.ra.degree < 0:
        image.center.ra.degree += 360.0
    ra_hms_L = Angle(image.center.ra - image.ra_radius()).hms
    ra_hms_H = Angle(image.center.ra + image.ra_radius()).hms
    ra_str = "fov-ra-lim=%02i-%02i-%05.2f,%02i-%02i-%05.2f" % (
        ra_hms_L[0],
        ra_hms_L[1],
        ra_hms_L[2],
        ra_hms_H[0],
        ra_hms_H[1],
        ra_hms_H[2],
    )

    # Format the Dec query including half width.
    dec_str = ""
    dec_dms_L = Angle(image.center.dec - image.dec_radius()).dms
    if dec_dms_L[0] >= 0:
        dec_str = "fov-dec-lim=%02i-%02i-%05.2f" % (dec_dms_L[0], dec_dms_L[1], dec_dms_L[2])
    else:
        dec_str = "fov-dec-lim=M%02i-%02i-%05.2f" % (-dec_dms_L[0], -dec_dms_L[1], -dec_dms_L[2])
    dec_dms_H = Angle(image.center.dec + image.dec_radius()).dms
    if dec_dms_H[0] >= 0:
        dec_str = "%s,02i-%02i-%05.2f" % (dec_str, dec_dms_H[0], dec_dms_H[1], dec_dms_H[2])
    else:
        dec_str = "%s,M%02i-%02i-%05.2f" % (dec_str, -dec_dms_H[0], -dec_dms_H[1], -dec_dms_H[2])

    # Only do the second (more accurate) pass.
    pass_str = "two-pass=true&suppress-first-pass=true"

    # Complete the full query.
    query = "%s&%s&%s&%s&%s&%s" % (base_url, obs_str, t_str, pass_str, ra_str, dec_str)

    return query


def jpl_search_frame(image):
    """
    Gets all known objects within the frame of a single FITS image from
        JPL Horizons API.

    Arguments:
        image : an ImageMetadata object.
    Returns:
        A list of objects, where each object row is in
            the form [Name, SkyCoord].
    """
    query_string = create_jpl_query_string(image)
    if not query_string:
        raise ValueError("WARNING: Insufficient data in image_metadata.")

    objects = []

    with libreq.urlopen(query_string) as url:
        feed = url.read().decode("utf-8")
        results = json.loads(feed)

        num_results = results["n_second_pass"]
        for item in results["data_second_pass"]:
            name = item[0]
            ra_str = item[1]
            dec_str = item[2].replace("'", " ").replace('"', "")
            sc = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
            objects.append([name, sc])

    return objects


def jpl_query_known_objects(potential_sources, image, tolerance=0.5):
    """
    Finds all known objects that should appear in an image
    given meta data from a FITS file in the form of an
    ImageMetadata and checks them against a list of PotentialSources
    using NASA JPL’s SSD (Solar System Dynamics) API service.

    Arguments:
        potential_sources : List of PotentialSource objects
        image : ImageMetadata object
           The metadata for the current image.
        tolerance : the allowed separation between objects in arcseconds
    Returns:
        matches : a matrix of matches, where row[0] is the index of the
            potential source, row[1] a list with the object name and a SkyCoord
            denoting the location of the given object.
    """
    query_string = create_jpl_query_string(image)
    if not query_string:
        raise ValueError("WARNING: Insufficient data in image_metadata.")

    observations = jpl_search_frame(image)

    matches = []
    for i in range(len(potential_sources)):
        ps = potential_sources[i]
        time = image.get_epoch().mjd
        ps_coord = SkyCoord(ps[time][0], ps[time][1], unit="deg")
        num_results = len(observations)
        for row in range(num_results):
            name = observations[row][0]
            row_coord = observations[row][1]
            sep = ps_coord.separation(row_coord)
            if sep.arcsecond < tolerance:
                matches.append([i, [name, row_coord]])
    return matches


def jpl_query_known_objects_stack(potential_sources, images, tolerance=0.5, min_observations=1):
    """
    Finds all known objects that should appear in a series of
    images given the meta data from the corresponding FITS files
    using NASA JPL’s SSD (Solar System Dynamics) API service.

    Arguments:
        potential_sources : List of PotentialSource objects
        images : an ImageMetadataSet holding the metadata for a stack of images.
        tolerance : the allowed separation between objects in arcseconds
        min_obeservations : minimum number of times a source has to be found throughout
            the frames to be returned in the results.
    Returns:
        matches : dictionary using the potential source id as a key that
            contains another dict of found objects with a corresponding count
            of observations (object name is the key).
    """
    matches = {}
    for i in range(len(potential_sources)):
        matches[i] = {}

    print("NOTE: JPL Horizons queries are rate limited and can take up to 5 minutes per query to complete.")

    for image in tqdm(images):
        frame_sources = jpl_query_known_objects(potential_sources, image, tolerance)
        for res in frame_sources:
            ps_id = res[0]
            obj_name = res[1][0]
            if obj_name in matches[ps_id].keys():
                matches[ps_id][obj_name] += 1
            else:
                matches[ps_id][obj_name] = 1

    for ps_id in matches.keys():
        bad_obs = []
        for obj in matches[ps_id].keys():
            if matches[ps_id][obj] < min_observations:
                bad_obs.append(obj)
        for rem in bad_obs:
            matches[ps_id].pop(rem)

    return matches
