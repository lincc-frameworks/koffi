# koffi
Known Objects From Fits Indices

A package that should hopefully be helpful to people working with shift-and-stack algorithms and other solar system science use cases. Given an input of possible solar system objects (PotentialSources) relating either to an x and y coordinate of a FITS image or an RA and Dec coordinate with a timestamp, we provide an easy API for converting x and y coordinates into RA and Dec and querying the most well known solar system dynamics api services. Uses the IMCCE's SkyBoT VO tool (Berthier et. al. 2006) and JPL’s SSD (Solar System Dynamics) [API service](https://ssd.jpl.nasa.gov/).

This package is based off of code developed by Jeremy Kubica for the [KBMOD](https://github.com/dirac-institute/kbmod) package.

## Setup
### Install from PyPI
koffi is now pip installable! To get the latest version, just run:
```bash
pip install koffi
```
### Install from source
In python virtual environment of your choice:
```bash
git clone https://github.com/lincc-frameworks/koffi.git
cd koffi
pip install .
python setup.py build
```

## Usage
At the most basic level, koffi can be used to find all the known objects based on FITS image metadata.
```python
import koffi

filename = 'path/to/image/data.fits'
image = koffi.ImageMetadata(filename)

# SkyBoT - get all possible objects in given image.
skybot_objects = koffi.skybot_search_frame(image)

# JPL Horizons = get all possible objects in given image.
jpl_objects = koffi.jpl_search_frame(image)
```

We can also use koffi to check against a list of provided sources; the search functions will return a list of matches.
```python
# a possible observation in ra dec + the observation time in mjd.
position = [200.501562, -14.166247]
time = image.get_epoch().mjd

# the PotentialSource class is our interface for checking possible discoveries against known objects.
# position at time can be access by ps[time]
ps = koffi.PotentialSource()
ps.build_from_times_and_known_positions([position], [time])

# return a list of possible matches, attached to the index of a potential source.
skybot_observations = koffi.skybot_query_known_objects([ps], image, tolerance = 0.25)

jpl_observations = koffi.jpl_query_known_objects([ps], image, tolerance = 0.25)
```

At the highest level, we can pass in a stack of images and potential sources and count the number of times they potentially appear in each image!
We do this by making use of the ImageMetadataStack class.
```python
filenames = [
	'/path/to/data/434593.fits', 
	'/path/to/data/434601.fits', 
	'/path/to/data/435478.fits', 
	'/path/to/data/435486.fits', 
	'/path/to/data/435494.fits'
]
# ImageMetadataStack will build the stack for you, just provide it with
# a list of filenames!
images = koffi.ImageMetadataStack(filenames)

# now let's make 2 PotentialSource objects, 2 different ways.
# we can build for a list of ra dec pairs...
positions1 = [
	[200.501433, -14.166194],
	[200.502433, -14.166194],
	[200.503433, -14.166194],
	[200.504433, -14.166194],
	[200.505433, -14.166194]
]
# or from a list of x and y coordinates in the image(s)!
positions2 = [
  [97, 200],
  [101, 204],
  [105, 208],
  [109, 212],
  [113, 216]
]

ps1 = koffi.PotentialSource()
ps2 = koffi.PotentialSource()

# it's highly recommended that you associated your potential source positions with times in mjd given
# from the ImageMetadataStack.get_mjds() method!
mjds = images.get_mjds()

ps1.build_from_times_and_known_positions(positions1, mjds)
ps2.build_from_images_and_xy_positions(positions2, mjds)
sources = [ps1, ps2]

# stack search functions will return a count of each object that it was able to be associated with a potential source.
# key for the dict is the index of potential sources.
possible_detections_skybot = koffi.skybot_query_known_objects_stack(sources, images)
possible_detections_jpl = koffi.jpl_query_known_objects_stack(sources, images)
```

