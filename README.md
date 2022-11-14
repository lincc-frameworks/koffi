# koffi
Known Objects From Fits Indices

A package that should hopefully be helpful to people working with stack-and-shift alogorithms and other solar system science use cases. Given an input of possible solar system objects (PotentialSources) relating either to an x and y coordinate of a FITS image or a ra dec coordinate with a timestamp, we provide an easy API for converting x and y coordinates into ra and dec and querying the most well known solar system dynamics api services, the IMCCE's SkyBoT VO tool and JPL’s SSD (Solar System Dynamics) API service.

This package is based off of code developed by Jeremy Kubica for the [KBMOD](https://github.com/dirac-institute/kbmod) package.
