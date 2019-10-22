[![DOI](https://zenodo.org/badge/215828311.svg)](https://zenodo.org/badge/latestdoi/215828311)

# gemini_reduction  -  Support for reducing Gemini-N and Gemini-S GMOS imaging data


### Environment creation and dependency installation:
    conda create -n geminiconda python=2.7 iraf-all pyraf-all stsci gemini
    conda install -y -c astropy     aplpy

### Example usage
    python gmos_reduce_image_example.py
    
This script implements:
- download data from Gemini archive given an object name
- associate and download corresponding calibrations (bias, flat)
- generate the master calibration files
- reduce the imaging data using the master calibrations
- generate pdf images
    
    
     


### Documentation

That's it, sorry.

### Contributing
Please open a new issue or new pull request for bugs, feedback, or new features you would like to see. If there is an issue you would like to work on, please leave a comment and we will be happy to assist. New contributions and contributors are very welcome!   
 

### References

Programmatic interface to Gemini archive implemented with the help of
https://archive.gemini.edu/help/api.html


### Author

Johannes Sahlmann, 2016 --

