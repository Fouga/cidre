CIDRE for tile stitching pipeline
=====

About
-----
This is a fork of the [CIDRE repository](https://github.com/smithk/cidre). The only modification of the code is concerned the way the images are loaded. It is adjusted to our in-house 2-photon microscopy data.   


Details of the [CIDRE](https://github.com/smithk/cidre) algorithm are described in 
<ol>
<li>
K. Smith, Y. Li, F. Ficcinini, G. Csucs, A. Bevilacqua, and P. Horvath<br>
<a href="http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3323.html">CIDRE: An Illumination Correction Method for Optical Microscopy</a>,
Nature Methods, <em>Early Online Access 16 March 2015</em>, doi:10.1038/NMETH.3323
</li>
</ol>.

Contents
--------

- ``matlab`` a modified Matlab implementation of CIDRE.


How to use
------------

In odre to run CIDRE correction in the [Stichit](https://github.com/BaselLaserMouse/StitchIt), you need to download this repository and add the folder to the path in the Stitchit pipeline(https://github.com/BaselLaserMouse/StitchIt). An example of use you can see on [Stichit fork](https://github.com/Fouga/StitchIt). 

Acknowledgements
============

- Rob Campbell for [Stichit](https://github.com/BaselLaserMouse/StitchIt) algorithm.
- Kevin Smith, 2015 for [CIDRE](https://github.com/smithk/cidre). 
- And Gus Brown for recurcive directory listing [rdir](https://uk.mathworks.com/matlabcentral/fileexchange/19550-recursive-directory-listing).
