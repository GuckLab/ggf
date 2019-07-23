Sungrid Engine files
====================
These scripts are used to compute the GGF look-up table (LUT) and
the data for various visualizations in the documentation.
Tested with jobmanager 0.1.0 and Python 3.6.


lut.py
------
server: `python lut.py server`
clients: `python lut.py client`

Distributed computation of LUTs defined in `lut.py`.

scan.py
-------
server: `python scan.py server`
clients: `python scan.py client`

Distributed computation of 1D line-scans defined in `scan.py`.

plot_scan.py
------------
Visualization of line scans (to get an idea about the necessary
spacing in `lut.py`).

compress_lut.py
---------------
Compress computed LUTs for smaller distribution sizes.
