# VERSION 0.9.5 (24 Oct 2026)
- Nmeneyi: New type of plot, with plottype="ranges"
- Nemenyi: Forced ordering of methods in plottype="matrix" for better readabilty
- Nemenyi: Added output of model groups in console and plots (with option)
- Nemenyi: Respects previous graphical settings, and better adjustment of plotting area to fit labels
- Nemenyi: Fixed a few typos in the help file and added some more error trapping.
- Sthief: Improved input handling and allow to skip levels of the hierarchy. 
- plotSThief: This is now depreciated and Will be removed in following versions of the package. Use the new plotThief instead.
- New: plotS: General function to plot hierarchies from an S-matrix. Calls plotSblock and plotStree.
- New: plotSblock: New function to plot hierarchies using the block representation. 
- New: plotStree: New function to plot hierarchies using the typical tree representation
- New: plotThief: Replaces plotSThief, relies on the other new functions. 

# VERSION 0.9.4 (12 Nov 2023)
- Fix for seasplot graphs
- Fix for tsutils.Rd
- Fix documentation formatting

# VERSION 0.9.3 (07 Jul 2022)
- Recompiled for new version of R
- Added more data input filters for nemenyi.

# VERSION 0.9.2 (06 Feb 2020)
- Added functions Sthief and plotSthief.
- Fixed class() check for abc.R, xyz.R and theta.R

# VERSION 0.9.1
- Fixed output list names for lambdaseq.
- Fixed example for lambdaseq.

# INITIAL RELEASE 0.9.0 (01 Jan 2018)
- Initial release on CRAN. 
