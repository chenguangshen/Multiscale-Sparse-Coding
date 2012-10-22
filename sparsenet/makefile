MEX = /Applications/MATLAB_R2012a.app/bin/mex 
NRFDIR = ./nrf
MEXT = mexmaci64

cgf.$(MEXT): cgf.c
	$(MEX) -I$(NRFDIR) -L$(NRFDIR) -lnrfopt cgf.c
