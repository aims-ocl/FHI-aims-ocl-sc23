How to run the script:

- First compile the fortran code "home_made_ft.f90". The output has to be called "home_made_ft.x"

- Copy the python script "auto-correlate-XX.py" and the binary "home_made_ft.x" to the directory you want to work on. 
  If you don't want to do so, you can edit the paths in the script accordingly.

- Use:

  python auto-correlate-XX.py -h

  to learn how to use the scripts.

The PI version should be used for path-integral simulations only.


- Four files will be generated, namely: 
   - autocorr.dat: containing 3 columns, the first being time in ps, the second being the autocorrelation function, 
     and the third being the autocorrelation function times a windowing function that makes it go to zero on the edges. 
   - raw_fourier_transform.dat: containing 2 columns, the first being wavenumbers in cm-1 and the second the intensities in arbitrary units
   - norm_raw_fourier_transform.dat: containing 2 columns, the first being wavenumbers in cm-1 and the second the intensities normalized by the total length of the FT
   - convoluted_fourier_transform.dat: containing 2 columns, the first being wavenumbers in cm-1 and the second the intensities in arbitrary 
     units convoluted with a gaussian curve.

Happy usage.


MR 2009
