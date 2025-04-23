# halflife.c

A program to fit a Gaussian convoluted with an exponential decay.

This is ideally suited to fitting time spectra for radioactive decay
(exponential decay) while at the same time taking into account detector
resolution (Gaussian).

For a derivation of the mathematical convolution function used see file
[**halflife_gauss_exp_conv.pdf**](halflife_gauss_exp_conv.pdf)

Program `halflife` accepts the following formats:

- ASCII: 1 (y), 2 (x,y) and 3 (x,y,dy) column formats;
- ORTEC (MAESTRO) .Spe (1-column ASCII) and .Chn (binary) formats with headers
and trailers.

On achieving a good fit `halflife` writes out a .fit file with the
parameters values as well as the fit function and individually all the
components of the fit.

The integral of the whole spectrum and of each component over the spectrum
range is also written to the output .fit file

The current version of **`halflife`** uses
[**xmgrace**](https://plasma-gate.weizmann.ac.il/Grace/)
for plotting the results of the fit.

The 'Auto-fill initial parameters' option can be disabled from appearing
in the menu list by setting the `AUTO` macro to 0 prior to compiling.
This is useful for teaching settings.

The full list of `halflife` options is:

1. Enter filename and read data;
2. Enter initial parameters;
3. Fix or free parameters;
4. Perform fit and write output to file;
5. Compress data by factor of 2;
6. Explore background value as function of chi^2;
7. Print current values of parameters to screen;
8. Write output with current values of parameters or different point spacing;
9. Write input to file as 2 (x y) or 3 (x y dy) column data;
0. Auto-fill initial parameters and perform fit;
1. Display spectrum;
2. Kill all or selected plot windows;
3. Set or adjust spectrum limits for fitting;
0. Quit.
