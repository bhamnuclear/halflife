# halflife.c and halflife_double.c

**Program `halflife`** to fit a Gaussian convoluted with an exponential decay.

**Program `halflife_double`** to fit back-to-back curves of a Gaussian
convoluted with an exponential decay and a prompt Gaussian component.

Ideally suited to fitting time spectra for radioactive decay
(exponential decay) while at the same time taking into account detector
resolution (Gaussian).

For a derivation of the mathematical convolution function used see file
[**halflife_gauss_exp_conv.pdf**](halflife_gauss_exp_conv.pdf)

Programs `halflife` and `halflife_double` accept the following formats:

- ASCII: 1 (y), 2 (x,y) and 3 (x,y,dy) column formats;
- ORTEC (MAESTRO) .Spe (1-column ASCII) and .Chn (binary) formats with headers
and trailers.

On achieving a good fit `halflife` and `halflife_double` write out
a .fit file with the parameter values as well as the fit function and,
individually, all the component functions of the fit.

The integrals of the whole spectrum and of each component over the spectrum
range are also written to the output .fit file

The current versions of **`halflife`** and **`halflife_double`** use
[**xmgrace**](https://plasma-gate.weizmann.ac.il/Grace/)
for plotting the results of the fit. The plot scales are set automatically
for each spectrum.

The 'Auto-fill initial parameters' option can be prevented from appearing
in the menu list by setting the `AUTO` macro to 0 prior to compiling, but is
always enabled and accessible via option *a* in the menu.
This is useful for teaching settings.

The full list of `halflife` and `halflife_double` options is:

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
