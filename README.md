# Autocorrelation Function of Brownian Motion

## Summary

A sample code for calculating autocorrelation function (ACF) of
Brownian motion. It is also a simple example on how to use
fftw3. The ACF of momenta will be computed by using DFT.
FFTW 3.0 or higher is required to build this sample.

## Usage

    $ make clean
    $ make graph

You will obtain a graph which shows a comparison between
ACF and the theoretical value.

## Files

* main.cc: A source file
* acf.dat: Data file of ACF.
* acf.plt: Plot file for gnuplot
* acf.png: Plot image by gnuplot
