# Purpose
This repo provides a C interface to the 2decomp&FFT library originally written in Fortran.

# Use

This interface should be integrated into the original library.
An extra library will be created called *2decomp_fft_iface*.

## Install
In order to install the interface, do:
```
make PREFIX=<2decomp_root_path> install
```

Or, if this repo is placed at the root of the 2decomp&FFT library, just do:
```
make install
```

Note: you also can create a make.inc that sets the variable PREFIX permanently.
If needed, do `make` to get a help.

