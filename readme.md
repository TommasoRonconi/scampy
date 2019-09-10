# Library for cosmological modelling and mock-catalogue building

## Installation guide

### Pre-requisites

- `meson` build system tool
- `ninja`

can both be installed either via `conda install` or with `pip install`

### Install

```bash
# from root directory
meson build_dir -Dprefix=/path/to/install_directory
cd build_dir

ninja install
```
If no `-Dprefix` is specified the library will be installed in `/path/to/root_lib_dir/install` directory.
