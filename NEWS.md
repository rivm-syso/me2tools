# me2tools 0.2.1

This release features a breaking change for the function `metcor_import`.

- **BREAKING CHANGE**: The function `metcor_import` has undergone an overhaul to support grid sizes for the Northern hemisphere, Southern hemisphere and the world through the new `extent` parameter in the function call. Additionally, MetCor originally classifies the longitudes in the range [0,360], which might cause issues when plotting the raster. A new option has been implemented to output the SpatRaster using the default longitudes in the range of [-180, 180]. This format is now also the default output format. The old format, featuring columns from [0,360], can be obtained by setting `type = "metcor"` within the function call. This change has no influence of the plotting of the raster using `metcor_plot()`, which continues to work with both formats.

For a complete overview of the changes, please go through the [list of commits](https://github.com/rivm-syso/me2tools/commits/main)

# me2tools 0.2.0

This release comes with a few breaking changes and major updates to some of the functions. One of the major additions to this version is the reading and plotting of the results from the error estimate analysis (BS, DISP, BSDISP). Also, for some functions the documentation has been updated. A short overview of the changes are listed below

- **BREAKING CHANGE**: The function `me2_DISP_read_minmax()` has been renamed to `me2_DISP_read_res()`. Internally the old function name now points to the new function and a deprecated message is shown when this is used. Please update your references accordingly.
- **BREAKING CHANGE**: In an attempt to standardize the output of various functions the output of `me2_DISP_read_res()` and `me_BS_read_F()` have been updated to a list object with additional items.
- **NEW FUNCTIONS**: A couple of new functions for reading the error estimates from BSDISP (`me2_BSDISP_read_res()`), plotting an error summary based on DISP, BS and BSDISP results (`epa_plot_errorsummary()`), the boxplot for the BS results (`epa_plot_BSboxplot()`), the contributions from the G matrix (`epa_plot_contribution`) and the plotting of scaled residuals (`epa_plot_residuals()`) have been added.
- The function `compare_obs_mod()` has been completely overhauled with more options (i.e., regression with intercept, showing the CI and/or PI intervals). A bug regarding a wrong formula used for linair regression has also been fixed.
- To accomodate the information available for the error estimates BS, DISP and BSDISP the function `epa_plot_profile()` has been updated and is now able to plot the errorbars for one of the BS, DISP or BSDISP results.

For a complete overview of the changes, please go through the [list of commits](https://github.com/rivm-syso/me2tools/commits/main)

# me2tools 0.1.0

- Added a `NEWS.md` file to track changes to the package.
