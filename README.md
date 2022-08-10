# NCpi0-Cross-Section

Using the MINERvA Analysis Toolkit (MAT) and data from the single photon low energy excess analysis for MicroBooNE, uncertainty systematics are evaluated for the calculation of neutral current neutral pion cross-section.

# USAGE

1. Take processed TH1Ds and propogate uncertainties to create MnvH1Ds of the signal, background, event rate, efficiency, flux, number of targets & POT. Define response matrix.
`python translateHists.py <output_dir> <test> <fakedata>`

**Optional Arguments**
`test` uses only 10 systematic universes
`fakedata` scales fake data to the correct POT

OUTPUT - `<date>_out.root`

2. Unfold event rate.
`root -l "unfold.C(\"<infile_dir>/<infile_date>_out\")"`

Must manually uncomment lines to change between Wiener-SVD and D'Agostini unfolding

OUTPUT - `<date>_out_unfolded.root`

3. Calculate cross section.
`python calculateXsection.py <inout_dir> <in_date>`

**Optional Arguments**
`in_date` specifies the creation date of the input hist (from unfold.C). Defaults to a file in inout_dir dated today if argument is not provided.

OUTPUT - `<date>_out_final.root`

4. Make cross section plots.
`python makePlots.py <in_dir> <in_date> <out_dir>`

**Optional Arguments**
`in_date` specifies the creation date of the input hist (from calculateXsection.py). Defaults to a file in in_dir dated today if argument is not provided.
`out_dir` specifies the location of the output folder. Defaults to the same as in_dir if argument is not provided.

OUTPUT - `<date>_xsec-plots`
