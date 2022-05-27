# GMDApp
This is a short implementation in Julia of the PEER Ground Motion Database Web Application Algorithm for Time Series Selection.
The app is an extension of the PEER procedure for user-defined databases.

Compiled binaries can be found [here](https://github.com/gaaraujo/GMDApp/releases/tag/v0.2.0).

## Authors:
- Diego Casas, M.Sc. student at Federal University of Paraná
- Gustavo Araujo, Ph.D. student at Oregon State University
- Alexander Arciniegas, M.Sc.

## Notes:
This app uses the Gaston.jl module to plot, which requires gnuplot (v5.2 is recommended). If gnuplot is not installed, GMDApp will still work but will not show any plots. Make sure to enable the option "Add application directory to your PATH environment variable" when installing gnuplot.