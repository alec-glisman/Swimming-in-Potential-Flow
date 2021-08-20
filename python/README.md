# Directory: python

Python scripts to generate GSD files to input to simulation as well as analyze GSD files output from simulation.

## Subdirectory: analysis

Directory contains Python scripts that analyze GSD files (output from C++ simulation).  

## Subdirectory: initial_configurations

`initial_configurations`: Directory contains Python scripts that create GSD files (input to C++ simulation).  

## Files

`GSDUtil.py`: Python class to create/load GSD files.

- Initializer will either create a GSD and output to input path, or load a GSD from a given path based on the value of `create_gsd`
- `setLogParameters()` outputs integrator, material, and other parameters to GSD.
- `setParticleParameters()` outputs configuration and particle type parameters to GSD.
- `setKinematics()` outputs kinematic data in both float and double data types to GSD.
**Note that the acceleration components are output to moment_inertia as GSD does not store acceleration by default.**

`plotStyling.py`: Python class to create plots in a consistent format using matplotlib.pyplot

- Initializer sets all class parameters
- `make_plot()` should be called after initializer and creates the matplotlib.pyplot object.
If running interactively (jupyter notebook), set `showPlot=True` to use backend that can show plot live.
Otherwise, use `showPlot=False` to use more efficient backend.
- `save_plot()` save plot to desired output path.
- Data plotting options
  - `curve()` plot data as a solid line between data points.  
  - `scatter()` plot data as points.  
  - `scatter_dashed()` plot data as points connected via dashed lines.  
- Plot elements
  - `legend()` add a legend to plot
  - `textbox()` add text to plot
- Axes elements
  - `set_xaxis_scientific()` and `set_yaxis_scientific()` convert specified axis labelling to use scientific notation
  - `set_xaxis_tick_scalar_formatter()` and `set_yaxis_tick_scalar_formatter()` convert specified axis labelling to use traditional decimal notation.
  - `set_major_minor_ticks()` specify tick spacing for major (x) and minor (y) axes.
  - `set_latex_axis_labels()` convert axis label to use LaTeX symbol.
