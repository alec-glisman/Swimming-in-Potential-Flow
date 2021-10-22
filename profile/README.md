# Directory: profile

Scripts to profile code on various platforms and find performance improvement areas.

## Subdirectory (submodule): [FlameGraph](https://github.com/brendangregg/FlameGraph)

Perl repository that converts outputs of Linux `perf` command into FlameGraphs.

## Subdirectory: Linux

Contains script to profile code using `perf` to record data and `FlameGraph` to analyze results.

## Subdirectory: MacOS

Contains script to profile code using the instruments app ([documentation](https://help.apple.com/instruments/mac/current/#/devb14ffaa5)).
Hard to pass command line arguments that are required. This is made easier by calling instruments via command line rather than GUI.
Sample bash script is in this directory.
