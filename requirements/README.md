# Directory: requirements

Files and scripts relevant for loading C++, Perl, and Python dependencies.

## Subdirectory: C++

Contains shell script to install all relevant dependencies.
The script makes two assumptions: (1) [homebrew](https://brew.sh/) is installed locally on the computer and (2) the apt package manager exists to install Intel MKL dependencies via Intel oneAPI.

## Subdirectory: Perl

Perl dependencies are installed via the [cpanm](https://metacpan.org/dist/App-cpanminus/view/bin/cpanm) package manager.

```[shell]
# Install Perl modules via cpanm
cd requirements/Perl
cpanm --installdeps .
```

## Subdirectory: Python

Python dependencies are installed via [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) package manager.

```[shell]
# Create conda environment from given file
conda env create -f requirements/Python/environment.yml

# Export conda environment to a given file
conda env export --from-history  > requirements/Python/environment.yml
```
