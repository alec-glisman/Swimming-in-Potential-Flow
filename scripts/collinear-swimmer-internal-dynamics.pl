#!/usr/bin/env perl

# SECTION: Dependencies

# External dependencies
use Cwd;                                     # current working directory
use Sys::Hostname;                           # get the hostname
use Sys::Info;                               # hardware information
use Sys::Info::Constants qw( :device_cpu );  # CPU information

use File::Path qw(make_path);      # make directories like 'mkdir -p'
use File::Copy;                    # copy files like 'cp'
use File::Slurp;                   # load data files
use File::Tee qw(tee);             # Have std-out and std-err print and save to file like 'tee'
use Path::Tiny qw(path);           # more data I/O

use English;                       # use nice English (or awk) names for ugly punctuation variables
use English qw( -no_match_vars );  # Avoids regex performance penalty in perl 5.18 and earlier
use POSIX;                         # datetime
use Term::ANSIColor;               # color in print statements


use Switch;                        # switch-case control flow statements
use strict;                        # conform to better style and conventions
use warnings;                      # give warnings

# !SECTION (Dependencies)


# SECTION: Input variables that user must specify before running script

# Compiler
my $enableTesting  = "True";                    # Unit tests,        OPTIONS: (False) OFF, (True) ON
my $enableCoverage = "False";                   # Testing coverage,  OPTIONS: (False) OFF, (True) ON
my $build          = "Release";                 # CMake built type,  OPTIONS: Release, Debug, Profile
my $generator      = "Unix Makefiles";          # CMake generator,   OPTIONS: "Unix Makefiles", "Ninja"
my $buildDir       = "build/" . lc $build;      # Build folder path

# C++ Simulation
my $simulationTag    = "collinear-swimmer-internal-dynamics";
my $projectName      = "bodies-in-potential-flow";
my $inputDir         = "input";
my $runSimulationSimulan = 1; # NOTE: 0 only runs one simulation at a time
my @inputData        = ( "varyZHeight", "varyEpsilon" );  # other options: [ "varyDt", "varyRelDisp", "varyPhaseAngle",  "varyVarEpsilon" ]
my $numSimulationTypes = scalar @inputData;

# Simulation parameters
my $numberBodies         = 2; # integer number of bodies to simulate
my $boolImageSystem      = 1; # if system has an image (wall) along z=0 plane
my $boolInternalDynamics = 0; # if swimmer has internal dynamics (0: no, 1: yes)

my $ti          = 0.00e+0;
my $tf          = 1.00e+8;

# Python Numerical Analysis
my $pythonGSDCreation = "python/initial_configurations/" . $simulationTag . "-configuration.py";
my $pythonAggregrateAnalysis   = "python/analysis/" . $simulationTag ."-aggregate-analysis.py";
my $pythonIndividualAnalysis   = "python/analysis/" . $simulationTag ."-individual-analysis.py";


# Host path
my $host = `uname -n`;
chomp(my $home = `echo ~`);
my $cwd           = cwd();

# Output path
my $curDate          = strftime('%Y-%m-%d', localtime);
my $analysisDir      = "figures";

# Host hardware
my $numThreads = 0;
if ((index($host, "MacBook-Pro") != -1) or (index($host, "mbp") != -1)) {
	$numThreads = 12;
} elsif ((index($host, "alec-glisman-PC-Ubuntu") != -1) or (index($host, "alec-glisman-PC-Windows") != -1)) {
	$numThreads  = 20;
} elsif ((index($host, "stokes") != -1 )) {
	$numThreads  = 6;
} elsif ((index($host, "shear") != -1 ) or ( index($host, "s") != -1 )) {
	$numThreads  = 1;
} else {
	my %options;
	my $info = Sys::Info->new;
	my $cpu  = $info->device( CPU => %options );
	my $numThreads = $cpu->count / 2;
}
if ($numThreads == 0) {
	$numThreads = 2;
}

# Terminal display outputs (@DOCS: https://perldoc.perl.org/Term::ANSIColor#Constant-Interface)
print "\n";
print colored ("WELCOME TO PERL SIMULATION SCRIPT", 'bold bright_yellow'), "\n";
print colored ("HOST: $host", 'bold green');
print colored ("NUMBER OF THREADS TO UTILIZE SIMULTANEOUSLY: $numThreads", 'bold green'), "\n\n\n";

# !SECTION (Input variables that user must specify before running script)


# SECTION: Build project and make output directory

# Make build directory and move into it
make_path($buildDir);
chdir $buildDir
  or die "Could not move to build directory: $!";

# Configure and build the project
$ENV{'CLICOLOR_FORCE'} = 0;  # vcpkg fails without this

system( " cmake \"${cwd}\" -G \"${generator}\" -DCMAKE_BUILD_TYPE=${build} -DENABLE_TESTING=${enableTesting} -DENABLE_COVERAGE=${enableCoverage} -DCMAKE_TOOLCHAIN_FILE=\"${cwd}/vcpkg/scripts/buildsystems/vcpkg.cmake\"")
  and die "Configuring project failed: $!";

if($generator eq "Unix Makefiles") {

	# Build project
	system( "make -j" )
	  and die "Building project failed: $!";

	# Run Catch2 unit tests
	if(${enableCoverage} eq "True") {
		system( "make coverage" )
		  and die "Code coverage failed: $!";

	}elsif (${enableTesting} eq "True") {
		system( "make test" )
		  and die "Unit tests failed: $!";
	}
}elsif($generator eq "Ninja") {

	# Build project
	system( "ninja -j ${numThreads}" )
	  and die "Building project failed: $!";

	# Run Catch2 unit tests
	if(${enableCoverage} eq "True") {
		system( "ninja coverage" )
		  and die "Code coverage failed: $!";

	}elsif (${enableTesting} eq "True") {
		system( "ninja test" )
		  and die "Unit tests failed: $!";
	}
}


# Change back to main directory
chdir $cwd
  or die "Could not move back to base directory: $!";

# Make data directory
make_path( "data" );


# Add newline characters
print "\n\n\n";

# !SECTION (Build project and make output directory)


# SECTION: Loop through simulation types
for (my $i = 0; $i < $numSimulationTypes; $i += 1 ){

	# Current datetime
	my $curDateTime = strftime('%Y-%m-%d.%H-%M-%S', localtime);

	# Make temp output directory
	my $tempOutputDir = "temp/" . $curDateTime . "_" . ${simulationTag} . "_" .  ${inputData[$i]};
	make_path( $tempOutputDir );

	# tee the command line outputs
	tee STDERR, '>>', "${tempOutputDir}/std_err.txt";

	# Print the number of simulation runs
	my $ii = $i + 1;
	print "** Simulation type ${ii} of ${numSimulationTypes}: ". ${inputData[$i]} ." **\n";

	# Data file to read from and modify input preferences
	my $fileToOpen = ${inputDir} . '/' . ${inputData[$i]} . '.dat';

	# Get number of simulations from line count in data input file
	my $wcCmd = `wc -l ${fileToOpen}`;
	my @wcLines = split /\s+/, $wcCmd;
	my $numSimulations;
	if ( index(getpwuid $UID, "alec" ) != -1) {
		$numSimulations = $wcLines[1] + 1;
	} else {
		$numSimulations = $wcLines[0] + 1;
	}

	# Load data files
	open( my $dataFile, '<', $fileToOpen )
	  or die "Could not open data file: $!";
	my @data = read_file($dataFile, chomp => 1);
	close( $dataFile)
	  or die "Could not close data file: $!";

	my $simulationIter = 1;

	# Loop through the number of simulations of each type
	for (my $j = 0; $j < $numSimulations; $j += 1 ){

		# Print the number of simulation runs
		my $jj = $j + 1;
		print "Run ${jj} of ${numSimulations} \n";

		# Create individual simulation directory
		my $simulation_dir = ${tempOutputDir} . "/" . "sim_" . ${simulationIter};
		make_path( $simulation_dir );
		my $gsd_path = ${simulation_dir} . "/" . "data.gsd";

		# Simulation variables
		my $dt          = 1.00e+2;

		my $R_avg       = 3.50e+0;
		my $Z_height    = 1.00e+1;

		my $phase_angle = -1.57079632679e+0;
		my $U0          = 1.40e+0;
		my $omega       = 1.00e+0;

		if (!$boolInternalDynamics){
			$U0 = 0.00e+0;
		}

		# Modify default preferences for each simulation run
		switch($inputData[$i]) {
			case ( "varyDt" ) {
				$dt = ${data[$j]};
			}
			case ( "varyRelDisp" ) {
				$R_avg = ${data[$j]};
			}
			case ( "varyZHeight" ) {
				$Z_height = ${data[$j]};
			}
			case ( "varyPhaseAngle" ) {
				$phase_angle = ${data[$j]};
			}
			case ( "varyEpsilon" ) {
				my $epsilon = ${data[$j]};
				$U0 = $epsilon * $R_avg * $omega
			}
			case ( "varyVarEpsilon" ) {  # Golestanian (2004, PRE) Swimmer analog
				 # continuous equivalent collinear swimmer
				my $varEpsilon = ${data[$j]};
				$U0 = $varEpsilon * $omega * 0.50;

				# Recalculate relDispEqbm so that max separation (D) is 10.0
				my $D = 10.0;
				$R_avg = $D - ( 0.50 * $varEpsilon );
			}
		}

		# Generate GSD file
		my $pythonConfigCommand = "python3 ${pythonGSDCreation} --GSD-path=${gsd_path} --dt=${dt} --ti=${ti} --tf=${tf} --R-avg=${R_avg} --Z-height=${Z_height} --phase-angle=${phase_angle} --U0=${U0} --omega=${omega} --number-bodies=${numberBodies} --image-system=${boolImageSystem}";

		system( "${pythonConfigCommand}" )
		  and die "Unable to generate GSD file: $?, $!";

		# Prepare for simulation
		make_path( "${simulation_dir}/${analysisDir}" );
		my $shellSimulationCmd = "${buildDir}/src/./" . ${projectName} . " " . ${gsd_path} . " " . ${simulation_dir};
		my $shellPythonCmd = "python3 ${pythonIndividualAnalysis} --relative-path=${simulation_dir} --output-dir=${analysisDir}";


		# ANCHOR: Run executable: [executable] [input gsd] [output directory]
		if ( ($jj < $numSimulations) and ($jj % int($numThreads) != 0) and ($runSimulationSimulan) ) {

			system( "${shellSimulationCmd} && ${shellPythonCmd} &" )
			  and die "Main project executable or Python individual analysis script failed: $?, $!";

			sleep(5);  # brief pause before next simulation

		} else {

			system( ${shellSimulationCmd} )
			  and die "Main project executable failed: $?, $!";

			system( ${shellPythonCmd} )
			  and warn "Python individual analysis script failed: $!";
		}

		++${simulationIter};
	}

	# NOTE: Pause for all simulations to finish
	sleep(60);

	# !SECTION


	# SECTION: Run analysis scripts and clean-up

	# Run the analysis scripts
	make_path( "${tempOutputDir}/${analysisDir}" );
	system( "python3 ${pythonAggregrateAnalysis} --relative-path=${tempOutputDir} --output-dir=${analysisDir}" )
	  and warn "Python aggregate analysis script failed: $!";

	# Move all output into the "data" directory
	my $outputDir = "data/${curDateTime}" . "_" . ${simulationTag} . "_" . ${inputData[$i]};
	make_path($outputDir);
	system( "mv ${tempOutputDir}/* ${outputDir}" )
	  and die "Moving temporary output to final output failed: $!";

	# Delete temporary output
	system( "rm -rf ${tempOutputDir}" )
	  and die "Removing temporary output directory failed: $!";

	# Add newline character on std-out
	print "\n";

	# !SECTION
}
