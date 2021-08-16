#!/usr/bin/env perl

# SECTION: Dependencies

# External dependencies
use Cwd;                          # get the current working directory
use Sys::Hostname;                # get the hostname
use English qw( −no_match_vars ); 
use File::Path qw(make_path);     # make directories like 'mkdir -p'
use File::Copy;                   # copy files like 'cp'
use File::Slurp;                  # load data files
use File::Tee qw(tee);            # Have std-out and std-err print and save to file like 'tee'
use Path::Tiny qw(path);          # more data I/O
use POSIX;                        # get the datetime
use Switch;                       # switch-case control flow statements
use strict;                       # conform to better style and conventions
use warnings;                     # give warnings

# !SECTION (Dependencies)



# SECTION: Input variables that user must specify before running script

# Input variables
my $simulationTag    = "collinear-swimmer";
my $projectName      = "bodies_in_potential_flow";
my $inputDir         = "input";
my @inputData        = ( "varyDt", "varyEpsilon", "varyPhaseAngle", "varyRelDisp" );
my $numSimulationTypes = scalar @inputData;
my $runSimulationSimulan = 1; # 0 only runs one simulation at a time

# Python variables
my $pythonGSDCreation = "python/initial_configurations/" . "collinear-swimmer-configuration.py";
my $pythonAnalysis   = "python/analysis/" . "collinear-swimmer-analysis.py";

# Output variables
my $curDate          = strftime('%Y-%m-%d', localtime);
my $analysisDir      = "analysis";

# Get the name of computer script is running on
my $host = `uname -n`;
chomp(my $home = `echo ~`);

# Compiler variables
my $compiler = "";
my $numThreads = "";
if ( (index($host, "MacBook-Pro") != -1) or (index($host, "MBP") != -1) ) {
    $compiler    = "macOS"; 
    $numThreads = "16";
} elsif (( index($host, "Alec-Glisman-PC-Ubuntu") != -1 ) or ( index($host, "Alec-Glisman-PC-Windows") != -1 ) or ( index($host, "s") != -1 )) {
    $compiler    = "PC-Ubuntu";
    $numThreads  = "24";
} elsif (( index($host, "shear") != -1 ) or ( index($host, "s") != -1 )) {
    $compiler    = "shear";
    $numThreads  = "";
}

my $build         = "Release";                # OPTIONS: Release, Debug, Profile
my $enableTesting = "True";                   # OPTIONS: (False) OFF, (True) ON
my $buildDir      = "build";                  # Title whatever you want build folder to be
my $generator     = "Unix Makefiles";         # ONLY TESTED WITH UNIX
my $cwd           = cwd();

# !SECTION (Input variables that user must specify before running script)



# SECTION: Build project and make output directory

# Make build directory and move into it
make_path($buildDir);
chdir $buildDir
    or die "Could not move to build directory: $!";

# Configure and build the project
system( "cmake \"${cwd}\" -G \"${generator}\" -DMY_COMPILER_OPTION=${compiler} -DCMAKE_BUILD_TYPE=${build} -DENABLE_TESTING=${enableTesting}" ) 
    and die "Configuring project failed: $!";
system( "make -j" ) 
    and die "Building project failed: $!";

# Change back to main directory
chdir $cwd 
    or die "Could not move back to base directory: $!";

# Make data directory
make_path( "data" );

# Run Catch2 unit tests
if (${enableTesting} eq "True") {
    system( "\"${buildDir}/tests/./tests\"" ) 
            and die "Unit test failed: $!";
}

# Add newline characters
print "\n\n\n";

# !SECTION (Build project and make output directory)



# SECTION: Loop through simulation types

for (my $i = 0; $i < $numSimulationTypes; $i += 1 )
{

    # Current datetime
    my $curDateTime = strftime('%Y-%m-%d.%H-%M-%S', localtime);

    # Make temp output directory
    my $tempOutputDir = "temp_" . $curDateTime . "_" . ${simulationTag} . "_" .  ${inputData[$i]};
    make_path( $tempOutputDir );

    # tee the command line outputs
    tee STDERR, '>>', "${tempOutputDir}/std_err.out";

    # Print the number of simulation runs
    my $ii = $i + 1;
    print "** Simulation type ${ii} of ${numSimulationTypes} **\n";

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
    

    # Loop through the number of simulations of each type
    for (my $j = 0; $j < $numSimulations; $j += 1 )
    {  
        # Print the number of simulation runs
        my $jj = $j + 1;
        print "Run ${jj} of ${numSimulations} \n";

        # Create individual simulation directory TODO: add system time and unique tag (increment from beginning)
        my $simulation_dir = ${tempOutputDir} . "/";
        make_path( $simulation_dir );
        my $gsd_path    = ${simulation_dir} . "/" . "data.gsd";

        # TODO: Create individual simulation directory and alter $gsd_path to reflect this
        # TODO: Then update simulation call to keep all files separate

        # Default parameters
        my $dt          = 1e-5;
        my $R_avg       = 10.0;
        my $phase_angle = 1.57079632679;
        my $epsilon     = 1e-4;

        # Modify default preferences for each simulation run
        switch($inputData[$i]) {
            case ( "varyDt" ) {    
                $dt = ${data[$j]};
            }
            case ( "varyRelDisp" ) {
                $R_avg = ${data[$j]};
            }
            case ( "varyPhaseAngle" ) {
                $phase_angle = ${data[$j]};
            }
            case ( "varyEpsilon" ) {
                $epsilon = ${data[$j]};
            }
        }

        # Generate GSD file
        system( "python " . ${pythonGSDCreation} . " --GSD-path=" . ${gsd_path} . " --dt=" . ${dt} . " --R_avg=" . ${R_avg} . " --phase-angle=" . ${phase_angle} . " --epsilon=" . ${epsilon} ) and die "Unable to generate GSD file: $?, $!";


        # ANCHOR: Run executable: [executable] [input gsd] [output directory]
        if ( ($jj < $numSimulations) and ($jj % int($numThreads / 2) != 0) and ($runSimulationSimulan) ) {

            system( "${buildDir}/src/./" . ${projectName} . " " . ${gsd_path} . " " . ${simulation_dir} . " &" ) 
                and die "Main project executable failed: $?, $!";
            sleep(5);  # brief pause before next simulation
        
        } else {
        
            system( "${buildDir}/src/./" . ${projectName} . " " . ${gsd_path} . " " . ${simulation_dir} ) 
                and die "Main project executable failed: $?, $!";
          
        }
    }

    # NOTE: Pause for all simulations to finish
    sleep(60);

# !SECTION
   

# SECTION: Run analysis scripts and clean-up

    # Run the analysis scripts
    make_path( "${tempOutputDir}/${analysisDir}" );
    system( "python3 ${pythonAnalysis} --relPath=${tempOutputDir} --outputDir=${analysisDir}" ) 
        and warn "Python analysis script failed: $!";

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
