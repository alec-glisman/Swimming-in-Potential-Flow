#!/usr/bin/env perl

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


# TODO: Write this script
# 1) specify variables
# 2) generate GSD
# 3) generate output data files
# 4) move GSD to output location
# 5) call simulation (configure and build first)
# 6) analyze trajectory file





# SECTION: Input variables that user must specify before running script

# Input variables
my $inputDir         = "input/collinear_constrained";
my $inputPreferences = "defaultPref.in";
my @inputData        = (  );  # TODO
my $numSimulationTypes = scalar @inputData;
my $runSimulationSimulan = 1; # 0 only runs one simulation at a time

# Output variables
my $curDate          = strftime('%Y-%m-%d', localtime);
my $analysisDir      = "analysis";

# Get the name of computer script is running on
my $host = `uname -n`;
chomp(my $home = `echo ~`);

# REVIEW: Python variables, check analysis script is correct for simulations run
my $pythonAnalysis   = "src/python_scripts/analysis/continuousForcedOscillation.py";

# Compiler variables
my $compiler = "";
my $numThreads = "";
if ( (index($host, "MacBook-Pro") != -1) or (index($host, "MBP") != -1) ) {
    $compiler    = "macOS"; 
    $numThreads = "16";
} elsif ( index($host, "stokes") != -1 ) {
    $compiler    = "stokes"; 
    $numThreads  = "16";
} elsif (( index($host, "Alec-Glisman-PC-Ubuntu") != -1 ) or ( index($host, "Alec-Glisman-PC-Windows") != -1 ) or ( index($host, "s") != -1 )) {
    $compiler    = "PC-Ubuntu";
    $numThreads  = "24";
} elsif (( index($host, "shear") != -1 ) or ( index($host, "s") != -1 )) {
    $compiler    = "shear";
    $numThreads  = "";
} else {
    $compiler    = "Docker"; 
    $numThreads  = "";
}
my $build         = "Release";                # OPTIONS: Release, Debug
my $enableTesting = "True";                   # OPTIONS: (False) OFF, (True) ON
my $buildDir      = "build";                  # Title whatever you want build folder to be
my $generator     = "Unix Makefiles";         # ONLY TESTED WITH UNIX
my $cwd           = cwd();

# !SECTION  


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

# !SECTION


# SECTION: Loop through simulation types

for (my $i = 0; $i < $numSimulationTypes; $i += 1 )
{

    # Current datetime
    my $curDateTime = strftime('%Y-%m-%d.%H-%M-%S', localtime);

    # Make temp output directory
    my $tempOutputDir = "temp_" . $curDateTime . "_collinear_" .  ${inputData[$i]};
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


        # SECTION: Modify default preferences for each simulation run

        # Make a copy of 'defaultPref.in' in current tempOutputDir
        my $uniquePrefs = "${tempOutputDir}/pref.in";
        copy("${inputDir}/defaultPref.in", $uniquePrefs) 
            or die "Could not copy default preferences: $!";

        # Modify input data as needed
        my $file = path( $uniquePrefs );
        my $dataChange = path( $file )->slurp_utf8;
        my $oldKey;
        my $newKey;

        # REVIEW: If input preferences changes, these regex strings must change
        switch($inputData[$i]) {

            case ( "varyDt" ) {
            
                $oldKey = '1e-6               \| ND: finite timestep';
                $newKey = ${data[$j]} . '               | ND: finite timestep';
            }

            case ( "varyEpsilon" ) {

                # Recalculate relDispEqbm so that max separation (D) is 10.0
                $oldKey = '10.0               \| ND: particle pair eqbm separation';
                $newKey = '4.0                | ND: particle pair eqbm separation';
                $dataChange =~ s/$oldKey/$newKey/g;
                $file->spew_utf8( $dataChange );

                # Recalculate the max velocity to for the given epsilon
                my $U0 = $data[$j] * $defaultOmegaCon * $defaultR0;
                my $U0Str = sprintf( "%1.6e", $U0 );
                $oldKey = '0.001              \| D \[length / time\]: particle pair max velocity';
                $newKey = $U0Str . '            | D [length / time]: particle pair max velocity';
            }

        	case ( "varyVarEpsilon" ) {  # Golestanian (2004, PRE) Swimmer analog

                # Run through both "fictitious" steps detailed in Eq. (9)
		        for ( my $k = 0; $k < 2; $k++ ) {

                    # Make a copy of 'defaultPref.in' in current tempOutputDir
                    copy("${inputDir}/defaultPref.in", $uniquePrefs) 
                        or die "Could not copy default preferences: $!";
                    # Modify input data as needed
                    $file = path( $uniquePrefs );
                    $dataChange = path( $file )->slurp_utf8;
             
                    # Make Golestanian Swimmer flag true
                    $oldKey = 'false               \| Golestanian Swimmer';
                    $newKey = 'true                | Golestanian Swimmer';
                    $dataChange =~ s/$oldKey/$newKey/g;
                    $file->spew_utf8( $dataChange );

                    # Make continuously deforming flag false
                    $oldKey = 'true                \| Continuously deforming swimmer';
                    $newKey = 'false               | Continuously deforming swimmer';
                    $dataChange =~ s/$oldKey/$newKey/g;
                    $file->spew_utf8( $dataChange );

                    my $new_delta_len;
                    my $delta_len_str;

                    # First step: \Delta_f (D) 
                    if ($k == 0) {

                        # Nothing to change, default compression value is false

  
                    # Second step: \Delta_f (D - \epsilon)
                    } else {
                        
                        # Change boolean to have the swimmer compressed on 2--3 side
                        $oldKey = 'false               \| If Golestanian Swimmer is compressed';
                        $newKey = 'true                | If Golestanian Swimmer is compressed';
                    }

                    $dataChange =~ s/$oldKey/$newKey/g;
                    $file->spew_utf8( $dataChange );

                    # Recalculate U0 to that we oscillate correct magnitude in position space
                    my $varEpsilon = ${data[$j]};
                    my $U0 = $varEpsilon * $defaultOmegaCon * 0.50;
                    my $U0Str = sprintf( "%1.4e", $U0 );
                    $oldKey = '0.001              \| D \[length / time\]: particle pair max velocity';
                    $newKey = $U0Str . '        | D [length / time]: particle pair max velocity';
                    $dataChange =~ s/$oldKey/$newKey/g;
                    $file->spew_utf8( $dataChange );
                    
                    # Recalculate relDispEqbm so that max separation (D) is 10.0
                    my $D = 10.0;
                    my $newRelDispEqbm = $D - ( 0.50 * $varEpsilon );
                    my $newRelDispEqbmStr = sprintf( "%1.4e", $newRelDispEqbm );
                    $oldKey = '10.0               \| ND: particle pair eqbm separation';
                    $newKey = $newRelDispEqbmStr . '              | ND: particle pair eqbm separation';

                    # ANCHOR: Run executable: [executable] [input preferences] [output directory]
                    if ( ($jj < $numSimulations) and ($jj % int($numThreads / 2) != 0) and ($runSimulationSimulan) ) {
                        system( "\"${buildDir}/src/./potential_swimmer_dynamics\" ${uniquePrefs} ${tempOutputDir}/ &" ) 
                            and die "Main project executable failed: $?, $!";
                        sleep(5);  # brief pause before next simulation
                    } else {
                        system( "\"${buildDir}/src/./potential_swimmer_dynamics\" ${uniquePrefs} ${tempOutputDir}/ " ) 
                            and die "Main project executable failed: $?, $!";            
                    }
                }  


                # Make a copy of 'defaultPref.in' in current tempOutputDir
                copy("${inputDir}/defaultPref.in", $uniquePrefs) 
                    or die "Could not copy default preferences: $!";
                # Modify input data as needed
                $file = path( $uniquePrefs );
                $dataChange = path( $file )->slurp_utf8;

                # Make "continuous" swimmer analog
                # Recalculate U0 to that we oscillate correct magnitude in position space
                my $varEpsilon = ${data[$j]};
                my $U0 = $varEpsilon * $defaultOmegaCon * 0.50;
                my $U0Str = sprintf( "%1.4e", $U0 );
                $oldKey = '0.001              \| D \[length / time\]: particle pair max velocity';
                $newKey = $U0Str . '         | D [length / time]: particle pair max velocity';
                $dataChange =~ s/$oldKey/$newKey/g;
                $file->spew_utf8( $dataChange );
                
                # Recalculate relDispEqbm so that max separation (D) is 10.0
                my $D = 10.0;
                my $newRelDispEqbm = $D - ( 0.50 * $varEpsilon );
                my $newRelDispEqbmStr = sprintf( "%1.4e", $newRelDispEqbm );
                $oldKey = '10.0               \| ND: particle pair eqbm separation';
                $newKey = $newRelDispEqbmStr . '              | ND: particle pair eqbm separation';

			}

            case ( "varyPhaseAngle" ) {

                $oldKey = '1.57079632679      \| ND \[radians\]: phase angle between oscillators';
                $newKey = ${data[$j]} . '                | ND [radians]: phase angle between oscillators';
            }

            case ( "varyRelDisp" ) {

                $oldKey = '10.0               \| ND: particle pair eqbm separation';
                $newKey = ${data[$j]} . '              | ND: particle pair eqbm separation';
            }

        }
        
        $dataChange =~ s/$oldKey/$newKey/g;
        $file->spew_utf8( $dataChange );

        # !SECTION


        # ANCHOR: Run executable: [executable] [input preferences] [output directory]
        if ( ($jj < $numSimulations) and ($jj % int($numThreads / 2) != 0) and ($runSimulationSimulan) ) {
            system( "\"${buildDir}/src/./potential_swimmer_dynamics\" ${uniquePrefs} ${tempOutputDir}/ &" ) 
                and die "Main project executable failed: $?, $!";
            sleep(5);  # brief pause before next simulation
        } else {
            system( "\"${buildDir}/src/./potential_swimmer_dynamics\" ${uniquePrefs} ${tempOutputDir}/ " ) 
                and die "Main project executable failed: $?, $!";            
        }
    }

# !SECTION
   
    # NOTE: Pause for python plots to finish (individual)
    if ( $compiler eq "Server" ) {
        sleep(40 * 60);
    } elsif ( $compiler eq "Docker" ) {
        sleep(20 * 60);
    } else {
        sleep(2 * 60);
    } 

    # SECTION: Run analysis scripts and clean-up

    # Change to other analysis scripts as needed
    if ($inputData[$i] eq "varyVarEpsilon" ) {
        $pythonAnalysis = "src/python_scripts/analysis/GolestanianOscillation.py";
    } else {
        $pythonAnalysis   = "src/python_scripts/analysis/continuousForcedOscillation.py";

    }

    # Run the analysis scripts
    make_path( "${tempOutputDir}/${analysisDir}" );
    system( "python3 ${pythonAnalysis} --relPath=${tempOutputDir}/ --outputDir=${analysisDir}" ) 
        and warn "Python analysis script failed: $!";

    # Move all output into the "data" directory
    my $outputDir = "data/${curDateTime}" . "_collinear_" .  ${inputData[$i]};
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
