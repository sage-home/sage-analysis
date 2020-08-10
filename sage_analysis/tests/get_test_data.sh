#!/bin/bash
cwd=$(pwd)
datadir=test_data/

# The bash way of figuring out the absolute path to this file
# (irrespective of cwd). parent_path should be $ROOT/tests
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# Create the directories to host all the data.
mkdir -p "$parent_path"/$datadir
if [[ $? != 0 ]]; then
    echo "Could not create directory $parent_path/$datadir"
    echo "Failed."
    echo "If the fix to this isn't obvious, please feel free to open an issue on our GitHub page."
    echo "https://github.com/sage-home/sage-analysis/issues/new"
    exit 1
fi

cd "$parent_path"/$datadir
if [[ $? != 0 ]]; then
    echo "Could not cd into directory $parent_path/$datadir"
    echo "Failed"
    exit 1
fi

# If there isn't a Mini-Millennium data, download it.
if [ ! -f correct-mini-millennium-output_z0.000_0 ]; then

    # To download the trees, we use either `wget` or `curl`. By default, we want to use `wget`.
    # However, if it isn't present, we will use `curl` with a few parameter flags.
    echo "First checking if either 'wget' or 'curl' are present in order to download trees."

    clear_alias=0
    command -v wget

    if [[ $? != 0 ]]; then
        echo "'wget' is not available. Checking if 'curl' can be used."
        command -v curl

        if [[ $? != 0 ]]; then
            echo "Neither 'wget' nor 'curl' are available to download the Mini-Millennium trees."
            echo "Please install one of these to download the trees."
            exit 1
        fi

        echo "Using 'curl' to download trees."

        # `curl` is available. Alias it to `wget` with some options.
        alias wget="curl -L -O -C -"

        # We will need to clear this alias up later.
        clear_alias=1
    else
        echo "'wget' is present. Using it!"
    fi

    # Get the data files.
    wget "https://www.dropbox.com/s/lud6k0m2zvcqmwj/mini-millennium-sage-correct-output.tar.gz?dl=0" -O "mini-millennium-sage-correct-output.tar"
    if [[ $? != 0 ]]; then
        echo "Could not download correct model output from the Manodeep Sinha's Dropbox...aborting."
        echo "Failed."
        echo "If the fix to this isn't obvious, please feel free to open an issue on our GitHub page."
        echo "https://github.com/sage-home/sage-analysis/issues/new"
        exit 1
    fi

    # If we used `curl`, remove the `wget` alias.
    if [[ $clear_alias == 1 ]]; then
        unalias wget
    fi

    tar -xvzf mini-millennium-sage-correct-output.tar
    if [[ $? != 0 ]]; then
        echo "Could not untar the correct model output...aborting."
        echo "Failed."
        echo "If the fix to this isn't obvious, please feel free to open an issue on our GitHub page."
        echo "https://github.com/sage-home/sage-analysis/issues/new"
        exit 1
    fi

fi

# The parameter file refers to SAGE output files as "test_sage" files. However, the data has been saved using 'correct-mini-millennium-output'.
tmpfile="$(mktemp)"
sed '/^FileNameGalaxies /s/.*$/FileNameGalaxies        correct-mini-millennium-output/' mini-millennium.par > ${tmpfile}
mv ${tmpfile} mini-millennium.par

# The parameter file refers to "./tests/test_data/" for the output and input files.  Instead, we change these to just be "./".
sed '/^OutputDir /s/.*$/OutputDir        .\//' mini-millennium.par > ${tmpfile}
mv ${tmpfile} mini-millennium.par

sed '/^SimulationDir /s/.*$/SimulationDir        .\//' mini-millennium.par > ${tmpfile}
mv ${tmpfile} mini-millennium.par

sed '/^FileWithSnapList /s/.*$/FileWithSnapList        .\/millennium.a_list/' mini-millennium.par > ${tmpfile}
mv ${tmpfile} mini-millennium.par

exit 0
