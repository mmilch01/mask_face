# Source this file to enable face masking script.

# prerequisites
#SCR=/nrgpackages/scripts
#source ${SCR}/fsl5_setup.sh
#source ${SCR}/nil-tools_setup.sh
#source ${SCR}/nrg-tools_setup.sh
#source ${SCR}/xnat-tools_setup.sh
#source ${SCR}/matlab_setup.sh

# custom scripts

RT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
R="$( cd $RT && cd ../ && pwd )"

BIN_HOME=$R/bin
SCRIPT_HOME=$R/bin
MASKFACE_MATLAB_ROOT=$R/bin/matlab

# update path
PATH=${BIN_HOME}:${SCRIPT_HOME}:${PATH}

# update FSL configuration
FSLOUTPUTTYPE=NIFTI_PAIR

MCR_HOME=/home/shared/NRG/mmilch/lib/MATLAB_Runtime.8.5
MASKFACE_MCR_HOME=/home/NRG/mmilch01/src/maskface/matlab_runtime/for_redistribution_files_only

# Finalize
export FSLOUTPUTTYPE PATH MASKFACE_MATLAB_ROOT MCR_HOME MASKFACE_MCR_HOME

