#!/bin/bash

#EOC
#------------------------------------------------------------------------------
#              Harvard Biogeochemistry of Global Pollutants Group             !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: make_new_pcb_rundir
#
# !DESCRIPTION: Creates an additional pcb run directory in your local
#  directory that replicates the settings used by H. Amos for January 2016 
#  runs on Svante, but adapted for Odyssey. It should be used after you have
#  created a local pcb run directory called MITgcm_code since it creates
#  symbolic links to input data already downloaded to your regal directory.
#  This save time and memory, and is particularly useful if you would like
#  different run directories for different congeners so that you do not need
#  to change source code between runs.
#\\
#\\
# !CALLING SEQUENCE:
#  ./make_additional_hamos_pcb_rundir
#
#  AUTHOR: Lizzie Lundgren, 19 Dec 2016
#  HISTORY:
#       - Feb 13, 2018 - upadated location of control file directory (ccw)
#EOP
#------------------------------------------------------------------------------
#BOC

# Define your MITgcm repository directory
ROOT_LOCAL=$(pwd -P)

# Check that you are running this on regal. Exit with a message if you are not.
if [[ ${ROOT_LOCAL} != "/n/regal/"* ]]; then
    echo "You must run this script from a clone of the MITgcm repository stored on regal. Exiting."
    exit 1
fi

#===========================
# Define local directories
#===========================
VERIF_LOCAL=${ROOT_LOCAL}/verification
OCE_INPUTS=${VERIF_LOCAL}/global_oce_input_fields
PCB_LOCAL=${VERIF_LOCAL}/global_pcb_llc90
RUN_LOCAL=${PCB_LOCAL}/run
CODE_LOCAL=${PCB_LOCAL}/code
BUILD_LOCAL=${PCB_LOCAL}/build

#===========================
# Define remote directories
#===========================
SUNDERLAND=/n/sunderland_lab/Lab
FORCING_REMOTE=${SUNDERLAND}/eccov4_input/era-interim
HAMOS_SVANTE=${SUNDERLAND}/svante/helen/MITgcm_ECCOv4/
MASTER_REMOTE=${SUNDERLAND}/MITgcm_code
VERIF_REMOTE=${MASTER_REMOTE}/verification
PCB_REMOTE=${VERIF_REMOTE}/global_pcb_llc90
HG_REMOTE=${VERIF_REMOTE}/global_hg_llc90

#=========================================================================
# Accept user input on whether to create a new run directory from scratch.
# NOTE: It is much faster to choose to skip copying input files!!!
#=========================================================================
KEEP_INPUTS=0
if [[ -d ${ROOT_LOCAL}/verification ]]; then
   echo "You are choosing to create a completely new run directory from scratch."
   echo "If you have already created the directory, you can simply clean it to save time."
   echo "Do this by running ./clean_pcb_rundir.sh within this directory."
   echo "Do you still want to proceed with creating a new rundir from scratch? (y/n)"
   read answer
   if echo "$answer" | grep -iq "^y"; then
      echo "Creating new ./verification directory."
   else
      echo "Exiting."
      exit 0
   fi
fi

# Create your local 'verification' directory which will contain your input
# files and your run directory
mkdir -p verification

# Create symbolic links to ocean inputs
ln -s ${ROOT_LOCAL}/../MITgcm_code/verification/global_oce_cs32 ${VERIF_LOCAL}/global_oce_cs32
ln -s ${ROOT_LOCAL}/../MITgcm_code/verification/global_oce_input_fields ${VERIF_LOCAL}/global_oce_input_fields

# Create the global_pcb_llc90 directory
echo "Creating ${VERIF_REMOTE}/global_pcb_llc90 directory"
mkdir -p ${VERIF_LOCAL}/global_pcb_llc90

# Copy remote global_pcb_llc90 directories 
echo "Copying ${VERIF_REMOTE}/global_pcb_llc90 subdirectories:"
echo "    (1) code"
echo "    (2) build"
echo "    (3) run"
cp -rfL ${PCB_REMOTE}/run ${PCB_LOCAL}
cp -rfL ${PCB_REMOTE}/build ${PCB_LOCAL}
cp -rfL ${PCB_REMOTE}/code ${PCB_LOCAL}

# Copy the control subdirectory from the eccov4_input/controls_new run directory
echo "Creating ${RUN_LOCAL}/control"
mkdir ${RUN_LOCAL}/control
echo "Copying ${FORCING_REMOTE}/../controls_new/* to ${RUN_LOCAL}/control/"
cp -rfL ${FORCING_REMOTE}/../controls_new/* ${RUN_LOCAL}/control/

# Copy the data configuration files to the run directory
echo "Copying data config files from the run directory template"
cp -f ${PCB_REMOTE}/config/* ${RUN_LOCAL}

# Create symbolic links to forcing data stored in MITgcm_code run/forcing/
ln -s ${ROOT_LOCAL}/../MITgcm_code/verification/global_pcb_llc90/run/forcing/* ${RUN_LOCAL}/forcing

# Create symbolic links to inputs within the PCB run directory after
# removing any existing ones
echo "Creating local symbolic links in global_pcb_llc90/run"
cd ${RUN_LOCAL}/
ln -s ../../global_oce_cs32/input_fields/CORE2*    ${RUN_LOCAL}/
ln -s ../../global_oce_input_fields/ecmwf/*        ${RUN_LOCAL}/
ln -s ../../global_oce_input_fields/ecco_v4/*      ${RUN_LOCAL}/
ln -s ../../global_oce_input_fields/llc90/*        ${RUN_LOCAL}/
ln -s ../../global_oce_input_fields/llc90_fields/* ${RUN_LOCAL}/
cd ${ROOT_LOCAL}

# **** TEMPORARY TO REPLICATE HELEN'S CODE (DIFFERENT THAN REPO CODE) ****
# Replace some of the default repository directories with Helen's svante
echo "WARNING: Replacing eesupp and model directories with Helen Amos' "
echo "svante directories"
rm -r ${ROOT_LOCAL}/eesupp
rm -r ${ROOT_LOCAL}/model
cp -r ${HAMOS_SVANTE}/eesupp ${ROOT_LOCAL}
cp -r ${HAMOS_SVANTE}/model  ${ROOT_LOCAL}

#========================================================
# Copy utility and run scripts for setting up directory,
# compiling, and running
#========================================================
# Remove setup scripts that came with the repo but are not for hamos_pcb
rm -f ${ROOT_LOCAL}/make_rundir
rm -f ${ROOT_LOCAL}/setcpus

echo "Copying PCB utility scripts to top-level directory"
cp -f ${PCB_REMOTE}/util_scripts/* ${ROOT_LOCAL}

echo "Copying template run scripts to top level directory (update congener as needed)"
cp -f ${PCB_REMOTE}/run_scripts/* ${ROOT_LOCAL}

#==================================================
# Copy H. Amos' pickup files for easy access
#==================================================
cp -rf ${PCB_REMOTE}/hamos_pickups ${ROOT_LOCAL}

##==========================================================================
## Copy select files from local pkg/pcb so that the initial versions match
## NOTE: Edit these within code/ rather than in pkg/pcb/ for running.
## Put your updates in /pkg/pcb only if you want to commit changes as the 
## default in the repository. This code is currently commented out
## since H. Amos' files in the code directory are used instead.
##==========================================================================
#pcb_options=${ROOT_LOCAL}/pkg/pcb/PCB_OPTIONS.h
#pcb_size=${ROOT_LOCAL}/pkg/pcb/PCB_SIZE.h
#if [[ -f pcb_options ]]; then
#    echo "Copying local versions pkg/pcb/PCB_OPTIONS.h and "
#    echo "pkg/pcb/PCB_SIZE.h to global_pcb_llc90/code"
#    cp pcb_options ${CODE_LOCAL}
#    cp pcb_size ${CODE_LOCAL}
#fi

   # Replace the PCB_OPTIONS.h and PCB_SIZE.h with Helen's svante files
   echo "Replacing code/PCB_OPTIONS.h and code/PCB_SIZE.h with H. Amos' svante files"
   cp ${HAMOS_SVANTE}/verification/global_pcb_llc90/code/PCB_OPTIONS.h ${CODE_LOCAL}
   cp ${HAMOS_SVANTE}/verification/global_pcb_llc90/code/PCB_SIZE.h ${CODE_LOCAL}

#=====================================================================
# Replace pkg/ with bitbucket repository. Prompt the user to do this
# this manually after the automatic deletion of the copied directory.
# This forces the user to always know what code is being compiled/run.
#=====================================================================
echo " "
echo "WARNING: Deleting pkg directory. You must replace it with Helen Amos' code: "
rm -rf ${ROOT_LOCAL}/pkg
echo "  $ git clone -b master https://bitbucket.org/harvard-bgc/hamos_mitgcm_pkg pkg"
echo " "
echo "Make sure you checkout the version of pkg you want to run and change to the congener you want to simulate prior to compiling!"







