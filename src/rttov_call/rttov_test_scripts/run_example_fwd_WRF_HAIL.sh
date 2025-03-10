#!/bin/sh
# Script to run the example_fwd.F90 for rttov hail case study
# Vito Galligani
# RTTOV-14
#------------------------------------------------------------------------

BIN=$(perl -e 'for(@ARGV){m/^BIN=(\S+)$/o && print "$1";}' $*)
if [ "x$BIN" = "x" ]
then
  BIN=bin
fi

# Profile dir
# '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RTTOVout/InputAll'

#---------- CONFIGURE THE FOLLOWING INFO ------------#

# Input data
COEF_FILENAME='rtcoef_gcom-w_1_amsr2.dat'                           # Location of this file is set below in $COEF_DIR
PROF_FILENAME='atm_em1.0rttov14_WRF-WSM6_2018-11-10_20_AMSR2.dat'   # Input profile(s), usually found in $TEST_DIR 

# Insstrument info
NCHAN=14                           # Number of channels to simulate for each profile
CHAN_LIST=$(seq -s ' ' $NCHAN)     # Space-separated channel-list
ZENITH=45

# Experiment naming
EXP='WRF_WSM6_20181110_20_AMSR2'

#----------------------------------------------------#
# Main rttov dir
rttov_main_dir="/home/vito.galligani/Work/Studies/HAILCASE_10112018/rttov14.0_beta"

# Confis
DO_SOLAR=0                         #  0 = solar off / 1 = solar on
NTHREADS=1                         # Number of threads to use (compile RTTOV with OpenMP to exploit this)
CHECK_REF=0                        # Set to 0 to omit check against test reference
EMISSIVITY_VAL=1                   # Set to 1 later explore emissivity atlas?

# Path relative to the rttov_test directory:
TEST_DIR=./../../../../../../datosmunin3/Work/HAILCASE_10112018_datos/RTTOVout/InputAll/

# Output folder
out_folder=/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RTTOVout/${EXP}
mkdir -p "${out_folder}"
 
# Paths relative to the rttov_test/${TEST_DIR} directory:
BIN_DIR=$rttov_main_dir/$BIN                               # BIN directory (may be set with BIN= argument)
COEF_DIR=$rttov_main_dir/'rtcoef_rttov14/rttov13pred54L'

#----------------------------------------------------#
CWD=$(pwd)
cd $TEST_DIR

# Coefficient file
COEF_FILENAME=$COEF_DIR/$COEF_FILENAME
if [ ! -f $COEF_FILENAME ]; then
  echo "Coef file $COEF_FILENAME not found, aborting..."
  exit 1
fi

$BIN_DIR/example_fwd.exe << EOF
"${COEF_FILENAME}", Coefficient filename
"${PROF_FILENAME}",   Input profile filename
${DO_SOLAR}       ,   Turn solar radiation on/off
${NCHAN}          ,   Number of channels
${CHAN_LIST}      ,   Channel numbers
${NTHREADS}       ,   Number of threads
${EMISSIVITY_VAL} ,   Channel numbers
${ZENITH} ,   Zenith angle
EOF

if [ $? -ne 0 ]; then
  echo " "
  echo "TEST FAILED"
  echo " "
  exit 1
fi

#----------------------------------------------------#
# copy outputs a outdir
mv output_tb_45 ${out_folder}/.
rm surface2space_transm_45
rm output_transm_45


exit
