#!/bin/bash


#source setup_env.sh
source /gapps/root/Linux_RHEL7-x86_64-gcc4.8.2/root-6.18.00/bin/thisroot.csh
# only root needed

RUNNUM=${1-none}
MAXEVT=${2-0}
FRSTEVT=${3-0}

if [[ ${RUNNUM} == "none" ]] ; then
    echo "================================="
    echo " Usage: ./$0 <RunNum> [Max_Events] [First_Event]"
    echo "================================="
    exit 0;
fi

echo "====>  Process RUN=$RUNNUM <=========="

root --web=off -l <<EOC
.L trdclass.C+g
trdclass t(${RUNNUM},${MAXEVT},${FRSTEVT})
t.Loop()
EOC
