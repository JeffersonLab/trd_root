#!/bin/bash

source setup_env.sh

RUNNUM=${1-none}
MAXEVT=${2-0}

if [[ ${RUNNUM} == "none" ]] ; then
    echo "================================="
    echo " Usage: ./$0 <RunNum> [Max_Events] "
    echo "================================="
    exit 0;
fi

echo "====>  Process RUN=$RUNNUM <=========="

root --web=off -l <<EOC
.L trdclass.C+
trdclass t(${RUNNUM},${MAXEVT})
t.Loop()
EOC
