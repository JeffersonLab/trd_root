#!/bin/bash


#source setup_env.sh
# only root needed

RUNNUM=${1-none}
MAXEVT=${2-0}
FRSTEVT=${3-0}
FILENUM=${4-0}

if [[ ${RUNNUM} == "none" ]] ; then
    echo "================================="
    echo " Usage: ./$0 <RunNum> [Max_Events] [First_Event] [FileNum] "
    echo "================================="
    exit 0;
fi

echo "====>  Process RUN=$RUNNUM FILE=$FILENUM <=========="

root --web=off -l <<EOC
.L trdclass.C+
trdclass t(${RUNNUM},${MAXEVT},${FRSTEVT})
t.Loop()
EOC
