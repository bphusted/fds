#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_101CR.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_102.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_102R1.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_102R.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_103.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_104.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_105.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_106.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_107.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_108.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_109.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_110.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_111.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_112A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_113A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_115A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_126B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_126BR1.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_128B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_202.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_203.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_205.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_207A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_208A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_210.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_212.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_214A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_215A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_216A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_217A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_218B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_223.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_226.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_227A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_229.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_230.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_231.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_233.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_235.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_236.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_238A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_239.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_244B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_245B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_246B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_247B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_248B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_249B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_250B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_251B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_252B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_301A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_302A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_303A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_305A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_306A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_309A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_312A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_313A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_314.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_315A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_316.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_317A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_318A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_319A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_320A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_321A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_338B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_339B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_340B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_341B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_342B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_343B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_344B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_345B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_346B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_401A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_403A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_404A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_407B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_408B.fds

$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Cold_Flow_Series_1.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_501.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_502.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_605.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_606A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_607.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_608.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_610.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_611.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_612B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_615B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_617A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_618A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_621A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_622B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_623B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_624B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_625B.fds

echo FDS cases submitted