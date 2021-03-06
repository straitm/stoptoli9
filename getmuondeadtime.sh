#!/bin/bash

# Finds the amount of deadtime in a run due to my muon cut. My cut for
# selecting a stopping muon (as distinct from the cut for selecting
# decay events *after* a stopping muon!) is for 0.5ms after each "muon",
# where a muon is defined as coinov || fido_qiv > 5000 || ctEvisID > 60.
# Since these 0.5ms windows can overlap AND muons can be correlated to
# each other, this is not as easy as just counting the number of events
# in a run and multiplying by the window size, nor as easy as a small
# correction on that. You have to step through the run and see what it
# really adds up to if you want an exact answer.

run=$1

runfile10=$(printf '/cp/s4/dchooz/cheetah/prod-08-05_p01_v2/reduced.Run%07d_Seq010.root' $run)
runfile11=$(printf '/cp/s4/dchooz/cheetah/seq11/reduced.Run%07d_Seq011.root' $run)

if [ -e $runfile10 ]; then
  runfile=$runfile10
elif [ -e $runfile11 ]; then
  runfile=$runfile11
else
  echo I don\'t have your run $run
  exit 1
fi

(root -l -b $runfile << EOF

data->SetScanField(0)
data->SetBranchStatus("*", 0)
data->SetBranchStatus("trgtime", 1)
data->SetBranchStatus("coinov", 1)
data->SetBranchStatus("fido_qiv", 1)
data->SetBranchStatus("ctEvisID", 1)
data->Scan("trgtime", "coinov || fido_qiv > 5000 || ctEvisID > 60", "colsize=20 col=.17f")
.q
EOF
) 2> /dev/null | grep -E '^\* +[0-9]' | awk \
 'BEGIN{
    window[500000] = 500000.;
    window[1000000]=1000000.;
    # We always veto the first windows-worth of a run in case there was
    # a muon right before it. In other words, we treat it as though
    # there is a muon at time=0.
    last = 0;
    vetoed[500000] = 0;
    vetoed[1000000] = 0;
  }
  {
    this = $4; # time of this muon
    since = this - last;

    # How long was the veto window for the *last* muon? It was either
    # a fullwindow if it has been longer than that, or the time between muons
    # if it has not been that long.
    #
    # This works great except for the last muon, which has no next muon
    # to come along and trigger us to find out how long its window was.
    # Without looking up the run length, we just have to assume the last
    # muon imposes its whole window. (It is no good looking at the time
    # of the last trigger, since triggers come in much less often than
    # once every window-length.)
    for(i in window){
      if(since > window[i]) vetoed[i] += window[i];
      else vetoed[i] += since;
      last = this;
    }
  }
  END{
    for(i in window)
      vetoed[i] += window[i]; # assume last muon has the whole window
    printf "%0.9f\t%0.9f\n", vetoed[500000]/1e9, vetoed[1000000]/1e9
  }'
