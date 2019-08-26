#!/bin/bash
export TMUX=""
export SESSION=$1
echo Spawning new session with following command
echo "tmux new-session -s $SESSION bash /cvmfs/icecube.opensciencegrid.org/users/steinrob/combo-realtime/build/skymap_scanner/resources/scripts/start_env.sh $2"
tmux new-session -d -s $SESSION "bash /cvmfs/icecube.opensciencegrid.org/users/steinrob/combo-realtime/build/skymap_scanner/resources/scripts/start_env.sh $2"
