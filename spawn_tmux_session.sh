#!/bin/bash
export TMUX=""
export PYTHONPATH=`which python`
echo $PYTHONPATH
export CODEDIR="$(cd "$(dirname "$3")"; pwd -P)/$(basename "$3")"
echo $CODEDIR

export CMD="tmux new-session $PYTHONPATH $CODEDIR/slack_bot.py -t $1 -c $2"
echo $CMD
#tmux new-session $PYTHONPATH $CODEDIR/slack_bot.py -t $1 -c $2
#echo Spawning new session with following command
#echo "tmux new-session -s $SESSION bash /cvmfs/icecube.opensciencegrid.org/users/steinrob/combo-realtime/build/skymap_scanner/resources/scripts/start_env.sh $2"
#tmux new-session -d -s $SESSION "bash /cvmfs/icecube.opensciencegrid.org/users/steinrob/combo-realtime/build/skymap_scanner/resources/scripts/start_env.sh $2"
