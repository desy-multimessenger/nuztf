#!/bin/bash
export TMUX=""
export PYTHONPATH=`which python`
echo $PYTHONPATH
export CODEDIR="$(cd "$(dirname "$4")"; pwd -P)/$(basename "$4")"
echo $CODEDIR

export CMD="tmux new-session -d -s $1 $PYTHONPATH $CODEDIR/slack_bot.py -t $2 -c $3"
echo $CMD
tmux new-session -d -s $1 "$PYTHONPATH $CODEDIR/slack_bot.py -t $2 -c $3"
