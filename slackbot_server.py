#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)

import os

from flask import Flask  # type: ignore
from slack import WebClient  # type: ignore
from slackeventsapi import SlackEventAdapter  # type: ignore

nuztf_slackbot = Flask(__name__)

slack_events_adapter = SlackEventAdapter(
    os.environ.get("SLACK_EVENTS_TOKEN"), "/slack/events", nuztf_slackbot
)
slack_web_client = WebClient(token=os.environ.get("SLACK_TOKEN"))
