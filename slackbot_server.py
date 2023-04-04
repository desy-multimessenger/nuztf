#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)

import datetime
import logging
import os

from flask import Flask  # type: ignore
from slack import WebClient  # type: ignore
from slackbot import Slackbot
from slackeventsapi import SlackEventAdapter  # type: ignore

nuztf_slackbot = Flask(__name__)

slack_events_adapter = SlackEventAdapter(
    os.environ.get("SLACK_EVENTS_TOKEN"), "/slack/events", nuztf_slackbot
)
slack_web_client = WebClient(token=os.environ.get("SLACK_TOKEN"))


def scan(channel, ts, name: str, event_type: str, do_gcn: bool):
    """ """
    slack_bot = Slackbot(
        channel=channel, name=name, event_type=event_type, do_gcn=do_gcn
    )

    if do_gcn:
        slack_web_client.chat_postMessage(
            channel=channel,
            text=slack_bot.gcn,
            thread_ts=ts,
        )


def fuzzy_parameters(param_list) -> list:
    """ """
    fuzzy_parameters = []
    for param in param_list:
        for character in ["", "-", "--", "–"]:
            fuzzy_parameters.append(f"{character}{param}")
    return fuzzy_parameters


def parse_name(name: str) -> str:
    """ """
    if name[0:2] == "IC":
        return "nu"
    else:
        raise ValueError()


def get_help_message(user: str) -> str:
    """
    Get the help message to display all commands for the user
    """
    blocks = [
        {
            "type": "section",
            "text": {
                "type": "mrkdwn",
                "text": f"Hi <@{user}>. This is a bot for skymap scanning for neutrino/GW/GRB events with AMPEL. Just type *Neutrino IceCube-Name*, or *GW GW-Event*",
            },
        }
    ]

    return blocks


ts_old = []


@slack_events_adapter.on("message")
def message(payload):
    """ """
    event = payload.get("event", {})
    text = event.get("text")
    user = event.get("user")
    ts = event.get("ts")
    if ts not in ts_old:
        ts_old.append(ts)

        text = text.replace("*", "")
        split_text = text.split()
        logging.info(split_text)

        do_scan = False
        do_gcn = False

        if len(split_text) == 0:
            return

        elif split_text[0] in ["Scan", "SCAN", "scan"]:
            channel_id = event.get("channel")

            if len(split_text) == 1:
                blocks = get_help_message(user)
                slack_web_client.chat_postMessage(
                    channel=channel_id,
                    text=" ",
                    blocks=blocks,
                    thread_ts=ts,
                )
                return

            for i, parameter in enumerate(split_text):
                if parameter in fuzzy_parameters(["gcn", "GCN"]):
                    do_gcn = True

            do_scan = True
            display_help = False
            name = split_text[1]
            event_type = parse_name(name)

            if event_type == "nu":
                message = (
                    f"Hi there; running neutrino scan for *{name}*. One moment please."
                )

            slack_web_client.chat_postMessage(
                channel=channel_id,
                text=message,
                thread_ts=ts,
            )

            if do_scan:
                scan(
                    channel=channel_id,
                    ts=ts,
                    name=name,
                    event_type=event_type,
                    do_gcn=do_gcn,
                )


# for running directly with Flask (for debugging)
if __name__ == "__main__":
    nuztf_slackbot.run(port=4000, debug=True)
