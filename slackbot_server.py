#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)

import datetime
import logging
import os

from flask import Flask  # type: ignore
from slack import WebClient  # type: ignore
from slackeventsapi import SlackEventAdapter  # type: ignore

from nuztf.utils import is_icecube_name, is_ligo_name
from slackbot import Slackbot

nuztf_slackbot = Flask(__name__)

slack_events_adapter = SlackEventAdapter(
    os.environ.get("SLACK_EVENTS_TOKEN"), "/slack/events", nuztf_slackbot
)
slack_web_client = WebClient(token=os.environ.get("SLACK_TOKEN"))


def scan(
    channel: str,
    ts: str,
    name: str,
    dl_results: bool,
    event_type: str,
    do_gcn: bool,
    time_window: int | None,
):
    """ """
    slack_bot = Slackbot(
        channel=channel,
        ts=ts,
        name=name,
        dl_results=dl_results,
        event_type=event_type,
        do_gcn=do_gcn,
        time_window=time_window,
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
    if is_icecube_name(name):
        return "nu"
    elif is_ligo_name(name):
        return "gw"
    else:
        return "invalid"


def get_help_message(user: str) -> list[dict]:
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

            time_window = None
            dl_results = True

            for i, parameter in enumerate(split_text):
                if parameter in fuzzy_parameters(["gcn", "GCN"]):
                    do_gcn = True
                elif parameter in fuzzy_parameters(["rerun", "nodl"]):
                    dl_results = False
                elif parameter in fuzzy_parameters(
                    ["window", "timewindow", "time-window"]
                ):
                    try:
                        time_window = int(split_text[i + 1])
                    except ValueError:
                        wc.chat_postMessage(
                            channel=channel_id,
                            text="Error: --window has to be an integer.",
                            thread_ts=ts,
                        )
                        return

            do_scan = True
            display_help = False
            name = split_text[1]
            event_type = parse_name(name)

            if event_type == "nu":
                message = (
                    f"Hi there; running a neutrino scan for *{name}*. One moment please"
                )
            elif event_type == "gw":
                message = f"Hi there; running a GW scan for *{name}*. One moment please"

            elif event_type == "invalid":
                message = f"Hi there; please enter either the name of a GW event (e.g. S190814bv) or a neutrino event (e.g. IC200620A)"
                do_scan = False

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
                    dl_results=dl_results,
                    event_type=event_type,
                    do_gcn=do_gcn,
                    time_window=time_window,
                )
            else:
                return


# for running directly with Flask (for debugging)
if __name__ == "__main__":
    nuztf_slackbot.run(port=4000, debug=True)
