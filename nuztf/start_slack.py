import numpy as np
import os
from slack import RTMClient
from slack_bot import bot_access_token, access_token

submit_file = os.path.join(
    os.path.dirname(os.path.abspath(__file__)) + "/spawn_tmux_session.sh"
)

ampel_bot_user = ["UMNJK00CU", "DMBKJG00K"]

keywords = ["<@{0}>".format(ampel_bot_user), "LIGO", "banana", "GRB", "Fermi", "Ampel"]


def run_on_event(data):
    print(data)
    ts = data["ts"]
    channel = data["channel"]
    sid = int(float(ts) * 1.0e6)
    cmd = "bash {0} {1} {2} {3}".format(submit_file, sid, ts, channel)
    print(cmd)
    os.system(cmd)


@RTMClient.run_on(event="message")
def say_hello(**payload):
    data = payload["data"]
    if "user" in data.keys():
        try:
            if not np.logical_and(
                np.sum([x in data["text"] for x in keywords]) == 0,
                "DMBKJG00K" not in data["user"],
            ):
                run_on_event(data)
        except KeyError:
            pass


if __name__ == "__main__":
    print("Running master client!")
    rtm_client = RTMClient(token=bot_access_token)
    rtm_client.start()
