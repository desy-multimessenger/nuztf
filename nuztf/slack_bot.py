from slack import WebClient
import getpass
import numpy as np
import logging
import io
import os
from nuztf.gw_multi_process import MultiGwProcessor
import traceback
import time

slack_token = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), ".slack_access_token.txt"
)

try:
    with open(slack_token, "r") as f:
        access_token = f.read()
except FileNotFoundError:
    access_token = getpass.getpass(prompt="Slack Access Token: ", stream=None)
    with open(slack_token, "wb") as f:
        f.write(access_token.encode())

bot_token = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), ".slack_bot_access_token.txt"
)

try:
    with open(bot_token, "r") as f:
        bot_access_token = f.read()
except FileNotFoundError:
    bot_access_token = getpass.getpass(prompt="Slack Bot Access Token: ", stream=None)
    with open(bot_token, "wb") as f:
        f.write(bot_access_token.encode())


def upload_fig(fig, data, filename, channel_id, thread_ts):
    imgdata = io.BytesIO()
    fig.savefig(imgdata, format="png", dpi=600, transparent=True)
    imgdata.seek(0)
    wc = WebClient(token=bot_access_token)
    wc.files_upload(
        file=imgdata.getvalue(),
        filename=filename,
        channels=channel_id,
        thread_ts=thread_ts,
        icon_emoji=":ampel-mm:",
        text="<@{0}>, here's the file {1} I've uploaded for you!".format(
            data["user"], filename
        ),
    )
    # fig.close()


def run_on_event(thread_ts, channel_id):

    web_client = WebClient(token=access_token)

    payload = web_client.conversations_history(
        channel=channel_id,
        oldest=str(float(thread_ts) - 1),
        latest=str(float(thread_ts) + 1),
    )

    data = payload["messages"][0]

    user = data["user"]

    web_client.chat_postMessage(
        channel=channel_id,
        text="Hi <@{0}>! You are interested in Ampel multi-messenger stuff, right? Let me get right on that for you.".format(
            user
        ),
        thread_ts=thread_ts,
        icon_emoji=":ampel-mm:",
    )

    message = data["text"].replace(u"\xa0", u" ")

    split_message = message.split(" ")

    gw_name = None
    gw_file = None
    rev_no = None
    prob_threshold = 0.95
    n_days = None
    try:
        for x in split_message:
            if len(x) > 0:
                if x[0] in ["s", "S"]:
                    if np.sum([y.isdigit() for y in x[1:7]]) == 6:
                        gw_name = x
                elif ".fit" in x:
                    # print(gw_file)
                    # gw_file = x[1:-1]
                    gw_file = x[1:-1]
                elif "rev" in x:
                    rev_no = int(x.split("=")[1])
                elif "prob_threshold" in x:
                    prob_threshold = float(x.split("=")[1])
                elif "n_days" in x:
                    n_days = float(x.split("=")[1])

        message = ""

        if gw_name is not None:
            message = "You are interested in LIGO event {0}. ".format(gw_name)

        if gw_file is not None:
            if gw_name is not None:
                message = (
                    "You have also specified a fits file. The fits file will be used.  "
                )
            else:
                message = "You are interested in the following fits file: {0}. ".format(
                    gw_file
                )

        if message == "":
            message = "No file was specified. I will just assume that you want the most recent LIGO event. "

        if gw_file is None:

            if rev_no is not None:
                if gw_name is not None:
                    message += "You have specified revision number {0}. ".format(rev_no)
            else:
                message += "No revision number has been specified. I will just take the most recent revision for this event. "

        message += "The Skymap will be scanned up to {0}% of the probability. ".format(
            100.0 * prob_threshold
        )

        if n_days is None:
            message += "No time range has been specified. I will scan from merger time to now. "
        else:
            message += "I will scan for objects first detected between merger time and {0} days after merger. ".format(
                n_days
            )

    except:
        message = "Sorry <@{0}>, I have run into a parsing error with your message. All your bases are belong to us. \n {1}".format(
            data["channel"], split_message
        )

    web_client.chat_postMessage(
        channel=channel_id, text=message, thread_ts=thread_ts, icon_emoji=":ampel-mm:"
    )

    logger = logging.getLogger("quiet_logger")
    logger.setLevel(logging.INFO)

    try:
        gw = MultiGwProcessor(
            n_days=n_days,
            gw_name=gw_name,
            gw_file=gw_file,
            rev=rev_no,
            logger=logger,
            prob_threshold=prob_threshold,
        )
        web_client.chat_postMessage(
            channel=channel_id,
            text="Scanning method: {0} \n Effective sky number: {1}".format(
                gw.scan_method, gw.n_sky
            ),
            thread_ts=thread_ts,
            icon_emoji=":ampel-mm:",
        )
        fig = gw.plot_skymap()
        upload_fig(fig, data, "LIGO_skymap.png", channel_id, thread_ts)
        gw.clean_cache()
        gw.fill_queue()
        gw.terminate()
        gw.combine_cache()
        gw.clean_cache()
        wc = WebClient(token=bot_access_token)
        wc.files_upload(
            file=gw.output_path,
            filename=os.path.basename(gw.output_path),
            channels=channel_id,
            thread_ts=thread_ts,
            icon_emoji=":ampel-mm:",
        )
        fig, overlap = gw.plot_overlap_with_observations(first_det_window_days=n_days)
        web_client.chat_postMessage(
            channel=channel_id,
            text=overlap,
            thread_ts=thread_ts,
            icon_emoji=":ampel-mm:",
        )
        upload_fig(fig, data, "LIGO_overlap.png", channel_id, thread_ts)
        web_client.chat_postMessage(
            channel=channel_id,
            thread_ts=thread_ts,
            icon_emoji=":ampel-mm:",
            text="Here's a draft GCN:",
        )
        web_client.chat_postMessage(
            channel=channel_id,
            thread_ts=thread_ts,
            icon_emoji=":ampel-mm:",
            text=gw.draft_gcn(),
        )
        web_client.chat_postMessage(
            channel=channel_id,
            thread_ts=thread_ts,
            icon_emoji=":ampel-mm:",
            text="<@{0}>, I'm all finished. Go find that counterpart!".format(
                data["user"]
            ),
        )
    except Exception as e:
        web_client.chat_postMessage(
            channel=channel_id,
            text="Sorry <@{0}>, we need to talk. It's not you, it's me. I have run into an error, and cannot process your request further. I wish you the best of luck in all your future endeavours. \n\n `{1}`. ".format(
                data["user"], e
            ),
            thread_ts=thread_ts,
            icon_emoji=":ampel-mm:",
        )
        traceback.print_exc()
        time.sleep(120)
        raise

    # Session will not die until multi-processes have been terminated

    try:
        gw.terminate()
    except:
        pass
    print("Done!")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--timestamp", type=str)
    parser.add_argument("-c", "--channel", type=str)

    cfg = parser.parse_args()

    run_on_event(cfg.timestamp, cfg.channel)
