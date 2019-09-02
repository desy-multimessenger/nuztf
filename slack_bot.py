from slack import RTMClient, WebClient
import getpass
import numpy as np
import logging
import matplotlib.pyplot as plt
import io
import os
from gw_scanner import GravWaveScanner
from gw_multi_process import MultiGwProcessor

try:
    with open(".slack_access_token.txt", "r") as f:
        access_token = f.read()
except FileNotFoundError:
    access_token = getpass.getpass(prompt='Slack Access Token: ', stream=None)
    with open(".slack_access_token.txt", "wb") as f:
        f.write(access_token.encode())


def upload_fig(fig, web_client, data, filename):
    imgdata = io.BytesIO()
    fig.savefig(imgdata, format='png', dpi=600, transparent=True)
    imgdata.seek(0)
    web_client.files_upload(
        file=imgdata.getvalue(),
        filename=filename,
        channels=data['channel'],
        thread_ts = data['ts'],
        icon_emoji=':ligo:',
        text="<@{0}>, here's the file {1} I've uploaded for you!".format(data["user"], filename)
    )
    #fig.close()

def run_on_event(data, web_client):
    channel_id = data['channel']
    thread_ts = data['ts']
    user = data['user']
    web_client.chat_postMessage(
        channel=channel_id,
        text=f"Hi <@{user}>! You are interested in LIGO stuff, right? Let me get right on that for you.",
        thread_ts=thread_ts,
        icon_emoji=':ligo:'
    )
    split_message = data['text'].split(" ")

    gw_name = None
    gw_file = None
    rev_no = None
    prob_threshold = 0.9
    try:
        for x in split_message:
            if len(x) > 0:
                if x[0] in ["s", "S"]:
                    if np.sum([y.isdigit() for y in x[1:7]]) == 6:
                        gw_name = x
                elif ".fits" in x:
                    print(gw_file)
                    gw_file = x[1:-1]
                elif "rev" in x:
                    rev_no = int(x.split("=")[1])
                elif "prob_threshold" in x:
                    prob_threshold = float(x.split("=")[1])     



        message = ""

        if gw_name is not None:
            message = "You are interested in LIGO event {0}. ".format(gw_name)

        if gw_file is not None:
            if gw_name is not None:
                message = "You have also specified a fits file. The fits file will be used.  "
            else:
                message = "You are interested in the following fits fille: {0}. ".format(gw_file)

        if message == "":
            message = "No file was specified. I will just assume that you want the most recent LIGO event. "

        if rev_no is not None:
            if gw_name is not None:
                message += "You have specified revision number {0}. ".format(rev_no)
        else:
            message += "No revision number has been specified. I will just take the most recent revision for this event. "
    
        message += "The LIGO Skymap will be scanned up to {0}% of the probability.".format(100. * prob_threshold)

    except:
        message = "Sorry <@{0}>, I have run into a parsing error with your message. All your bases are belong to us. \n {1}".format(data["channel"], split_message)

    web_client.chat_postMessage(
        channel=channel_id,
        text=message,
        thread_ts=thread_ts,
        icon_emoji=':ligo:'
    )

    logger = logging.getLogger("quiet_logger")
    logger.setLevel(logging.ERROR)

    try:
        gw = MultiGwProcessor(gw_name=gw_name, gw_file=gw_file, rev=rev_no, logger=logger)
        web_client.chat_postMessage(
            channel=channel_id,
            text="Scanning method: {0} \n Effective sky number: {1}".format(gw.scan_method, gw.n_sky),
            thread_ts=thread_ts,
            icon_emoji=':ligo:'
        )
        fig = gw.plot_skymap()
        upload_fig(fig, web_client, data, "LIGO_skymap.png")
        gw.clean_cache()
        gw.fill_queue()
        gw.terminate()
        gw.combine_cache()
        gw.clean_cache()
        web_client.files_upload(
            file=gw.output_path,
            filename=os.path.basename(gw.output_path),
            channels=data['channel'],
            thread_ts=data['ts'],
            icon_emoji=':ligo:'
        )
        fig, overlap = gw.plot_overlap_with_observations()
        web_client.chat_postMessage(
            channel=channel_id,
            text=overlap,
            thread_ts=thread_ts,
            icon_emoji=':ligo:'
        )
        upload_fig(fig, web_client, data, "LIGO_overlap.png")
        web_client.chat_postMessage(
            channel=data['channel'],
            thread_ts=data['ts'],
            icon_emoji=':ligo:',
            text="<@{0}>, I'm all finished. Go find that kilonova!".format(data["user"])
        )
    except Exception as e:
        web_client.chat_postMessage(
            channel=channel_id,
            text="Sorry <@{0}>, we need to talk. It's not you, it's me. I have run into an error, and cannot process your request further. I wish you the best of luck in all your future endeavours. \n\n `{1}`. ".format(data["user"], e),
            thread_ts=thread_ts,
            icon_emoji=':ligo:'
        )

ampel_bot_user = ["UMNJK00CU", "DMBKJG00K"]

keywords = ["<@{0}>".format(ampel_bot_user), "LIGO", "banana"]

@RTMClient.run_on(event="message")
def say_hello(**payload):
    data = payload['data']
    web_client = payload['web_client']
    if "user" in data.keys():
        try:
            if not np.logical_and(np.sum([x in data['text'] for x in keywords]) == 0, "DMBKJG00K" not in data["user"]):
                run_on_event(data, web_client)
        except KeyError:
            pass 
# slack_token = os.environ["SLACK_API_TOKEN"]
rtm_client = RTMClient(token=access_token)
rtm_client.start()
