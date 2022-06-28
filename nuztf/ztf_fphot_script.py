# ztf_ForcedPhotometryRequestExample.py Updated: 2022-06-02
# Written for Python 3.9.7, likely compatible with other python3 versions

# requests is used to query the Forced Photometry Service
import requests

# bs4 (Beauiful Soup 4) is used to parse the html output of request's methods
from bs4 import BeautifulSoup

# Python packages
from time import sleep

# Example to place a single request
# Easily iterable with for/while loops, though the service limits users to 100 requests at a time

# URL used to place a request
SendReqURL = "https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi"

# IPAC Auth (Don't change these)
IUN = "ztffps"
IPW = "dontgocrazy!"

# Request Parameters to send (USER DETERMINED)
SendRA = "280.8058788"
SendDec = "45.2077645"
SendJDStart = "2458231.891227"
SendJDEnd = "2458345.025359"
Email = "xxxx"
UserPass = "xxxx"

# Send request and return unformatted HTML output of get() method
Request = requests.get(
    SendReqURL,
    auth=(IUN, IPW),
    params={
        "ra": SendRA,
        "dec": SendDec,
        "jdstart": SendJDStart,
        "jdend": SendJDEnd,
        "email": Email,
        "userpass": UserPass,
    },
)

# Formatted (Parseable) HTML of request
ReqSoup = BeautifulSoup(Request.text, "html.parser")

# For one pending job (ie: one request), check forced photometry job tables

# Grab the table section of the HTML
ReqTable = ReqSoup.find("table")
# Limit to rows in said table
for row in ReqTable.find_all("tr"):
    # Find the items in the row (omit header cells)
    cols = row.find_all("td")
    if len(cols) > 0:
        # Find RA, Dec, and JD values of request as recorded by the service
        # Use a nested for loop here if submitting multiple requests, append values to list/arr/etc
        ReqRA = float(cols[0].text.strip())
        ReqDec = float(cols[1].text.strip())
        ReqJDS = float(cols[2].text.strip())
        ReqJDE = float(cols[3].text.strip())

# Check if request is fulfilled using a request status check and parse with beautifulsoup4
# Note: Requests older than 30 days will not show up here
# Iterable, as beautifulsoup can parse html columns and output lists

TrueIfPending = True
slp = True

# URL for checking the status of jobs sent by user
# Note: This webpage only updates once an hour, on the hour
StatusURL = "https://ztfweb.ipac.caltech.edu/cgi-bin/getForcedPhotometryRequests.cgi"

# Open loop to periodically check the Forced Photometry job status page
while TrueIfPending:

    # Query service for jobs sent in past 30 days
    OutputStatus = requests.get(
        StatusURL,
        auth=(IUN, IPW),
        params={
            "email": Email,
            "userpass": UserPass,
            "option": "All recent jobs",
            "action": "Query Database",
        },
    )

    # Check if job has ended and lightcurve file was created
    # Note: If an exotic error occurs and the ended field is not populated, this will go on forever
    # Table has 11 cols, reqid=0, ended=7, lc=10

    # Format HTML
    OutputSoup = BeautifulSoup(OutputStatus.text, "html.parser")
    # Get Job information in HTML table
    OutputTable = OutputSoup.find("table")
    OutputEnded = ""
    # Verify Table exists (If no requests have been sent in past 30 days or if service is updating, this will catch it)
    if OutputTable != "" and OutputTable is not None:
        # Parse Table rows
        for row in OutputTable.find_all("tr"):
            # Parse Table entries
            cols = row.find_all("td")
            if len(cols) > 0:
                # Check if values contained in a given row coorespond to the current request
                # Use a nested for loop here to check all requests submitted if there are more than one pending
                OutputRA = float(cols[1].text.strip())
                OutputDec = float(cols[2].text.strip())
                OutputJDS = float(cols[3].text.strip())
                OutputJDE = float(cols[4].text.strip())
                OutputEnded = cols[7].text.strip()
                # Check if job is finished (OutputEnded)
                # Check for equality between recorded request params and what is known to be in the Job Status table, accounting for rounding
                if (
                    OutputEnded != ""
                    and abs(OutputRA - ReqRA) <= 0.00001
                    and abs(OutputDec - ReqDec) <= 0.00001
                    and abs(OutputJDS - ReqJDS) <= 0.00001
                    and abs(OutputJDE - ReqJDE) <= 0.00001
                ):
                    # Get end time of job and lightcurve path
                    OutputLC = cols[10].text.strip()
                    OutputEnded = cols[7].text.strip()
                    # Set open loop hooks
                    TrueIfPending = False
                    slp = False
    # Pace the script so it doesn't spam the service
    if slp:
        sleep(30)

# If Job Status table values are not null, set the path for file to be downloaded to the path recorded in Job Status Table
if OutputEnded != "":
    if OutputLC != "":
        ReqPath = OutputLC
    else:
        # Catch errors in processing and alert user
        print("Job is done but no lightcurve was produced, quitting")
        quit()

# Set link to download lightcurve from
ForcedPhotometryReqURL = f"https://ztfweb.ipac.caltech.edu{ReqPath}"

# Set local path and filename for lightcurve (currently same as what is assigned by the service)
LocalFilename = ForcedPhotometryReqURL.split("/")[-1]

# Download the file with a get request
with requests.get(ForcedPhotometryReqURL, stream=True, auth=(IUN, IPW)) as r:
    # Write to the local file in chunks in case file is large
    with open(LocalFilename, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
