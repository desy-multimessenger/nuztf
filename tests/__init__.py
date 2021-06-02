import os
from ztfquery import io
from ampel.log.AmpelLogger import AmpelLogger


logger = AmpelLogger()

try:
    io.set_account("ampel_api",
                   username=os.environ["AMPEL_API_USER"],
                   password=os.environ["AMPEL_API_PASSWORD"])
    logger.info('Set up AMPEL credentials')

except KeyError:
    logger.info('No Credentials for AMPEL API found in environment'
                'Assume they are set.')

if not os.path.exists(io.LOCALSOURCE):
    os.mkdir(io.LOCALSOURCE)