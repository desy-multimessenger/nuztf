"""
Shared configuration for Kowalski
"""

from penquins import Kowalski

from nuztf.credentials import load_credentials

fp_mapping = {"mag": "magpsf", "magerr": "sigmapsf"}

kwargs = {
    "protocol": "https",
    "host": "kowalski.caltech.edu",
    "port": 443,
    "verbose": False,
    "timeout": 300.0,
}


def get_kowalski():
    """
    Get Kowalski object
    :return: Kowalski object
    """
    kowalski_token = load_credentials("kowalski", token_based=True)
    return Kowalski(token=kowalski_token, **kwargs)
