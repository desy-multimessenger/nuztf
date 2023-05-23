import re

import astropy.units as u
import requests
from astropy.coordinates import Angle
from astropy.time import Time
from bs4 import BeautifulSoup

from nuztf.ampel_api import ampel_api_name


def check_mpc_coord(
    ra_deg: float,
    dec_deg: float,
    jd: float,
    rad_arcmin: float = 1.0,
    obs_code: str = "I41",
    lim_mag=24.0,
) -> bool:
    """
    Function to search for matches in the MPC database

    :param ra_deg: Right ascension (deg)
    :param dec_deg: Declination (deg)
    :param jd: Julian date
    :param rad_arcmin: Search radius (arcmin)
    :param obs_code: Observatory code (ZTF: I41)
    :param lim_mag: Limiting magnitude
    :return:
    """

    time = Time(jd, format="jd")
    year = time.strftime("%Y")
    month = time.strftime("%m")
    hour = int(time.strftime("%H"))
    min = int(time.strftime("%M"))
    sec = float(time.strftime("%S.%f"))
    day = f"{int(time.strftime('%d'))+(hour+min/60+sec/3600)/24:.2f}"
    ra = Angle(ra_deg, unit="degree")
    ra_hms = ra.to_string(unit="hourangle", sep=" ", precision=2, pad=True).split(" ")
    dec = Angle(dec_deg, unit="degree")
    dec_dms = dec.to_string(unit=u.deg, sep=" ", precision=1, pad=True).split(" ")
    url = (
        f"https://minorplanetcenter.net/cgi-bin/mpcheck.cgi?"
        f"year={year}&month={month}&day={day}&which="
        f"pos&ra={ra_hms[0]}+{ra_hms[1]}+{ra_hms[2]}&"
        f"decl={dec_dms[0]}+{dec_dms[1]}+{dec_dms[2]}&"
        f"TextArea=&radius={rad_arcmin}&oc={obs_code}&limit={lim_mag}&"
        f"sort=d&mot=h&tmot=s&pdes=u&needed=f&ps=n&type=p"
    )

    session = requests.Session()
    response = session.get(url)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, "html.parser")

    match = False

    text = soup.find("pre").text

    print(len(re.findall("The following objects,.*", response.text)))
    #     print([0])
    #     print(text)
    # except Exception:
    #     try:
    #         print(re.findall("No known minor planets,.*", response.text)[0])
    #     except Exception:
    #         print(response.text)

    raise

    return match


def check_mpc_alert(alert, **kwargs):
    print(alert)

    return check_mpc_coord(
        ra_deg=alert["ra"], dec_deg=alert["dec"], jd=alert["jd"], **kwargs
    )


def check_mpc_name(name: str, **kwargs):
    alert = ampel_api_name(name, with_history=False)[0]["candidate"]
    return check_mpc_alert(alert, **kwargs)


check_mpc_name("ZTF23aajquxz")
