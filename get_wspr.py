# see http://wspr.rocks/liveembed/ for more config options

import urllib.request
import json
import datetime as dt
import nc_utils
import pandas as pd
import xarray


def main(
    out_fn_fmt = '~/data/wspr/%Y%b%d-wspr.nc',
    stime = dt.datetime(2024, 4, 5),
    etime = dt.datetime(2024, 4, 10),
):
    """ 
    bands = def_bands()
    cmd = gen_cmd(band=bands[160], latlim=[25, 50], lonlim=[-125, -66], mindist=500, maxlinks=10)
    print(cmd)
    print(wsprlive_get(cmd))
    #print(wsprlive_get("SELECT * FROM wspr.rx where band = 10 and (time between '2024-05-08' and '2024-05-09') LIMIT 1"))
    """

    multiband_dl(
        out_fn_fmt, stime, etime, 
        latlim=[25, 50], lonlim=[-125, -66], mindist=500, min_band=2,
    )


def multiband_dl(
        out_fn_fmt, stime, etime, 
        latlim=[25, 50], lonlim=[-125, -66], mindist=500, min_band=10,
):
    bands = def_bands()
    time = stime 
    while time < etime:
        out_fn = time.strftime(out_fn_fmt)
        frames = []
        for band, idnum in bands.items():
            if band < min_band:
                continue
            cmd = gen_cmd(band=idnum, latlim=latlim, lonlim=lonlim, mindist=mindist)
            frames.append(reformat(wsprlive_get(cmd), band))
       
        df = pd.concat(frames) 
        xarray.Dataset.from_dataframe(df).to_netcdf(out_fn)
        print(f'wrote to {out_fn}')
        time += dt.timedelta(days=1)


def wsprlive_get(query):
    # put together the request url
    url = "https://db1.wspr.live/?query=" + urllib.parse.quote_plus(query + " FORMAT JSON")

    # download contents from wspr.live
    contents = urllib.request.urlopen(url).read()

    # return the json decoded data
    return json.loads(contents.decode("UTF-8"))["data"]


def def_bands():
    """ WSPRnet band numbers in {wavelen (m): idx} """
    return {
        160: 1,
        80: 3,
        60: 5,
        40: 7,
        30: 10,
        20: 14,
        17: 18,
        15: 21,
        12: 24,
        10: 28,
        6: 50,
        4: 70,
        2: 144,
        0.70: 432,
        0.23: 1296,
}

def gen_cmd(
        band=10, 
        times=[dt.datetime.now() - dt.timedelta(days=2), dt.datetime.now() - dt.timedelta(days=1)],
        latlim=None,
        lonlim=None,
        mindist=None,
        maxlinks=1E6,
):
    #timestr = '%Y-%m-%d %H:%M' 
    timestr = '%Y-%m-%d' 
    cmd = [
        f"band = {band}",
        f"(time between '{times[0].strftime(timestr)}' and '{times[1].strftime(timestr)}')",
    ]

    if latlim:
        cmd += [
        f"(tx_lat between {latlim[0]} and {latlim[1]}) ",
        f"(rx_lat between {latlim[0]} and {latlim[1]}) ",
        ]
    if lonlim:
        cmd += [
        f"(tx_lon between {lonlim[0]} and {lonlim[1]}) ",
        f"(rx_lon between {lonlim[0]} and {lonlim[1]}) ",
        ]
    if mindist:
        cmd += [f"distance > {mindist}",]
    cmd = "SELECT * FROM wspr.rx where " + " and ".join(cmd)
    maxlinks = int(maxlinks)
    cmd += f" LIMIT {maxlinks}"
    return cmd
    

def reformat(data, wlen):
    fields = 'time', 'frequency', 'distance', 'tx_lat', 'tx_lon', 'rx_lat', 'rx_lon',
    d2 = []
    for entry in data:
        e2 = {}
        for k, v in entry.items():
            if k in fields:
                e2[k] = v
        e2['wavelength'] = wlen
        d2.append(e2)

    df = pd.DataFrame.from_records(d2)
    if not df.empty:
        df = df.set_index('time')
    
    return df


if __name__ == "__main__":
    main()








