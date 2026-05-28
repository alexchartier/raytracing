"""
see http://wspr.rocks/liveembed/ for more config options
Note < 10m propagation is contaminated by repeaters
"""
import urllib.request
import json
import datetime as dt
import nc_utils
import pandas as pd
import xarray


def main(
    out_fn_fmt = '~/data/wspr/wspr_%s_%i.nc',
    datestr_fmt = '%Y%m%d',
    stime = dt.datetime(2024, 4, 5),
    etime = dt.datetime(2024, 4, 6),
    latlim = [-90, 90],
    lonlim = [-180, 180],
    bandlim = [10, 160],
    mindist = 250,
    maxlinks = 10000,
):
    """ 
    bands = def_bands()
    query = gen_query(band=bands[160], latlim=[25, 50], lonlim=[-125, -66], mindist=500, maxlinks=10)
    print(cmd)
    print(wsprlive_get(cmd))
    #print(wsprlive_get("SELECT * FROM wspr.rx where band = 10 and (time between '2024-05-08' and '2024-05-09') LIMIT 1"))
    """

    multiband_dl(
        out_fn_fmt, stime, etime, datestr_fmt,  
        latlim=latlim, lonlim=lonlim, mindist=mindist, bandlim=bandlim, maxlinks=maxlinks,
    )


def multiband_dl(
        out_fn_fmt, stime, etime, datestr_fmt,
        latlim=[25, 50], lonlim=[-125, -66], mindist=500, bandlim=[2, 160],
        maxlinks = 1E6,
        timeinc = dt.timedelta(days=1),
):
    bands = def_bands()
    time = stime 

    while time < etime:
        datestr = time.strftime(datestr_fmt)

        for band, idnum in bands.items():
            frames = []
            if band < min(bandlim) or band > max(bandlim):
                continue
            query = gen_query(times=[time, time + timeinc], band=idnum, 
                latlim=latlim, lonlim=lonlim, mindist=mindist, maxlinks=maxlinks)
            frames.append(reformat(wsprlive_get(query), band))
       
            df = pd.concat(frames) 
            out_fn = out_fn_fmt % (datestr, band)
            xarray.Dataset.from_dataframe(df).to_netcdf(out_fn)
            print(f'wrote to {out_fn}')

        time += timeinc


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


def gen_query(
        band=10, 
        times=[dt.datetime.now() - dt.timedelta(days=2), dt.datetime.now() - dt.timedelta(days=1)],
        latlim=None,
        lonlim=None,
        mindist=None,
        maxlinks=1E8,
):
    #timestr = '%Y-%m-%d %H:%M' 
    timestr = '%Y-%m-%d' 
    query = [
        f"band = {band}",
        f"(time between '{times[0].strftime(timestr)}' and '{times[1].strftime(timestr)}')",
    ]

    if latlim:
        query += [
        f"(tx_lat between {latlim[0]} and {latlim[1]}) ",
        f"(rx_lat between {latlim[0]} and {latlim[1]}) ",
        ]
    if lonlim:
        query += [
        f"(tx_lon between {lonlim[0]} and {lonlim[1]}) ",
        f"(rx_lon between {lonlim[0]} and {lonlim[1]}) ",
        ]
    if mindist:
        query += [f"distance > {mindist}",]
    query = "SELECT * FROM wspr.rx where " + " and ".join(query)
    maxlinks = int(maxlinks)
    query += f" LIMIT {maxlinks}"

    return query
    

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








