import numpy as np
import pandas as pd
import os
import datetime as dt
import pickle
import bisect
import scipy.stats
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import nc_utils
import coord_utils

"""
ionosonde_val.py
validation of SAMI3 model output using ML-cleaned ionosonde data

authors: Alex T Chartier & Glenn Sugar
Johns Hopkins University Applied Physics Laboratory
"""

def main(
    sday, eday, step, slist_fn, iono_data_fn_fmt, mod_fn_fmt, class_fn, 
    out_fn_fmt, dl_iono_data,
):
    """ 
    Validate the SAMI3 model against autoscaled, ML-cleaned ionosonde NmF2 data
    """
   
    #  runthrough, download and preprocess the ionosonde data 
    slist = get_station_list(slist_fn)

    day = sday - step
    while day < eday:
        # increment 
        day += step

        # filename handling 
        iono_data_fn = os.path.expanduser(day.strftime(iono_data_fn_fmt))

        # optional data download 
        if dl_iono_data:
            dl_ionosonde_data(day, day + step, slist, iono_data_fn)

    """ Clean up the ionosonde data using ML. Use the full set because the classifier requires a window of data.  """ 
    # combine all the files together into one
    iono_data = combine_preproc_files(sday, eday, step, iono_data_fn_fmt) 

    # kick out data from ionosondes that are not in the sitelist
    iono_sites = list(iono_data)
    for k in iono_sites:
        if k not in slist:
            iono_data.pop(k)

    # Load in the scaler and classifier
    [scaler, classifiers, betas] = nc_utils.unpickle(class_fn)

    # choose a single classifier 
    classifier = classifiers[3]

    for station, data in iono_data.items():
        print('Cleaning %s' % station)
        iono_data[station] = clean_ionosonde_data(data, scaler, classifier)

    """ Produce model expectation values of the good data """
    for station, data in iono_data.items():
        # Create a column for the model expectation values
        iono_data[station]['nmF2'] = 1.24E10 * data['foF2']  ** 2
        iono_data[station]['nmF2_sami'] = np.ones(len(data)) * np.nan

        # throw out the bad data
        print('Removing %i bad points ' % np.sum(data['qualflag'] == 1))
        iono_data[station] = data[data['qualflag'] == 0]

    day = sday - step
    while day < eday:
        # increment 
        day += step
        print(day)
        mod_fn = os.path.expanduser(day.strftime(mod_fn_fmt))
        try:
            sami = nc_utils.ncread_vars(mod_fn)
        except:
            print('%s not found' % mod_fn)
            continue

        for station, data in iono_data.items():
            print(station)
            iono_data[station] = calc_nmf2_errs(sami, data)

    """ calculate the LTs for each point and save out"""
    global_atts = global_atts_define()
    metadata = metadata_define()
    for station, data in iono_data.items():
        iono_data[station]['lt'] = coord_utils.calc_lt(data.index, data['lon'])
        nc_fname = out_fn_fmt % station
        nc_utils.write_netcdf_from_df(iono_data[station], metadata, global_atts, nc_fname)


def global_atts_define():
    return {
        'Conventions':'CF-1.6', 
        'title':'Ionosonde validation of SAMI3 with ML data classification', 
        'summary':'Data generated',
    }   


def metadata_define():
    return {
        'POSIXtime': {'units':'Seconds', 'long_name':'POSIX time (seconds since 1/1/1970)'},
        'station_id' : {'units':'None', 'long_name':'Station identifier'},
        'lat': {'units':'degrees', 'long_name':'Latitude'},
        'lon': {'units':'degrees', 'long_name':'Longitude'},
        'cs': {'units':'None', 'long_name':'Confidence score'},
        'foF2': {'units':'MHz', 'long_name':'F2-layer Critical frequency'},
        'hmF2': {'units':'km', 'long_name':'F2-layer peak height'},
        'foE'    : {'units':'MHz', 'long_name':'E-layer critical frequency'},
        'foEs': {'units':'MHz', 'long_name':'Sporadic E-layer critical frequency'},
        'qualflag': {'units':'None', 'long_name':'Machine-learning quality flag'},
        'nmF2': {'units':'El. m-3', 'long_name':'F2-layer peak density'},
        'nmF2_sami': {'units':'El. m-3', 'long_name':'F2-layer peak density from SAMI3'},
        'lt': {'units':'hours', 'long_name':'Local time'},
    }


def calc_nmf2_errs(sami, data):
    """ 
    Calculate the NmF2 error of SAMI3 compared to ionosonde data 
    Visualize using plot_station 
    """ 

    sami['dt'] = np.array([dt.datetime.utcfromtimestamp(t) for t in sami['time']])
    
    # loop over times in the ionosonde dataframe
    for time, row in data.iterrows():
        # find the closest SAMI3 time to measurement
        sami_tidx = find_closest_tidx(sami['dt'], time)
        if sami_tidx == None:
            continue

        # 2D interpolate the NmF2 field to measurement location
        intf = interp2d(sami['lon'], sami['lat'], np.squeeze(sami['nmf2'][sami_tidx, :, :]))
        nmF2_sami = intf(row['lon'], row['lat'])[0] * 1E6
        if np.isnan(nmF2_sami):
            breakpoint()
        data.at[time, 'nmF2_sami'] = nmF2_sami
       
    return data



def clean_ionosonde_data(iono_data, scaler, classifier,
    window_size = 9,
    tstep_min = 15,
):
    """ Use the ML classifier to identify good autoscaled ionosonde data points """

    # sort out the timestamps
    ts_rd = pd.to_datetime(iono_data.index).round('1s')
    times_julian = pd.to_datetime(iono_data.index).to_julian_date()


    # create some qual flags (0 good, 1 bad)
    iono_data['qualflag'] = np.ones(len(iono_data['foF2']))
    
    # sliding window filter
    assert window_size % 2, 'Need odd-sized window'
    halfwin = int(window_size / 2)
    data_idx = halfwin

    #  loop
    while data_idx < (len(ts_rd) - 1):

        # Increment the counter
        data_idx += 1

        # Set up the indexing (ID a single window of data)
        start_idx = data_idx - int(np.ceil(halfwin))
        mid_idx = data_idx
        end_idx = data_idx + int(np.ceil(halfwin)) + 1
        iono_data_t = iono_data[start_idx:end_idx]

        # Check feature vec is long enough and full of data
        if len(iono_data_t) < window_size:
            continue
        if iono_data_t['foF2'].isna().any():
            continue
        if iono_data_t['hmF2'].isna().any():
            continue

        # Generate the feature vector for the data point
        feature_vector = make_feature_vector(
            times_julian[start_idx:end_idx], 
            np.array(iono_data_t['foF2'], dtype=float), 
            slope=True,
            hmf_artist=np.array(iono_data_t['hmF2'], dtype=float), 
            confidence_scores=np.array(iono_data_t['cs'], dtype=float), 
            window_stats=True,
            mid_index=mid_idx - start_idx,
        ).reshape(1, -1) 

        # Scale and classify feature vector
        feature_vector_scaled = scaler.transform(feature_vector)
        is_outlier = classifier.predict(feature_vector_scaled)

        # Get the output
        if not is_outlier:
            iono_data.iat[data_idx, iono_data.columns.get_loc('qualflag')] = 0

    return iono_data


def combine_preproc_files(sday, eday, step, iono_data_fn_fmt):
    """ combine the full set of data """ 

    # load all the files
    iono_errs = []
    day = sday - step
    while day < eday:
        # increment 
        day += step

        # filenames
        iono_data_fname = os.path.expanduser(day.strftime(iono_data_fn_fmt))
        if not os.path.isfile(iono_data_fname):
            continue

        # load
        with open(iono_data_fname, 'rb') as f:
            iono_errs.append(pickle.load(f))

    # find set of sites
    sites = set()
    sites.update(*iono_errs)

    # concatenate each dataframe for all the stations available
    out = {}
    for site in sites:
        for iono_data in iono_errs:
            if site not in iono_data.keys():
                continue
            if iono_data[site].empty:
                continue
            if site in out.keys():
                out[site] = pd.concat([out[site], iono_data[site]])
            else:
                out[site] = iono_data[site]

    for k, v in out.items():
        out[k] = v.set_index('datetime')

    return out


def plot_station(nmf2_data, station):
    """ plot the output of calc_nmf2_errs """
    nmf2_data[station]['nmf2_obs'].plot()
    nmf2_data[station]['nmf2_sami'].plot()
    plt.ylabel("NmF2 (el./m3)")
    plt.grid()
    plt.legend()
    plt.title(station)
    plt.show()


def find_closest_tidx(tlist, time, dtmin=dt.timedelta(minutes=10)):
    """ return the index of the closest time in tlist to time """
    idx1 = bisect.bisect_left(tlist, time) - 1
    idx2 = bisect.bisect_right(tlist, time)

    if idx2 >= (len(tlist) - 1):
        idx1 = len(tlist) - 1
        idx2 = len(tlist) - 1

    t1 = abs(tlist[idx1] - time)
    t2 = abs(tlist[idx2] - time)

    # check we're within the dt specified
    if (t1 > dtmin) & (t2 > dtmin):
        return None

    # return the closest one
    if t1 < t2:
        idx = idx1
    else:
        idx = idx2

    return idx


def dl_ionosonde_data(
        stime, etime, station_ids, out_fn,
        params=['foF2', 'hmF2', 'foE', 'foEs'],
):
    """ get the ionosonde data from the DIDBase """
    out_dir = os.path.abspath(os.path.expanduser(os.path.dirname(out_fn)))
    os.makedirs(out_dir, exist_ok=True)
    out_fn = os.path.join(out_dir, os.path.basename(out_fn))
    data = {}

    if type(station_ids) is str:
        station_ids = [station_ids]
    for station_id in station_ids:
        data[station_id] = get_iono_params(station_id, params, stime, etime) 

    with open(out_fn, 'wb') as f:
        pickle.dump(data, f)
    print('Written to %s' % out_fn)



def get_iono_params(station_id, params, start_datetime, end_datetime):
    """
    Download the ionospheric parameters from the DIDBase
     To access this servlet use format: DIDBGetValues?ursiCode=CODE&charName=NAME&fromDate=yyyy.mm.dd%20(...)%20hh:mm:ss&toDate=yyyy.mm.dd%20(...)%20hh:mm:ss
     For example: DIDBGetValues?ursiCode=DB049&charName=foF2&fromDate=2007.06.25%20(...)%2000:00:00&toDate=2007.06.26%20(...)%2000:00:00
     or: DIDBGetValues?ursiCode=DB049&charName=foF2&fromDate=2007.06.25&toDate=2007.06.26

    """
    # List of acceptable params, see https://giro.uml.edu/didbase/scaled.php for more info
    params_ok = ['foF2', 'foF1', 'foE', 'foEs', 'fbEs', 'foEa', 'foP', 'fxI', 'MUFD', 'MD', 'hF2', 'hF', 'hE', 'hEs',
                 'hEa', 'hP', 'TypeEs', 'hmF2', 'hmF1', 'zhalfNm', 'yF2', 'yF1', 'yE', 'scaleF2', 'B0', 'B1', 'D1',
                 'TEC', 'FF', 'FE', 'QF', 'QE', 'fmin', 'fminF', 'fminE', 'fminEs', 'foF2p']
    # Create the params string
    if type(params) is list:
        params_str = ''
        for param in params:
            params_str += param + ','
        params_list = params
    else:
        params_str = params
        params_list = [params]

    # Validate each parameter
    for param in params_list:
        if param not in params_ok:
            raise(ValueError('Invalid parameter: ' + param))

    url = 'https://lgdc.uml.edu/common/DIDBGetValues?ursiCode={}&charName={}&fromDate={}&toDate={}'.\
        format(station_id, params_str, start_datetime.strftime('%Y.%m.%d+%H:%M:%S'), end_datetime.strftime('%Y.%m.%d+%H:%M:%S'))

    # Download the data
    print('Downloading data from: {}'.format(url))
    page = requests.get(url)
    page_lines = page.text.splitlines()
    past_header = False

    # Create a pandas table to hold data
    results = pd.DataFrame(columns=['station_id', 'lat', 'lon', 'datetime', 'cs', *params_list])
    for line in page_lines:
        if not past_header:
            if line[0:5] == '#Time':
                line_header = line.split()
                past_header = True
            elif line[0:10] == '# Location':
                lat_str = line.split()[3]
                station_lat = float(lat_str[0:-1])
                if lat_str[-1] == 'S':
                    station_lat = -station_lat
                lon_str = line.split()[4][0:-1]
                station_lon = float(lon_str[0:-1])
                if lon_str[-1] == 'W':
                    station_lon = -station_lon
        else:
            # Extract the data
            line_split = line.split()
            line_dict = {'station_id': station_id,
                         'lat': station_lat,
                         'lon': station_lon,
                         'datetime': datetime.datetime.strptime(line_split[0], '%Y-%m-%dT%H:%M:%S.%fZ'),
                         'cs': float(line_split[1])}
            param_num = 0
            for i, line_el in enumerate(line_split):
                # Skip the first and last elements
                if i <= 1 or i == (len(line_split)-1) or line_header[i] == 'QD':
                    continue
                # Parse the value
                try:
                    line_dict[line_header[i]] = float(line_el)
                except ValueError:
                    line_dict[line_header[i]] = np.nan
                param_num += 1
            results = pd.concat([results,  pd.DataFrame.from_dict(line_dict, orient='index').T], ignore_index=True)

    results.set_index('datetime')

    return results


def get_station_list(slist_fn):
    """ Lowell provides a site list online. Download it manually and read here"""

    # colnames = ['station_id', 'station_name', 'LAT', 'LONG']   # the list contains these, but station names can split in two
    slist = []
    with open(slist_fn, 'r') as f:
        txt = f.readlines()
    for line in txt:
        try:
            int(line[0])
            slist.append(line.split()[1])
            # Sites have a number at the beginning
        except:
            continue
            
    return slist


def make_feature_vector(times_julianday, fof2_artist, fof2_filtered=None, slope=True, hmf_artist=None,
                        confidence_scores=None, window_stats=True, mid_index=None):
    """ This will build a feature vector used for the classifier """
    if mid_index is None:
        mid_index = int(len(times_julianday) / 2)

    # Check if a filtered fof2 was passed
    if fof2_filtered is None:
        # See if we want to calculate window stats or just use the results
        if window_stats:
            features = np.hstack([fof2_artist[mid_index],
                                  np.nanmean(fof2_artist[0:mid_index]),
                                  np.nanstd(fof2_artist[0:mid_index]),
                                  np.nanmean(fof2_artist[mid_index+1:]),
                                  np.nanstd(fof2_artist[mid_index+1:])])
        else:
            features = fof2_artist
    else:
        fof2_filtered_minus_artist = fof2_filtered - fof2_artist
        if window_stats:
            try:
                features = np.hstack([fof2_artist[mid_index],
                                      np.nanmean(fof2_filtered_minus_artist[0:mid_index]),
                                      np.nanstd(fof2_filtered_minus_artist[0:mid_index]),
                                      np.nanmean(fof2_filtered_minus_artist[mid_index + 1:]),
                                      np.nanstd(fof2_filtered_minus_artist[mid_index + 1:])])
            except:
                print('error with stacking features')
        else:
            features = np.hstack([fof2_filtered_minus_artist,
                                  fof2_artist[mid_index]])

    if slope:
        finite_inds = np.isfinite(fof2_artist)
        finite_inds_start = finite_inds[0:mid_index+1]
        finite_inds_end = finite_inds[mid_index:]
        if sum(finite_inds) >= 2:
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(times_julianday[finite_inds],
                                                                                 fof2_artist[finite_inds])
        else:
            slope = np.nan
            r_value = np.nan
            p_value = np.nan
            std_err = np.nan
        if sum(finite_inds_start) >= 2:
            slope_start, intercept_start, r_value_start, p_value_start, std_err_start = scipy.stats.linregress(
                times_julianday[0:mid_index + 1][finite_inds_start], fof2_artist[0:mid_index + 1][finite_inds_start])
        else:
            slope_start = np.nan
            r_value_start = np.nan
            p_value_start = np.nan
            std_err_start = np.nan
        if sum(finite_inds_end) >= 2:
            slope_end, intercept_end, r_value_end, p_value_end, std_err_end = scipy.stats.linregress(
                times_julianday[mid_index:][finite_inds_end], fof2_artist[mid_index:][finite_inds_end])
        else:
            slope_end = np.nan
            r_value_end = np.nan
            p_value_end = np.nan
            std_err_end = np.nan
        features = np.hstack([features, slope, r_value, p_value, std_err,
                              slope_start, r_value_start, p_value_start, std_err_start,
                              slope_end, r_value_end, p_value_end, std_err_end])

    if hmf_artist is not None:
        if window_stats:
            # Calculate the hmf stats
            finite_inds = np.isfinite(hmf_artist)
            if sum(finite_inds) >= 2:
                slope_hmf, intercept_hmf, r_value_hmf, p_value_hmf, std_err_hmf = scipy.stats.linregress(
                    times_julianday[finite_inds], hmf_artist[finite_inds])
            else:
                slope_hmf = np.nan
                r_value_hmf = np.nan
                p_value_hmf = np.nan
                std_err_hmf = np.nan
            features = np.hstack([features, hmf_artist[mid_index], np.nanstd(hmf_artist), np.nanmean(hmf_artist),
                                  slope_hmf, r_value_hmf, p_value_hmf, std_err_hmf])
        else:
            slope_hmf, intercept_hmf, r_value_hmf, p_value_hmf, std_err_hmf = scipy.stats.linregress(times_julianday,
                                                                                                     hmf_artist)
            features = np.hstack([features, hmf_artist, slope_hmf, r_value_hmf, p_value_hmf, std_err_hmf])

    if confidence_scores is not None:
        # Add confidence scores
        if window_stats:
            features = np.hstack([features, confidence_scores[mid_index], np.nanmin(confidence_scores[0:mid_index]),
                                  np.nanmin(confidence_scores[mid_index+1:])])
        else:
            features = np.hstack([features, confidence_scores])

    return features




if __name__ == '__main__':
    
    sday = dt.datetime(2019, 3, 1)
    eday = dt.datetime(2019, 4, 1)
    step = dt.timedelta(days=1)
    slist_fn = 'ionosonde_classifier/lowell_iono_list.txt'
    iono_data_fn_fmt = '~/data/ionosonde/iono_data_%Y%b%d.pkl'
    out_fn_fmt = 'data/ionosonde_val/iono_%s_val.nc'
    mod_fn_fmt = '~/data/sami3/2019/sami3_regulargrid_ancillary_%Y%b%d.nc'
    class_fn = 'ionosonde_classifier/classifier_final.sav'
    dl_iono_data = False

    main(
        sday, eday, step, slist_fn, iono_data_fn_fmt, mod_fn_fmt, class_fn, 
        out_fn_fmt, dl_iono_data,
    )
