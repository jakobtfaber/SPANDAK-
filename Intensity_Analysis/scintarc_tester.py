import copy
import datetime
import glob
import numpy as np
import os
import re
import time
import pickle
import chime_frb_api
import requests
import matplotlib.pyplot as plt

plt.switch_backend("agg")

from iautils.scripts.sburst_classifier import find_subbursts


response = requests.get("http://frb-vsop.chime:8001/v1/candidates")
response.raise_for_status()
candidates = response.json()
#eventids = [i["events"][0] for i in candidates]
eventids = [30465284, 41643446, 44961057]

results = dict()

for eventid in eventids:

    # interact with CHIME/FRB's API
    master = chime_frb_api.frb_master.FRBMaster()

    #try:
    header = master.events.get_event(event_number=eventid, full_header=True)

    # TODO extract beam number from header

    detected_beams = []
    snrs = []
    for beam_header in header["event_beam_header"]:
        detected_beams.append("{:04d}".format(beam_header["beam_no"]))
        snrs.append(beam_header["snr"])

    max_idx = np.argmax(snrs)
    max_beam_no = detected_beams[max_idx]
    # print('max beam no', max_beam_no)

    event_time = datetime.datetime.strptime(
        header["event_beam_header"][max_idx]["timestamp_utc"], "%Y%m%d%H%M%S.%f"
    )
    # print('event time', event_time)

    # TODO convert event_time into date path

    path = os.path.join(*event_time.strftime("%Y %m %d").split(" "))
    path = os.path.join(path, "astro_{}".format(eventid), "intensity", "processed")
    path = os.path.join(path, max_beam_no)
    # print('path', path)

    filepath = os.path.join("/data/frb-archiver", path)
    # print('filepath', filepath)

    # TODO construct full path
    fname = glob.glob(os.path.join(filepath, "*snapshot.npz"))
    # print('fname', fname)

    # TODO check if file exists

    if fname:  # TODO load file
        fname_ = os.path.join(filepath, fname[0])
        # print('fname_', fname_)

        waterfall = np.load(fname_)
        #intensity_array = waterfall['intensity']
        try:
            intensity_array, freq, time_fpga,\
            event_time_fpga, centroid_coord_x_array, centroid_coord_x_units, \
            centroid_coord_y_array, centroid_coord_y_units, \
            centroid_prominence, centroid_widths, component_number,\
            centroid_drift_rate, sidelobe, frequency_widths = find_subbursts(waterfall, eventid, plot_dyn = True, plot_mask = True)
            #print(centroid_location_time)
            results[eventid] = dict()
            results[eventid]['centroid_coord_x'] = centroid_coord_x_units
            results[eventid]['centroid_coord_y'] = centroid_coord_y_units
            results[eventid]['centroid_prominence'] = centroid_prominence
            results[eventid]['centroid_widths'] = centroid_widths
            results[eventid]['component_number'] = component_number
            results[eventid]['centroid_drift_rate'] = centroid_drift_rate
            results[eventid]['sidelobe_detection'] = sidelobe
            results[eventid]['frequency widths'] = frequency_widths
            #results[eventid]['freq'] = freq
            #results[eventid]['time_fpga'] = time_fpga
            #results[eventid]['event_time_fpga'] = event_time_fpga
        except ValueError:
            print("No Hint of Burst Detected: ", eventid)
            results[eventid] = dict()
            results[eventid]['centroid_coord_x'] = np.nan
            results[eventid]['centroid_coord_y'] = np.nan
            results[eventid]['centroid_prominence'] = np.nan
            results[eventid]['centroid_widths'] = np.nan
            results[eventid]['component_number'] = np.nan
            results[eventid]['centroid_drift_rate'] = np.nan
            results[eventid]['sidelobe_detection'] = np.nan
            results[eventid]['frequency widths'] = np.nan
            #results[eventid]['freq'] = np.nan
            #results[eventid]['time_fpga'] = np.nan
            #results[eventid]['event_time_fpga'] = np.nan
            pass
    else:
        continue
    #except:
        #print('HTTPError: 500 Server Error', eventid)
        #pass
    # TODO save dictionary to disk as pickle file

    # Store data (serialize)
    with open('sbclass_results_mat.pickle', 'wb') as handle:
        pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Load data (deserialize)
    with open('sbclass_results_mat.pickle', 'rb') as handle:
        unserialized_data = pickle.load(handle)

    print(results == unserialized_data)
    print("event done: ", eventid)
