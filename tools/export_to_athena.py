import argparse
import csv
from pathlib import Path

import h5py

CHANNEL_NNLS_NAMES_PATH = "/MAPS/XRF_Analyzed/NNLS/Channel_Names"
CHANNEL_NNLS_VALUES_PATH = "/MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec"

CHANNEL_FITTED_NAMES_PATH = "/MAPS/XRF_Analyzed/Fitted/Channel_Names"
CHANNEL_FITTED_VALUES_PATH = "/MAPS/XRF_Analyzed/Fitted/Counts_Per_Sec"

SCALER_NAMES_PATH = "/MAPS/Scalers/Names"
SCALER_VALUES_PATH = "/MAPS/Scalers/Values"

X_AXIS_PATH = "/MAPS/Scan/x_axis"

QUANT_NNLS_LABELS_PATH = '/MAPS/Quantification/Calibration/NNLS/Calibration_Curve_Labels'
QUANT_NNLS_US_PATH = '/MAPS/Quantification/Calibration/NNLS/Calibration_Curve_US_IC'
QUANT_NNLS_DS_PATH = '/MAPS/Quantification/Calibration/NNLS/Calibration_Curve_DS_IC'

QUANT_FITTED_LABELS_PATH = '/MAPS/Quantification/Calibration/Fitted/Calibration_Curve_Labels'
QUANT_FITTED_US_PATH = '/MAPS/Quantification/Calibration/Fitted/Calibration_Curve_US_IC'
QUANT_FITTED_DS_PATH = '/MAPS/Quantification/Calibration/Fitted/Calibration_Curve_DS_IC'

I_PATH = 0
I_C_NAME = 1
I_C_VAL = 2
I_N_NAME = 3
I_N_LABEL = 4
I_N_VAL = 5

def find_quant_values(h5f, to_export, labels_path, values_path):
    labels = [
        n.decode() if isinstance(n, bytes) else n
        for n in h5f[labels_path][()].flatten()
    ]
    values = h5f[values_path][()]
    results = {}
    for wanted in to_export:
        if wanted in labels:
            idx = labels.index(wanted)
            results[wanted] = values[:, idx]
    return results

def find_matching_entries(h5f, to_export, CHANNEL_NAMES_PATH):
    matches = {}
    if CHANNEL_NAMES_PATH is not None:
        matches = {CHANNEL_NAMES_PATH: [], SCALER_NAMES_PATH: []}
    else:
        matches = {SCALER_NAMES_PATH: []}

    for ds_path in matches:
        names = [
            n.decode() if isinstance(n, bytes) else n
            for n in h5f[ds_path][()]
        ]
        for wanted in to_export:
            if wanted in names:
                matches[ds_path].append((wanted, names.index(wanted)))
    return matches

def get_entrie_values(h5f, to_export):
    matches = {SCALER_NAMES_PATH: []}
    values = {}
    for ds_path in matches:
        names = [
            n.decode() if isinstance(n, bytes) else n
            for n in h5f[ds_path][()]
        ]
        for wanted in to_export:
            if wanted in names:
                #matches[ds_path].append((wanted, names.index(wanted)))
                #print(wanted, ds_path, names.index(wanted))
                values[wanted] = h5f[SCALER_VALUES_PATH][names.index(wanted), :, :]
    return values


def write_header(writer, matches):
    header = ["##", "axis_pos"]
    extra_to_add = 10 - len(matches)
    for entries in matches.values():
        for name, _ in entries:
            header.append(f"{name}[]")
    for e in range(extra_to_add):
        header.append(f"[]")
    writer.writerow(header)

def write_rows(writer, x_axis, h5f, matches, CHANNEL_NAMES_PATH, CHANNEL_VALUES_PATH, normalizer_arr, quat_dict):
    extra_to_add = 10 - len(matches)
    data_for = { CHANNEL_NAMES_PATH: h5f[CHANNEL_VALUES_PATH][()],
                 SCALER_NAMES_PATH: h5f[SCALER_VALUES_PATH][()] }   

    x_flat = x_axis.flatten()
    for i, x_val in enumerate(x_flat):
        row = [x_val, x_val / 1000.0]
        for ds_path, entries in matches.items():
            data = data_for[ds_path]
            for el_name, idx in entries:
                values = data[idx].flatten()
                if normalizer_arr is not None and quat_dict is not None:
                    if el_name in quat_dict:
                        idx = 0
                        elen = len(el_name)
                        if elen > 2:
                            if el_name[elen -2] == '_':
                                if el_name[elen -1] == 'L':
                                    idx = 1
                                elif el_name[elen -1] == 'M':
                                    idx = 2
                        values[i] = values[i] / normalizer_arr[0][i] / quat_dict[el_name][0]
                row.append(values[i])
        for e in range(extra_to_add):
            row.append(0)
        writer.writerow(row)

def export_to_csv(hdf5_path, to_export):
    hdf5_path = Path(hdf5_path)
    new_name = f"{hdf5_path.stem}_nnls_cts.csv"
    hdf5_nnls_cts_path = hdf5_path.with_name(new_name)
    new_name = f"{hdf5_path.stem}_fitted_cts.csv"
    hdf5_fitted_cts_path = hdf5_path.with_name(new_name)
    new_name = f"{hdf5_path.stem}_nnls_us_ic_ugcm2.csv"
    hdf5_nnls_us_path = hdf5_path.with_name(new_name)
    new_name = f"{hdf5_path.stem}_fitted_us_ic_ugcm2.csv"
    hdf5_fitted_us_ugcm2_path = hdf5_path.with_name(new_name)
    new_name = f"{hdf5_path.stem}_nnls_ds_ic_ugcm2.csv"
    hdf5_nnls_ds_ugcm2_path = hdf5_path.with_name(new_name)
    new_name = f"{hdf5_path.stem}_fitted_ds_ic_ugcm2.csv"
    hdf5_fitted_ds_ugcm2_path = hdf5_path.with_name(new_name)

    csv_nnls_cts_path = (hdf5_nnls_cts_path, CHANNEL_NNLS_NAMES_PATH, CHANNEL_NNLS_VALUES_PATH, None, None, None)
    csv_fitted_cts_path = (hdf5_fitted_cts_path, CHANNEL_FITTED_NAMES_PATH, CHANNEL_FITTED_VALUES_PATH, None, None, None)
    csv_nnls_us_ugcm2_path = (hdf5_nnls_us_path, CHANNEL_NNLS_NAMES_PATH, CHANNEL_NNLS_VALUES_PATH, ['US_IC',], QUANT_NNLS_LABELS_PATH, QUANT_NNLS_US_PATH)
    csv_fitted_us_ugcm2_path = (hdf5_fitted_us_ugcm2_path, CHANNEL_FITTED_NAMES_PATH, CHANNEL_FITTED_VALUES_PATH, ['US_IC',], QUANT_FITTED_LABELS_PATH, QUANT_FITTED_US_PATH)
    csv_nnls_ds_ugcm2_path = (hdf5_nnls_ds_ugcm2_path, CHANNEL_NNLS_NAMES_PATH, CHANNEL_NNLS_VALUES_PATH, ['DS_IC',], QUANT_NNLS_LABELS_PATH, QUANT_NNLS_DS_PATH)
    csv_fitted_ds_ugcm2_path = (hdf5_fitted_ds_ugcm2_path, CHANNEL_FITTED_NAMES_PATH, CHANNEL_FITTED_VALUES_PATH, ['DS_IC',], QUANT_FITTED_LABELS_PATH, QUANT_FITTED_DS_PATH)
    all_csv_paths = [csv_nnls_cts_path, csv_fitted_cts_path, csv_nnls_us_ugcm2_path, csv_fitted_us_ugcm2_path, csv_nnls_ds_ugcm2_path, csv_fitted_ds_ugcm2_path] 

    with h5py.File(hdf5_path, "r") as h5f:
        x_axis = h5f[X_AXIS_PATH][()]
        for csv_obj in all_csv_paths:
            normalizer_arr = None
            norm_vals = None
            if csv_obj[I_N_NAME] is not None:
                normalizer_arr = get_entrie_values(h5f, csv_obj[I_N_NAME])[csv_obj[I_N_NAME][0]]
                #print(normalizer_arr)
                norm_vals = find_quant_values(h5f, to_export, csv_obj[I_N_LABEL], csv_obj[I_N_VAL])
            matches = find_matching_entries(h5f, to_export, csv_obj[I_C_NAME])
            print(csv_obj[I_PATH])
            csv_f = open(csv_obj[I_PATH], "w", newline="")
            writer = csv.writer(csv_f, delimiter="\t")
            write_header(writer, matches)
            write_rows(writer, x_axis, h5f, matches, csv_obj[I_C_NAME], csv_obj[I_C_VAL], normalizer_arr, norm_vals)
            csv_f.close()
            

def main():
    parser = argparse.ArgumentParser(description="Export XRF data to Athena format.")
    parser.add_argument("hdf5_path", type=str, help="Path to the HDF5 file.")
    parser.add_argument(
        "to_export",
        type=str,
        nargs="+",
        help="List of strings.",
    )
    args = parser.parse_args()

    print (args)
    export_to_csv(args.hdf5_path, args.to_export)
    return args


if __name__ == "__main__":
    main()
