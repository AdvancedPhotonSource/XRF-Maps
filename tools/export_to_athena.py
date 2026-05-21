import argparse
import csv
from pathlib import Path

import h5py


#CHANNEL_NAMES_PATH = ["/MAPS/XRF_Analyzed/NNLS/Channel_Names", "/MAPS/XRF_Analyzed/Fitted/Channel_Names"]
#CHANNEL_VALUES_PATH = ["/MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec", "/MAPS/XRF_Analyzed/Fitted/Counts_Per_Sec"]
CHANNEL_NAMES_PATH = "/MAPS/XRF_Analyzed/NNLS/Channel_Names"
CHANNEL_VALUES_PATH = "/MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec"
SCALER_NAMES_PATH = "/MAPS/Scalers/Names"
SCALER_VALUES_PATH = "/MAPS/Scalers/Values"
X_AXIS_PATH = "/MAPS/Scan/x_axis"
#QUANT_PATH


def find_matching_entries(h5f, to_export):
    matches = {CHANNEL_NAMES_PATH: [], SCALER_NAMES_PATH: []}

    for ds_path in matches:
        names = [
            n.decode() if isinstance(n, bytes) else n
            for n in h5f[ds_path][()]
        ]
        for wanted in to_export:
            if wanted in names:
                matches[ds_path].append((wanted, names.index(wanted)))

    return matches

def write_header(writer, matches):
    header = ["##", "axis_pos"]
    extra_to_add = 10 - len(matches)
    for entries in matches.values():
        for name, _ in entries:
            header.append(f"{name}[]")
    for e in range(extra_to_add):
        header.append(f"[]")
    writer.writerow(header)

def write_rows(writer, x_axis, h5f, matches):
    extra_to_add = 10 - len(matches)
    data_for = {
        CHANNEL_NAMES_PATH: h5f[CHANNEL_VALUES_PATH][()],
        SCALER_NAMES_PATH: h5f[SCALER_VALUES_PATH][()],
    }

    x_flat = x_axis.flatten()
    for i, x_val in enumerate(x_flat):
        row = [x_val, x_val / 1000.0]
        for ds_path, entries in matches.items():
            data = data_for[ds_path]
            for _, idx in entries:
                values = data[idx].flatten()
                row.append(values[i])
        for e in range(extra_to_add):
            row.append(0)
        writer.writerow(row)

def export_to_csv(hdf5_path, to_export):
    hdf5_path = Path(hdf5_path)
    csv_path = hdf5_path.with_suffix(".csv")

    with h5py.File(hdf5_path, "r") as h5f, open(csv_path, "w", newline="") as csv_f:
        writer = csv.writer(csv_f, delimiter="\t")
        x_axis = h5f[X_AXIS_PATH][()]
        matches = find_matching_entries(h5f, to_export)        
        write_header(writer, matches)
        write_rows(writer, x_axis, h5f, matches)

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
