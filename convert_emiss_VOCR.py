import os
import xarray as xr
import numpy as np
import glob

KOH = {
    "PAR":   8.1e-13,
    "ETH":   2.5e-13,
    "OLE":   6.6e-12,
    "IOLE":  7.0e-12,
    "TOL":   5.6e-12,
    "ARO1":  3.4e-12,
    "ARO2":  2.0e-11,
    "ISOP":  1.0e-10,
    "TERP":  5.0e-11,
    "FORM":  8.5e-12,
    "ALD2":  1.5e-11,
    "ALDX":  2.7e-11,
    "KET":   1.7e-12,
    "MEOH":  9.0e-13,
    "ETOH":  3.3e-12
}

NEI_TO_CB06 = {
    # Paraffins
    "ETHA": "PAR",
    "HEXANE": "PAR",
    "NVOL": "PAR",
    "PRPA": "PAR",
    "UNR": "PAR",
    "IVOC": "PAR",

    # Olefins
    "ETHY": "OLE",
    "BUTADIENE13": "OLE",
    "CHLOROPRENE": "OLE",
    "STYRENE": "OLE",
    "IOLE": "IOLE",

    # Aromatics
    "BENZ": "ARO1",
    "ETHYLBENZ": "ARO2",
    "XYLENES": "ARO2",
    "XYLMN": "ARO2",
    "TOL": "TOL",
    "TOLU": "TOL",

    # Biogenics
    "ISOP": "ISOP",
    "APIN": "TERP",
    "TERP": "TERP",

    # Oxygenated VOCs
    "FORM": "FORM",
    "FORM_PRIMARY": "FORM",
    "ALD2": "ALD2",
    "ALD2_PRIMARY": "ALD2",
    "ALDX": "ALDX",
    "ACET": "KET",
    "ETOH": "ETOH",
    "MEOH": "MEOH",
    "GLY": "ALDX",
    "GLYD": "ALDX",
    "MGLY": "ALDX",

    # Optional slow VOC
    "ETH": "ETH"
}

# Where your daily files live
input_pattern = "./EMIS_NEI_2022_EPA_All_Anthro_OneLayer_*"
output_dir    = "./VOCRE_output"
os.makedirs(output_dir, exist_ok=True)

files = sorted(glob.glob(input_pattern))


for infile in files:
    fname = os.path.basename(infile)
    outfile = os.path.join(output_dir, fname.replace("EMIS_NEI_2022_EPA_All_Anthro_OneLayer_", "VOCRE_"))

    print(f"Processing {fname}")

    ds = xr.open_dataset(infile)

    vocre = None

    for nei_sp, cb06_sp in NEI_TO_CB06.items():
        if nei_sp not in ds.data_vars:
            continue
        if cb06_sp not in KOH:
            continue

        emis = ds[nei_sp]
        weighted = emis * KOH[cb06_sp]

        if vocre is None:
            vocre = weighted
        else:
            vocre = vocre + weighted

    if vocre is None:
        raise RuntimeError("No VOC species found for VOCRE calculation")

    vocre.name = "VOCRE"
    vocre.attrs = {
        "long_name": "OH-reactivity-weighted anthropogenic VOC emissions",
        "description": "NEI VOCs mapped to CB06 surrogates and weighted by KOH",
        "units": "mol s-1 * cm3 molecule-1 s-1"
    }

    out_ds = xr.Dataset(
        data_vars={"VOCRE": vocre},
        coords=ds.coords,
        attrs=ds.attrs
    )

    if "TFLAG" in ds.variables:
        out_ds["TFLAG"] = ds["TFLAG"]

    out_ds.to_netcdf(
        outfile,
        format="NETCDF4",
        encoding={"VOCRE": {"zlib": True, "complevel": 4}}
    )

    ds.close()
    out_ds.close()

print("All VOCRE files created.")
