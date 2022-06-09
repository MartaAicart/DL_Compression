#!/usr/bin/env python
import os
import sys
import xarray as xr
import numpy as np
import s3fs
import caterva as cat
def open_zarr(year, month, datestart, dateend):
    fs = s3fs.S3FileSystem(anon=True)
    datestring = "era5-pds/zarr/{year}/{month:02d}/data/".format(year=year, month=month)
    s3map = s3fs.S3Map(datestring + "air_temperature_at_2_metres_1hour_Minimum.zarr/", s3=fs)
    precip_zarr = xr.open_dataset(s3map, engine="zarr")
    precip_zarr = precip_zarr.sel(time1=slice(np.datetime64(datestart), np.datetime64(dateend)))
    return precip_zarr.air_temperature_at_2_metres_1hour_Minimum
print("Fetching data from S3 (era5-pds)...")
precip_m0 = open_zarr(1987, 10, "1987-10-01", "1987-10-30 23:59")
if os.path.exists("temp1.cat"):
    cat.remove(path)
# ia.set_config_defaults(favor=ia.Favor.SPEED)
m_shape = precip_m0.shape
m_chunks = (128, 128, 256)
m_blocks = (32, 32, 32)
cat_precip0 = cat.empty(m_shape, itemsize=4, chunks=m_chunks, blocks=m_blocks,
                        urlpath="temp1.cat", contiguous=True)
print("Fetching and storing 1st month...")
values = precip_m0.values
cat_precip0[:] = values
