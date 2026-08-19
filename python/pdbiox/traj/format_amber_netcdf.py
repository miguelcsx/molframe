"""AMBER NetCDF trajectory bindings."""

from .._trajectory import (
    AmberNetcdfError,
    AmberNetcdfMetadata,
    AmberNetcdfPrecision,
    AmberNetcdfTrajectory,
    AmberNetcdfWriteOptions,
    parse_amber_netcdf,
    parse_amber_netcdf_record,
    write_amber_netcdf,
)

__all__ = [
    "AmberNetcdfError",
    "AmberNetcdfMetadata",
    "AmberNetcdfPrecision",
    "AmberNetcdfTrajectory",
    "AmberNetcdfWriteOptions",
    "parse_amber_netcdf",
    "parse_amber_netcdf_record",
    "write_amber_netcdf",
]
