# PROCESS — reusable functions for the verified lat/lon grid
import geopandas as gpd
import numpy as np
from shapely.geometry import box


def read_one_adm1(path, layer, label_column, label):
    """Read exactly one ADM1 polygon and return it in lon/lat coordinates."""
    field = '"' + label_column.replace('"', '""') + '"'
    value = "'" + label.replace("'", "''") + "'"

    region = gpd.read_file(
        path, layer=layer, columns=[label_column], where=f"{field} = {value}"
    )
    if len(region) != 1:
        raise ValueError(f"Expected one feature for {label!r}; found {len(region)}.")
    return region.to_crs("EPSG:4326")


def build_overlap_weights(region, grid, cell_size, area_crs):
    """Return the grid indices and overlap area in m² for each intersecting cell."""
    west, south, east, north = region.total_bounds
    lon = grid["lon"].values
    lat = grid["lat"].values
    half = cell_size / 2

    # Keep a margin so cells crossing the boundary are included.
    cols = np.flatnonzero((lon >= west - cell_size) & (lon <= east + cell_size))
    rows = np.flatnonzero((lat >= south - cell_size) & (lat <= north + cell_size))

    records = []
    for row in rows:
        for col in cols:
            # Snap edges to this product's regular grid to remove float32 noise.
            left = round((float(lon[col]) - half) / cell_size) * cell_size
            bottom = round((float(lat[row]) - half) / cell_size) * cell_size
            records.append((row, col, box(left, bottom, left + cell_size, bottom + cell_size)))

    cells = gpd.GeoDataFrame(
        records, columns=["row", "col", "geometry"], geometry="geometry", crs=region.crs
    )
    pieces = gpd.overlay(cells, region[["geometry"]], how="intersection",
                         keep_geom_type=True)
    pieces["area_m2"] = pieces.to_crs(area_crs).geometry.area
    return pieces.loc[pieces["area_m2"] > 0, ["row", "col", "area_m2"]].copy()


def clean_pm25(grid, missing_value):
    """Mask the documented workflow sentinel and flag unexpected negatives."""
    unexpected = (grid < 0) & (grid != missing_value)
    if unexpected.any().item():
        raise ValueError("Found a negative PM2.5 value other than the chosen sentinel.")
    return grid.where(grid != missing_value)


def summarize_one_year(clean_grid, weights):
    """Return the area-weighted mean and the areas used to calculate coverage."""
    values = clean_grid.transpose("lat", "lon").values[
        weights["row"].to_numpy(dtype=int), weights["col"].to_numpy(dtype=int)
    ]
    areas = weights["area_m2"].to_numpy()
    valid = np.isfinite(values)

    footprint_area = areas.sum()
    valid_area = areas[valid].sum()
    mean = np.average(values[valid], weights=areas[valid]) if valid_area > 0 else np.nan

    return {
        "mean_pm25": mean,
        "footprint_area_km2": footprint_area / 1e6,
        "valid_area_km2": valid_area / 1e6,
        "valid_fraction_of_footprint": valid_area / footprint_area if footprint_area > 0 else np.nan,
    }

# PROCESS — prepare boundary coordinates and verify annual grid alignment
import numpy as np
from pyproj import CRS


def prepare_boundary_crs(boundaries, grid_crs):
    """Return boundaries in the grid CRS. Require a known source CRS."""
    if boundaries.crs is None:
        raise ValueError("The ADM1 layer has no recorded CRS. Check its source metadata.")

    target = CRS.from_user_input(grid_crs)
    if boundaries.crs.equals(target, ignore_axis_order=True):
        return boundaries

    return boundaries.to_crs(target)


def check_annual_grid(ds, reference, *, year, variable, grid_crs):
    """Verify that saved row/column weights can be reused for this annual grid."""
    grid = ds[variable]

    if set(grid.dims) != {"lat", "lon"}:
        raise ValueError(f"{year}: expected lat/lon dimensions, found {grid.dims}.")

    # rioxarray reads the CRS if the NetCDF declares one.
    import rioxarray  # Registers the .rio accessor
    declared_crs = grid.rio.crs
    if declared_crs is not None:
        expected = CRS.from_user_input(grid_crs)
        if not expected.equals(declared_crs, ignore_axis_order=True):
            raise ValueError(f"{year}: declared CRS {declared_crs} differs from {grid_crs}.")

    for axis in ("lat", "lon"):
        if not np.array_equal(grid[axis].values, reference[axis]):
            raise ValueError(f"{year}: {axis} coordinates differ from the reference grid.")

    return "file metadata" if declared_crs is not None else "provider documentation"