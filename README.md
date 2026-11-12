# From a PM2.5 grid to a California state average

California’s **2024 area-weighted mean PM2.5 concentration was 6.1139 µg/m³**. This notebook shows how I calculated that value from a grid of annual fine particulate matter concentrations. It uses California as one example of a first-level administrative region (ADM1).

[Open the California notebook and see the map](demo/california_pm25_2024.ipynb)

## How the average works

A grid cell on the state border may lie partly outside California. I clip each cell to the state boundary and measure the area of the piece inside the state. That area is the cell’s weight.

I multiply each valid concentration by its overlap area, add those products, and divide by the total valid overlap area. Cells with missing concentrations contribute to neither sum. In the 2024 calculation, **409,583.9 km²** contributed to the mean. The result describes concentrations across California’s area. Estimating average exposure per resident requires population weights.

## Run the California example

You’ll need two local inputs:

- The **2024 global annual PM2.5 NetCDF** from the CNN V6.GL.03 product on its 0.1° grid. The notebook reads the `PM25` variable.
- A **GADM 4.1 ADM1 GeoPackage** with the California boundary. The example selects California using `VARNAME_1 = CA|Calif.`. If your GeoPackage has multiple layers, set `boundary_layer` in the notebook’s CONFIG cell.

From the repository folder, run these commands if you need a new Conda environment:

```bash
conda create -n adm1-demo -c conda-forge python=3.11 jupyterlab xarray netcdf4 geopandas pyogrio shapely pyproj numpy matplotlib
conda activate adm1-demo
jupyter lab demo/california_pm25_2024.ipynb
```

Run the cells from top to bottom. Enter the two file paths when prompted. The notebook builds the overlap weights, checks the result against the California benchmark, and draws the map. The calculation functions are in [aggregation.py](demo/aggregation.py).

