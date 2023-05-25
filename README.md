# napari-hdf5-netcdf-viewer

A napari plugin for visualizing geoscience HDF5 and netCDF data

![](./img/napari_geo_demo.gif)


Notes on this current version:
1. Xarray is used to visualize netCDF, and h5py is used for HDF5. The dimensions stored in HDF5 doesn't work easily with xarray, so h5py is used instead. 

2. For netCDF, it displays world map as a napari shape layer. This only works well for entire globe domain, not for a regional domain, which will be added in future. 

3. Map projections are not available yet. 




