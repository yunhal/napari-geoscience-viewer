import napari
import h5py
import netCDF4 as nc
import numpy as np
import pathlib
from pathlib import Path
import xarray as xr
from magicgui import magicgui, magic_factory
import sys
from napari_geojson import napari_get_reader


# potential file extensions: 
H5_FILENAME_LIST = ['h5', 'hdf5', 'hdf']
NETCDF_FILENAME_LIST = ['nc', 'netcdf', 'ncf', 'nc4']  

class FileReader:

    def __init__(self, viewer):
        self.data_list = []
        self.viewer = viewer

    def print_h5_objs(self, name, obj):
        """
        Parameters
        ----------
        name : hdf5 data path 
                (e.g., HDFEOS/GRIDS/OMI Column Amount O3/Data Fields/ColumnAmountO3)
            
        obj: hdf5 data
        """

        # save the data name into data_list 
        self.data_list.append(name)
        if isinstance(obj, h5py.Group): # test for dataset
            print("Group ", name)
        if isinstance(obj, h5py.Dataset): # test for dataset
            print('    Dateset {} : shape = {} and dtype = {}'.format(name, str(obj.shape), str(obj.dtype) ) )

            if (h5py.check_string_dtype(obj.dtype) is not None):
                (encod, ln) = h5py.check_string_dtype(obj.dtype)

                if (ln is None): 
                    ## This is variable length string 
                    print(' =============== {} contains the information below ========= '.format(name))
                    if (len(obj) > 30):
                        print("Only first 30 elements are printed here. It contains the following unique elements", list(set(obj)))
                        print(obj[:30])
                    else:
                        print(obj[:])
                    print(' ===============  end of {} data =========================== '.format(name))       
                    print("                                               ")                          
                
                else:
                    print(' ================ {} contains the information below ======== '.format(name))
                    print(np.array(obj).item().decode(encod))
                    print(' ================  end of {} data ========================== '.format(name))       
                    print("                                               ")


    def print_netcdf(self, fn):

        """
        Parameters
        ----------
        fn : netcdf file path 
        """

        with nc.Dataset(fn) as f:
            print("                                               ")
            print(" =======", fn, " file structure start =========")
            for var in f.variables.keys():

                # save data name into data_list 
                self.data_list.append(var)
                print('Variables  {} :   shape = {}, dimesion = {}'.format(var, f[var].shape, f[var].dimensions))

            for d in f.dimensions.keys():
                print('Dimensions  {} :   size = {} '.format(f.dimensions[d].name, f.dimensions[d].size))

            for attr in f.ncattrs():
                print("Global attrs  {} = {}".format(attr, getattr(f, attr)))

            print(" =======", fn, " file structure end============")
            print("                                               ")

            # For netcdf, print(f) will do the same thing as above but it is not pretty. 


    def all_h5_to_napari(self, fn, h5_varname):

        """
        Parameters
        ----------
        fn : netcdf file path 
        h5_varname : variable name
        """

        label_name = 'label'
        point_name = 'annotation'

        with h5py.File(fn,'r') as f:

            var_path = str(h5_varname)
            data = f[var_path]
            np_data = np.array(data)

            if hasattr(data, 'shape'):

                if data.ndim >= 2: 

                    if (label_name.casefold() in var_path.casefold() ):
                        # read as labels
                        print("WARNING: label layer is determined, based on matching {} in the data path".format(label_name))
                        layer = self.viewer.add_labels(np_data, name= var_path)

                    elif (point_name.casefold() in var_path.casefold() ):
                        # read as points
                        print("WARNING: point layer is determined, based on matching {} in the data path".format(point_name))
                        layer = self.viewer.add_points(np_data, name= var_path)
                    else: 
                        # read as image
                        for key, val in data.attrs.items():
                            if (key == "_FillValue" ): 
                                print("_FillValue exits", len(np_data [np_data == data.attrs["_FillValue"][0]]))
                                np_data [np_data == data.attrs["_FillValue"][0]] = np.nan
                                
                            elif (key == "MissingValue" ): 
                                print("MissingValue exits", len(np_data [np_data == data.attrs["MissingValue"][0]]))
                                np_data [np_data == data.attrs["MissingValue"][0]] = np.nan

                        layer = self.viewer.add_image(np_data, name= var_path)


    def all_netcdf_to_napari(self, fn):

        """
        Parameters
        ----------
        fn : netcdf file path 
        """

        with nc.Dataset(fn) as f:

            for d in self.data_list:

                var_path = str(d)
                data = f[var_path]
                np_data = np.array(data)
                    
                if hasattr(data, 'shape'):

                    if data.ndim >= 2: 

                        if '_FillValue' in data.ncattrs() :
                            fvalue = getattr(data, '_FillValue')
                            print ("_Fill_Value was {} and is replaced with np.nan".format(fvalue))
                            print(np_data.shape, fvalue.shape)
                            np_data [np_data == fvalue] = np.nan

                        if 'missing_value' in data.ncattrs() :
                            mvalue = getattr(data, 'missing_value')
                            print ("missing_value was {} and is replaced with np.nan".format(mvalue))
                            np_data [np_data == mvalue] = np.nan

                        # read as image
                        layer = self.viewer.add_image(np_data, name= var_path)
    

    def get_geo_dataset(self, fn, xr_data):

        """
        Parameters
        ----------
        fn : h5netcdf file path 
        xr_data : data variable name you want to visualize (e.g., "O3", "PM25")
        """

        file =xr.open_dataset(fn, engine="h5netcdf")

        df = file[xr_data]
        print("data has the following coords : ", df.coords)
        # boolean variables to check the dimensions for NAPARI
        LONG_EXIST = False
        LAT_EXIST = False
        VERT_EXIST = False

        # get longitude coords
        long_list =['lon', 'longitude', 'long', 'x','LON','LONGITUDE', 'LONG', 'X','XLONG']

        for x in long_list:
            
            if(x in df.coords):

                LONG_EXIST=True
                longitude = df.coords[x].values
                lon_atrib = df.coords[x].attrs
                longitude_scale = (np.amax(longitude) - np.amin(longitude))/len(longitude)

                print('original long is ', longitude[:10])
                print(np.amin(longitude),np.amax(longitude)) 

                if(np.amin(longitude) >= 0 and np.amax(longitude) > 180):
                    print("longitude is shifting to 0 to -180")
                    # This condition won't work well for regional cases

                    # longitude is shifting from 0~360 to -180~180
                    df.coords[x] = (df.coords[x] + 180) % 360 - 180
                    df = df.sortby(df.lon)
                    df.coords[x].attrs = lon_atrib
                
                    # overwrite with shifted longitude
                    longitude = df.coords[x].values
                
                    print('shifted long is ', longitude[:10])
                
        # Xarray converts 'missing_values' and '_FillValue' to NaN as a default. 
        if hasattr(df.attrs, 'fmissing_value'): 
            fvalue = df.attrs['fmissing_value']
            print ("fmissing_value was {} and is replaced with np.nan".format(fvalue))
            print(np_data.shape, fvalue.shape)
            np_data [np_data == fvalue] = np.nan

        # get latitude coords
        lat_list =['lat', 'latitude', 'y', 'LAT', 'LATITUDE', 'Y','XLAT']

        for x in lat_list:
            if(x in df.coords):
                LAT_EXIST=True
                latitude = df.coords[x].values
                latitude_scale = (np.amax(latitude) - np.amin(latitude))/len(latitude)
        #print('lat is ', latitude[:10]) 

        # get vertical coords
        lev_list =['lev', 'levels', 'level','vertical', 'height','altitude', 'pressure','ver']

        for x in lev_list:
            if(x in df.coords):
                VERT_EXIST=True
                vertical = df.coords[x].values
                #print('lev is ', df.coords[x].values[:10])

                vertical_scale = 5 # hard-coded number  

                # first vertical layer to be near the groun
                if(vertical[0] < vertical[-1]):
                    df = df.sortby(x, ascending=False)
        
        # for WRF or CMAQ, it may fail to get lat and long information, so try this
        if (LONG_EXIST == False) or (LAT_EXIST == False):
            print("Either Longitute or Latitude is available: LONG_EXIST {}; LAT_EXIST {}".format(LONG_EXIST, LAT_EXIST) )
            print("Current version may not work for WRF, WRF-Chem, and CMAQ model outputs")
            sys.exit("geo coords is not available")
        

        np_data = np.asarray(df.values)
        print(xr_data, df.attrs)

        from napari.layers import Shapes

        my_shapes = [layer for layer in self.viewer.layers if isinstance(layer, Shapes)]
        print(my_shapes)

        if "Country borders" not in str(my_shapes): 
            worldmap, world_shape_type = get_world_geojson()
            shape_layer = self.viewer.add_shapes(worldmap, name = "Country borders", shape_type=world_shape_type, edge_width = 0)
            if (LONG_EXIST == True) and (LAT_EXIST == True) and (VERT_EXIST == False ):
                print("latitude array has been inversed for napari visualization")
                set_scale_at_axis(shape_layer, axis=-2, value=-1)

        if (LONG_EXIST == True) and (LAT_EXIST == True) and (VERT_EXIST == False ):
            
            print("latitude array has been inversed for napari visualization")

            layer = self.viewer.add_image(np_data, name= xr_data)
            
            self.viewer.dims.axis_labels = ("lat", "lon") 
            layer.colormap = 'hsv'
            layer.opacity = 0.5
            layer.scale = (1, latitude_scale*-1, longitude_scale)
            layer.translate  = (0, 90, 0)
            

        if (LONG_EXIST == True) and (LAT_EXIST == True) and (VERT_EXIST == True ):
            layer = self.viewer.add_image(np_data, name= xr_data)
            self.viewer.dims.axis_labels = ("height", "lat", "lon") 
            layer.colormap = 'hsv'
            layer.opacity = 0.5
            self.viewer.dims.ndisplay = 3
            self.viewer.camera.angles = (-1.7571935971733401, -26.823353526707475, -77.73048528666025)
            layer.translate  = (0, 0, -90, 0)
            layer.scale = (1, vertical_scale, latitude_scale, longitude_scale)

def get_world_geojson():

    # Set the domain for defining the plot region.
    fname= "/Users/yunhalee/Desktop/Napari/napari-geo-worldmap/world_0_360.geojson"
    reader = napari_get_reader(fname)
    layer_data = reader(fname)

    return layer_data[0][0], layer_data[0][1]["shape_type"]

def set_scale_at_axis(layer, axis=0, value=1):

    """
    Parameters
    ----------
    layer : napari image layer 
    axis : the axis in the napari image layer that needs to be modified 
            (e.g., for data = (Z,Y,X), axis = 2 to modify Z)
    value : value assigned to the axis (e.g., -1 for flipping along the axis; 2 to increase the axis' physical scale)
    """

    curr_scale = layer.scale
    curr_scale[axis] = value
    layer.scale = curr_scale

def load_path(viewer, path: str):

    """print the file structure and create a list of data to pass to the napari.
    Parameters
    ----------
    viewer - napari viewer
    path: str
        hdf5 or netcdf file path
    -------
    """

    plugin = FileReader(viewer)

    # print the file structure and create the data list for napari
    if path.endswith(tuple(H5_FILENAME_LIST)) :
        with h5py.File(path,'r') as f:
            print("                                               ")
            print(" =======", path, " file structure start =========")

            f.visititems(plugin.print_h5_objs)

            print(" =======", path, " file structure end============")
            print("                                               ")
        # turn off plugin.all_h5_to_napari(path)
    # check if a netcdf file
    elif path.endswith(tuple(NETCDF_FILENAME_LIST)) :
        print("                                               ")
        print(" =======", path, " file structure start =========")   

        plugin.print_netcdf(path)

        print(" =======", path, " file structure end============")
        print("                                               ")
        # turn off plugin.all_netcdf_to_napari(path)
    # file is not either h5 or netcdf
    else:
        print("Failure: {} is recognized as neither hdf5 nor netcdf".format(path))
        sys.exit("Not acceptable files")

    return plugin.data_list

def data_to_napari(viewer, path: str):
    """read the data arrays and visualized with napari.
    Parameters
    ----------
    viewer: napari viewer
    path: str
        hdf5 or netcdf file path
    -------
    """
    plugin = FileReader(viewer)

    data_list = load_path(viewer, path)

    @magic_factory (dropdown={"choices": data_list})
    def make_widget_dropdown(dropdown=data_list[0]):

        varname = str(dropdown)
        print(f"varname {varname}")

        # HDF seems not to work well with Xarray (to get dimensions), so for now only simple visualization is used 
        if path.endswith(tuple(H5_FILENAME_LIST)) :
            plugin.all_h5_to_napari(path, varname)
        else: 
            plugin.get_geo_dataset(path, varname)

    viewer.window.add_dock_widget(make_widget_dropdown(), area="right") 

@magic_factory ()
def make_widget(viewer: "napari.viewer.Viewer", file_path: "pathlib.Path" = Path()):
    filename = str(file_path)
    load_path(viewer, filename)
    data_to_napari(viewer, filename)

if __name__ == "__main__":

    viewer = napari.Viewer()
    viewer.window.add_dock_widget(make_widget(), area="right")  

    napari.run()
