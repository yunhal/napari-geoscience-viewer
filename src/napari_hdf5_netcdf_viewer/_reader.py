"""
This module is an example of a barebones numpy reader plugin for napari.

It implements the Reader specification, but your plugin may choose to
implement multiple readers or even other plugin contributions. see:
https://napari.org/stable/plugins/guides.html?#readers
"""
import napari
import h5py
import numpy as np
import pathlib
from pathlib import Path

from magicgui import magicgui, magic_factory
import netCDF4 as nc

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


    def h5_to_napari(self, fn):

        """
        Parameters
        ----------
        fn : netcdf file path 
        """

        label_name = 'label'
        point_name = 'annotation'

        with h5py.File(fn,'r') as f:

            for d in self.data_list:

                var_path = str(d)
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


    def netcdf_to_napari(self, fn):

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


def napari_get_reader(path):
    """read file structure and open the images with napari if the path file format is hdf or netcdf.
    Parameters
    ----------
    path: str
        hdf5 or netcdf file path
    -------
    """

    viewer = napari.Viewer()
    viewer.axes.visible = True

    plugin = FileReader(viewer)

    if path.endswith(tuple(H5_FILENAME_LIST)) :
        with h5py.File(path,'r') as f:
            print("                                               ")
            print(" =======", path, " file structure start =========")

            f.visititems(plugin.print_h5_objs)

            print(" =======", path, " file structure end============")
            print("                                               ")
        plugin.h5_to_napari(path)
    
    # check if a netcdf file
    elif path.endswith(tuple(NETCDF_FILENAME_LIST)) :
        print("                                               ")
        print(" =======", path, " file structure start =========")   

        plugin.print_netcdf(path)

        print(" =======", path, " file structure end============")
        print("                                               ")

        plugin.netcdf_to_napari(path)
    
    # file is not either h5 or netcdf
    else:
        print("Failure: {} is recognized as neither hdf5 nor netcdf".format(path))

@magic_factory
def make_widget(file_path: "pathlib.Path" = Path()):
    filename = str(file_path)
    load_path(filename)

