import napari
from napari_hdf5_netcdf_viewer import make_widget


def test_hdf5_netcdf_viewer():
    """
    test to read file structure and open the images with napari if the path file format is hdf or netcdf.

    Parameters
    ----------
    fpath: str
        hdf5 or netcdf file path
    -------
    """

    viewer = napari.Viewer()
    viewer.window.add_dock_widget(make_widget(), area="right")  


