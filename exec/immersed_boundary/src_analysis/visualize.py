# used to visualize data (and load plotfiles)
import yt
from yt.frontends.boxlib.data_structures import AMReXDataset
# used pretty much everywhere
import numpy as np
# defines data module
from data_model import SoA, Particle, AoS



# control yt's log level
YT_DEFAULT_LOG_LEVEL=50

def yt_set_verbosity(log_level=YT_DEFAULT_LOG_LEVEL):
    yt.funcs.mylog.setLevel(log_level)

# turn verbosity to low by default
yt_set_verbosity()




class PlotProperies(object):
    
    def __init__(self, **kwargs):

        # dict.get's second argument is the default value
        
        # immersed boundary marker presentation
        self.n_marker_skip     = kwargs.get("n_marker_skip", 1)
        self.marker_rad        = kwargs.get("marker_rad", 2e-5)
        self.inner_marker_f    = kwargs.get("inner_marker_f", 0.7)
        self.inner_marker_prop = kwargs.get("inner_marker_prop",
                                            {'circle_args': 
                                             {'fill':True, 'color':'red'}
                                             }
                                            )
        self.outer_marker_prop = kwargs.get("outer_marker_prop",
                                            {'circle_args':
                                             {'fill':True, 'color':'white'}
                                             }
                                            )

        # vector data fields
        self.x_name = kwargs.get("x_name", "avg_velx")
        self.y_name = kwargs.get("y_name", "avg_vely")
        self.z_name = kwargs.get("z_name", "avg_velz")

        # axes properties
        self.x_log   = False
        self.y_log   = False
        self.z_log   = False
        self.slc_max = kwargs.get("slc_max", 1e-4)

        # plot slicing
        self.slc_name = kwargs.get("slc_name", "avg_velz")
        self.slc_axis = kwargs.get("slc_axis", "z")
        self.quiver_x = kwargs.get("quiver_x", self.x_name)
        self.quiver_y = kwargs.get("quiver_y", self.y_name)


    def set_slc_axis(self, slc_axis: int):
        if slc_axis == 0:
            self.slc_axis = "x"
            self.quiver_x = self.y_name
            self.quiver_y = self.z_name
        elif slc_axis == 1:
            self.slc_axis = "y"
            self.quiver_x = self.z_name
            self.quiver_y = self.x_name
        elif slc_axis == 2:
            self.slc_axis = "z"
            self.quiver_x = self.x_name
            self.quiver_y = self.y_name
        else:
            RuntimeError("Not Implemented!")


    def set_log(self, plt):
        plt.set_log(self.x_name, self.x_log)
        plt.set_log(self.y_name, self.y_log)
        plt.set_log(self.z_name, self.z_log)


    def set_slc_lim(self, plt):
        if self.slc_axis == "x":
            plt.set_xlim(self.z_name, -self.slc_max, self.slc_max)
        elif self.slc_axis == "y":
            plt.set_ylim(self.z_name, -self.slc_max, self.slc_max)
        elif self.slc_axis == "z":
            plt.set_zlim(self.z_name, -self.slc_max, self.slc_max)
        else:
            RuntimeError("Not Implemented!")


    def set_axes(self, plt):
        self.set_log(plt)
        self.set_slc_lim(plt)




def slice_plt(file_name: str, prop: PlotProperies):

    # load data
    ds = yt.load(file_name)

    # plot data
    slc = yt.SlicePlot(ds, prop.slc_axis, prop.slc_name)
    slc.annotate_quiver(prop.quiver_x, prop.quiver_y)
    prop.set_axes(slc)

    return slc



def add_ibm(plt, file_name: str, prop: PlotProperies):

    # extract immersed boundary data
    amrex_ds = AMReXDataset(file_name)
    ad = amrex_ds.all_data()
    aos = AoS(ad)

    for j, part in enumerate(aos.particles[::prop.n_marker_skip]):
        p = part.pos
        p.tolist()
        plt.annotate_sphere(p, radius=prop.marker_rad, **prop.outer_marker_prop)
        inner_rad = prop.marker_rad*prop.inner_marker_f
        plt.annotate_sphere(p, radius=inner_rad, **prop.inner_marker_prop)
        #plt.annotate_marker(p, plot_args={'color':'black'})

    return plt
