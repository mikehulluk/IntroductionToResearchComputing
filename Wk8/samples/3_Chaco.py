from enthought.traits.api \
    import HasTraits, Array, Range, Float, Enum, on_trait_change, Property
from enthought.traits.ui.api import View, Item
from enthought.chaco.chaco_plot_editor import ChacoPlotItem
from numpy import arange

class Data(HasTraits):
    volume = Array
    pressure = Property(Array, depends_on=['temperature', 'attraction',
                                       'totVolume'])
    attraction = Range(low=-50.0,high=50.0,value=0.0)
    totVolume = Range(low=.01,high=100.0,value=0.01)
    temperature = Range(low=-50.0,high=50.0,value=50.0)
    r_constant= Float(8.314472)
    plot_type = Enum("line", "scatter")

    traits_view = View(ChacoPlotItem("volume", "pressure",
                               type_trait="plot_type",
                               resizable=True,
                               x_label="Volume",
                               y_label="Pressure",
                               x_bounds=(-10,120),
                               x_auto=False,
                               y_bounds=(-2000,4000),
                               y_auto=False,
                               color="blue",
                               bgcolor="white",
                               border_visible=True,
                               border_width=1,
                               title='Pressure vs. Volume',
                               padding_bg_color="lightgray"),
                       Item(name='attraction'),
                       Item(name='totVolume'),
                       Item(name='temperature'),
                       Item(name='r_constant', style='readonly'),
                       Item(name='plot_type'),
                       resizable = True,
                       buttons = ["OK"],
                       title='Van der Waal Equation',
                       width=900, height=800)


    def _volume_default(self):
        """ Default handler for volume Trait Array. """
        return arange(.1, 100, 0.1)

    def _get_pressure(self):
        """Recalculate when one a trait the property depends on changes."""
        return ((self.r_constant*self.temperature)
              /(self.volume - self.totVolume)
             -(self.attraction/(self.volume*self.volume)))

if __name__ == '__main__':
    viewer = Data()
    viewer.configure_traits()
