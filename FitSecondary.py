import HelperFunctions
import FittingUtilities
from astropy import units as u, constants
import StellarModel
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import sys


def LM_Model(x, vsini, rv, temperature, metal, scale, model_getter=None, **mgargs):
    if model_getter is None:
        raise KeyError("Must give model_getter keyword!")
    c = constants.c.cgs.to(u.km / u.s).value

    # Get the model from the ModelGetter instance
    mgargs['vsini'] = vsini
    model = model_getter(temperature, metal, **mgargs)
    model.y = (model.y - 1.0) * scale + 1.0

    # Finally, interpolate to the x values
    modelfcn = spline(model.x, model.y)
    return modelfcn(x * (1.0 - rv / c))


def Fit(arguments, mg=None):
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/Sorted/Stellar/Vband/"
    T_min = 3500
    T_max = 4500
    metal_min = -0.5
    metal_max = 0.5
    good_orders = range(41)
    good_orders.pop(33)

    # Eventually, override these with command line arguments

    if mg is None:
        mg = StellarModel.PhoenixGetter(modeldir,
                                        T_min=T_min,
                                        T_max=T_max,
                                        metal_min=metal_min,
                                        metal_max=metal_max,
                                        wavemin=350.0)

    fitter = HelperFunctions.ListModel(LM_Model, independent_vars=['x'], model_getter=mg)
    fitter.set_param_hint("rv", value=0.0, min=-50, max=50, vary=True)
    fitter.set_param_hint('vsini', value=20.0, min=0.0, max=500.0, vary=True)
    fitter.set_param_hint('temperature', value=4100, min=T_min, max=T_max, vary=True)
    fitter.set_param_hint('metal', value=0.1, min=metal_min, max=metal_max, vary=True)
    fitter.set_param_hint('scale', value=0.1, min=0.0, max=1.0, vary=True)


    #Eventually, make a loop
    if 1:
        filename = "test_obs_t4000.fits"
        filename = "GeneratedObservations/HIP_78106_unscaled_t4300_m0.0_rv-100_vsini10.0"
        all_orders = HelperFunctions.ReadExtensionFits(filename)
        orders = [o[1] for o in enumerate(all_orders) if o[0] in good_orders]

        optdict = {"epsfcn": 1e-2}
        params = fitter.make_params()
        result = fitter.fit(orders, fitcont=False, fit_kws=optdict, params=params)
        result.fit_report()
        print result.fit_report()


if __name__ == "__main__":
    Fit(sys.argv[1:])