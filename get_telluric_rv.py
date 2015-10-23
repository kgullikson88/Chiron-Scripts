import h5py
import Telluric_Wavecal
import os
import logging
import telfit

logging.basicConfig(level=logging.INFO)

home = os.environ['HOME']
RV_HDF5_FILE='{}/School/Research/BinaryInference/RV_data.h5'.format(home)

if __name__ == '__main__':
    # Make a generic telluric model
    modeler = telfit.Modeler()
    lowfreq = 1e7/(775+2)
    highfreq = 1e7/(680-2)
    telluric_model = modeler.MakeModel(lowfreq=lowfreq, highfreq=highfreq)

    with h5py.File(RV_HDF5_FILE, 'a') as f:
        inst_grp = f['CHIRON']
        for starname, star_grp in inst_grp.iteritems():
            for date, date_ds in star_grp.iteritems():
                print(date_ds.name)
                if 'telluric_rv_o2' in date_ds.attrs:
                    continue
                filename = date_ds.attrs['path'].replace('~', home)

                o2_fitter = Telluric_Wavecal.VelocityFitter(filename, tell_orders=(690, 765), telluric_model=telluric_model)
                h2o_fitter = Telluric_Wavecal.VelocityFitter(filename, tell_orders=(700., 715., 725., 735.), telluric_model=telluric_model)

                if not (o2_fitter.successful_init and h2o_fitter.successful_init):
                    logging.warn('Skipping dataset {}'.format(date_ds.name))
                    continue

                logging.info('O2 Fit:\n================')
                o2_rv, o2_rv_err = o2_fitter.fit()
                logging.info('\n\nH2O Fit:\n===============')
                h2o_rv, h2o_rv_err = h2o_fitter.fit()

                date_ds.attrs['telluric_rv_o2'] = o2_rv
                date_ds.attrs['telluric_rv_err_o2'] = o2_rv_err
                date_ds.attrs['telluric_rv_h2o'] = h2o_rv
                date_ds.attrs['telluric_rv_err_h2o'] = h2o_rv_err

                f.flush()
                