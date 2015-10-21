import h5py
import Telluric_Wavecal
import os
import logging

logging.basicConfig(level=logging.INFO)

home = os.environ['HOME']
RV_HDF5_FILE='{}/School/Research/BinaryInference/RV_data.h5'.format(home)

if __name__ == '__main__':
	with h5py.File(RV_HDF5_FILE, 'a') as f:
		inst_grp = f['CHIRON']
		for starname, star_grp in inst_grp.iteritems():
			for date, date_grp in star_grp.iteritems():
				filename = date.attrs['path'].replace('~', home)

				o2_fitter = Telluric_Wavecal.VelocityFitter(filename, tell_orders=(690, 765))
                h2o_fitter = Telluric_Wavecal.VelocityFitter(filename, tell_orders=(700., 715., 725., 735.))

                logging.info('O2 Fit:\n================')
				o2_rv, o2_rv_err = o2_fitter.fit()
				logging.info('\n\nH2O Fit:\n===============')
				h2o_rv, h2o_rv_err = h2o_fitter.fit()

				date.attrs['telluric_rv_o2'] = o2_rv
				date.attrs['telluric_rv_err_o2'] = o2_rv_err
				date.attrs['telluric_rv_h2o'] = h2o_rv
				date.attrs['telluric_rv_err_h2o'] = h2o_rv_err

				f.flush()