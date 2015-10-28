""" Measure the primary star RV from CCFS.
"""
import sys
from collections import defaultdict

import h5py
import pandas as pd


ADDMODE = 'ml'


def measure_rv(hdf5_file, output_log=None):
    """ Measure the RV of each star/date combo in the HDF5 file
    """
    with h5py.File(hdf5_file, 'r') as f:
        df_list = []
        i = 0
        for starname, star_grp in f.iteritems():
            for date, date_grp in star_grp.iteritems():
                print(date_grp.name)
                summary = defaultdict(list)
                for ds_name, dataset in date_grp.iteritems():
                    if dataset.attrs['addmode'] != ADDMODE:
                        continue
                    summary['teff'].append(dataset.attrs['T'])
                    summary['vsini'].append(dataset.attrs['vsini'])
                    summary['logg'].append(dataset.attrs['logg'])
                    summary['feh'].append(dataset.attrs['[Fe/H]'])
                    summary['addmode'].append(dataset.attrs['addmode'])
                    summary['RV'].append(dataset.attrs['vel_max'])
                    summary['CCF'].append(dataset.attrs['ccf_max'])
                df = pd.DataFrame(data=summary)
                df['star'] = starname
                df['date'] = date
                df_list.append(df)

                # Save the maximum value, if requested
                if output_log is not None:
                    best = df.sort_values(by='CCF').tail(1)
                    best = {k: v for k, v in zip(best.columns, best.values[0])}
                    with open(output_log, 'a') as log:
                        log.write('{star},{date},{teff},{logg},{feh},{vsini},{addmode},{RV},{CCF}\n'.format(**best))

    return pd.concat(df_list, ignore_index=True)


if __name__ == '__main__':
    for fname in sys.argv[1:]:
        output = fname.replace('hdf5', 'rv.txt')
        with open(output, 'w'):
            output.write('star,date,teff,logg,feh,vsini,addmode,rv,ccf\n')
        summary = measure_rv(fname, output_log=output)
        summary.to_csv(output.replace('.txt', '_summary.txt'), index=False)