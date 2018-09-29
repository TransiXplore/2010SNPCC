#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import gzip
import pickle
import numpy as np

from pathlib import Path
from collections import defaultdict

from astropy.table import Table

DES_FILTERS = 'griz'

def parse_des_snana_file(path):
    """
    Parameters
    ----------
    path : str
        Path to DES light curve file

    Returns
    -------
    result : dict
        Dictionary of header and filter data for the light curve
    istarget : bool
        boolean flag to distinguish between train/test sample

    """
    # set types
    Ibc = [1,5,6,7,8,9,10,11,13,14,16,18,22,23,29,45,28]
    II = [2,3,4,12,15,17,19,20,21,24,25,26,27,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44]


    with open(path, 'r') as f:
        lin1 = f.read().splitlines()

    data1 = [elem.split() for elem in lin1]

    # get header data
    raw_data = dict([[line[0], line[1:]] for line in data1
                     if len(line) > 1 and line[0] != 'OBS:'])

    header = {}
    header['snid'] = int(raw_data['SNID:'][0])
    header['z'] = float(raw_data['SIM_REDSHIFT:'][0])

    sntype  =  int(raw_data['SIM_NON1a:'][0])
    if sntype in Ibc:   
        header['type'] = 'Ibc'
    elif sntype in II:
        header['type'] = 'II'
    elif sntype == 0:
        header['type'] = 'Ia'
    

    header['pkmjd'] = float(raw_data['SIM_PEAKMJD:'][0])
    for j, filt in enumerate(DES_FILTERS):
        header['pkmag_%s' % filt] = float(raw_data['SIM_PEAKMAG:'][j])
    header['comment'] = ' '.join(raw_data['SIM_COMMENT:'])

    istarget = int(raw_data['SNTYPE:'][0]) == -9

    light_curve = {}
    light_curve['header'] = header

    # get flux measurements
    for filt in DES_FILTERS:
        light_curve[filt] = {}
        flist = []
        for line in data1:
            if len(line) > 2 and line[0] == 'OBS:' and line[2] == filt:
                flist.append([float(line[1]), float(line[4]),
                              float(line[5]), float(line[6])])

        # extract columns from the table
        table = np.array(flist).T
        columns = ['mjd', 'fluxcal', 'fluxcalerr']
        for idx, name in enumerate(columns):
            light_curve[filt][name] = table[idx].tolist()

    return light_curve, istarget


def serialize_des_snana(path, filename):
    p = Path.cwd() / path
    lc_files = p.glob('DES_SN*.DAT')

    train = {}
    target = {}

    for lc in lc_files:
        data, istarget = parse_des_snana_file(lc)
        sn_id = data['header']['snid']
        if istarget:
            target[sn_id] = data
        else:
            train[sn_id] = data

    train_filename = filename.replace('.pkl', '_train.pkl')
    with gzip.open(train_filename, 'wb') as output:
        pickle.dump(train, output, protocol=2)
    target_filename = filename.replace('.pkl', '_target.pkl')
    with gzip.open(target_filename, 'wb') as output:
        pickle.dump(target, output, protocol=2)



if __name__ == '__main__':
    import time
    start = time.time()
    serialize_des_snana('SIMGEN_PUBLIC_DES/', 'des.pkl')
    end1 = time.time()

    end2 = time.time()
    print("DES light curves serialized in {:.1f} seconds".format(end1 - start))

