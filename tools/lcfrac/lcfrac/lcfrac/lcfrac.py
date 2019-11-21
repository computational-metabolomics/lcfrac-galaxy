#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
import configparser
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import numpy as np
import shutil
import sys
import os
import re
import sqlite3
import csv
from scipy.stats import rankdata
from operator import itemgetter


from db import update_db
from spectra import process_frac_spectra
from annotation import process_frac_annotation, combine_annotations, summarise_annotations

def pths2dict(pths, galaxy):
    pthl = pths.split(',')
    pthl = [i.strip() for i in pthl]
    pthl = [i for i in pthl if i]
    d = {}

    if galaxy:
        # if galaxy was used we have the name within Galaxy and the identifier (which has the well number)
        for c, i in enumerate(pthl):
            print(c, i, c % 2)
            if (c % 2) == 0:
                bn = os.path.basename(i)
                fileid = os.path.splitext(bn)[0]
                well = os.path.splitext(bn)[0].split('_')[0]

                if well in d:
                    d[well]['full_name'] = i
                else:
                    d[well] = {'full_name': i}
            else:
                d[well]['pth'] = i

    else:
        for i in enumerate(pthl):
            well = os.path.splitext(os.path.basename(i))[0].split('_')
            if well in d:
                d[well]['full_name'] = i
                d[well]['pth'] = i
            else:
                d[well] = {'full_name': i, 'pth': i}
    return d



def extract_pth(spths, tech, dimsn_name):
    return spths[tech]['pth'][spths[tech]['full_name'] == dimsn_name]


def process_all_fracs(dims_pls,
                      dimsn_trees,
                      lcms_sqlite_pth,
                      comp_db_pth,
                      library_db_pth,
                      beams,
                      metfrag,
                      sirius,
                      spectral_matching,
                      additional,
                      frac_times_pth,
                      time_tolerance,
                      ppm_lcms,
                      ppm_dims,
                      weights,
                      ms1_lookup_source,
                      out_sqlite,
                      out_tsv,
                      rank_limit,
                      galaxy,
                      ms1_lookup_keepAdducts,
                      ms1_lookup_checkAdducts):
    # get file name in dims_pl_pth
    dims_pl_pths = pths2dict(dims_pls, galaxy)

    # get filenames in dimsn_tree_pth
    dimsn_tree_pths = pths2dict(dimsn_trees, galaxy)

    # Get BEAMS file paths
    beams_pths = pths2dict(beams, galaxy)

    # get sirius file
    sirius_pths = pths2dict(sirius, galaxy)

    # get metrag file paths
    metfrag_pths = pths2dict(metfrag, galaxy)

    # get spectral matching paths
    spectral_matching_pths = pths2dict(spectral_matching, galaxy)

    # get additional peak information
    additional_peak_info_pths = pths2dict(additional, galaxy)

    # Loop through files
    frac_spectra = {}
    for well in dims_pl_pths.keys():
        frac_spectra[well] = {}
        frac_spectra[well]['dims_data_file'] = dims_pl_pths[well]['full_name']
        frac_spectra[well]['dims_pl_pth'] = dims_pl_pths[well]['pth']
        if well in dimsn_tree_pths:
            frac_spectra[well]['dimsn_tree_pth'] = dimsn_tree_pths[well]['pth']
        else:
            frac_spectra[well]['dimsn_tree_pth'] = ''

        if well in beams_pths:
            frac_spectra[well]['beams_pth'] = beams_pths[well]['pth']
        else:
            frac_spectra[well]['beams_pth'] = ''

        if well in sirius_pths:
            frac_spectra[well]['sirius_pth'] = sirius_pths[well]['pth']
        else:
            frac_spectra[well]['sirius_pth'] = ''

        if well in metfrag_pths:
            frac_spectra[well]['metfrag_pth'] = metfrag_pths[well]['pth']
        else:
            frac_spectra[well]['metfrag_pth'] = ''

        if well in spectral_matching_pths:
            frac_spectra[well]['spectral_matching_pth'] = \
                spectral_matching_pths[well]['pth']
        else:
            frac_spectra[well]['spectral_matching_pth'] = ''

        if well in additional_peak_info_pths:
            frac_spectra[well]['additional_peak_info_pth'] = \
                additional_peak_info_pths[well]['pth']
        else:
            frac_spectra[well]['additional_peak_info_pth'] = ''

    if os.path.exists(lcms_sqlite_pth):
        if out_sqlite:
            lcms_sqlite_pth = shutil.copy(lcms_sqlite_pth, out_sqlite)
        conn, cur = update_db(lcms_sqlite_pth)
    else:
        if out_sqlite:
            # conn, cur = create_db(out_sqlite)
            print("Requires LC-MS(/MS) sqlite database")
            exit()
        else:
            conn, cur = create_db('frac_exp.sqlite')

    # connect to compound database
    comp_conn = sqlite3.connect(comp_db_pth)
    frac_times_d = {}

    with open(frac_times_pth, 'r') as ft:
        dr = csv.DictReader(ft)
        for row in dr:
            frac_times_d[row['well']] = {'frac_num': int(row['frac_num']),
                                         'frac_start_minutes':
                                             float(row['frac_start_minutes']),
                                         'frac_end_minutes':
                                             float(row['frac_end_minutes']),
                                         }

    for well, spths in frac_spectra.items():

        pid_d = process_frac_spectra(spths['dimsn_tree_pth'],
                                     spths['dims_pl_pth'],
                                     well,
                                     conn,
                                     frac_times_d,
                                     time_tolerance,
                                     ppm_lcms,
                                     ppm_dims
                                     )

        # find relevant sirius, metfrag and spectral matching for the dimsn file
        process_frac_annotation(beams_pth=spths['beams_pth'],
                                sirius_pth=spths['sirius_pth'],
                                metfrag_pth=spths['metfrag_pth'],
                                spectral_matching_pth=spths['spectral_matching_pth'],
                                additional_peak_info_pth=spths[
                                    'additional_peak_info_pth'],
                                conn=conn,
                                pid_d=pid_d,
                                ms1_lookup_source=ms1_lookup_source,
                                comp_conn=comp_conn,
                                well=well,
                                ms1_lookup_keepAdducts=ms1_lookup_keepAdducts,
                                ms1_lookup_checkAdducts=ms1_lookup_checkAdducts
                                )

        # process the MS1 annotation separately (as only 1 DIMS file)


    # combine annotations
    combine_annotations(conn, comp_conn, weights)

    # add any library spectra information (from spectral matching results
    if os.path.exists(library_db_pth):
        library_conn = sqlite3.connect(library_db_pth)
        add_library_spectra(conn, library_conn)
    else:
        print('Spectral library path does not exist so can not add library spectra to database')
    # create summary table
    # Include both the LC-MS and DI-MS annotations
    summarise_annotations(conn, out_tsv, rank_limit)




