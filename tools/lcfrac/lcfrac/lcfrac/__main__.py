#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
import configparser
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import os
import sqlite3
import shutil

from lcfrac import LcFracExp

def config_strip(pth):
    pth = pth.strip('"')
    return pth.strip("'")


if __name__ == "__main__":
    print('LC fractionation tool')
    parser = ArgumentParser(description='Combine spectra and '
                                        'annotation results '
                                        'from and LC-MS/MS '
                                        'fractionation '
                                        'experiment',
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-l', '--lcms_sqlite', type=str, required=True,
                        help="Path to the LC-MS sqlite results")
    parser.add_argument('-d', '--dims_pls', type=str, required=True,
                        help="Path to the directory of dims peaklists (hdf5)")
    parser.add_argument('-n', '--dimsn_merged', type=str, required=True,
                        help="Path to the directory of merged dimsn peaklists (hdf5)")
    parser.add_argument('-k', '--dimsn_non_merged', type=str, required=True,
                        help="Path to the directory of non merged dimsn peaklists (hdf5)")
    parser.add_argument('-p', '--dimsn_ms1_precursors', type=str, required=True,
                        help="Path to the directory of dimsn MS1 precursors peaklists (hdf5)")
    parser.add_argument('-b', '--beams', type=str, required=True,
                        help="Path to the directory of BEAMS results")
    parser.add_argument('-m', '--metfrag', type=str, required=True,
                        help="Path to the directory of MetFrag results")
    parser.add_argument('-s', '--sirius', type=str, required=True,
                        help="Path to the directory of Sirius CSI:FingerID results")
    parser.add_argument('-x', '--spectral_matching', type=str, required=True,
                        help="Path to the directory of the spectral "
                             "matching results")
    parser.add_argument('-g', '--mf_annotation', type=str, required=True,
                        help="Molecular formula annotation from MSnPy")
    parser.add_argument('-f', '--frac_times_pth', type=str, required=False,
                        help="Path to the mapping of the fraction to LC times")
    parser.add_argument('-a', '--additional_info', default=100,
                        help="Additional information for each peak (e.g. precursor ion purity)")
    parser.add_argument('-o', '--out_sqlite', type=str, required=False,
                        help="Out path for sqlite database")
    parser.add_argument('--time_tolerance', type=float, default=10,
                        help="+/- this time in seconds to LC-MS peaks when searching for associated dims spectra")
    parser.add_argument('--lcms_ppm', type=float, default=5,
                        help="ppm accuracy of LC-MS experiment")
    parser.add_argument('--dims_ppm', type=float, default=5,
                        help="ppm accuracy of DI-MS experiment")
    parser.add_argument('--weight_sirius_csifingerid', type=float, default=0.25,
                        help="Weight for sirius_csifingerid annotation")
    parser.add_argument('--weight_metfrag', type=float, default=0.15,
                        help="Weight for metfrag annotation")
    parser.add_argument('--weight_beams', type=float, default=0.1,
                        help="Weight for beams annotation")
    parser.add_argument('--weight_spectral_matching', type=float, default=0.45,
                        help="Weight for spectral matching")
    parser.add_argument('--weight_biosim', type=float, default=0.05,
                        help="Weight for biosim score")
    parser.add_argument('--ms1_lookup_source', type=str, default='hmdb',
                        help="Database used for MS1 lookup source (e.g. hmdb, kegg or pubchem)")
    parser.add_argument('--galaxy', action='store_true', help="Flag if running from Galaxy")
    parser.add_argument('--out_tsv', type=str, default='lcfrac_results.tsv',
                        help="Out path for summary of results (tsv)")
    parser.add_argument('--rank_limit', default=0,
                        help="Limit the number of annotation by rank shown in the summary tsv file")
    parser.add_argument('--sirius_rank_limit', default=25, type=int,
                        help="Limit the annotations results from SIRIUS CSI:FingerID")
    parser.add_argument('--ms1_lookup_keepAdducts', default='',
                        help="Provide a list of adducts that should be used from the MS1 lookup (e.g. [M+H]+, [M+Na]+)")
    parser.add_argument('--ms1_lookup_checkAdducts', action='store_true',
                        help="""Check if adducts match to those found in CAMERA (adducts to be present in the additional
                                 peak information file)""")
    parser.add_argument('--comp_db_pth', default="", type=str, required=False,
                        help="Path to the compound sqlite database")
    parser.add_argument('--library_db_pth', default="", type=str, required=False,
                        help="Path to the spectral library sqlite database")

    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.ini'))

    if args.comp_db_pth:
        comp_db_pth = args.comp_db_pth
    else:
        comp_db_pth = config['db']['comp_db_pth']
        comp_db_pth = config_strip(comp_db_pth)

    if args.library_db_pth:
        library_db_pth = args.library_db_pth
    else:
        library_db_pth = config['db']['library_db_pth']
        library_db_pth = config_strip(library_db_pth)

    print(comp_db_pth)
    print(library_db_pth)

    ms1_lookup_keepAdducts = args.ms1_lookup_keepAdducts.replace('__cb__', ']')
    ms1_lookup_keepAdducts = ms1_lookup_keepAdducts.replace('__ob__', '[')

    print(ms1_lookup_keepAdducts)

    weights = {'sirius_csifingerid': args.weight_sirius_csifingerid,
               'metfrag': args.weight_metfrag,
               'biosim': args.weight_biosim,
               'spectral_matching': args.weight_spectral_matching,
               'beams': args.weight_beams
            }
    print(weights)
    if args.out_sqlite:
        shutil.copy(args.lcms_sqlite, args.out_sqlite)
        out_sqlite = args.out_sqlite
    else:
        out_sqlite = lcms_sqlite

    lc_frac_exp = LcFracExp(time_tolerance=args.time_tolerance,
                            dims_ppm=args.dims_ppm,
                            lcms_ppm=args.lcms_ppm,
                            rank_limit=int(args.rank_limit),
                            ms1_lookup_source=args.ms1_lookup_source,
                            ms1_lookup_keepAdducts=ms1_lookup_keepAdducts,
                            ms1_lookup_checkAdducts=args.ms1_lookup_checkAdducts,
                            weights=weights,
                            galaxy=args.galaxy,
                            lcms_db_pth=out_sqlite,
                            metab_compound_db_pth=comp_db_pth,
                            library_spectra_db_pth=library_db_pth,
                            frac_times_pth=args.frac_times_pth,
                            sirius_rank_limit=int(args.sirius_rank_limit))

    # Get connection to relevant databases
    lc_frac_exp.connect2dbs()

    # Assign data paths to wells
    lc_frac_exp.assign_data_pths_to_wells(args.dims_pls,
                                          args.dimsn_merged,
                                          args.dimsn_non_merged,
                                          args.dimsn_ms1_precursors,
                                          args.spectral_matching,
                                          args.metfrag,
                                          args.sirius,
                                          args.beams,
                                          args.additional_info,
                                          args.mf_annotation)

    # Update the lcms database for the fractionation results
    print('########### Update database schema ###############################')
    lc_frac_exp.update_lcms_db()

    # Add the data for each wells to database
    print('########### Add wells to db        ###############################')
    lc_frac_exp.add_wells_to_db()

    # Combine annotations
    print('########### Combine annotations     ###############################')
    lc_frac_exp.combine_annotations()

    # Add library spectra
    if os.path.exists(args.library_db_pth):
        lc_frac_exp.add_library_spectra()

    # summarise spectra
    lc_frac_exp.summarise_annotation(out_tsv_pth=args.out_tsv)

