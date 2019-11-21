#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
import configparser
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import os
import sqlite3

from lcfrac import process_all_fracs



def config_strip(pth):
    pth = pth.strip('"')
    return pth.strip("'")


if __name__ == "__main__":
    parser = ArgumentParser(description='Combine spectra and '
                                        'annotation results '
                                        'from and LC-MS/MS '
                                        'fractionation '
                                        'experiment',
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-l', '--lcms_sqlite',
                        type=str, required=True,
                        help="Path to the LC-MS sqlite results")
    parser.add_argument('-d', '--dims_pls',
                        type=str, required=True,
                        help="Path to the directory of dims peaklists (hdf5)")
    parser.add_argument('-n', '--dimsn_trees',
                        type=str, required=True,
                        help="Path to the directory of dimsn peaklists (hdf5)")
    parser.add_argument('-b', '--beams',
                        type=str, required=True,
                        help="Path to the directory of BEAMS results")
    parser.add_argument('-m', '--metfrag',
                        type=str, required=True,
                        help="Path to the directory of MetFrag results")
    parser.add_argument('-s', '--sirius',
                        type=str, required=True,
                        help="Path to the directory of Sirius " \
                             "CSI:FingerID results")
    parser.add_argument('-x', '--spectral_matching',
                        type=str, required=True,
                        help="Path to the directory of the spectral "
                             "matching results")
    parser.add_argument('-f', '--frac_times_pth',
                        type=str, required=False,
                        help="Path to the mapping of the fraction to LC"
                             "times")
    parser.add_argument('--time_tolerance',
                        type=float, default=10,
                        help="+/- this time in seconds to LC-MS peaks"
                             "when searching for associated dims spectra")
    parser.add_argument('--ppm_lcms',
                        type=float, default=5,
                        help="ppm accuracy of LC-MS experiment")
    parser.add_argument('--ppm_dims',
                        type=float, default=5,
                        help="ppm accuracy of DI-MS experiment")
    parser.add_argument('--weight_sirius_csifingerid',
                        type=float, default=0.2,
                        help="Weight for sirius_csifingerid annotation")
    parser.add_argument('--weight_metfrag',
                        type=float, default=0.2,
                        help="Weight for metfrag annotation")
    parser.add_argument('--weight_beams',
                        type=float, default=0.05,
                        help="Weight for beams annotation")
    parser.add_argument('--weight_spectral_matching',
                        type=float, default=0.3,
                        help="Weight for spectral matching")
    parser.add_argument('--weight_biosim',
                        type=float, default=0.25,
                        help="Weight for biosim score")
    parser.add_argument('--ms1_lookup_source',
                        type=str, default='hmdb',
                        help="Database used for MS1 lookup source (e.g."
                             "hmdb, kegg or pubchem)")
    parser.add_argument('-o', '--out_sqlite',
                        type=str, required=False,
                        help="Out path for sqlite database")
    parser.add_argument('-g', '--galaxy', action='store_true',
                        help="Flag if running from Galaxy")
    parser.add_argument('-y', '--out_tsv',
                        type=str, default='lcfrac_results.tsv',
                        help="Out path for summary of results (tsv)")
    parser.add_argument('-r', '--rank_limit',
                        default=100,
                        help="Limit the number of annotation by rank shown in "
                             "the summary tsv file")
    parser.add_argument('-a', '--additional',
                        default=100,
                        help="Additional information for each peak (e.g. precursor ion purity)")
    parser.add_argument('--ms1_lookup_keepAdducts',
                        default='',
                        help="Provide a list of adducts that should be used from the MS1 lookup (e.g. [M+H]+, [M+Na]+)")
    parser.add_argument('--ms1_lookup_checkAdducts', action='store_true',
                        help="""Check if adducts match to those found in CAMERA (adducts to be present in the additional
                     peak information file)""")
    parser.add_argument('--comp_db_pth', default="",
                        type=str, required=False,
                        help="Path to the compound sqlite database")
    parser.add_argument('--library_db_pth', default="",
                        type=str, required=False,
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

    weights = {'sirius_csifingerid': args.weight_sirius_csifingerid,
               'metfrag': args.weight_metfrag,
               'biosim': args.weight_biosim,
               'spectral_matching': args.weight_spectral_matching,
               'beams': args.weight_beams}

    process_all_fracs(args.dims_pls,
                      args.dimsn_trees,
                      args.lcms_sqlite,
                      comp_db_pth,
                      library_db_pth,
                      args.beams,
                      args.metfrag,
                      args.sirius,
                      args.spectral_matching,
                      args.additional,
                      args.frac_times_pth,
                      args.time_tolerance,
                      args.ppm_lcms,
                      args.ppm_dims,
                      weights,
                      args.ms1_lookup_source,
                      args.out_sqlite,
                      args.out_tsv,
                      args.rank_limit,
                      args.galaxy,
                      ms1_lookup_keepAdducts,
                      args.ms1_lookup_checkAdducts)