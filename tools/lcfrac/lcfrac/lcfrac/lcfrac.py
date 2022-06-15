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

from db import update_lcms_db
from spectra import process_ms1_spectra, process_frag_spectra, add_library_spectra, lcms_dims_link
from annotation import process_ms1_annotation, process_frag_annotation, combine_annotations, summarise_annotations
from collections import defaultdict

def dims_well_check(msnpths, dims_wells):
    return(' '.join([i for i in msnpths.split(' ') if i.strip()[0:3] in dims_wells]))

class LcFracExp(object):
    """Class to handle LC-MS/MS fractionation spectra and annotations
    """
    def __init__(self,
                 time_tolerance=10,
                 dims_ppm=5,
                 lcms_ppm=5,
                 rank_limit=0,
                 ms1_lookup_source="hmdb",
                 ms1_lookup_keepAdducts="[M+H]+",
                 ms1_lookup_checkAdducts=True,
                 weights={'sirius_csifingerid': 0.2,
                          'metfrag': 0.2,
                          'biosim': 0.25,
                          'spectral_matching': 0.3,
                          'beams': 0.05},
                galaxy = False,
                lcms_db_pth='',
                metab_compound_db_pth='',
                library_spectra_db_pth='',
                frac_times_pth='',
                sirius_rank_limit=25
                 ):
        self.wells = {}
        self.lcms_db_pth = lcms_db_pth
        self.metab_compound_db_pth = metab_compound_db_pth
        self.library_spectra_db_pth = library_spectra_db_pth
        self.galaxy = galaxy
        self.time_tolerance = time_tolerance
        self.dims_ppm = dims_ppm
        self.lcms_ppm = lcms_ppm
        self.rank_limit = rank_limit
        self.ms1_lookup_source = ms1_lookup_source
        self.ms1_lookup_keepAdducts = ms1_lookup_keepAdducts
        self.ms1_lookup_checkAdducts = ms1_lookup_checkAdducts
        self.weights = weights
        self.sirius_rank_limit = sirius_rank_limit

        if frac_times_pth:
            self.frac_times_pth = frac_times_pth
        else:
            self.frac_times_pth = 'frac_times.csv'
        # create multiple well objects for the LC fractionation experiment
        self.init_wells()

    def init_wells(self):
        with open(self.frac_times_pth, 'r') as ft:
            dr = csv.DictReader(ft)
            for row in dr:
                rtmin = float(row['frac_start_minutes'])
                rtmax = float(row['frac_end_minutes'])
                self.wells[row['well']] = WellInfo(well_number=row['well'],
                                                   frac_number=int(row['frac_num']),
                                                   rtmin=rtmin,
                                                   rtmax=rtmax)

    def connect2dbs(self):
        if self.lcms_db_pth:
            print('connecting to lcms sqlite database {}'.format(self.lcms_db_pth))
            self.lcms_conn = sqlite3.connect(self.lcms_db_pth)
        if self.metab_compound_db_pth:
            print('connecting to metab compound sqlite database {}'.format(self.metab_compound_db_pth))
            self.comp_conn = sqlite3.connect(self.metab_compound_db_pth)
        if self.library_spectra_db_pth:
            print('connecting to library spectra sqlite database {}'.format(self.library_spectra_db_pth))
            self.library_conn = sqlite3.connect(self.library_spectra_db_pth)


    def update_lcms_db(self):
        update_lcms_db(self.lcms_conn)


    def pths2wells(self, pths, dims_msn, tech):
        pthl = pths.split(',')
        pthl = [i.strip() for i in pthl]
        pthl = [i.strip('\n') for i in pthl if i]

        if self.galaxy:
            # if galaxy was used we have the name within Galaxy and the identifier (which has the well number)
            for c, i in enumerate(pthl):

                if (c % 2) == 0:
                    bn = os.path.basename(i)
                    element = os.path.splitext(bn)[0]
                    well = os.path.splitext(bn)[0].split('_')[0]
                    if well in self.wells:

                        if dims_msn == 'dims':
                            if not self.wells[well].dims_pths.element:
                                setattr(self.wells[well].dims_pths, 'element', element)
                            elif self.wells[well].dims_pths.element != element:
                                if os.path.splitext(element)[0].split('_')[0] == self.wells[well].dims_pths.element:
                                    setattr(self.wells[well].dims_pths, 'element', os.path.splitext(element)[0].split('_')[0])
                                else:
                                    #print(well, 'ELEMENT1', element, 'ELEMENT2')
                                    print(self.wells.keys())
                                    print('Element identifiers for DIMS data do not match')
                        elif dims_msn == 'msn':
                            if not self.wells[well].msn_pths_c[element]:
                                self.wells[well].msn_pths_c[element] = MsnPths(element=element)
                else:
                    if well in self.wells:

                        if dims_msn == 'dims':
                            setattr(self.wells[well].dims_pths, tech, i)
                        elif dims_msn == 'msn':
                            setattr(self.wells[well].msn_pths_c[element], tech, i)
                        setattr(self.wells[well], 'data_assigned', True)
        else:
            for i in enumerate(pthl):
                well = os.path.splitext(os.path.basename(i))[0].split('_')
                if well in self.wells:
                    if dims_msn == 'dims':
                        setattr(self.wells[well].dims_pths, 'element', i)
                        setattr(self.wells[well].dims_pths, tech, i)
                    elif dims_msn == 'msn':
                        setattr(self.wells[well].msn_pths_c[element], tech, i)
                    setattr(self.wells[well], 'data_assigned', True)

    def assign_data_pths_to_wells(self,
                                  dims_pls_pths,
                                  dimsn_merged_pths,
                                  dimsn_non_merged_pths,
                                  dimsn_ms1_precursors_pths,
                                  spectral_matching_pths,
                                  metfrag_pths,
                                  sirius_pths,
                                  beams_pths,
                                  additional_info_pths,
                                  mf_annotation_pths
                                  ):
        dims_wells = [i.split(',')[0] for i in dims_pls_pths.split(' ')]

        #print('DIMS spectra')
        self.pths2wells(dims_pls_pths, 'dims', 'spectra')
        #print('MSn spectra (merged)')
        self.pths2wells(dims_well_check(dimsn_merged_pths, dims_wells), 'msn', 'merged_spectra')
        #print('MSn spectra (non merged)')
        self.pths2wells(dims_well_check(dimsn_non_merged_pths, dims_wells), 'msn', 'non_merged_spectra')
        #print('MSn spectra (MS1 precursors)')
        self.pths2wells(dims_well_check(dimsn_ms1_precursors_pths, dims_wells), 'msn', 'ms1_precursor_spectra')
        #print('MS1 annotation')
        self.pths2wells(beams_pths, 'dims', 'beams')
        #print('MS1 additional info')
        self.pths2wells(additional_info_pths, 'dims', 'additional_info')
        #print('MSn spectral matching')
        self.pths2wells(dims_well_check(spectral_matching_pths, dims_wells), 'msn', 'spectral_matching')
        #print('MSn metfrag')
        self.pths2wells(dims_well_check(metfrag_pths, dims_wells), 'msn', 'metfrag')
        #print('MSn sirius')
        self.pths2wells(dims_well_check(sirius_pths, dims_wells), 'msn', 'sirius')
        # print('MSn MF msnpy')
        self.pths2wells(dims_well_check(mf_annotation_pths, dims_wells), 'msn', 'mf_annotation')

    def add_well_to_db(self, wellinfo):
        print(wellinfo.well_number)
        #print('process ms1 spectra')
        orig_dims_d, orig_dims_pid = process_ms1_spectra(wellinfo, self.lcms_conn)
        #print('process ms1 annotation')
        process_ms1_annotation(wellinfo,
                               self.lcms_conn,
                               self.comp_conn,
                               self.ms1_lookup_source,
                               self.ms1_lookup_keepAdducts,
                               self.ms1_lookup_checkAdducts)

        #print('process frag spectra')
        for element_i, msn_pths in wellinfo.msn_pths_c.items():
            #print(element_i, msn_pths)
            pid_d, msn_prec_d = process_frag_spectra(wellinfo,
                                         element_i,
                                         self.lcms_conn,
                                         self.dims_ppm,
                                         orig_dims_d)

        # add LC-MS to DI-MS link
        #print('Link LC-MS to DI-MS')
        lcms_dims_link(self.lcms_conn, orig_dims_pid, self.time_tolerance,
                       self.lcms_ppm, self.dims_ppm, wellinfo.rtmin,
                       wellinfo.rtmax)

        #print('Process frag annotation')
        for element_i, msn_pths in wellinfo.msn_pths_c.items():
            process_frag_annotation(msn_pths,
                                    lcms_conn=self.lcms_conn,
                                    pid_d=pid_d,
                                    msn_prec_d=msn_prec_d,
                                    sirius_rank_limit=self.sirius_rank_limit)



    def add_wells_to_db(self):
        for wellid, wellinfo in self.wells.items():
            if wellinfo.data_assigned:
                self.add_well_to_db(wellinfo)

    def combine_annotations(self):
        combine_annotations(self.lcms_conn, self.comp_conn, self.weights)

    def add_library_spectra(self):
        if self.library_spectra_db_pth:
            add_library_spectra(self.lcms_conn, self.library_conn)
        else:
            print('Library spectra database path needs to be defined')

    def summarise_annotation(self, out_tsv_pth):
        # create summary table
        # Include both the LC-MS and DI-MS annotations
        summarise_annotations(self.lcms_conn, out_tsv_pth, self.rank_limit)


class WellInfo(object):
    """Class to handle LC-MS/MS wells of LcFrac object
    """
    def __init__(self, well_number, frac_number, rtmin, rtmax):
        self.well_number = well_number
        self.frac_number = frac_number
        self.rt = np.mean([rtmin, rtmax]),
        self.rtmin = rtmin
        self.rtmax = rtmax
        self.dims_pths = DimsPths()
        self.msn_pths_c = defaultdict(str)
        self.data_assigned = False


class MsnPths(object):
    def __init__(self, element=''):
        self.element = element
        self.non_merged_spectra = ''
        self.merged_spectra = ''
        self.ms1_precursor_spectra = ''
        self.sirius = ''
        self.metfrag = ''
        self.spectral_matching = ''
        self.mf_annotation = ''


class DimsPths(object):
    def __init__(self, element=''):
        self.element = element
        self.spectra = ''
        self.beams = ''
        self.additional_info = ''


# if __name__ == "__main__":
#
#     # dims_pls_pths = """
#     # A01,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/dims_pls/A01.h5,
#     # A02,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/dims_pls/A02.h5
#     # """
#     #
#     # additional_info_pths = """
#     # A01,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/additional/A01.tsv,
#     # A02,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/additional/A02.tsv
#     # """
#     #
#     # beams_pths = """
#     # A01,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/beams/A01.tsv,
#     # A02,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/beams/A02.tsv
#     # """
#     #
#     # dimsn_tree_pths = """
#     # A01,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/dimsn_trees/A01.json,
#     # A02_1,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/dimsn_trees/A02_1.json,
#     # A02_2,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/dimsn_trees/A02_2.json
#     # """
#     #
#     # metfrag_pths = """
#     # A01,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/metfrag/A01.tabular,
#     # A02_1,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/metfrag/A02_1.tabular
#     # """
#     #
#     # sirius_pths = """
#     # A01,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/sirius/A01.tsv,
#     # A02_1,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/sirius/A02_1.tsv,
#     # A02_2,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/sirius/A02_2.tsv
#     # """
#     #
#     # spectral_matching_pths = """
#     # A01,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/spectral_matching/A01.tsv,
#     # A02_1,/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/spectral_matching/A02_1.tsv
#     # """
#     #
#     # lcms_sqlite = "/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/lcmsms.sqlite"
#     # metab_compound_sqlite = "/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data" \
#     #                         "/metab_compound_subset.sqlite"
#     #
#     # lcms_sqlite_test = "/home/tomnl/Dropbox/code/post-doc-brum/lcfrac_new/tools/lcfrac/test-data/lcmsms_TEST.sqlite"
#     # shutil.copy(lcms_sqlite, lcms_sqlite_test)
#     # # Setup the LCFractionation experiment object
#     # lc_frac_exp = LcFracExp(galaxy=True, lcms_db_pth=lcms_sqlite_test, metab_compound_db_pth=metab_compound_sqlite)
#     #
#     # # Get connection to relevant databases
#     # lc_frac_exp.connect2dbs()
#     #
#     # # Assign data paths to wells
#     # lc_frac_exp.assign_data_pths_to_wells(dims_pls_pths,
#     #                                       dimsn_tree_pths,
#     #                                       spectral_matching_pths,
#     #                                       metfrag_pths,
#     #                                       sirius_pths,
#     #                                       beams_pths,
#     #                                       additional_info_pths)
#     #
#     # # Update the lcms database for the fractionation results
#     # lc_frac_exp.update_lcms_db()
#     #
#     # # Add the data for each wells to database
#     # lc_frac_exp.add_wells_to_db()
#     #
#     # # Combine annotations
#     # lc_frac_exp.combine_annotations()
#     #
#     # # Add library spectra
#     # lc_frac_exp.add_library_spectra()
#     #
#     # # summarise spectra
#     # lc_frac_exp.summarise_annotation(out_tsv_pth='test.tsv')
