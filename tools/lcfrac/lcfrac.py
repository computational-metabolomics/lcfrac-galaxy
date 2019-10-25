#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import numpy as np
import shutil
import sys
import os
import sqlite3
import csv
from scipy.stats import rankdata
from operator import itemgetter
from msnpy.convert import tree2peaklist
from dimspy.portals.hdf5_portal import load_peak_matrix_from_hdf5
from dimspy.process.peak_alignment import align_peaks


def update_db(file_pth):
    """
    """
    conn = sqlite3.connect(file_pth)
    c = conn.cursor()
    r = c.execute('PRAGMA table_info(s_peak_meta)')
    cols = [i[1] for i in r.fetchall()]

    if 'msnpy_convert_id' not in cols:
        # this means that we are missing the dims specific columns
        c.execute('ALTER TABLE s_peak_meta ADD msnpy_convert_id integer')
        c.execute('ALTER TABLE s_peak_meta ADD well text')
        c.execute('ALTER TABLE s_peak_meta ADD well_rtmin real')
        c.execute('ALTER TABLE s_peak_meta ADD well_rtmax real')
        c.execute('ALTER TABLE s_peak_meta ADD well_rt real')

        c.execute('''CREATE TABLE c_peak_groups_X_s_peaks(
                                        gxsid integer PRIMARY KEY,
                                        grpid integer,
                                        sid integer,
                                        mzdiff real,
                                        FOREIGN KEY(grpid) REFERENCES 
                                        c_peak_groups(grpid),
                                        FOREIGN KEY(sid) REFERENCES 
                                        s_peaks(sid) 
                                   )'''
                  )

        c.execute('''CREATE TABLE s_peaks_X_s_peaks(
                                        sxsid integer PRIMARY KEY,
                                        sid1 integer,
                                        sid2 integer,
                                        mzdiff real,
                                        link_type text,
                                        FOREIGN KEY(sid1) REFERENCES
                                        s_peaks(sid), 
                                        FOREIGN KEY(sid2) REFERENCES 
                                        s_peaks(sid)
                                   )'''
                  )

        c.execute('''CREATE TABLE s_peak_meta_X_s_peaks(
                                        smxsid integer PRIMARY KEY,
                                        sid integer,
                                        pid integer,
                                        mzdiff real,
                                        link_type text,
                                        FOREIGN KEY(sid) REFERENCES 
                                        s_peaks(sid), 
                                        FOREIGN KEY(pid) REFERENCES 
                                        s_peak_meta(pid)
                                   )'''
                  )

        r = c.execute('PRAGMA table_info(ms1_lookup_results)')
        cols = r.fetchall()
        if not cols:
            c.execute('''CREATE TABLE ms1_lookup_results(
                                            id integer PRIMARY KEY,
                                            grpid integer,
                                            pid integer,
                                            msnpy_convert_id	integer,
                                            well text,
                                            name text,	
                                            mz real,
                                            rt	real,
                                            intensity real,
                                            exact_mass real,
                                            ppm_error real,	
                                            adduct text,
                                            C integer,
                                            H integer,
                                            N integer,
                                            O integer,
                                            P integer,	
                                            S integer,
                                            molecular_formula text,	
                                            compound_name text,
                                            compound_id text,
                                            inchikey text,
                                            score real,
                                            FOREIGN KEY(grpid) REFERENCES 
                                                c_peak_groups(grpid), 
                                            FOREIGN KEY(pid) REFERENCES 
                                                s_peak_meta(pid)
                                       )'''
                      )
        elif 'msnpy_convert_id' not in cols:
            c.execute('ALTER TABLE ms1_lookup_results ADD id integer')
            c.execute(
                'ALTER TABLE ms1_lookup_results ADD msnpy_convert_id integer')
            c.execute('ALTER TABLE ms1_lookup_results ADD pid integer')
            c.execute('ALTER TABLE ms1_lookup_results ADD well text')
            c.execute('ALTER TABLE ms1_lookup_results ADD inchikey text')

    r = c.execute('PRAGMA table_info(sirius_csifingerid_results)')
    cols = r.fetchall()
    if not cols:
        c.execute('''CREATE TABLE sirius_csifingerid_results(
                                        id integer PRIMARY KEY,
                                        grpid integer,
                                        pid integer,
                                        msnpy_convert_id	integer,
                                        adduct text,
                                        inchikey2D text,
                                        InChI text,
                                        molecularFormula text,
                                        Rank integer,	
                                        Score real,
                                        Name text,
                                        smiles text,
                                        xlogp real,	
                                        pubchemids text,
                                        links text,
                                        bounded_score real,
                                        FOREIGN KEY(grpid) REFERENCES 
                                            c_peak_groups(grpid), 
                                        FOREIGN KEY(pid) REFERENCES 
                                            s_peak_meta(pid)
                                   )'''
                  )
    elif 'msnpy_convert_id' not in cols:
        c.execute('ALTER TABLE sirius_csifingerid_results ADD id integer')
        c.execute('ALTER TABLE sirius_csifingerid_results ADD '
                  'msnpy_convert_id '
                  'integer')
        c.execute('ALTER TABLE sirius_csifingerid_results ADD pid integer')
        c.execute('ALTER TABLE sirius_csifingerid_results ADD well text')
        c.execute('ALTER TABLE sirius_csifingerid_results ADD mz real')

    r = c.execute('PRAGMA table_info(metfrag_results)')
    cols = r.fetchall()
    if not cols:
        c.execute('''CREATE TABLE metfrag_results(
                                        id integer PRIMARY KEY,
                                        grpid integer,
                                        pid integer,
                                        msnpy_convert_id	integer,
                                        adduct text,
                                        CompoundName text,	
                                        ExplPeaks text,
                                        FormulasOfExplPeaks text,	
                                        FragmenterScore	text,
                                        FragmenterScore_Values text,	
                                        Identifier text,
                                        InChI text,	
                                        InChIKey text,	
                                        InChIKey1 text,
                                        InChIKey2 text,	
                                        InChIKey3 text,
                                        MaximumTreeDepth text,	
                                        MolecularFormula text,	
                                        MonoisotopicMass real,
                                        Name text,	
                                        NoExplPeaks real,
                                        NumberPeaksUsed real,	
                                        SMILES text,
                                        Score real,
                                        FOREIGN KEY(grpid) REFERENCES 
                                            c_peak_groups(grpid), 
                                        FOREIGN KEY(pid) REFERENCES 
                                            s_peak_meta(pid)
                                   )'''
                  )
    elif 'msnpy_convert_id' not in cols:
        c.execute('ALTER TABLE metfrag_results ADD id integer')
        c.execute('ALTER TABLE metfrag_results ADD msnpy_convert_id integer')
        c.execute('ALTER TABLE metfrag_results ADD pid integer')
        c.execute('ALTER TABLE metfrag_results ADD well text')
        c.execute('ALTER TABLE metfrag_results ADD mz real')

    r = c.execute('PRAGMA table_info(sm_matches)')
    cols = r.fetchall()
    if not cols:
        c.execute('''CREATE TABLE sm_matches(
                                        msnpy_convert_id integer,
                                        mid integer,
                                        lpid integer,
                                        qpid integer,
                                        dpc real,	
                                        rdpc real,
                                        cdpc real,
                                        mcount integer,
                                        allcount integer,	
                                        mpercent real,
                                        library_rt real,
                                        query_rt real,	
                                        rtdiff real,
                                        library_precursor_mz real,	
                                        query_precursor_mz real,	
                                        library_precursor_ion_purity real,	
                                        query_precursor_ion_purity real,	
                                        library_accession text,	
                                        library_precursor_type text,	
                                        library_entry_name text,
                                        inchikey text,	
                                        library_source_name text,	
                                        library_compound_name text

                                   )'''
                  )
    elif 'msnpy_convert_id' not in cols:
        c.execute('ALTER TABLE sm_matches ADD msnpy_convert_id integer')

    r = c.execute('PRAGMA table_info(combined_annotations)')
    cols = r.fetchall()

    if not cols:
        c.execute('''CREATE TABLE combined_annotations(
                                        sid integer,
                                        grpid integer,
                                        inchikey text,
                                        sirius_id integer,	
                                        sirius_score real,
                                        sirius_adduct text,	
                                        sirius_wscore real,
                                        metfrag_id integer,	
                                        metfrag_score real,
                                        metfrag_adduct text,	
                                        metfrag_wscore	real,
                                        sm_lpid integer,	
                                        sm_adduct text,
                                        sm_mid integer,
                                        sm_score real,	
                                        sm_wscore real,
                                        probmetab_id integer,	
                                        probmetab_score real,	
                                        probmetab_wscore real,	
                                        ms1_lookup_score real,	
                                        ms1_lookup_wscore real,	
                                        ms1_lookup_adduct text,
                                        ms1_lookup_id integer,		
                                        biosim_max_score real,	
                                        biosim_wscore real,
                                       	wscore real,	
                                        adduct_overall text,
                                        rank integer)
                                   ''')
        # Don't bother with all the foreign key references - handle this later
        # in mogi if needed.

    elif 'sid' not in cols:
        c.execute('ALTER TABLE combined_annotations ADD sid integer')
        c.execute('ALTER TABLE combined_annotations ADD ms1_lookup_id integer')

    conn.commit()

    c.execute('CREATE UNIQUE INDEX inchikey ON metab_compound(inchikey);')

    return conn, c


def insert_query_m(data, table, conn, columns=None, db_type='mysql',
                   ignore_flag=False):
    """ Insert python list of tuples into SQL table
    Args:
        data (list): List of tuples
        table (str): Name of database table
        conn (connection object): database connection object
        columns (str): String of column names to use if not assigned then
        all columns are presumed to be used [Optional]
        db_type (str): If "sqlite" or "mysql"
    """
    # if length of data is very large we need to break into chunks the
    # insert_query_m is then used recursively untill
    # all data has been inserted
    if len(data) > 10000:
        _chunk_query(data, 10000, columns, conn, table, db_type)
    else:
        # sqlite and mysql have type string (? or %s) reference to use
        if db_type == 'sqlite':
            type_sign = '?'
        else:
            type_sign = '%s'
        # create a string of types for the insertion string (e.g. ?,?,? if
        # inserting 3 columns of data)
        type_com = type_sign + ", "
        type = type_com * (len(data[0]) - 1)
        type = type + type_sign

        if ignore_flag:
            ignore_str = "OR IGNORE"
        else:
            ignore_str = ""

        # if using specific columns to insert data
        if columns:
            stmt = "INSERT " + ignore_str + " INTO " + table + \
                   "( " + columns + ") VALUES (" + type + ")"
        else:
            stmt = "INSERT " + ignore_str + " INTO " + table + " VALUES (" + \
                   type + ")"

        # execute query
        cursor = conn.cursor()
        cursor.executemany(stmt, data)
        conn.commit()


def _chunk_query(l, n, cn, conn, table, db_type):
    """ Call for inserting SQL query in chunks based on n rows
    Args:
        l (list): List of tuples
        n (int): Number of rows
        cn (str): Column names
        conn (connection object): Database connection object
        table (str): Table name
        db_type (str): If "sqlite" or "mysql"
    """
    # For item i in a range that is a length of l,
    [insert_query_m(l[i:i + n], table, conn, cn, db_type) for i in
     range(0, len(l), n)]


def _make_sql_compatible(ll):
    """ Convert any python list of lists (or tuples) so that the strings are
    formatted correctly for insertion into
    Args:
        ll (list): List of lists (or tuples)
    """

    new_ll = []
    for l in ll:
        new_l = ()
        for i in l:
            if not i:
                new_l = new_l + (None,)
            else:

                if isinstance(i, str):
                    if sys.version_info < (3, 0):

                        val = i.decode('utf8').encode('ascii', errors='ignore')
                    else:
                        # in py3 strings should be ok...
                        val = i
                else:
                    val = i
                new_l = new_l + (val,)
        new_ll.append(new_l)

    return new_ll


def max_row(conn, table, id):
    c = conn.cursor()
    c.execute('SELECT max({}) FROM {}'.format(id, table))
    pidr = c.fetchone()[0]
    if pidr:
        pid = pidr + 1
    else:
        pid = 1
    return pid


def pm_sqlite(pm, conn, name, well, well_rtmin, well_rtmax, plid=0):
    cpid = max_row(conn, 's_peak_meta', 'pid')
    s_peak_meta = [(cpid, name, 1, 'dimspy{}'.format(plid), well, well_rtmin,
                    well_rtmax,
                    np.mean([well_rtmin, well_rtmax]))]
    insert_query_m(s_peak_meta,
                   's_peak_meta',
                   conn=conn,
                   columns='pid, name, ms_level, '
                           'spectrum_type,'
                           'well, well_rtmin, well_rtmax, well_rt',
                   db_type='sqlite')

    csid = max_row(conn, 's_peaks', 'sid')

    mz_i = list(zip(pm.mz_matrix[plid], pm.intensity_matrix[plid]))
    sids = range(csid, csid + len(mz_i))
    s_peaks = [(sids[i],) + tuple(mz_i[i]) + (cpid,) for i in range(0,
                                                                    len(mz_i))]
    insert_query_m(s_peaks,
                   's_peaks',
                   conn=conn,
                   columns='sid, mz, i, pid',
                   db_type='sqlite')

    return {row[1]: row[0] for row in s_peaks}, cpid


def pm_link_sqlite(pm, d1, d2, conn):
    mz1 = pm.mz_matrix[0][pm.mz_matrix[1] > 0]
    mz2 = pm.mz_matrix[1][pm.mz_matrix[1] > 0]
    sxs = []
    c = max_row(conn, 's_peaks_X_s_peaks', 'sxsid')
    for i in range(0, len(mz2)):
        if mz1[i] == 0:
            continue
        sxs.append((c, d1[mz1[i]], d2[mz2[i]], mz1[i] - mz2[i], 'chosen_ms1'))
        c += 1

    insert_query_m(sxs,
                   's_peaks_X_s_peaks',
                   conn=conn,
                   db_type='sqlite')


def dimsn_sqlite(pls, msn_prec_d, msn_prec_pid, conn, well):
    # get dictionary of precursors
    mzd = {0: msn_prec_d}
    smd = {0: msn_prec_pid}

    # add all peaks from pl_non_merged (Checking links as we go)
    cpid = max_row(conn, 's_peak_meta', 'pid')
    csid = max_row(conn, 's_peaks', 'sid')
    csmxsid = max_row(conn, 's_peak_meta_X_s_peaks', 'smxsid')

    sm_rows = []
    s_rows_all = []

    for pli in pls:
        ms_level = max(pli.metadata['parent'].keys()) + 1
        precursor_mz = pli.metadata['parent'][ms_level - 1]['mz']
        msnpy_convert_id = pli.metadata['convert_id']
        if 'header' in pli.metadata:
            header = pli.metadata['header']
        else:
            header = pli.ID
        sm_row = (cpid, header, ms_level, precursor_mz, 'msnpy', well,
                  msnpy_convert_id)
        sm_rows.append(sm_row)

        mz = pli.peaks['mz']
        intensity = pli.peaks['intensity']

        s_rows = list(zip(range(csid, csid + len(mz)), mz, intensity, [cpid]
                          * len(mz)))

        s_rows_all.extend(s_rows)

        csid = csid + len(mz) + 1
        cpid += 1
        mzd[pli.ID] = {r[1]: r[0] for r in s_rows}
        smd[pli.ID] = sm_row[0]

    insert_query_m(s_rows_all,
                   's_peaks',
                   columns="sid, mz, i, pid",
                   conn=conn,
                   db_type='sqlite')
    insert_query_m(sm_rows,
                   's_peak_meta',
                   columns="pid, name, ms_level, precursor_mz, "
                           "spectrum_type, well, msnpy_convert_id",
                   conn=conn,
                   db_type='sqlite')

    sm_x_s_rows = []
    for pli in pls:

        mid = max(pli.metadata['parent'].keys())
        precd = pli.metadata['parent'][mid]

        if mid == 1:
            sid = mzd[0][precd['mz']]
        else:
            sid = mzd[precd['ID']][precd['mz']]

        smid = smd[pli.ID]
        sm_x_s_rows.append((csmxsid, smid, sid, 0, 'chosen precursor'))
        csmxsid += 1

    insert_query_m(sm_x_s_rows,
                   's_peak_meta_X_s_peaks',
                   columns="smxsid, pid, sid, mzdiff, link_type",
                   conn=conn,
                   db_type='sqlite')


def ppm_tol_range(mz, ppm):
    mz_low = mz - ((mz * 0.000001) * ppm)
    mz_high = mz + ((mz * 0.000001) * ppm)
    return mz_low, mz_high


def lcms_dims_link(conn, time_tolerance, ppm_lcms,
                   ppm_dims, orig_dims_pid, frac_times_d, well):
    cursor = conn.cursor()
    lcms_r = cursor.execute("""SELECT grpid, mz, rtmin, rtmax 
                                    FROM c_peak_groups
                            """)

    lcms_rows = lcms_r.fetchall()
    dims_r = cursor.execute("SELECT sid, mz  FROM s_peaks WHERE pid={}".
                            format(orig_dims_pid))

    dims_rows = dims_r.fetchall()

    well_seconds_low = frac_times_d[well]['frac_start_minutes'] * 60
    well_seconds_high = frac_times_d[well]['frac_end_minutes'] * 60

    matches = []

    for lcms_row in lcms_rows:

        # check if rtmin rtmax within range
        grpid = lcms_row[0]
        lcms_mz = lcms_row[1]

        lcms_mz_low, lcms_mz_high = ppm_tol_range(lcms_mz, ppm_lcms)
        lcms_seconds_low = lcms_row[2] - time_tolerance
        lcms_seconds_high = lcms_row[3] + time_tolerance
        if (lcms_seconds_high >= well_seconds_low) and \
                (well_seconds_high >= lcms_seconds_low):
            for dims_row in dims_rows:

                sid = dims_row[0]
                dims_mz = dims_row[1]
                dims_mz_low, dims_mz_high = ppm_tol_range(dims_mz,
                                                          ppm_dims)

                if (lcms_mz_high >= dims_mz_low) and \
                        (dims_mz_high >= lcms_mz_low):
                    matches.append((grpid, sid, lcms_mz - dims_mz))

    if matches:
        insert_query_m(matches,
                       'c_peak_groups_X_s_peaks',
                       conn,
                       'grpid,sid,mzdiff',
                       'sqlite')


def process_frac_spectra(tree_pth,
                         dims_pm_pth,
                         well,
                         conn,
                         frac_times_d,
                         time_tolerance,
                         ppm_lcms,
                         ppm_dims):
    dimsn_non_merged_pls, \
    dimsn_merged_pls, \
    dimsn_precursors_pl = tree2peaklist(tree_pth, adjust_mz=False)

    dims_pm = load_peak_matrix_from_hdf5(dims_pm_pth, compatibility_mode=True)

    orig_align_dims_pl = dims_pm.to_peaklist('orig_align_peaklist')

    dims_dimsn_pls = [orig_align_dims_pl, dimsn_precursors_pl[0]]

    dimsn_aligned_pm = align_peaks(dims_dimsn_pls, ppm=2)

    well_rtmin = frac_times_d[well]['frac_start_minutes'] * 60
    well_rtmax = frac_times_d[well]['frac_end_minutes'] * 60

    # Add all peaks from aligned_pm (as 2 scan events):
    #  1) The originally aligned peaks
    #  2) the MS1 precursors
    orig_dims_d, orig_dims_pid = pm_sqlite(dimsn_aligned_pm, conn,
                                           'original_aligned_peaklist',
                                           well, well_rtmin, well_rtmax, 0)

    # msn aligned peaklist
    msn_prec_d, msn_prec_pid = pm_sqlite(dimsn_aligned_pm, conn,
                                         'dimsn_ms1_precursors',
                                         well, well_rtmin, well_rtmax, 1)

    # add the links
    pm_link_sqlite(dimsn_aligned_pm, orig_dims_d, msn_prec_d, conn)

    # add dimsn non merged peaklists
    dimsn_sqlite(dimsn_non_merged_pls, msn_prec_d, msn_prec_pid, conn, well, )

    # add dimsn merged peaklist
    dimsn_sqlite(dimsn_merged_pls, msn_prec_d, msn_prec_pid, conn, well)

    # add LC-MS to DI-MS link
    lcms_dims_link(conn,
                   time_tolerance,
                   ppm_lcms,
                   ppm_dims,
                   orig_dims_pid,
                   frac_times_d,
                   well)

    cursor = conn.cursor()
    r = cursor.execute("SELECT pid, msnpy_convert_id FROM s_peak_meta WHERE "
                       "well='{}'".format(well))
    return {spm[1]: spm[0] for spm in r.fetchall()}


def add_beams(beams_pth, conn, comp_conn, pid_d, ms1_lookup_source='hmdb'):
    '''
    '''
    if ms1_lookup_source == 'hmdb':
        table_nm = 'hmdb'
        column_nm = 'hmdb_id'
    elif ms1_lookup_source == 'kegg':
        table_nm = 'kegg'
        column_nm = 'kegg_cid'
    elif ms1_lookup_source == 'pubchem':
        table_nm = 'pubchem'
        column_nm = 'pubchem_cid'
    else:
        return 0

    c = max_row(conn, 'ms1_lookup_results', 'id')

    with open(beams_pth, 'r') as bf:
        dr = csv.DictReader(bf, delimiter='\t')
        rows = []

        for drow in dr:
            curs = comp_conn.cursor()
            r = curs.execute("""SELECT inchikey FROM {} WHERE
                                               {}='{}'
                                         """.format(table_nm,
                                                    column_nm,
                                                    drow['compound_id']))
            inchikey_r = r.fetchall()
            if not inchikey_r:
                inchikey = 'unknown'
            else:
                inchikey = inchikey_r[0][0]
            pid = pid_d[int(drow['msnpy_convert_id'])]
            rows.append((int(c),
                         int(pid),
                         int(drow['msnpy_convert_id']),
                         drow['well'],
                         drow['name'],
                         float(drow['mz']),
                         float(drow['exact_mass']),
                         float(drow['ppm_error']),
                         drow['adduct'],
                         int(drow['C']),
                         int(drow['H']),
                         int(drow['N']),
                         int(drow['O']),
                         int(drow['P']),
                         int(drow['S']),
                         drow['molecular_formula'],
                         drow['compound_name'],
                         drow['compound_id'],
                         inchikey,
                         1  # score is always 1 for beams
                         ))
            c += 1

    cols = "id, pid, msnpy_convert_id,well, name, mz, exact_mass, " \
           "ppm_error," \
           "adduct,C,H,N,O,P,S,molecular_formula,compound_name,compound_id," \
           " inchikey, score"

    insert_query_m(rows,
                   'ms1_lookup_results',
                   conn,
                   columns=cols,
                   db_type='sqlite')


def add_metfrag(metfrag_pth, conn, pid_d):
    '''
    '''

    c = max_row(conn, 'metfrag_results', 'id')

    with open(metfrag_pth, 'r') as bf:
        dr = csv.DictReader(bf, delimiter='\t')
        rows = []

        for drow in dr:
            pid = pid_d[int(drow['msnpy_convert_id'])]
            rows.append((int(c),
                         int(pid),
                         int(drow['msnpy_convert_id']),
                         drow['adduct'],
                         drow['CompoundName'],
                         drow['ExplPeaks'],
                         drow['FormulasOfExplPeaks'],
                         drow['FragmenterScore'],
                         drow['FragmenterScore_Values'],
                         drow['Identifier'],
                         drow['InChI'],
                         drow['InChIKey'],
                         drow['InChIKey1'],
                         drow['InChIKey2'],
                         drow['InChIKey3'],
                         drow['MaximumTreeDepth'],
                         drow['MolecularFormula'],
                         drow['MonoisotopicMass'],
                         drow['Name'],
                         drow['NoExplPeaks'],
                         drow['NumberPeaksUsed'],
                         drow['SMILES'],
                         drow['Score']
                         ))
            c += 1

    cols = "id, pid, msnpy_convert_id, adduct, CompoundName, ExplPeaks," \
           "FormulasOfExplPeaks, FragmenterScore, FragmenterScore_Values, " \
           "Identifier, InChI, InChIKey, InChIKey1, InChIKey2, InChIKey3," \
           "MaximumTreeDepth, MolecularFormula,MonoisotopicMass,Name," \
           "NoExplPeaks,NumberPeaksUsed, SMILES, Score"

    insert_query_m(rows,
                   'metfrag_results',
                   conn,
                   columns=cols,
                   db_type='sqlite')


def add_spectral_matching(spectral_matching_pth, conn, pid_d):
    '''
    '''
    c = max_row(conn, 'sm_matches', 'mid')
    with open(spectral_matching_pth, 'r') as bf:
        dr = csv.DictReader(bf, delimiter='\t')
        rows = []

        for drow in dr:
            pid = pid_d[int(drow['msnpy_convert_id'])]
            rows.append((int(c),
                         int(pid),
                         int(drow['msnpy_convert_id']),
                         drow['lpid'],
                         drow['dpc'],
                         drow['rdpc'],
                         drow['cdpc'],
                         drow['mcount'],
                         drow['allcount'],
                         drow['mpercent'],
                         drow['library_rt'],
                         drow['query_rt'],
                         drow['rtdiff'],
                         drow['library_precursor_mz'],
                         drow['query_precursor_mz'],
                         drow['library_precursor_ion_purity'],
                         drow['query_precursor_ion_purity'],
                         drow['library_accession'],
                         drow['library_precursor_type'],
                         drow['library_entry_name'],
                         drow['inchikey'],
                         drow['library_source_name'],
                         drow['library_compound_name']
                         ))
            c += 1

    cols = "mid, qpid, msnpy_convert_id, lpid, dpc, rdpc, cdpc,	" \
           "mcount, allcount, mpercent, library_rt,	query_rt,	" \
           "rtdiff,	library_precursor_mz,	query_precursor_mz,	" \
           "library_precursor_ion_purity,	query_precursor_ion_purity,	" \
           "library_accession,	library_precursor_type,	library_entry_name,	" \
           "inchikey,	library_source_name,	library_compound_name"

    insert_query_m(rows,
                   'sm_matches',
                   conn,
                   columns=cols,
                   db_type='sqlite')


def neg_min_max(x):
    if np.equal.reduce(x):
        return [1] * len(x)
    else:
        x = np.array(x)
        xn = (x - min(x)) / (max(x) - min(x))
        return abs(xn - 1)


def add_sirius(sirius_pth, conn, pid_d):
    '''
    '''

    with open(sirius_pth, 'r') as bf:
        dr = csv.DictReader(bf, delimiter='\t')
        rows = []
        c = max_row(conn, 'sirius_csifingerid_results', 'id')
        annotation_group = {}
        for drow in dr:
            pid = pid_d[int(drow['msnpy_convert_id'])]
            rows.append((c,
                         pid,
                         drow['msnpy_convert_id'],
                         drow['well'],
                         drow['adduct'],
                         drow['inchikey2D'],
                         drow['InChI'],
                         drow['molecularFormula'],
                         drow['Rank'],
                         drow['Score'],
                         drow['Name'],
                         drow['smiles'],
                         drow['xlogp'],
                         drow['pubchemids'],
                         drow['links']
                         ))
            if pid not in annotation_group:
                annotation_group[pid] = {}
            annotation_group[pid][c] = float(drow['Score'])
            c += 1

    ##################
    # calculate bound sirius score
    ##################
    bounded_score_d = {}
    for pid, scores in annotation_group.items():
        score_l = []
        cid_l = []
        for cid, score in scores.items():
            score_l.append(abs(float(score)))
            cid_l.append(cid)
        bounded_score = neg_min_max(score_l)
        for i in range(0, len(bounded_score)):
            bounded_score_d[cid_l[i]] = bounded_score[i]

    # add bounded score to rows
    rows = [row + (bounded_score_d[row[0]],) for row in rows]

    cols = "id, pid, msnpy_convert_id, well, adduct, inchikey2D, InChI, " \
           "molecularFormula," \
           "Rank, Score, Name, smiles, xlogp, pubchemids, links, bounded_score"

    insert_query_m(rows, 'sirius_csifingerid_results', conn, columns=cols,
                   db_type='sqlite')


def process_frac_annotation(beams_pth,
                            sirius_pth,
                            metfrag_pth,
                            spectral_matching_pth,
                            conn,
                            pid_d,
                            ms1_lookup_source,
                            comp_conn):
    # write beams to sqlite database
    add_beams(beams_pth, conn, comp_conn, pid_d, ms1_lookup_source)

    # write sirius to sqlite database
    add_sirius(sirius_pth, conn, pid_d)

    # write metfrag to sqlite database
    add_metfrag(metfrag_pth, conn, pid_d)

    # write spectral matching to sqlite database
    add_spectral_matching(spectral_matching_pth, conn, pid_d)


def col_multiple_check(r, cols):
    if "," in cols or "*" == cols:
        return r
    else:
        return [i[0] for i in r]


def get_metab_compound_rows(comp_conn, inchikeys, inchi_level='inchikey',
                            cols='*'):
    c = comp_conn.cursor()
    inchi_str = "','".join(inchikeys)

    sql_stmt = "SELECT {} FROM metab_compound WHERE {} IN ('{" \
               "}')".format(cols, inchi_level, inchi_str)
    c.execute(sql_stmt)
    r = c.fetchall()
    return col_multiple_check(r, cols)


def sql_simple_select(conn, table_nm, cols="*"):
    c = conn.cursor()
    c.execute("SELECT {} FROM {}".format(cols, table_nm))
    r = c.fetchall()

    return col_multiple_check(r, cols)


def get_inchikey_sid(conn, table_nm, inchi_sid_d, pid='pid',
                     mid="id", score='score',
                     adduct='adduct', weight=1.0):
    # Need to update to only keep the best score (in case doing scan by scan)

    c = conn.cursor()
    sql_stmt = """SELECT m.inchikey, spXsp.sid1, m.{},
                        m.{}, m.{}*{}, m.{} 
                    FROM {} AS m 
                    LEFT JOIN s_peak_meta AS spm ON spm.pid=m.{}
                    LEFT JOIN s_peak_meta_X_s_peaks AS sxs ON sxs.pid=spm.pid  
                    LEFT JOIN s_peaks AS sp ON sp.sid=sxs.sid
                    LEFT JOIN s_peaks_X_s_peaks AS spXsp ON spXsp.sid2=sp.sid
                    WHERE m.msnpy_convert_id NOT NULL 
    """.format(mid, score, score, weight, adduct, table_nm, pid)

    c.execute(sql_stmt)
    r = c.fetchall()
    return inchi_sid_d_update(r, inchi_sid_d, table_nm)


def inchi_sid_d_update(r, inchi_sid_d, table_nm):
    # NOTE we keep the best score for inchikey-mz-tablenm - so if we are using
    # multiple scans or energies that give the same inchikey annotation we
    # only keep the one with the best score (highest)
    for i in r:
        inchi_sid = "{}_{}".format(i[0], i[1])
        rowd = {'mid': i[2], 'score': i[3], 'wscore': i[4], 'adduct': i[5]}
        if inchi_sid in inchi_sid_d:
            if table_nm in inchi_sid_d[inchi_sid]:
                if float(rowd['wscore']) > float(inchi_sid_d[inchi_sid][
                                                     table_nm]['wscore']):
                    inchi_sid_d[inchi_sid][table_nm] = rowd
            else:
                inchi_sid_d[inchi_sid][table_nm] = rowd
        else:
            inchi_sid_d[inchi_sid] = {table_nm: rowd}
    return inchi_sid_d


def get_inchikey_sid_sirius(conn, weight):
    # Note we want to get the sid for the original peak list
    c = conn.cursor()
    c.execute("""SELECT mc.inchikey, spXsp.sid1, c.id,
                        c.bounded_score, c.bounded_score*{} AS wscore,
                        c.adduct   
                    FROM 
                    sirius_csifingerid_results AS c 
                  LEFT JOIN metab_compound AS mc ON mc.inchikey1=c.inchikey2D
                  LEFT JOIN s_peak_meta AS spm ON spm.pid=c.pid
                  LEFT JOIN s_peak_meta_X_s_peaks AS sxs ON sxs.pid=spm.pid  
                  LEFT JOIN s_peaks AS sp ON sp.sid=sxs.sid
                  LEFT JOIN s_peaks_X_s_peaks AS spXsp ON spXsp.sid2=sp.sid
                  WHERE c.pid NOT NULL    
    """.format(weight))
    r = c.fetchall()
    inchi_sid_d = {}
    table_nm = 'sirius_csifingerid_results'
    return inchi_sid_d_update(r, inchi_sid_d, table_nm)


def add_to_row(results, name, row):
    if name in results:
        sr = results[name]
        row.extend([sr['mid'], sr['score'], sr['wscore'], sr['adduct']])
    else:
        row.extend([0, 0, 0, 0])
        sr = {'mid': 0, 'score': 0, 'wscore': 0, 'adduct': 0}
    return row, sr


def combine_annotations(conn, comp_conn, weights):
    # First get all the compound information for each inchikey from each
    # result (metfrag, sirius, spectral matching)
    # First get inchikeys from sirius result
    inchikeys1 = sql_simple_select(conn, 'sirius_csifingerid_results',
                                   'inchikey2D')

    inchikeys = get_metab_compound_rows(comp_conn, inchikeys1, 'inchikey1',
                                        'inchikey')

    # get inchikeys from beams
    inchikeys.extend(sql_simple_select(conn, 'ms1_lookup_results', 'inchikey'))

    # get inchikeys from metfrag
    inchikeys.extend(sql_simple_select(conn, 'metfrag_results', 'inchikey'))

    # get inchikeys from spectral matching
    inchikeys.extend(sql_simple_select(conn, 'sm_matches', 'inchikey'))

    # remove duplicates
    inchikeys = list(set(inchikeys))

    # add to metab_compound in conn
    metab_rows = get_metab_compound_rows(comp_conn, inchikeys, 'inchikey',
                                         '*')
    # biosim dict
    biosim_d = {row[0]: row[len(row) - 3] for row in metab_rows}

    # Add rows to sqlite results database
    insert_query_m(metab_rows, 'metab_compound', conn, columns=None,
                   db_type="sqlite", ignore_flag=True)

    # for each pid - get the (best) results for each inchikey and annotations
    # approach
    # make combined dict
    # get all inchikey_pid combinations from all approaches.
    # Then loop through and check each approach and make a row for each!
    inchi_sid_d = get_inchikey_sid_sirius(conn,
                                          weight=weights['sirius_csifingerid'])
    inchi_sid_d = get_inchikey_sid(conn,
                                   'ms1_lookup_results',
                                   inchi_sid_d,
                                   weight=weights['beams'])
    inchi_sid_d = get_inchikey_sid(conn,
                                   'metfrag_results',
                                   inchi_sid_d,
                                   weight=weights['metfrag'])
    inchi_sid_d = get_inchikey_sid(conn,
                                   'sm_matches',
                                   inchi_sid_d,
                                   pid='qpid',
                                   score='dpc',
                                   adduct='library_precursor_type',
                                   mid="mid",
                                   weight=weights['spectral_matching'])
    # Creat rows
    sid_d = {}
    for inchi_sid, results in inchi_sid_d.items():
        row = inchi_sid.split("_")
        sid = row[1]
        row, sirius_d = add_to_row(results, 'sirius_csifingerid_results', row)
        row, metfrag_d = add_to_row(results, 'metfrag_results', row)
        row, beams_d = add_to_row(results, 'ms1_lookup_results', row)
        row, sm_d = add_to_row(results, 'sm_matches', row)

        # Get biosim score
        biosim_score = biosim_d[row[0]]
        biosim_wscore = biosim_score * 0.25

        # calculated weighted score
        wscore = sirius_d['wscore'] + \
                 metfrag_d['wscore'] + \
                 beams_d['wscore'] + \
                 sm_d['wscore'] + \
                 biosim_wscore
        # get overall adduct column
        adducts = set([str(sirius_d['adduct']),
                       str(metfrag_d['adduct']),
                       str(beams_d['adduct']),
                       str(sm_d['adduct'])])
        adducts = [a for a in adducts if a]
        adducts.remove('0')
        adduct_str = ",".join(adducts)

        row.extend([biosim_score, biosim_wscore, adduct_str, wscore])

        if sid in sid_d:
            sid_d[sid].append(row)
        else:
            sid_d[sid] = [row]

    final_rows = []
    for sid, rows in sid_d.items():
        # Get all scores
        scores = [r[len(r) - 1] for r in rows]

        ranks = rankdata(np.array(scores), 'dense')

        for i in range(0, len(ranks)):
            rows[i].append(int(ranks[i]))
        rows = sorted(rows, key=itemgetter(len(rows[0]) - 1))
        final_rows.extend(rows)

    columns = """
    inchikey,sid,
    sirius_id,sirius_score,sirius_wscore,sirius_adduct,
    metfrag_id,metfrag_score,metfrag_wscore,metfrag_adduct,
    ms1_lookup_id,ms1_lookup_score,ms1_lookup_wscore,ms1_lookup_adduct,
    sm_mid,sm_score,sm_wscore,sm_adduct,
    biosim_max_score, biosim_wscore,
    adduct_overall, wscore, rank
    """

    insert_query_m(final_rows, "combined_annotations", conn,
                   db_type='sqlite', columns=columns)


def pth2dict(dr):
    return {os.path.splitext(os.path.basename(f))[0]:
                os.path.join(dr, f) for f in os.listdir(dr)}


def process_all_fracs(dims_pl_dir,
                      dimsn_tree_dir,
                      lcms_sqlite_pth,
                      comp_pth,
                      beams_dir,
                      metfrag_dir,
                      sirius_dir,
                      spectral_matching_dir,
                      frac_times_pth,
                      time_tolerance,
                      ppm_lcms,
                      ppm_dims,
                      weights,
                      ms1_lookup_source,
                      out_sqlite,
                      out_csv):
    # get file name in dims_pl_pth
    dims_pl_pths = pth2dict(dims_pl_dir)

    # get filenames in dimsn_tree_pth
    dimsn_tree_pths = pth2dict(dimsn_tree_dir)

    # Get BEAMS file paths
    beams_pths = pth2dict(beams_dir)

    # get sirius file
    sirius_pths = pth2dict(sirius_dir)

    # get metrag file paths
    metfrag_pths = pth2dict(metfrag_dir)

    # get spectral matching paths
    spectral_matching_pths = pth2dict(spectral_matching_dir)

    # Loop through files
    frac_spectra = {}
    for well in dims_pl_pths.keys():
        frac_spectra[well] = {}
        frac_spectra[well]['dims_pl_pth'] = dims_pl_pths[well]
        if well in dimsn_tree_pths:
            frac_spectra[well]['dimsn_tree_pth'] = dimsn_tree_pths[well]
        else:
            frac_spectra[well]['dimsn_tree_pth'] = ''

        if well in beams_pths:
            frac_spectra[well]['beams_pth'] = beams_pths[well]
        else:
            frac_spectra[well]['beams_pth'] = ''

        if well in sirius_pths:
            frac_spectra[well]['sirius_pth'] = sirius_pths[well]
        else:
            frac_spectra[well]['sirius_pth'] = ''

        if well in metfrag_pths:
            frac_spectra[well]['metfrag_pth'] = metfrag_pths[well]
        else:
            frac_spectra[well]['metfrag_pth'] = ''

        if well in spectral_matching_pths:
            frac_spectra[well]['spectral_matching_pth'] = \
                spectral_matching_pths[well]
        else:
            frac_spectra[well]['spectral_matching_pth'] = ''

    # os.remove('frac_exp.sqlite')
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

    comp_conn = sqlite3.connect(comp_pth)

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

        process_frac_annotation(beams_pth=spths['beams_pth'],
                                sirius_pth=spths['sirius_pth'],
                                metfrag_pth=spths['metfrag_pth'],
                                spectral_matching_pth=spths[
                                    'spectral_matching_pth'],
                                conn=conn,
                                pid_d=pid_d,
                                ms1_lookup_source=ms1_lookup_source,
                                comp_conn=comp_conn
                                )

    # combine annotations
    combine_annotations(conn, comp_conn, weights)

    # create summary table
    # Include both the LC-MS and DI-MS annotations
    summarise_annotations(conn, 'test.tsv')


def summarise_annotations(conn, out_file, rank_limit=100):
    # Creat two summary tables (one for LC and one for DIMS and then output
    # to single csv file
    c = conn.cursor()
    sql_stmt = """
       SELECT
       'dims' AS ms_type, 
       sp.sid,
       '' AS grpid,
       '' AS grp_name,
       round(sp.mz, 6) AS mz,
       round(sp.i, 2) AS i,
       round(spm.well_rt, 3) AS rt,
       round(spm.well_rtmin,3) AS rtmin,
       round(spm.well_rtmax,3) AS rtmax,
       GROUP_CONCAT(DISTINCT (CAST (cpgXsp.grpid AS INTEGER) ) ) AS lc_grpid_mtchs,
       '' AS dims_sid_mtchs,
       spm.well,  
       mc.inchi,
       mc.inchikey,
       mc.inchikey1,
       mc.inchikey2,
       mc.inchikey3,
       mc.name,
       mc.exact_mass,
       mc.molecular_formula,
       mc.pubchem_cids,
       mc.kegg_cids,
       mc.kegg_brite,
       mc.kegg_drugs,
       mc.hmdb_ids,
       mc.hmdb_bio_custom_flag,
       mc.hmdb_drug_flag,
       mc.biosim_max_count,
       mc.biosim_hmdb_ids,
       '' AS fragmentation_acquistion_num,
       '' AS mean_precursor_ion_purity,
       l.accession,
       l.id AS lpid,
       ca.sirius_score,
       ca.sirius_wscore,
       ca.metfrag_score,
       ca.metfrag_wscore,
       ca.sm_score,
       ca.sm_wscore,
       ca.probmetab_score,
       ca.probmetab_wscore,
       ca.ms1_lookup_score,
       ca.ms1_lookup_wscore,
       mc.biosim_max_score,
       ca.biosim_wscore,
       ca.wscore,
       ca.rank,
       ca.adduct_overall
    FROM s_peaks AS sp
       LEFT JOIN
       combined_annotations AS ca ON sp.sid = ca.sid
       LEFT JOIN
       metab_compound AS mc ON ca.inchikey = mc.inchikey
       LEFT JOIN
       l_s_peak_meta AS l ON l.id = ca.sm_lpid
       LEFT JOIN
       s_peak_meta AS spm ON spm.pid = sp.pid
       LEFT JOIN
       c_peak_groups_X_s_peaks AS cpgXsp ON cpgXsp.sid = sp.sid

    WHERE (sp.sid IS NOT NULL) AND (IFNULL(ca.rank<={}, 1)) AND  (
    spm.spectrum_type IS 'dimspy0')
    GROUP BY sp.sid,
             IFNULL(ca.inchikey, sp.sid)
    ORDER BY sp.sid,
             IFNULL(ca.rank, sp.sid)
    """.format(rank_limit)

    r = c.execute(sql_stmt)
    dims_annotations = r.fetchall()

    sql_stmt = """    SELECT 
       'lcms' AS ms_type,
       '' AS sid,
       cpg.grpid,
       cpg.grp_name,
       round(cpg.mz, 6) AS mz,
       ROUND(AVG(cp._into),3) AS i,
       round(cpg.rt, 3) AS rt,
       round(cpg.rtmin,3) AS rtmin,
       round(cpg.rtmax,3) AS rtmax,
       '' AS lcms_grpid_mtchs,
       GROUP_CONCAT(DISTINCT (CAST (cpgXsp.sid AS INTEGER) ) ) AS dims_sid_mtchs,
       spm.well,  
       mc.inchi,
       mc.inchikey,
       mc.inchikey1,
       mc.inchikey2,
       mc.inchikey3,
       mc.name,
       mc.exact_mass,
       mc.molecular_formula,
       mc.pubchem_cids,
       mc.kegg_cids,
       mc.kegg_brite,
       mc.kegg_drugs,
       mc.hmdb_ids,
       mc.hmdb_bio_custom_flag,
       mc.hmdb_drug_flag,
       mc.biosim_max_count,
       mc.biosim_hmdb_ids,
       GROUP_CONCAT(DISTINCT(cast(spm.acquisitionNum  as INTEGER) )) AS fragmentation_acquistion_num,
       ROUND(AVG(spm.inPurity),3) AS mean_precursor_ion_purity,
       l.accession,
       l.id AS lpid,
       ca.sirius_score,
       ca.sirius_wscore,
       ca.metfrag_score,
       ca.metfrag_wscore,
       ca.sm_score,
       ca.sm_wscore,
       ca.probmetab_score,
       ca.probmetab_wscore,
       ca.ms1_lookup_score,
       ca.ms1_lookup_wscore,
       mc.biosim_max_score,
       ca.biosim_wscore,
       ca.wscore,
       ca.rank,
       ca.adduct_overall
    FROM c_peak_groups AS cpg
    LEFT JOIN
    combined_annotations AS ca ON ca.grpid=cpg.grpid
    LEFT JOIN
    metab_compound AS mc ON ca.inchikey = mc.inchikey
    LEFT JOIN
    l_s_peak_meta AS l ON l.id = ca.sm_lpid
    LEFT JOIN
    c_peak_groups_X_s_peaks AS cpgXsp ON cpgXsp.grpid = cpg.grpid
    LEFT JOIN
    c_peak_X_c_peak_group AS cpXcpg ON cpXcpg.grpid = cpg.grpid
    LEFT JOIN
    c_peaks AS cp ON cp.cid = cpXcpg.cid
    LEFT JOIN
    c_peak_X_s_peak_meta AS cpXspm ON cpXspm.cid = cp.cid
    LEFT JOIN
    s_peak_meta AS spm ON spm.pid = cpXspm.pid
    LEFT JOIN
    fileinfo AS fi ON fi.fileid = spm.fileid

    WHERE IFNULL(ca.rank<={}, 1) AND fi.class NOT LIKE '%blank%'
    GROUP BY
     cpg.grpid,
     IFNULL(ca.inchikey, cpg.grpid)
    ORDER BY 
     cpg.grpid, 
     IFNULL(ca.rank, cpg.grpid)
    """.format(rank_limit)

    r = c.execute(sql_stmt)
    lcms_annotations = r.fetchall()
    colnames = list(map(lambda x: x[0], c.description))

    with open(out_file, 'w') as summary_f:

        sw = csv.writer(summary_f, delimiter='\t')
        sw.writerow(colnames)
        sw.writerows(lcms_annotations)
        sw.writerows(dims_annotations)


if __name__ == "__main__":
    parser = ArgumentParser(description='Combine spectra and '
                                        'annotation results '
                                        'from and LC-MS/MS '
                                        'fractionation '
                                        'experiment',
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-d', '--dims_pl_dir',
                        type=str, required=True,
                        help="Path to the directory of dims peaklists (hdf5)")
    parser.add_argument('-n', '--dimsn_tree_dir',
                        type=str, required=True,
                        help="Path to the directory of dimsn peaklists (hdf5)")
    parser.add_argument('-l', '--lcms_sqlite_pth',
                        type=str, required=True,
                        help="Path to the LC-MS sqlite results")
    parser.add_argument('-c', '--comp_pth',
                        type=str, required=True,
                        help="Path to the compound sqlite database")
    parser.add_argument('-b', '--beams_dir',
                        type=str, required=True,
                        help="Path to the directory of BEAMS results")
    parser.add_argument('-m', '--metfrag_dir',
                        type=str, required=True,
                        help="Path to the directory of MetFrag results")
    parser.add_argument('-s', '--sirius_dir',
                        type=str, required=True,
                        help="Path to the directory of Sirius " \
                             "CSI:FingerID results")
    parser.add_argument('-x', '--spectral_matching_dir',
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
                        help="Database used for MS1 lookup source (e.g."
                             "hmdb, kegg or pubchem)")

    parser.add_argument('-y', '--out_csv',
                        type=str, default='lcfrac_results.csv',
                        help="Database used for MS1 lookup source (e.g."
                             "hmdb, kegg or pubchem)")

    args = parser.parse_args()

    weights = {'sirius_csifingerid': args.weight_sirius_csifingerid,
               'metfrag': args.weight_metfrag,
               'biosim': args.weight_biosim,
               'spectral_matching': args.weight_spectral_matching,
               'beams': args.weight_beams}

    process_all_fracs(args.dims_pl_dir,
                      args.dimsn_tree_dir,
                      args.lcms_sqlite_pth,
                      args.comp_pth,
                      args.beams_dir,
                      args.metfrag_dir,
                      args.sirius_dir,
                      args.spectral_matching_dir,
                      args.frac_times_pth,
                      args.time_tolerance,
                      args.ppm_lcms,
                      args.ppm_dims,
                      weights,
                      args.ms1_lookup_source,
                      args.out_sqlite,
                      args.out_csv)
