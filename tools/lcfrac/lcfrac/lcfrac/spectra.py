import numpy as np
from msnpy.convert import tree2peaklist
from dimspy.portals.hdf5_portal import load_peak_matrix_from_hdf5, load_peaklists_from_hdf5
from dimspy.process.peak_alignment import align_peaks
from db import max_row, insert_query_m


def dimsn_sqlite(pls, msn_prec_d, msn_prec_pid, conn, well, ppm=5):
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
        mzd[pli.ID] = {round(r[1], 8): r[0] for r in s_rows}
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
            # need to match within ppm range because the precursor could have been recalibrated based on MF!
            sid, mzdiff = match_mzd(mzd[0], precd['mz'], ppm)

        else:
            sid, mzdiff = match_mzd(mzd[precd['ID']], precd['mz'], ppm)

        if sid:
            smid = smd[pli.ID]
            sm_x_s_rows.append((csmxsid, smid, sid, mzdiff, 'chosen precursor'))
            csmxsid += 1

    insert_query_m(sm_x_s_rows,
                   's_peak_meta_X_s_peaks',
                   columns="smxsid, pid, sid, mzdiff, link_type",
                   conn=conn,
                   db_type='sqlite')

def match_mzd(mzd, mz, ppm):
    mz_low, mz_high = ppm_tol_range(mz, ppm)
    match_sid = ''
    match_mzdiff = ''
    for mzi, sid in mzd.items():
        mzi_low, mzi_high = ppm_tol_range(mzi, ppm)

        if (mz_high >= mzi_low) and \
                (mzi_high >= mz_low):
            mzdiff = mz-mzi
            if not match_sid or abs(mzdiff) < abs(match_mzdiff):
                match_sid = sid
                match_mzdiff = mzdiff

    return match_sid, match_mzdiff


def ppm_tol_range(mz, ppm):
    mz_low = mz - ((mz * 0.000001) * ppm)
    mz_high = mz + ((mz * 0.000001) * ppm)
    return mz_low, mz_high


def lcms_dims_link(conn, orig_dims_pid, time_tolerance, ppm_lcms, ppm_dims, well_seconds_low, well_seconds_high):
    cursor = conn.cursor()
    lcms_r = cursor.execute("""SELECT grpid, mz, rtmin, rtmax 
                                    FROM c_peak_groups
                            """)

    lcms_rows = lcms_r.fetchall()
    dims_r = cursor.execute("SELECT sid, mz  FROM s_peaks WHERE pid={}".
                            format(orig_dims_pid))

    dims_rows = dims_r.fetchall()

    matches = []

    for lcms_row in lcms_rows:

        # check if rtmin rtmax within range
        grpid = lcms_row[0]
        lcms_mz = lcms_row[1]

        lcms_mz_low, lcms_mz_high = ppm_tol_range(lcms_mz, ppm_lcms)
        lcms_seconds_low = lcms_row[2] - time_tolerance
        lcms_seconds_high = lcms_row[3] + time_tolerance

        if (lcms_seconds_high >= well_seconds_low) and (well_seconds_high >= lcms_seconds_low):
            for dims_row in dims_rows:
                sid = dims_row[0]
                dims_mz = dims_row[1]
                dims_mz_low, dims_mz_high = ppm_tol_range(dims_mz, ppm_dims)
                if (lcms_mz_high >= dims_mz_low) and (dims_mz_high >= lcms_mz_low):
                    matches.append((grpid, sid, lcms_mz - dims_mz))
                    print('MATCH!!!!', (grpid, sid, lcms_mz - dims_mz))

    if matches:
        insert_query_m(matches,
                       'c_peak_groups_X_s_peaks',
                       conn,
                       'grpid,sid,mzdiff',
                       'sqlite')


def peaklist_to_sqlite(pl, conn, name, well, well_rtmin, well_rtmax):
    # Set the meta data for the original peak list
    cpid = max_row(conn, 's_peak_meta', 'pid')
    s_peak_meta = [(cpid, name, 1, 'dimspy', well, well_rtmin,
                    well_rtmax,
                    np.mean([well_rtmin, well_rtmax]))]
    insert_query_m(s_peak_meta,
                   's_peak_meta',
                   conn=conn,
                   columns='pid, name, ms_level, '
                           'spectrum_type,'
                           'well, well_rtmin, well_rtmax, well_rt',
                   db_type='sqlite')

    # save the peaks for the original peaklist
    csid = max_row(conn, 's_peaks', 'sid')

    mz_i = list(zip(pl.peaks['mz'], pl.peaks['intensity']))

    sids = range(csid, csid + len(mz_i))
    s_peaks = [(sids[i],) + tuple(mz_i[i]) + (cpid,) for i in range(0,
                                                                    len(mz_i))]
    insert_query_m(s_peaks,
                   's_peaks',
                   conn=conn,
                   columns='sid, mz, i, pid',
                   db_type='sqlite')

    return {round(row[1], 8): row[0] for row in s_peaks}, cpid


def process_ms1_spectra(wellinfo, conn):

    # first step is add the MS1 DIMS data to the database and get ids
    print(wellinfo.dims_pths.spectra,)
    orig_align_dims_pl = load_peaklists_from_hdf5(wellinfo.dims_pths.spectra, compatibility_mode=False)[0]

    # save peaklist to sqlite
    orig_dims_d, orig_dims_pid = peaklist_to_sqlite(orig_align_dims_pl,
                                                    conn,
                                                    'original_aligned_peaklist',
                                                    wellinfo.well_number,
                                                    wellinfo.rtmin,
                                                    wellinfo.rtmax)

    return orig_dims_d, orig_dims_pid


def process_frag_spectra(wellinfo,
                         msn_element,
                         lcms_conn,
                         ppm_dims,
                         orig_dims_d):

    msn_pths = wellinfo.msn_pths_c[msn_element]
    print(wellinfo, msn_element, msn_pths.non_merged_spectra)
    # Now get dimsn precursors and align
    dimsn_non_merged_pls = load_peaklists_from_hdf5(msn_pths.non_merged_spectra, compatibility_mode=False)
    dimsn_merged_pls = load_peaklists_from_hdf5(msn_pths.merged_spectra, compatibility_mode=False)
    dimsn_precursors_pl = load_peaklists_from_hdf5(msn_pths.ms1_precursor_spectra, compatibility_mode=False)

    msn_prec_d, msn_prec_pid = peaklist_to_sqlite(dimsn_precursors_pl[0],
                                                  lcms_conn,
                                                  'dimsn_ms1_precursors_{}'.format(msn_element),
                                                  wellinfo.well_number,
                                                  wellinfo.rtmin,
                                                  wellinfo.rtmax)
    sxs = []
    csxsid = max_row(lcms_conn, 's_peaks_X_s_peaks', 'sxsid')
    for mz, sid_msn in msn_prec_d.items():

        sid_dims, mzdiff = match_mzd(orig_dims_d, mz, ppm_dims)
        if sid_dims:
            sxs.append((csxsid, sid_dims, sid_msn, mzdiff, 'chosen_ms1'))

            csxsid += 1
    if sxs:
        insert_query_m(sxs,
                   's_peaks_X_s_peaks',
                   conn=lcms_conn,
                   db_type='sqlite')

    # # add dimsn non merged peaklists
    dimsn_sqlite(dimsn_non_merged_pls, msn_prec_d, msn_prec_pid, lcms_conn, wellinfo.well_number, ppm_dims)

    # add dimsn merged peaklist
    dimsn_sqlite(dimsn_merged_pls, msn_prec_d, msn_prec_pid, lcms_conn, wellinfo.well_number, ppm_dims)

    cursor = lcms_conn.cursor()
    msnpy_convert_r = cursor.execute("SELECT pid, msnpy_convert_id FROM s_peak_meta WHERE "
                       "well='{}'".format(wellinfo.well_number))

    return {spm[1]: spm[0] for spm in msnpy_convert_r.fetchall()}, msn_prec_d


def add_library_spectra(conn, library_conn):
    print("Get all the library pids that do not have any relevant library spectra yet")
    lpid_r = conn.execute("""SELECT lpid FROM sm_matches AS sm LEFT JOIN l_s_peak_meta AS lspm ON lspm.id=sm.lpid
                                WHERE lspm.id IS NULL""")
    lpids = [str(int(i[0])) for i in lpid_r.fetchall()]

    lpids_str = ",".join(lpids)

    print('Get the meta data from the library')
    new_lsm_r = library_conn.execute("SELECT * FROM library_spectra_meta WHERE id IN ({})".format(lpids_str))
    insert_query_m([i for i in new_lsm_r.fetchall()],
                   'l_s_peak_meta',
                   conn=conn,
                   db_type='sqlite')

    print('Get relevant library spectra')
    new_ls_r = library_conn.execute("SELECT * FROM library_spectra WHERE library_spectra_meta_id IN ({})".format(
        lpids_str))

    insert_query_m([i for i in new_ls_r.fetchall()],
                   'l_s_peaks',
                   conn=conn,
                   db_type='sqlite')

    print('Get relevant library source information')
    lcms_source_id_r = conn.execute("SELECT id FROM l_source")
    lcms_source_ids = [str(int(i[0])) for i in lcms_source_id_r.fetchall()]
    lcms_source_str = ",".join(lcms_source_ids)

    new_source_r = library_conn.execute("SELECT * FROM library_spectra_source WHERE id NOT IN ({})".format(
        lcms_source_str))

    if new_source_r:
        insert_query_m([i for i in new_source_r.fetchall()],
                   'l_source',
                   conn=conn,
                   db_type='sqlite')
