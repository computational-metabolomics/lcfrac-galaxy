import sys
import csv
import re
import numpy as np
from scipy.stats import rankdata
from operator import itemgetter
from db import max_row, insert_query_m

def process_ms1_annotation(wellinfo, lcms_conn, comp_conn, ms1_lookup_source, ms1_lookup_keepAdducts,
                           ms1_lookup_checkAdducts):

    # write additional peak info (camera isotope/adduct annotations and precursor ion purity)
    print('Add additional ms1 info to database')
    add_additional_peak_info(wellinfo.dims_pths.additional_info, lcms_conn, wellinfo.well_number)

    # write beams to sqlite database
    print('Add BEAMS MS1 annotation info to database')
    add_beams(wellinfo.dims_pths.beams, lcms_conn, comp_conn, ms1_lookup_keepAdducts, ms1_lookup_checkAdducts,
              ms1_lookup_source,
              wellinfo.well_number)


def process_frag_annotation(msn_pths, lcms_conn, pid_d, msn_prec_d, sirius_rank_limit):

    # write sirius to sqlite database
    if msn_pths.sirius:
        print('Process sirius annotation')
        add_sirius(msn_pths.sirius, lcms_conn, pid_d, sirius_rank_limit)

    # write metfrag to sqlite database
    if msn_pths.metfrag:
        print('Process metfrag annotation')
        add_metfrag(msn_pths.metfrag, lcms_conn, pid_d)

    # write spectral matching to sqlite database
    if msn_pths.spectral_matching:
        print('Process spectral matching annotation')
        add_spectral_matching(msn_pths.spectral_matching, lcms_conn, pid_d)

    # write spectral matching to sqlite database
    if msn_pths.mf_annotation:
        print('Process MF annotation from MSnPy')
        add_mf_annotation(msn_pths.mf_annotation, lcms_conn, msn_prec_d)

def add_beams(beams_pth, conn, comp_conn, ms1_lookup_keepAdducts,  ms1_lookup_checkAdducts, ms1_lookup_source, well):
    '''
    '''
    if not beams_pth:
        return
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

    select ="""SELECT sp.mz, sp.sid, sp.adduct, sp.isotopes
                    FROM s_peaks AS sp
                    LEFT JOIN s_peak_meta AS spm ON spm.pid=sp.pid
                    WHERE spm.name = 'original_aligned_peaklist' AND spm.well='{}'
            """.format(well)
    cursor = conn.cursor()
    r = cursor.execute(select)
    peaks_to_align = r.fetchall()
    mz_d = {round(float(i[0]), 8): {'sid': int(i[1]), 'adduct': i[2],'isotopes': i[3]} for i in peaks_to_align}


    c = max_row(conn, 'ms1_lookup_results', 'id')

    c_curs = comp_conn.cursor()
    csv.field_size_limit(sys.maxsize)
    rows = []

    with open(beams_pth, 'r') as bf:
        dr = csv.DictReader(bf, delimiter='\t')


        for drow in dr:

            # get pid from mz and well
            sidi = mz_d[round(float(drow['mz']), 8)]

            kal = ms1_lookup_keepAdducts.split(',')



            if (ms1_lookup_keepAdducts and drow['adduct'] in kal) or (ms1_lookup_checkAdducts and drow['adduct'] in
                                                                      sidi['adduct']):
                valid_adduct = True
            elif not ms1_lookup_keepAdducts and not ms1_lookup_checkAdducts:
                valid_adduct = True
            else:
                valid_adduct = False

            if not valid_adduct:
                continue


            r = c_curs.execute("""SELECT inchikey FROM {} WHERE
                                               {}='{}'
                                         """.format(table_nm,
                                                    column_nm,
                                                    drow['compound_id']))
            inchikey_r = r.fetchall()
            if not inchikey_r:
                inchikey = 'unknown'
            else:
                inchikey = inchikey_r[0][0]

            rows.append((c,
                         sidi['sid'],
                         drow['name'],
                         float(drow['mz']),
                         drow['exact_mass'],
                         drow['ppm_error'],
                         drow['adduct'],
                         drow['C'],
                         drow['H'],
                         drow['N'],
                         drow['O'],
                         drow['P'],
                         drow['S'],
                         drow['molecular_formula'],
                         drow['compound_id'],
                         inchikey,
                         1  # score is always 1 for beams
                         ))


            c += 1

    cols = "id, sid, name, mz, exact_mass, " \
           "ppm_error," \
           "adduct,C,H,N,O,P,S,molecular_formula,compound_id," \
           " inchikey, score"
    if rows:
        insert_query_m(rows,
                   'ms1_lookup_results',
                   conn,
                   columns=cols,
                   db_type='sqlite')


def add_additional_peak_info(additional_peak_info_pth, conn, well):
    '''
    '''
    select ="""SELECT sp.mz, sp.sid
                    FROM s_peaks AS sp
                    LEFT JOIN s_peak_meta AS spm ON spm.pid=sp.pid
                    WHERE spm.name = 'original_aligned_peaklist' AND spm.well='{}'
            """.format(well)
    cursor = conn.cursor()
    sid_mz_d = {round(float(i[0]), 8): int(i[1]) for i in cursor.execute(select)}

    cursor = conn.cursor()
    if not additional_peak_info_pth:
        return
    with open(additional_peak_info_pth, 'r') as bf:
        dr = csv.DictReader(bf, delimiter='\t')
        for drow in dr:
            # get pid from mz and well
            sid = sid_mz_d[round(float(drow['mz']), 8)]
            stmt = """UPDATE s_peaks SET adduct='{}',isotopes='{}',dims_predicted_precursor_ion_purity={}
                    WHERE sid={}""".format(
                                    drow['adduct'] if 'adduct' in drow else '""',
                                    drow['isotopes'] if 'isotopes' in drow else '""',
                                    drow['medianPurity'] if 'medianPurity' in drow else '""',                                    
                                    sid)
            #print(stmt)
            cursor.execute(stmt)
            conn.commit()


def add_metfrag(metfrag_pth, conn, pid_d):
    '''
    '''

    c = max_row(conn, 'metfrag_results', 'id')
    if not metfrag_pth:
        return

    with open(metfrag_pth, 'r') as bf:
        dr = csv.DictReader(bf, delimiter='\t')
        rows = []

        for drow in dr:

            pn = parse_name(drow['name'], pid_d)
            if not pn:
                continue

            rows.append((int(c),
                         pn['pid'],
                         pn['msnpy_convert_id'],
                         drow['precursor_type'] if 'precursor_type' in drow else None,
                         drow['CompoundName'] if 'CompoundName' in drow else None,
                         drow['ExplPeaks'],
                         drow['FormulasOfExplPeaks'],
                         drow['FragmenterScore'],
                         drow['FragmenterScore_Values'],
                         drow['Identifier'],
                         drow['InChI'],
                         drow['InChIKey'],
                         drow['InChIKey1'],
                         drow['InChIKey2'],
                         drow['InChIKey3'] if 'InChiKey3'in drow else None,
                         drow['MaximumTreeDepth'],
                         drow['MolecularFormula'],
                         drow['MonoisotopicMass'],
                         drow['NoExplPeaks'],
                         drow['NumberPeaksUsed'],
                         drow['SMILES'],
                         drow['Score']
                         ))
            c += 1

    cols = "id, pid, msnpy_convert_id, adduct, CompoundName, ExplPeaks," \
           "FormulasOfExplPeaks, FragmenterScore, FragmenterScore_Values, " \
           "Identifier, InChI, InChIKey, InChIKey1, InChIKey2, InChIKey3," \
           "MaximumTreeDepth, MolecularFormula,MonoisotopicMass," \
           "NoExplPeaks,NumberPeaksUsed, SMILES, Score"

    if rows:
        insert_query_m(rows,
                   'metfrag_results',
                   conn,
                   columns=cols,
                   db_type='sqlite')


def add_spectral_matching(spectral_matching_pth, conn, pid_d):
    '''
    '''
    c = max_row(conn, 'sm_matches', 'mid')
    if not spectral_matching_pth:
        return

    with open(spectral_matching_pth, 'r') as bf:
        dr = csv.DictReader(bf, delimiter='\t')
        rows = []

        for drow in dr:
            pn = parse_name(drow['query_entry_name'], pid_d)
            if not pn:
                continue

            rows.append((int(c),
                         pn['pid'],
                         pn['msnpy_convert_id'],
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
    if rows:
        insert_query_m(rows,
                   'sm_matches',
                   conn,
                   columns=cols,
                   db_type='sqlite')



def add_mf_annotation(mf_annotation_pth, conn, msn_prec_d):
    '''
    '''

    c = max_row(conn, 'mf_annotations', 'mfaid')
    if not mf_annotation_pth:
        return

    # get the DIMSn peaks fo this file?

    with open(mf_annotation_pth, 'r') as bf:
        dr = csv.DictReader(bf, delimiter='\t')
        rows = []


        for drow in dr:

            if drow['mz'] and round(float(drow['mz']), 8) in msn_prec_d:
                sid = msn_prec_d[round(float(drow['mz']), 8)]
            else:
                sid = ''

            rows.append((int(c),
                         sid,
                         drow['mz'],
                         drow['tree_id'],
                         drow['group_id'],
                         drow['scan_events'],
                         drow['max_mslevel'],
                         drow['mf_id'],
                         drow['molecular_formula'],
                         drow['adduct'],
                         drow['mass'],
                         drow['ppm_error'],
                         drow['rank'],
                         drow['total_ranks'],
                         drow['ranked_equal'],
                         drow['trees'],
                         drow['neutral_losses_explained']

                         ))
            c += 1


    cols = """mfaid,sid,mz,tree_id,scan_events,max_mslevel,mf_id,molecular_formula,
              adduct,mass,ppm_error,rank,total_ranks,ranked_equal,trees,neutral_losses_explained
              ,filename"""
    if rows:
        insert_query_m(rows,
                       'mf_annotations',
                       conn,
                       columns=cols,
                       db_type='sqlite')


def neg_min_max(x):
    if np.all(x == x[0]):
        return [1] * len(x)
    elif len(x)==1:
        return [1]
    else:
        x = np.array(x)
        xn = (x - min(x)) / (max(x) - min(x))
        return abs(xn - 1)

def parse_name(name, pid_d):
    mtch = re.search('^header (.*)\|.*msnpy_convert_id (\d+).*', name)
    if not mtch:
        return None
    else:
        header = mtch.group(1)
        msnpy_convert_id = int(mtch.group(2))
        pid = pid_d[msnpy_convert_id]
    return {'h':header, 'msnpy_convert_id':msnpy_convert_id, 'pid':pid}


def add_sirius(sirius_pth, conn, pid_d, rank_limit=25):
    '''
    '''
    if not sirius_pth:
        return

    with open(sirius_pth, 'r') as bf:
        dr = csv.DictReader(bf, delimiter='\t')
        rows = []
        c = max_row(conn, 'sirius_csifingerid_results', 'id')
        annotation_group = {}
        for drow in dr:
            pn = parse_name(drow['name'], pid_d)
            if not pn:
                continue
            else:
                pid = pn['pid']
                msnpy_convert_id = pn['msnpy_convert_id']

            if rank_limit != 0 and int(drow['rank']) > rank_limit:
                continue
            drow = {k.lower():v for k,v in drow.items()}
            
            rows.append((c,
                         pid,
                         msnpy_convert_id,
                         drow['adduct'] if 'adduct' in drow else None,
                         drow['inchikey2d'],
                         drow['inchi'],
                         drow['molecularformula'],
                         drow['rank'],
                         drow['score'],
                         drow['name'],
                         drow['smiles'],
                         drow['xlogp'],
                         drow['pubchemids'],
                         drow['links']
                         ))
            if pid not in annotation_group:
                annotation_group[pid] = {}
            annotation_group[pid][c] = float(drow['score'])
            c += 1

    ##################
    # calculate bound sirius score
    ##################
    print('CALCULATE BOUND SIRIUS SCORE')
    bounded_score_d = {}
    for pid, scores in annotation_group.items():
        score_l = []
        cid_l = []
        for cid, score in scores.items():
            # Make it absolute value (as that is what we do for the msPurity equiv.
            # but essentially we could have just used a min max function that is isn't 
            # negative... the outcome is the same though
            score_l.append(abs(float(score)))
            cid_l.append(cid)
        bounded_score = neg_min_max(score_l)
        
        for i in range(0, len(bounded_score)):
            bounded_score_d[cid_l[i]] = bounded_score[i]
            if not bounded_score[i]:
                print('bounded_score none', i, cid_l[i])
            else:
                bounded_score_d[cid_l[i]] = bounded_score[i]


    # add bounded score to rows
    rows = [row + (bounded_score_d[row[0]],) for row in rows]

    cols = "id, pid, msnpy_convert_id, adduct, inchikey2D, InChI, " \
           "molecularFormula," \
           "Rank, Score, Name, smiles, xlogp, pubchemids, links, bounded_score"
    
    if rows:
      insert_query_m(rows, 'sirius_csifingerid_results', conn, columns=cols,
                   db_type='sqlite')





def col_multiple_check(r, cols):
    if "," in cols or "*" == cols:
        return r
    else:
        return [i[0] for i in r]


def get_metab_compound_rows(comp_conn, inchikeys, inchi_level='inchikey',
                            cols='*'):
    c = comp_conn.cursor()
    
    inchi_str = "','".join(i for i in inchikeys if i)

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
    # perhaps update to only keep the best score (in case doing scan by scan)

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
    # multiple scans or energies for the same precurosr
    # that give the same inchikey annotation we
    # only keep the result one with the best score (highest) for the combined
    # table
    for i in r:
        if not i[0]:
            # no annotation (with inchikey  so skip)
            continue
        inchi_sid = "{}_{}".format(i[0], i[1])
        rowd = {'mid': i[2], 'score': i[3], 'wscore': i[4], 'adduct': i[5]}
        if inchi_sid in inchi_sid_d:
            if table_nm in inchi_sid_d[inchi_sid]:
                if rowd['wscore']:
                    #print(inchi_sid_d[inchi_sid][table_nm], file=sys.stderr)
                    #print(rowd, file=sys.stderr)
                    #print(table_nm, file=sys.stderr)
                    if not inchi_sid_d[inchi_sid][table_nm]['wscore'] or float(rowd['wscore']) > float(inchi_sid_d[inchi_sid][table_nm]['wscore']):
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


def get_inchikey_sid_beams(conn, weight, inchi_sid_d):
    # Note we want to get the sid for the original peak list
    c = conn.cursor()
    c.execute("""SELECT inchikey, sid, id, score, score*{} AS wscore, adduct
                    FROM ms1_lookup_results   WHERE sid NOT NULL  
    """.format(weight))
    r = c.fetchall()
    table_nm = 'ms1_lookup_results'
    return inchi_sid_d_update(r, inchi_sid_d, table_nm)



def add_to_row(results, name, row):
    if name in results:
        sr = results[name]
        row.extend([sr['mid'], sr['score'], sr['wscore'], sr['adduct']])
    else:
        row.extend([0, 0, 0, 0])
        sr = {'mid': 0, 'score': 0, 'wscore': 0, 'adduct': 0}
    return row, sr


def get_metab_compound_all(conn):
    c = conn.cursor()
    c.execute("SELECT * FROM metab_compound")
    metab_rows = c.fetchall()
    return metab_rows


def combine_annotations(conn, comp_conn, weights):
    # First get all the compound information for each inchikey from each
    # result (metfrag, sirius, spectral matching)
    # First get inchikeys from sirius result
    print('get inchikeys')
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
    
    
    # Remove inchikeys we have already
    old_inchikeys = sql_simple_select(conn, 'metab_compound', 'inchikey')

    print(old_inchikeys)
    inchikeys = [i for i in inchikeys if i not in old_inchikeys]
    
    # add to metab_compound in conn
    print('Add new inchikeys to compound')
    metab_rows = get_metab_compound_rows(comp_conn, inchikeys, 'inchikey', '*')

    # Add rows to sqlite results database
    insert_query_m(metab_rows, 'metab_compound', conn, columns=None,
                   db_type="sqlite", ignore_flag=True)
    
    # Get all metab_rows
    metab_rows = get_metab_compound_all(conn)
    
    # biosim dict
    biosim_d = {row[0]: row[len(row) - 3] for row in metab_rows}
    #print(metab_rows)
    
    # for each sid - get the (best) results for each inchikey and annotations
    # approach
    # make combined dict
    # get all inchikey_sid combinations from all approaches.
    # Then loop through and check each approach and make a row for each!
    print('Get sids')
    inchi_sid_d = get_inchikey_sid_sirius(conn,
                                          weight=weights['sirius_csifingerid'])
    inchi_sid_d = get_inchikey_sid_beams(conn,
                                   weight=weights['beams'],
                                   inchi_sid_d=inchi_sid_d)
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
    print('Combine sids')
    sid_d = {}
    for inchi_sid, results in inchi_sid_d.items():
        row = inchi_sid.split("_")
        if row[0].lower() == 'unknown' or row[1] == 'None':
            continue
        sid = row[1]
        row, sirius_d = add_to_row(results, 'sirius_csifingerid_results', row)
        row, metfrag_d = add_to_row(results, 'metfrag_results', row)
        row, beams_d = add_to_row(results, 'ms1_lookup_results', row)
        row, sm_d = add_to_row(results, 'sm_matches', row)

        # Get biosim score
        biosim_score = biosim_d[row[0]] if row[0] in biosim_d else 0
        biosim_score = biosim_score if biosim_score else 0
        print(biosim_score, weights['biosim'])
        biosim_wscore = biosim_score * weights['biosim']

        # calculated weighted score
        print(sirius_d['wscore'], metfrag_d['wscore'], beams_d['wscore'], sm_d['wscore'], biosim_wscore, file=sys.stderr)
        wscore = (sirius_d['wscore'] if sirius_d['wscore'] else 0) + \
                 (metfrag_d['wscore'] if metfrag_d['wscore'] else 0) + \
                 (beams_d['wscore'] if beams_d['wscore'] else 0) + \
                 (sm_d['wscore'] if sm_d['wscore'] else 0) + \
                 biosim_wscore

        # get overall adduct column
        adducts = set([str(sirius_d['adduct']),
                       str(metfrag_d['adduct']),
                       str(beams_d['adduct']),
                       str(sm_d['adduct'])])
        adducts = [a for a in adducts if a]
        adducts.sort()

        if '0' in adducts:
            adducts.remove('0')

        if 'None' in adducts:
            adducts.remove('None')

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

        ranks = rankdata(1-np.array(scores), 'dense')

        for i in range(0, len(ranks)):
            rows[i].append(int(ranks[i]))
        rows = sorted(rows, key=itemgetter(len(rows[0]) - 1))
        print(rows)
        final_rows.extend(rows)

    columns = """
    inchikey,sid,
    sirius_id,sirius_score,sirius_wscore,sirius_adduct,
    metfrag_id,metfrag_score,metfrag_wscore,metfrag_adduct,
    ms1_lookup_id,ms1_lookup_score,ms1_lookup_wscore,ms1_lookup_adduct,
    sm_mid,sm_score,sm_wscore,sm_adduct,
    biosim_max_score, biosim_wscore,
    adduct_overall,wscore,rank
    """
    print('insert combined results')
    print(len(final_rows[0]))

    insert_query_m(final_rows, "combined_annotations", conn,
                   db_type='sqlite', columns=columns)


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
       round(sp.mz, 8) AS mz,
       round(sp.i, 2) AS i,
       round(spm.well_rt, 3) AS rt,
       round(spm.well_rtmin,3) AS rtmin,
       round(spm.well_rtmax,3) AS rtmax,
       sp.adduct AS camera_adduct,
       sp.isotopes AS camera_isotopes,
       GROUP_CONCAT(DISTINCT (CAST (cpgXsp.grpid AS INTEGER) ) ) AS lc_grpid_mtchs,
       '' AS dims_sid_mtchs,
       spm.well, 
       spm.name AS spm_name, 
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
       round(sp.dims_predicted_precursor_ion_purity, 3) AS precursor_ion_purity,
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

    WHERE (sp.sid IS NOT NULL) AND  (
    spm.name IS 'original_aligned_peaklist') AND (0=={} OR IFNULL(ca.rank<={}, 1))
    GROUP BY sp.sid,
             IFNULL(ca.inchikey, sp.sid)
    ORDER BY spm.well,
             sp.mz,
             IFNULL(ca.rank, sp.mz)
             
    """.format(rank_limit,rank_limit)

    r = c.execute(sql_stmt)
    dims_annotations = r.fetchall()

    sql_stmt = """    SELECT 
       'lcms' AS ms_type,
       '' AS sid,
       cpg.grpid,
       cpg.grp_name,
       round(cpg.mz, 8) AS mz,
       ROUND(AVG(cp._into),3) AS i,
       round(cpg.rt, 3) AS rt,
       round(cpg.rtmin,3) AS rtmin,
       round(cpg.rtmax,3) AS rtmax,
       cpg.adduct AS camera_adduct,
       cpg.isotopes AS camera_isotopes,

       '' AS lcms_grpid_mtchs,
       GROUP_CONCAT(DISTINCT (CAST (cpgXsp.sid AS INTEGER) ) ) AS dims_sid_mtchs,
       spm.well,
       spm.name AS spm_name,  
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
       ROUND(AVG(spm.inPurity),3) AS precursor_ion_purity,
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

    WHERE (cpg.grpid IS NOT NULL) AND IFNULL((fi.class NOT LIKE '%blank%'), 1)
        AND (0=={} OR IFNULL(ca.rank<={}, 1))
    GROUP BY
     cpg.grpid,
     IFNULL(ca.inchikey, cpg.grpid)
    ORDER BY 
     cpg.grpid, 
     IFNULL(ca.rank, cpg.grpid)
    """.format(rank_limit, rank_limit)

    r = c.execute(sql_stmt)
    lcms_annotations = r.fetchall()
    colnames = list(map(lambda x: x[0], c.description))

    with open(out_file, 'w') as summary_f:

        sw = csv.writer(summary_f, delimiter='\t')
        sw.writerow(colnames)
        sw.writerows(lcms_annotations)
        sw.writerows(dims_annotations)
