import sqlite3

def update_lcms_db(conn):
    """
    """
    c = conn.cursor()
    # this means that we are missing the dims specific columns
    c.execute('ALTER TABLE s_peak_meta ADD msnpy_convert_id integer')
    c.execute('ALTER TABLE s_peak_meta ADD well text')
    c.execute('ALTER TABLE s_peak_meta ADD well_rtmin real')
    c.execute('ALTER TABLE s_peak_meta ADD well_rtmax real')
    c.execute('ALTER TABLE s_peak_meta ADD well_rt real')

    c.execute('ALTER TABLE s_peaks ADD adduct real')
    c.execute('ALTER TABLE s_peaks ADD isotopes real')
    c.execute('ALTER TABLE s_peaks ADD dims_predicted_precursor_ion_purity real')

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


    c.execute('''CREATE TABLE mf_annotations(
                                     mfaid integer PRIMARY KEY,
                                     mz real,
                                     sid integer,
                                     tree_id integer,
                                     scan_events text,
                                     max_mslevel integer,
                                     mf_id integer,
                                     molecular_formula text,
                                     adduct text,
                                     mass real,
                                     ppm_error real,
                                     rank integer,
                                     total_ranks integer,
                                     ranked_equal integer,
                                     trees integer,
                                     neutral_losses_explained integer,
                                     filename text,
                                     FOREIGN KEY(sid) REFERENCES
                                     s_peaks(sid)
                                )'''
              )


    r = c.execute('PRAGMA table_info(ms1_lookup_results)')
    cols = r.fetchall()
    if not cols:
        c.execute('''CREATE TABLE ms1_lookup_results(
                                                id integer PRIMARY KEY,
                                                grpid integer,
                                                sid integer,
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
        c.execute('ALTER TABLE ms1_lookup_results ADD sid integer')
        c.execute('ALTER TABLE ms1_lookup_results ADD score real')
        # c.execute('ALTER TABLE ms1_lookup_results ADD inchikey text')

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
                                        inchi text,
                                        molecularFormula text,
                                        rank integer,	
                                        score real,
                                        name text,
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
                                        name text,	
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
        #c.execute('ALTER TABLE metfrag_results ADD mz real')

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
        #c.execute('ALTER TABLE combined_annotations ADD ms1_lookup_id integer')

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
