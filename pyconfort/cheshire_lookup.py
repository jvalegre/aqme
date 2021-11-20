import datetime
import sys
import pandas as pd

## convert HTML table to Pandas dataframe
def parse_html_table(table):
    """
    takes in an html table from the cheshire database and transforms it into 
    a pandas.Dataframe 

    Parameters
    ----------
    table : bs4.BeautifulSoup
        A data structure obtained from parsing an html file containing a table 
        of scaling factors. 

    Returns
    -------
    pandas.Dataframe
        A Dataframe version of the table contained in the html
    """
    n_rows = 0
    n_columns = 0
    column_names = []

    ## the tables contain different numbers of columns; get their names
    rows = table.find_all('tr')
    tds = rows[1].find_all('td')

    # avoid having two column names the same -
    for i in range(0,len(tds)):
        column_name = tds[i].get_text().split('\n')[0]
        if i == 2: 
            column_names.append(f'scale_{column_name}')
        elif i == 3 and column_name.upper() == '13C': 
            column_names.append(f'scale_{column_name}')
        else: 
            column_names.append(column_name)

    n_columns = len(column_names)

    # Determine the number of rows in the table
    i = 1
    for row in rows[2:]:
        td_tags = row.find_all('td')
        if len(td_tags) == n_columns:
            n_rows+=1

    columns = column_names
    df = pd.DataFrame(columns=columns,index=list(range(n_rows)))

    row_marker = 0
    for row in rows[2:]:
        column_marker = 0
        columns = row.find_all('td')
        for column in columns:
            #print row_marker, column_marker, ' '.join(column.get_text().split())
            df.iat[row_marker,column_marker] = ' '.join(column.get_text().split())
            column_marker += 1
        if len(columns) > 0:
            row_marker += 1
    return df

## Parse from local offline version of the CHESHIRE scaling factors html or directly from the web
def cheshire(online, nucleus, opt_method, opt_basis, opt_solv, nmr_method, nmr_basis, nmr_solv, nmr_aos,log):
    """
    Checks the cheshire database (either downloading or from a pre-downloaded 
    file and returns the scaling factor for the nucleus specified if the theory
    provided and the ones in the database match.

    Parameters
    ----------
    online : bool
        Control the download or not of the cheshire database.
    nucleus : str
        [description]
    opt_method : str
        method of optimization
    opt_basis : str
        basis of optimization
    opt_solv : str or None
        solvent used in the optimization
    nmr_method : str
        method 
    nmr_basis : str
        basis used for nmr prediction
    nmr_solv : str
        solvent used for nmr prediction
    nmr_aos : str
        aos used for nmr prediction
    log : pyconfort.Logger
        Log instance.

    Returns
    -------
    float ?
        scaling factor
    """
    try:
        from bs4 import BeautifulSoup
    except (ModuleNotFoundError,AttributeError):
        log.write('\nThe bs4 module is not installed correctly - CHESHIRE search is not available')
        sys.exit()

    ## current time for printing
    now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')

    if online:
        import requests
        url = 'http://cheshirenmr.info/ScalingFactors.htm'
        log.write(f"   READING FROM {url} {now}")
        response = requests.get(url)
        html = BeautifulSoup(response.text, "lxml")
    else:
        log.write(f"   READING FROM LOCAL VERSION OF CHESHIRE {now}")
        with open('./scaling_factors.html','r') as F: 
            txt = F.read()
        html = BeautifulSoup(txt, "lxml")

    calc_opt = f'{opt_method.upper()}/{opt_basis}'
    calc_nmr = f'{nmr_method.upper()}/{nmr_basis}'
    
    # Log the theory
    mid_name= f'{nmr_aos.upper()}-{calc_nmr}'
    if nmr_solv == None:
        top_name = ''
        end_name = calc_opt
        log.write(f'  {mid_name}//{calc_opt}')
    else:
        top_name = f'{nmr_solv[0].upper()}({nmr_solv[1]})-'
        if opt_solv == None:
            end_name = calc_opt
        else:
            end_name = f'{opt_solv[0].upper()}({opt_solv[1]})-{calc_opt}'
    log.write(f'  {top_name}{mid_name}//{end_name}')

    for table in html.find_all('table'):

        id = table['id']

        scaling_table = parse_html_table(table)

        # solvent details for the CHESHIRE database
        # manually entered would be better to parse from HTML - will add in due course
        if id == 'table1a': scrf = ['pcm', 'acetone']
        elif id == 'table1b': scrf = ['smd', 'chloroform']
        elif id == 'table1c': scrf = ['cpcm', 'chloroform'] #UAKS radii, nosymmcav
        elif id == 'table1d': scrf = ['smd', 'chloroform']
        elif id == 'table2': scrf = ['pcm', 'chloroform']
        elif id == 'table3a': scrf = ['pcm', 'toluene']
        elif id == 'table5-acetone': scrf = ['pcm', 'acetone']
        elif id == 'table5-acetonitrile': scrf = ['pcm', 'acetonitrile']
        elif id == 'table5-benzene': scrf = ['pcm', 'benzene']
        elif id == 'table5-chloroform': scrf = ['pcm', 'chloroform']
        elif id == 'table5-dichloromethane': scrf = ['pcm', 'dichloromethane']
        elif id == 'table5-dimethylsulfoxide': scrf = ['pcm', 'dimethylsulfoxide']
        elif id == 'table5-methanol': scrf = ['pcm', 'methanol']
        elif id == 'table5-tetrahydrofuran': scrf = ['pcm', 'tetrahydrofuran']
        elif id == 'table5-toluene': scrf = ['pcm', 'toluene']
        elif id == 'table5-water': scrf = ['pcm', 'water']
        elif id == 'table7': scrf = ['smd', 'chloroform']
        else: scrf = None

        # Look for a match between calculation and database (case insensitive)
        # Returns the first match and then breaks
        for index, row in scaling_table.iterrows():
            # Get database's NMR theory
            db_nmr = row['NMR'].lower().split()[0].split("/")
            if 'scrf' in row['NMR'].lower():
                db_nmr_solv = scrf
            else:
                db_nmr_solv = None
            same_nmr_solvent = db_nmr_solv == nmr_solv

            if 'cgst' in row['NMR'].lower():
                db_nmr_aos = 'CGST'
            else:
                db_nmr_aos = 'GIAO'
            same_nmr_aos = db_nmr_aos.lower() == nmr_aos.lower()

            try: 
                [db_nmr_method, db_nmr_basis] = db_nmr
            except ValueError: 
                same_nmr_theory = False
            else:
                if db_nmr_method[0] == '#': 
                        db_nmr_method = db_nmr_method[1:]
                same_nmr_method = db_nmr_method.lower() == nmr_method.lower()
                same_nmr_basis = db_nmr_basis.lower() == nmr_basis.lower()
                same_nmr_theory = same_nmr_method and same_nmr_basis and same_nmr_aos
            
            # Get database's OPT theory
            db_opt = row['Geometry'].lower().split()[0].split("/")
            if 'scrf' in row['Geometry'].lower(): 
                db_opt_solv = scrf
            else:
                db_opt_solv = None
            same_opt_solvent = db_opt_solv == opt_solv

            try: 
                [db_opt_method, db_opt_basis] = db_opt
            except ValueError: 
                pass
            else:
                if db_opt_method[0] == '#':
                    db_opt_method = db_opt_method[1:]
                same_opt_method = db_opt_method.lower() == opt_method.lower()
                same_opt_basis = db_opt_basis.lower() == opt_basis.lower()
                same_opt_theory = same_opt_method and same_opt_basis

            if same_nmr_theory and same_opt_theory and same_nmr_solvent and same_opt_solvent:
                log.write(f'   --- MATCH --- {id.upper()}')
                return row['scale_'+nucleus]
