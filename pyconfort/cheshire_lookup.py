#!/usr/bin/env python
from __future__ import print_function

import datetime
import sys
import unicodedata
import pandas as pd

## convert HTML table to Pandas dataframe
def parse_html_table(table):
    n_rows=0; n_columns = 0; column_names = []

    ## the tables contain different numbers of columns; get their names
    rows = table.find_all('tr')
    tds = rows[1].find_all('td')

    # avoid having two column names the same -
    for i in range(0,len(tds)):
        if i == 2: column_names.append('scale_'+tds[i].get_text().split('\n')[0])
        elif i == 3 and tds[i].get_text().split('\n')[0].upper() == '13C': column_names.append('scale_'+tds[i].get_text().split('\n')[0])
        else: column_names.append(tds[i].get_text().split('\n')[0])

    n_columns = len(column_names)

    # Determine the number of rows in the table
    i = 1
    for row in rows[2:]:
        td_tags = row.find_all('td')
        if len(td_tags) == n_columns:
            n_rows+=1

    columns = column_names
    df = pd.DataFrame(columns = columns,
                      index= list(range(0,n_rows)))

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
    try:
        from bs4 import BeautifulSoup
    except (ModuleNotFoundError,AttributeError):
        log.write('\nThe bs4 module is not installed correctly - CHESHIRE search is not available')
        sys.exit()

    ## current time for printing
    now = datetime.datetime.now()

    if online == False:
        log.write("   READING FROM LOCAL VERSION OF CHESHIRE {0}".format( now.strftime("%Y-%m-%d %H:%M")))
        html = BeautifulSoup(open('./scaling_factors.html'), "lxml")
    else:
        import requests
        log.write("   READING FROM http://cheshirenmr.info/ScalingFactors.htm {0}".format(now.strftime("%Y-%m-%d %H:%M")))
        url = 'http://cheshirenmr.info/ScalingFactors.htm'
        response = requests.get(url)
        html = BeautifulSoup(response.text, "lxml")

    calc_opt = opt_method.upper()+'/'+opt_basis
    calc_nmr = nmr_method.upper()+'/'+nmr_basis

    if nmr_solv == None:
        log.write("  ", nmr_aos.upper()+'-'+calc_nmr+'//'+calc_opt)
    else:
        if opt_solv == None: log.write("  ", nmr_solv[0].upper()+'('+nmr_solv[1]+')-'+nmr_aos.upper()+'-'+calc_nmr+'//'+calc_opt)
        else: log.write("  ", nmr_solv[0].upper()+'('+nmr_solv[1]+')-'+nmr_aos.upper()+'-'+calc_nmr+'//'+opt_solv[0].upper()+'('+opt_solv[1]+')-'+calc_opt)

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
            db_nmr_solv = None; db_opt_solv = None; db_nmr_aos = 'GIAO'
            try:
                db_nmr =  row['NMR'].lower().split()[0].split("/")
                if row['NMR'].lower().find('scrf') >-1: db_nmr_solv = scrf
                if row['NMR'].lower().find('cgst') >-1: db_nmr_aos = 'CGST'
                try: [db_nmr_method, db_nmr_basis] = db_nmr
                except ValueError: pass
                if db_nmr_method[0] == '#': db_nmr_method = db_nmr_method[1:]
                db_opt =  row['Geometry'].lower().split()[0].split("/")
                if row['Geometry'].lower().find('scrf') >-1: db_opt_solv = scrf
                try: [db_opt_method, db_opt_basis] = db_opt
                except ValueError: pass
                if db_opt_method[0] == '#': db_opt_method = db_opt_method[1:]

                if db_nmr_method.lower() == nmr_method.lower() and db_nmr_basis.lower() == nmr_basis.lower() and db_nmr_aos.lower() == nmr_aos.lower():
                    if db_opt_method.lower() == opt_method.lower() and db_opt_basis.lower() == opt_basis.lower():
                        #print "matched levels of theory"
                        #print db_nmr_solv, nmr_solv, db_opt_solv, opt_solv
                        if db_nmr_solv == nmr_solv and db_opt_solv == opt_solv:
                            log.write("   --- MATCH ---", id.upper()); return row['scale_'+nucleus]
            except: pass
