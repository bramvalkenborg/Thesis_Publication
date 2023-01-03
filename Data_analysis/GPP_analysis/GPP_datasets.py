import pandas as pd
import numpy as np
def load_dataset(input_dir, input_file):
    if input_file == 'FLX_SE-Deg_FLUXNET2015_FULLSET_HH_2001-2020_beta-3.csv':
        time_string = 'TIMESTAMP_END'
        format_time = '%Y%m%d%H%M'
        WTD_string = ''
        wtd_file = 'SE_DEG.csv'
        PAR_string = 'SW_IN_F'
        GPP_string = 'GPP_DT_CUT_REF'
        cm_to_m = False
        mean_elev_ref_correction = 0.0
    elif input_file == 'CA_MER_GPP_analysis.csv':
        time_string = 'TIMESTAMP_END'
        format_time = ''
        WTD_string = 'WTD'
        PAR_string = 'PAR'
        GPP_string = 'GPP'
        cm_to_m = True
        mean_elev_ref_correction = 0.0
    elif input_file == 'US-Los_HH_2000-2022.csv':
        time_string = 'TIMESTAMP_END'
        format_time = '%Y%m%d%H%M'
        WTD_string = 'WTD'
        PAR_string = 'SW_IN_1_1_1'
        GPP_string = 'GPP_F'
        cm_to_m = False
        mean_elev_ref_correction = 0.0
    elif input_file == 'mb_met_flux_data_1998_2018.txt':
        time_string = 'x_tv'
        format_time = ''
        WTD_string = 'wt'
        PAR_string = 'solar_f'
        GPP_string = 'gep_f'
        cm_to_m = True
        # Mer Bleue reference for water table depth is the elevation of hummocks
        # correction by 0.125 cm as in literature the mean difference between hollow and hummocks is 0.25 m
        mean_elev_ref_correction = 0.125
    elif input_file == 'SFZ_S_2010_2020_hh_Michl.csv':
        time_string = 'date [end]'
        format_time = '%d.%m.%Y %H:%M'
        WTD_string = ''
        wtd_file = 'mes_logger.txt'
        PAR_string = 'PAR'
        GPP_string = 'GPP'
        cm_to_m = False
        mean_elev_ref_correction = 0.0
    elif input_file == 'FI_LOM_GPP.csv':
        time_string = 'Date UTC (end time)'
        format_time = '%d/%m/%Y %H:%M'
        WTD_string = ''
        wtd_file = 'WTD_Lompolo_2006-2019.csv'
        PAR_string = 'PPFD'
        GPP_string = 'GPP modelled (mgCO2/m2s)'
        cm_to_m = False
        mean_elev_ref_correction = 0.0


    # Load the data
    Df = pd.read_csv(input_dir+input_file)
    if format_time == '':
        Df['Time'] = pd.to_datetime(Df[time_string])
    else:
        Df['Time'] = pd.to_datetime(Df[time_string], format=format_time)
    del Df[time_string]
    Df.replace(-9999, np.nan, inplace=True)
    Df = Df.set_index(['Time'])
    if WTD_string!='':
        WTD = pd.Series(Df[WTD_string])
    elif wtd_file=='SE_DEG.csv':
        Df_wtd = pd.read_csv(input_dir+wtd_file)
        Df_wtd['Time'] = pd.to_datetime(Df_wtd['Date']+' 12:00:00', format='%d/%m/%Y %H:%M:%S')
        Df_wtd = Df_wtd.set_index('Time')
        Df_wtd = Df_wtd.drop(columns=('Date'))
        Df = pd.concat([Df,Df_wtd],axis=1)
        Df_int = Df.interpolate(method='nearest',limit=25,limit_direction='both')
        WTD = Df_int['WTD']
    elif wtd_file=='mes_logger.txt':
        Df_wtd = pd.read_csv(input_dir+wtd_file)
        Df_wtd['Time'] = pd.to_datetime(Df_wtd['DateDaytime']+' 12:00:00', format='%Y-%m-%d %H:%M:%S')
        Df_wtd = Df_wtd.set_index('Time')
        Df_wtd = Df_wtd.drop(columns=('DateDaytime'))
        Df = pd.concat([Df,Df_wtd],axis=1)
        Df_int = Df.interpolate(method='nearest',limit=25,limit_direction='both')
        WTD = Df_int['wl_m']
    elif wtd_file=='WTD_Lompolo_2006-2019.csv':
        Df_wtd = pd.read_csv(input_dir+wtd_file)
        Df_wtd['Time'] = pd.to_datetime(Df_wtd['Date']+' 12:00:00', format='%d.%m.%Y %H:%M:%S')
        Df_wtd = Df_wtd.set_index('Time')
        Df_wtd = Df_wtd.drop(columns=('Date'))
        Df = pd.concat([Df,Df_wtd],axis=1)
        Df_int = Df.interpolate(method='nearest',limit=25,limit_direction='both')
        WTD = Df_int['WTD']


    GPP = pd.Series(Df[GPP_string])
    PAR = pd.Series(Df[PAR_string])
    if cm_to_m:
        WTD = WTD/100
    WTD = WTD + mean_elev_ref_correction

    # convert from PPFD to Watt/sec (typically used for PAR?)
    # µmol m-2 s-1  --> J m-2 s-1, which is W m-2
    if PAR_string=='PPFD':
        PAR = 0.327 * PAR


    return WTD, PAR, GPP
