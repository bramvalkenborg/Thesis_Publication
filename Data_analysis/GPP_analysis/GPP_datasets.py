def load_dataset(file):
    if file == 'FLX_SE-Deg_FLUXNET2015_FULLSET_HH_2001-2020_beta-3_WTD.csv':
        time_string = 'TIMESTAMP_END'
        format_time = '%Y%m%d%H%M'
        WTD_string = 'WTD'
        PAR_string = 'SW_IN_F'
        GPP_string = 'GPP_DT_CUT_REF'
    elif file == 'CA_MER_GPP_analysis.csv':
        time_string = 'TIMESTAMP_END'
        format_time = ''
        WTD_string = 'WTD'
        PAR_string = 'PAR'
        GPP_string = 'GPP'
    elif file == 'US-Los_HH_2000-2022.csv':
        time_string = 'TIMESTAMP_END'
        format_time = '%Y%m%d%H%M'
        WTD_string = 'WTD'
        PAR_string = 'SW_IN_1_1_1'
        GPP_string = 'GPP_F'
    elif file == 'mb_met_flux_data_1998_2018.txt':
        time_string = 'x_tv'
        format_time = ''
        WTD_string = 'wt'
        PAR_string = 'solar_f'
        GPP_string = ''
    return time_string, format_time, WTD_string, PAR_string, GPP_string
