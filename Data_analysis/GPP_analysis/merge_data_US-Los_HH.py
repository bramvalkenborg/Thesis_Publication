import pandas as pd
import glob

input_dir = '/data/leuven/336/vsc33653/Data/GPP_FluxTower/'
list_data = glob.glob(input_dir+'US-Los_HH_*.csv')

data = pd.read_csv(list_data[0])

for file in range(1, len(list_data)):
    df = pd.read_csv(list_data[file])
    data = data.append(df)

data.to_csv(input_dir+'US-Los_HH_2000-2022.csv')

