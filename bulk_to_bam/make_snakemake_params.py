import pandas as pd
import os
import subprocess
import yaml


sample_sheet = pd.read_csv()
sample_sheet['sample_name'] = sample_sheet['sample_name'].str.replace(' ', '_')
sample_names = sample_sheet['sample_name'].tolist()


data = {
    'SAMPLES' : sample_names,
    'working_dir' : ,
    'data_dir' : ,
    'star_genome' : ,
    'sample_sheet'
}

with open(file_path, 'w') as yaml_file:
    yaml.dump(data, yaml_file)