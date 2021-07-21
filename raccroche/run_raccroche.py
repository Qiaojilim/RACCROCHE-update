#!/usr/bin/env python3

import os
import subprocess
import sys
import yaml


def main(conf='config.yaml'):
    cwd = os.path.dirname(os.path.realpath(__file__)) 

    with open(conf, 'r') as stream:
        try:
            conf_dict = yaml.safe_load(stream)
            print(conf_dict['data.path'])
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit(1)
   
    # module 3
    
    
    
    print('===================================================') 
    print('==== Analyzing contig co-occurence ====')
    print('==== This step is very slow ... ')
    proc = subprocess.Popen(['Rscript', 'module3/analyze_autoC.R'], cwd=cwd)
    proc.wait()
    print('==== Check out results/clustering folder for clustering results')
    
    
    

if __name__ == '__main__':
    main()
