#!/bin/bash

root_dir=$1
parallel -a <(find "${root_dir}/unique/" -name "*.unique") --jobs 24 naapam_run_parse.py
