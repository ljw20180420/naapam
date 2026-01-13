#!/bin/bash

parallel -a <(find /home/ljw/sdb1/naapam/unique/ -name "*.unique") --jobs 24 ./run_parse.py
