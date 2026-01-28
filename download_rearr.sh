#!/bin/bash

if ! [ -f v1.0.11.tar.gz ]
then
    echo "do"
    # wget https://github.com/ljw20180420/rearr/archive/refs/tags/v1.0.11.tar.gz
fi

tar zxf v1.0.11.tar.gz
