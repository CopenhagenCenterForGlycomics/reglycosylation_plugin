#!/bin/bash

filename=${1:-reglycosylation_plugin.zip}

mkdir -p dist
bsdtar -c -a -f "dist/$filename" -s ",^,reglycosylation_plugin/," *.py *.md web