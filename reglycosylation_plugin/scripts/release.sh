#!/bin/bash

filename=${1:-reglycosylation_plugin.zip}

mkdir -p dist/reglycosylation_plugin
cp -r *.py *.md web dist/reglycosylation_plugin
(cd dist; zip -r "$filename" reglycosylation_plugin; rm -rf reglycosylation_plugin;)