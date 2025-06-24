#!/bin/bash

filename=${1:-reglycosylation_plugin.zip}

mkdir -p dist
zip -r "dist/$filename" *.py *.md web