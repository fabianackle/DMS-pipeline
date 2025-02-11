#!/usr/bin/env bash

nextflow run main.nf \
    -profile local \
    -with-report \
    -with-timeline \
    -params-file LmrCD.json \
    --subsample True
