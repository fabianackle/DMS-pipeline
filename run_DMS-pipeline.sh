#!/usr/bin/env bash

nextflow run main.nf \
    -resume \
    -profile local \
    -with-report \
    -with-timeline \
    -params-file LmrCD.json \
    --subsample False
