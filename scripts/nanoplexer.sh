#!/bin/bash

source activate nanopore
nanoplexer \
	-b barcode.fa \
	-t 16 \
	-p raw_data \
	pass.fastq.gz
