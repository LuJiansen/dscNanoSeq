#!/bin/bash
nanoplexer \
	-b barcode.fa \
	-t 16 \
	-p raw_data \
	pass.fastq.gz

pigz raw_data/*fastq