#!/bin/sh

awk '{sub(/;size=.*;/,"")}1' ASVs.fasta > ASVs_clean.fasta
