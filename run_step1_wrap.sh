#!/bin/bash
python $(realpath "$(dirname "${BASH_SOURCE[0]}")")/step1_hit_perLumi_analysis.py "$@"
