#!/bin/sh
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -V
$1
# Usage:
#           qsub /path/to/sge_run.sh  command

