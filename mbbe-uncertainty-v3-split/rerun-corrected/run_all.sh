#!/bin/bash
# Runs 01 then 02 in this folder; log -> output/run.log. Re-running resumes 02 from the last finished cell.
cd "$(dirname "$0")"
J=${JULIA:-"julia -t 2 --project=../.."}   # set JULIA to your Pumas invocation
[ -f output/parameter_sets.jls ] || $J 01_uncertainty.jl
$J 02_grid.jl
