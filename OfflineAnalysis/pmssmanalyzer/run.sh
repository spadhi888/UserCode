#!/bin/sh

make program=pmssmanalyzerosleptons
./pmssmanalyzerosleptons datafile.list > tt.log 2>/dev/null

