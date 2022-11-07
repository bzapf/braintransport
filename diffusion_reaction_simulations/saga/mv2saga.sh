#!/bin/bash
set -e

echo "setting pwd to cd /home/basti/Dropbox\ \(UiO\)/Sleep/"
cd /home/basti/Dropbox\ \(UiO\)/Sleep/

for PAT in 176 218 230 199 235 236 105 175 178 183 190 191 228 240 241 227 205 215 002 078 091 127 172 249
do
	scp -r ./${PAT}/FIGURES/* saga:/cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer/${PAT}/FIGURES/
	echo $PAT
done