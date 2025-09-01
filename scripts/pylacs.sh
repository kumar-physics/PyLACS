#!/bin/bash
tar -xzf my_venv.tar.gz
source my_venv/bin/active
python3 lacs.py $1 --data-id $2 --out /home/nmrbox/kbaskaran/PycharmProjects/PyLACS/scratch/out --method bayes
