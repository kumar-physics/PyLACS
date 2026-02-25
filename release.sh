#!/bin/sh
cd pylacs
rm -rf dist
python -m build
twine check dist/*
twine upload dist/*
