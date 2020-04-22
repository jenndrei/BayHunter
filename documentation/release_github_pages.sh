#!/bin/bash
make clean
make html
rsync -arv build/html/ ../docs/
