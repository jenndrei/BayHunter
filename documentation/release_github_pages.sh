#!/bin/bash
make clean
make html
rsync -arv build/html/ ../docs/
touch ../docs/.nojekyll