#!/bin/sh

sphinx-apidoc -f -o ./ ../src
make html
