#!/bin/bash
tar xzf HWWFitting.tar.gz
echo "$*"
eval "$*"
gzip HWW*.txt
gzip HWW*.log
