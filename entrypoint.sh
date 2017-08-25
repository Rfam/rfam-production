#!/bin/sh
source $LOC/virtualenvs/rfam-production/bin/activate
cd $PYTHONPATH
exec "$@"
