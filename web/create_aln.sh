#!/bin/sh
egrep -v '^>' $1 | sed 's/[a-z]//g'
