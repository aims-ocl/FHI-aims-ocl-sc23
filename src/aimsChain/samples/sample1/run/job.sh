#!/bin/bash
rm -r climbing_forces.log forces.log growing_forces.log iterations paths optimized grownstring

runchain.py

../getpath.sh

