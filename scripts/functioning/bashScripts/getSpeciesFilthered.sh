#!/bin/bash

grep 'species' filthering.xml > list.xml
grep -v 'subspecies' list.xml > loading.xml

