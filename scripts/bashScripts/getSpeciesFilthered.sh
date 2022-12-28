#!/bin/bash

#Made by Valente
grep 'species' test.xml > list.xml
grep -v 'subspecies' list.xml > loading.xml

