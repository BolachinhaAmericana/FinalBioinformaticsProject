#!/bin/bash
grep 'species' test.xml > list.xml
grep -v 'subspecies' list.xml > loading.xml

