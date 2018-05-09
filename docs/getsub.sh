#!/bin/bash
file=$1
gawk '/^subroutine / {print $1" "$2}' $file
gawk '/\) function / {print $1" "$2}' $file
