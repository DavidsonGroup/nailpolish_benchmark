#!/usr/bin/env bash

awk 'NR%4==1 {
    if ($0 ~ /-1$/) tag = "GG:Z:simplex"
    else tag = "GG:Z:duplex"
    header = ">" substr($0, 2) "\t" tag
    next
}

NR%4==2 { print header; print; next }
NR%4==3 || NR%4==0 { next }
' "$1" > "$2"