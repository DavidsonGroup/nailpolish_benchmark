#!/bin/sh
awk 'NR % 4 == 1 {
  if (match($0, /^@([ATCGNX]+)_#([ATCGNX]+)_#([ATCGNX]+)_([ATCGNX]+)#([^_ \t]+)/, m))
    $0 = "@" m[5] "\tCB:Z:" m[1] m[2] m[3] "\tUB:Z:" m[4]
} 1' SRR13948564.DMX.fastq > lr-split-seq-tagged.fastq
