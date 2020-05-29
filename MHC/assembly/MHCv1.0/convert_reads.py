from __future__ import print_function
from FastaReader import FastaReader


f = FastaReader("tmp.fa")
count = 0
for r in f:
    rlen = len(r.sequence)
    print(">ccs/{}/{}_{}".format(count, 0, rlen))
    for s in xrange(0, rlen, 60):
        print(r.sequence[s:s+60])
    count += 1

