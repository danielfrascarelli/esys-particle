#!/bin/sh

((XDIM=$1*40))
((YDIM=$2*40))
((GOUGE=$2*12))
((ROUGH=$2*5))

gengeo -x $XDIM -cx -y $YDIM -gouge $GOUGE -gougeMinRadius 0.2 -fault $ROUGH -faultMinRadius 0.2 -minDrivingTag 3 -maxDrivingTag 4 -distanceFromBBoxEdge 1.5 -o bench_gouge.geo -t 10000 -rs 42
