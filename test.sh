#!/usr/bin/sh

status=0

if ./poisson -n 7 -i 300 -t 2 --test  | cmp reference/7.txt; then
    echo "n=7 i=300 t=2 correct"
else
    echo "n=7 i=300 t=2 failed!"
    status=1
fi

if ./poisson -n 15 -i 300 -t 2 --test  | cmp reference/15.txt; then
    echo "n=15 i=300 t=2 correct"
else
    echo "n=15 i=300 t=2 failed!"
    status=1
fi

if ./poisson -n 51 -i 300 -t 1 --test  | cmp reference/51.txt; then
    echo "n=51 i=300 t=1 correct"
else
    echo "n=51 i=300 t=1 failed!"
    status=1
fi

if ./poisson -n 51 -i 300 -t 2 --test  | cmp reference/51.txt; then
    echo "n=51 i=300 t=2 correct"
else
    echo "n=51 i=300 t=2 failed!"
    status=1
fi

if ./poisson -n 51 -i 300 -t 4 --test  | cmp reference/51.txt; then
    echo "n=51 i=300 t=4 correct"
else
    echo "n=51 i=300 t=4 failed!"
    status=1
fi

if ./poisson -n 51 -i 300 -t 8 --test  | cmp reference/51.txt; then
    echo "n=51 i=300 t=8 correct"
else
    echo "n=51 i=300 t=8 failed!"
    status=1
fi

exit $status
