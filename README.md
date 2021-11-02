
# ENCE464 Assignment 2: Optimisation and Multithreading


**Objective:** Implement Jacobi Relaxation to solve Poisson's equation, while making the best use of cores and cashes. See assignment instructions [here](doc/instructions/instructions.pdf).

## Contents
 - `doc/` - assignment instructions, lab notes, report template.
 - `reference/` - correct output for test comparison.
 - `poisson_basic.c` - basic poisson solver. Non-optimised.
 - `poisson_threads.c` - Multithreaded poisson solver.
 - `threads.c` - example on how to use POSIX thread library.
 - `test.sh` - automatic testing script.
 - `makefile` - makefile to run program.

## Building

`Make` - builds **poisson_threads.c** into a binary **poisson_threads** file.
`Make run` - runs **poisson_threads.c**
`Make test` - tests **poisson_threads.c** against the referance files to check if output is correct.
`g++ -o poisson_threads poisson_threads.c -lpthread && ./poisson_threads`
Other make options can be found in Makefile.