cFDAP
=====

 cFDAP is a command-line Linux tool for fitting FRAP/FDAP data.

 The overall goal of cFDAP is to provide a fast and robust way to extract
kinetic parameters from FRAP/FDAP data gathered from living neuronal cells.
However, the code can be relatively simply adapted to any experiment outline,
which, in turn, requires theoretical expressions for FRAP(t)/FDAP(t) signals be
hardcoded explicitly (for more details please Igaev et al. (2015) Biophys. J.
107: 2567-2578).

Compilation
===========

 ```
 cc cFDAP.c -o cFDAP -lgsl -lgslcblas -lm
 ```

Usage
=====

 To get quick start just run
 ```
 ./cFDAP
 ```

Requirements
============

 To successfully compile cFDAP the GNU Scientific Library
(http://www.gnu.org/software/gsl/) must be installed.

Coping
======

 cFDAP is licensed under a GPL license. That means that distributing the code is only
allowed under the same terms. 
