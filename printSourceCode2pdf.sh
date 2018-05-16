#!/bin/sh
/usr/bin/enscript -f Courier8 --highlight --line-numbers -o - `ls -1 *.h *.cpp` | ps2pdf - CombEng_source_book.pdf
