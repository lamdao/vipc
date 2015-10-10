#!/bin/sh
for f in *.def; do
	x86_64-w64-mingw32-dlltool -d $f -l ${f/.def/.a};
done
