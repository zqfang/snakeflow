#!bin/sh
set -e
# if mac os
md5 $1 >> md5.txt
# if linux
md5sum $1 >> md5.txt
