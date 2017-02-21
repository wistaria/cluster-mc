#!/bin/sh

wget https://github.com/todo-group/standards/archive/develop.tar.gz
DIRS="config math stat"
for d in $DIRS; do
  tar zxf develop.tar.gz standards-develop/$d
  mkdir -p $d
  mv standards-develop/$d/* $d/
done
rm -rf develop.tar.gz standards-develop
