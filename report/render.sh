#!/bin/bash
for filename in *.gv
do
    name=$(echo "$filename" | cut -f 1 -d '.')
    dot -Tpng $name.gv -o $name.png #provided by graphviz
done
