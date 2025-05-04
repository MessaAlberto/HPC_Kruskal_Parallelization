#!/bin/bash

# Test graph of 1024 to 65536 nodes
v_start=1024
v_end=65536

v=$v_start
while [ $v -le $v_end ]; do
  SCRIPT_NAME="script/qsub_dir.sh"

  if [ ! -f $SCRIPT_NAME ]; then
    echo "Missing file: $SCRIPT_NAME"
    exit 1
  fi

  echo "Send job: Kruskal on graph of $v nodes"
  qsub -v V="$v" -o "./output/pivotSort_${v}_nodes.out" -e "./output_err/pivotSort_${v}_nodes.err" $SCRIPT_NAME
  # qsub -v V="$v" -o "./output/mergeSort_${v}_nodes.out" -e "./output_err/mergeSort_${v}_nodes.err" $SCRIPT_NAME

  # v is doubled for the next iteration -> it's always a power of 2
  v=$((v*2))
done