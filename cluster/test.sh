#!/bin/bash

SCRIPT_DIR="$(dirname "$(realpath "$0")")"

# Test graph of 1024 to 32768 nodes
v_start=1024
v_end=32768

v=$v_start
while [ $v -le $v_end ]; do
  SCRIPT_NAME="$SCRIPT_DIR/qsub_dir.sh"

  if [ ! -f $SCRIPT_NAME ]; then
    echo "Missing file: $SCRIPT_NAME"
    exit 1
  fi

  echo "Send job: Kruskal on graph of $v nodes"
  echo "qsub -v V="$v" -o "$SCRIPT_DIR/../output/test_${v}_nodes.out" -e "$SCRIPT_DIR/../output_err/test_${v}_nodes.err" $SCRIPT_NAME"

  v=$((v*2))
done