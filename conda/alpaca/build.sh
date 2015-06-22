#!/bin/bash

cargo build --release
mkdir -p $PREFIX/bin
cp target/release/alpaca $PREFIX/bin
