#!/bin/bash

echo "input ..."
cat ../input/*.tab > collected/input.tab


echo "nextclade ..."
cat ../output/*/nextclade/*.tab > collected/nextclade.tab


echo "pangolin ..."
cat ../output/*/pangolin/*.csv > collected/pangolin.csv


