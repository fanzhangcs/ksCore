# ksCore
The code is for computing the (k,s)-core which is published in the paper "Discovering Strong Communities with User Engagement and Tie Strength", Fan Zhang, Long Yuan, Ying Zhang, Lu Qin, Xuemin Lin, Alexander Zhou, DASFAA 2018


# files
khash.h - a hash table from https://github.com/attractivechaos/klib

ksCore.cpp - source code 

dataset.txt - toy friendship data with 5403 vertices and 20368 edges - data structure: vid \t nid \n... - note that each edge is stored twice and ordered here

the data file is a part of the Gowalla dataset from SNAP: https://snap.stanford.edu/data/


# compile and run
complie with g++ and -O3

run and input the values of 'k' and 's', such as 5 4 for k=5 and s=10

the program ouputs in result.txt

# note
If you have any question, please contact me by fanzhang.cs@gmail.com.

If you used this code, please kindly cite the paper.

