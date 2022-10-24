#!/bin/bash

#
# this is example for execution dockerized paean program
#
# prerequsites 
#     docker-ce
#     nvidia-docker2
#     nvidia-driver>=470 ( for support cuda11.3.1)
# 

DATAPATH=/SOME/DATA/FOLDER
FILENAME=TARGET_BAM_FILENAME.bam
THREAD=48

# nvidia-docker1
# nvidia-docker run -it --rm \
# ...
 
# nvidia-docker2
docker run --gpus all -it --rm \
 -v "${DATAPATH}":/data \
  -v `pwd`/output:/output bio-acc/paean:v1.1.0 \
    -b /usr/local/share/paean/input/gene.annotation.gff3 \
    -l /usr/local/share/paean/input/length_table.csv \
    -x SE,A3SS \
    -y /usr/local/share/paean/input/csv/SE.annotation.csv,/usr/local/share/paean/input/csv/A3SS.annotation.csv \
    -r /data/${FILENAME} \
    -o /output/ -t ${THREAD}  \
    -m 2

