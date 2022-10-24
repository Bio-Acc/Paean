#!/bin/bash

#
# this is example for building dockerized paean program
#
# prerequsites 
#     docker-ce
# 
# 1a Dockerfile.paean-ubuntu2004.dockerfile    
#
# 1b Dockerfile.paean-ubuntu2004-chn.dockerfile   (force mirror redirection)
#   * manually download src-pkg before building
#     https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2
#     https://github.com/Bio-Acc/Paean/archive/refs/tags/v1.1.0.tar.gz
#

# build 1a
docker build --network=host -f Dockerfile.paean-ubuntu2004.dockerfile -t bio-acc/paean:v1.1.0 .

# # build 1b
# docker build --network=host -f Dockerfile.paean-ubuntu2004-chn.dockerfile -t bio-acc/paean:v1.1.0 .

# test build result
docker run -it --rm bio-acc/paean:v1.1.0

