FROM nvidia/cuda:11.3.1-cudnn8-devel-ubuntu20.04 as builder

# as multi-stage builder, change to fast mirror
RUN     sed -i 's/archive.ubuntu.com/opentuna.cn/g' /etc/apt/sources.list && \
        sed -i "s+security.ubuntu.com+opentuna.cn+g" /etc/apt/sources.list

RUN     apt-get update && \
        DEBIAN_FRONTEND=nointeractive apt-get -y install build-essential zlib1g-dev libbz2-dev liblzma-dev cmake \
        autoconf automake make gcc perl libcurl4-gnutls-dev libssl-dev git && \
        rm -rf /var/lib/apt/lists/*

# manually download first!
# https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2
ADD     ./htslib-1.14.tar.bz2 /opt/src/
RUN     mkdir /opt/src -p && cd /opt/src && \
        cd htslib-1.14 && \
        autoheader && \
        autoconf && \
        ./configure && \
        make -j && \
        make install 

# manually download first!
# https://github.com/Bio-Acc/Paean/archive/refs/tags/v1.1.0.tar.gz
ADD     ./v1.1.0.tar.gz /opt/src/
RUN     mkdir /opt/src -p && cd /opt/src && \
        cd Paean-1.1.0 && \
        mkdir build && cd build && \
        cmake .. && \
        make


FROM nvidia/cuda:11.3.1-base-ubuntu20.04
LABEL "com.opencontainer.image.authors"="yaoruijie@picb.ac.cn"
LABEL version="1.1.0"
LABEL description="Paean: A unified and efficient transcriptome quantification system using heterogeneous computing"
LABEL "maintainer"="Ruijie YAO <yaoruijie@picb.ac.cn>"

COPY    --from=builder /opt/src/Paean-1.1.0/build/paean /usr/local/bin/paean
COPY    --from=builder /opt/src/Paean-1.1.0/input /usr/local/share/paean/input
COPY    --from=builder /usr/local/lib/libhts* /usr/local/lib/

# change mirror for fast mirror
RUN     mv /etc/apt/sources.list /etc/apt/sources.list.orig && cp /etc/apt/sources.list.orig /etc/apt/sources.list && \
        sed -i 's/archive.ubuntu.com/opentuna.cn/g' /etc/apt/sources.list && \
        sed -i "s+security.ubuntu.com+opentuna.cn+g" /etc/apt/sources.list && \
        apt-get update && \
        DEBIAN_FRONTEND=nointeractive apt-get -y install zlib1g libcurl4-gnutls-dev && \
        rm -rf /var/lib/apt/lists/* && rm /etc/apt/sources.list && mv /etc/apt/sources.list.orig /etc/apt/sources.list

RUN     mkdir -p /opt/paean

WORKDIR /opt/paean

ENTRYPOINT [ "/usr/local/bin/paean" ]
CMD     [ "--help" ]



