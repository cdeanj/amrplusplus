FROM debian:jessie-slim
MAINTAINER Chris Dean "cdean11@colostate.edu"

LABEL org.label-schema.vendor="Microbial Ecology Group" \
      org.label-schema.url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.name="Burrows Wheeler Aligner (BWA)" \
      org.label-description="A Docker image for BWA"\
      org.label-schema.version="1.0.0" \
      org.label-schema.vcs-url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.build-date="2017-07-25T08:36:44Z" \
      org.label-schema.schema-version="1.0.0-rc1"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    wget \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

ENV BWA_SOURCE="bwa-0.7.15.tar.bz2" \
    BWA_VERSION="0.7.15" \
    BWA_BIN="bwa" \
    BWA_DEST="/usr/local/bin/bwa"

RUN wget https://github.com/lh3/bwa/releases/download/v$BWA_VERSION/$BWA_SOURCE -O /opt/$BWA_SOURCE \
    && tar -xvf /opt/$BWA_SOURCE -C /opt \
    && cd /opt/bwa-$BWA_VERSION \
    && make \
    && ln -s /opt/bwa-$BWA_VERSION/$BWA_BIN $BWA_DEST \
    && rm /opt/$BWA_SOURCE
