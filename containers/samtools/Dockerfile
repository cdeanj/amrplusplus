FROM debian:jessie-slim
MAINTAINER Chris Dean "cdean11@colostate.edu"

LABEL org.label-schema.vendor="Microbial Ecology Group" \
      org.label-schema.url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.name="Samtools" \
      org.label-description="A Docker image for Samtools"\
      org.label-schema.version="1.0.0" \
      org.label-schema.vcs-url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.build-date="2017-07-25T08:36:44Z" \
      org.label-schema.schema-version="1.0.0-rc1"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    libbz2-dev \
    liblzma-dev \
    libncurses-dev \
    wget \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

ENV SAMTOOLS_SOURCE="samtools-1.4.1.tar.bz2" \
    SAMTOOLS_VERSION="1.4.1" \
    SAMTOOLS_BIN="samtools" \
    SAMTOOLS_DEST="/usr/local/bin/samtools"

RUN wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/$SAMTOOLS_SOURCE -O /opt/$SAMTOOLS_SOURCE \
    && tar -xvf /opt/$SAMTOOLS_SOURCE -C /opt \
    && cd /opt/samtools-$SAMTOOLS_VERSION \
    && ./configure \
    && make \
    && ln -s /opt/samtools-$SAMTOOLS_VERSION/$SAMTOOLS_BIN $SAMTOOLS_DEST \
    && rm /opt/$SAMTOOLS_SOURCE
