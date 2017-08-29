FROM debian:jessie-slim
MAINTAINER Chris Dean "cdean11@colostate.edu"

LABEL org.label-schema.vendor="Microbial Ecology Group" \
      org.label-schema.url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.name="SNPFinder" \
      org.label-description="A Docker image for SNPFinder"\
      org.label-schema.version="1.0.0" \
      org.label-schema.vcs-url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.build-date="2017-07-25T08:36:44Z" \
      org.label-schema.schema-version="1.0.0-rc1"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    wget \
 && rm -rf /var/lib/apt/lists/*

ENV SNPFINDER_SOURCE="v1.0.0.tar.gz" \
    SNPFINDER_VERSION="1.0.0" \
    SNPFINDER_BIN="snpfinder" \
    SNPFINDER_DEST="/usr/local/bin/snpfinder"

RUN wget https://github.com/cdeanj/snpfinder/archive/$SNPFINDER_SOURCE -O /opt/$SNPFINDER_SOURCE \
    && tar -xvf /opt/$SNPFINDER_SOURCE -C /opt \
    && cd /opt/snpfinder-$SNPFINDER_VERSION \
    && make \
    && ln -s /opt/snpfinder-$SNPFINDER_VERSION/$SNPFINDER_BIN $SNPFINDER_DEST \
    && rm /opt/$SNPFINDER_SOURCE
