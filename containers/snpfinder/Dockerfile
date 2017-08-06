# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Christopher Dean" \
      description="SNPFinder(v1.0.0) Docker image" \
      maintainer="cdean11@colostate.edu"

# Install core libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    wget \
 && rm -rf /var/lib/apt/lists/*

# Setup environment variables
ENV SNPFINDER_SOURCE="v1.0.0.tar.gz" \
    SNPFINDER_VERSION="1.0.0" \
    SNPFINDER_BIN="snpfinder" \
    SNPFINDER_DEST="/usr/local/bin/snpfinder"

# Install SNPFinder
RUN wget https://github.com/cdeanj/snpfinder/archive/$SNPFINDER_SOURCE -O /opt/$SNPFINDER_SOURCE \
    && tar -xvf /opt/$SNPFINDER_SOURCE -C /opt \
    && cd /opt/snpfinder-$SNPFINDER_VERSION \
    && make \
    && ln -s /opt/snpfinder-$SNPFINDER_VERSION/$SNPFINDER_BIN $SNPFINDER_DEST \
    && rm /opt/$SNPFINDER_SOURCE
