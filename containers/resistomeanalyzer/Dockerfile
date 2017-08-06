# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Christopher Dean" \
      description="ResistomeAnalyzer(v1.0.0) Docker image" \
      maintainer="cdean11@colostate.edu"

# Install core libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    git \
 && rm -rf /var/lib/apt/lists/*

# Install ResistomeAnalyzer
RUN cd /opt \
    && git clone https://github.com/cdeanj/resistomeanalyzer.git \
    && cd resistomeanalyzer \
    && make \
    && ln -s /opt/resistomeanalyzer/resistome /usr/local/bin/resistome
