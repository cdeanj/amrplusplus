FROM debian:jessie-slim
MAINTAINER Chris Dean "cdean11@colostate.edu"

LABEL org.label-schema.vendor="Microbial Ecology Group" \
      org.label-schema.url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.name="RarefactionAnalyzer" \
      org.label-description="A Docker image for RarefactionAnalyzer"\
      org.label-schema.version="1.0.0" \
      org.label-schema.vcs-url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.build-date="2017-07-25T08:36:44Z" \
      org.label-schema.schema-version="1.0.0-rc1"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    git \
 && rm -rf /var/lib/apt/lists/*

RUN cd /opt \
    && git clone https://github.com/cdeanj/rarefactionanalyzer.git \
    && cd rarefactionanalyzer \
    && make \
    && ln -s /opt/rarefactionanalyzer/rarefaction /usr/local/bin/rarefaction
