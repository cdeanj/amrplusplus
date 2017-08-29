FROM debian:jessie-slim
MAINTAINER Chris Dean "cdean11@colostate.edu"

LABEL org.label-schema.vendor="Microbial Ecology Group" \
      org.label-schema.url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.name="MultiQC" \
      org.label-description="A Docker image for MultiQC"\
      org.label-schema.version="1.0.0" \
      org.label-schema.vcs-url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.build-date="2017-07-25T08:36:44Z" \
      org.label-schema.schema-version="1.0.0-rc1"

# Thanks Maxime Garcia
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl  \
    gcc \
    python-dev \
    && curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /opt/get-pip.py \
    && python /opt/get-pip.py \
    && rm -rf /var/lib/apt/lists/* /opt/get-pip.py

ENV MULTIQC_VERSION="1.2"

RUN pip install multiqc==$MULTIQC_VERSION
