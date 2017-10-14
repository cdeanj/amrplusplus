FROM debian:jessie-slim
MAINTAINER Chris Dean "cdean11@colostate.edu"

LABEL org.label-schema.vendor="Microbial Ecology Group" \
      org.label-schema.url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.name="Kraken" \
      org.label-description="A Docker image for Kraken"\
      org.label-schema.version="1.0.0" \
      org.label-schema.vcs-url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.build-date="2017-07-25T08:36:44Z" \
      org.label-schema.schema-version="1.0.0-rc1"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    git \ 
    wget \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

ENV KRAKEN_DIR="/opt/kraken"

ENV MINIKRAKEN_SOURCE="minikraken.tgz" \
    MINIKRAKEN_DEST="minikraken" \
    MINIKRAKEN_DIR="minikraken_20141208"

RUN cd /opt \
    && git clone https://github.com/DerrickWood/kraken.git \
    && cd kraken \
    && ./install_kraken.sh $KRAKEN_DIR \
    && cp $KRAKEN_DIR/kraken /usr/local/bin \
    && cp $KRAKEN_DIR/kraken-build /usr/local/bin \
    && cp $KRAKEN_DIR/kraken-report /usr/local/bin

RUN cd /opt \
    && curl https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz -o $MINIKRAKEN_SOURCE \
    && tar -xvf $MINIKRAKEN_SOURCE \
    && mv $MINIKRAKEN_DIR $MINIKRAKEN_DEST \
    && rm $MINIKRAKEN_SOURCE
