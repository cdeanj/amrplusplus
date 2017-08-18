# Set the base image
FROM linuxbrew/linuxbrew

# Set image metadata
LABEL author="Christopher Dean" \
      description="Kraken (v0.10.5) Docker image" \
      maintainer="cdean11@colostate.edu"

# Install Kraken
RUN brew tap homebrew/science \
    && brew install kraken

# Set up environment variables
ENV MINIKRAKEN_SOURCE="minikraken.tgz" \
    MINIKRAKEN_DEST="minikraken" \
    MINIKRAKEN_DIR="minikraken_20141208"

# Install MiniKraken database
RUN cd /opt \
    && curl https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz -o $MINIKRAKEN_SOURCE \
    && tar -xvf $MINIKRAKEN_SOURCE \
    && mv $MINIKRAKEN_DIR $MINIKRAKEN_DEST \
    && rm $MINIKRAKEN_SOURCE
