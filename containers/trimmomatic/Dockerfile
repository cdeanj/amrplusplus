FROM openjdk:alpine
MAINTAINER Chris Dean "cdean11@colostate.edu"

LABEL org.label-schema.vendor="Microbial Ecology Group" \
      org.label-schema.url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.name="Trimmomatic" \
      org.label-description="A Docker image for Trimmomatic"\
      org.label-schema.version="1.0.0" \
      org.label-schema.vcs-url="https://github.com/colostatemeg/amrplusplus" \
      org.label-schema.build-date="2017-07-25T08:36:44Z" \
      org.label-schema.schema-version="1.0.0-rc1"

RUN apk add --update bash && rm -rf /var/cache/apk/*

ENV TRIMMOMATIC_SOURCE="Trimmomatic-0.36.zip" \
    TRIMMOMATIC_HOME="/opt/trimmomatic"

COPY $TRIMMOMATIC_SOURCE $TRIMMOMATIC_HOME/$TRIMMOMATIC_SOURCE

RUN cd $TRIMMOMATIC_HOME \
    && unzip $TRIMMOMATIC_SOURCE \
    && rm $TRIMMOMATIC_SOURCE
