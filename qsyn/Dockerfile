FROM ubuntu:20.04 
LABEL maintainer = "Kang, Chan Gu <changukang@korea.ac.kr>"


USER root
RUN chmod 777 /tmp
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        python3.8-dev \
        python3-pip \
        sudo \
        zsh && \
    apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/*

WORKDIR /workdir

COPY ./qsyn /workdir/qsyn
COPY ./benchmarks /workdir/benchmarks
COPY ./requirements.txt /workdir/
COPY ./transpile_qc.py /workdir/

RUN pip3 install -r requirements.txt