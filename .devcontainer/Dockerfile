FROM ubuntu:noble-20250404

RUN apt update

RUN apt install -y --no-install-recommends \
    git \
    python3 python3-matplotlib python3-scipy \
    ssh

ENV DISPLAY=:0

RUN useradd -m -s /bin/bash stauffertm

USER stauffertm
