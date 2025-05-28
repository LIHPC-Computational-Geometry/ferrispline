FROM ubuntu:noble-20250404 AS base

RUN apt update

FROM base AS python

RUN apt install -y --no-install-recommends \
    python3 python3-matplotlib python3-pytest python3-scipy

FROM python AS devcontainer

RUN apt install -y --no-install-recommends \
    ca-certificates \
    git \
    ssh \
    vim

ENV DISPLAY=:0

ARG USER

RUN useradd -m -s /bin/bash $USER

USER $USER
