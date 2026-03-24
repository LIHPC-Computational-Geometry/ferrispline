FROM ubuntu:noble-20260113 AS base

FROM base AS latest

RUN apt update && \
    apt install -y --no-install-recommends \
    git \
    pip \
    pre-commit \
    python3 \
    python3-matplotlib \
    python3-pytest \
    python3-scipy \
    libgl1 && \
    rm -rf /var/lib/apt/lists/* && \
    pip install pyvista --break-system-packages

FROM latest AS dev

RUN apt update && \
    apt install -y --no-install-recommends \
    ca-certificates \
    less \
    ssh \
    vim && \
    rm -rf /var/lib/apt/lists/*

ENV DISPLAY=:0

ARG USER

RUN useradd -m -s /bin/bash ${USER}

USER ${USER}

WORKDIR /home/${USER}/
