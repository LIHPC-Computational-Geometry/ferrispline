FROM ubuntu:noble-20250404 AS base

RUN apt update

FROM base AS gitlabci

RUN apt install -y --no-install-recommends \
    git \
    pre-commit \
    python3 python3-matplotlib python3-pytest python3-scipy

FROM gitlabci AS devcontainer

RUN apt install -y --no-install-recommends \
    ca-certificates \
    less \
    ssh \
    vim

RUN groupadd --gid 1001 vscode

RUN useradd --uid 1001 --gid 1001 -m -s /bin/bash vscode

USER vscode
