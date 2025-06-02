FROM ubuntu:noble-20250404

RUN apt update

RUN apt install -y --no-install-recommends \
    ca-certificates \
    git \
    python3 python3-matplotlib python3-scipy \
    ssh \
    vim

RUN groupadd --gid 1001 vscode

RUN useradd --uid 1001 --gid 1001 -m -s /bin/bash vscode

USER vscode
