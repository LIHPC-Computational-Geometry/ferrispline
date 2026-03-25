FROM ubuntu:noble-20260113 AS base

FROM base AS rust

RUN apt update && \
    apt install -y --no-install-recommends \
    ca-certificates \
    curl && \
    rm -rf /var/lib/apt/lists/*

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

FROM base AS dep

RUN apt update && \
    apt install -y --no-install-recommends \
    build-essential \
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

FROM dep AS latest

COPY --from=rust --link /root/.cargo/ /root/.cargo/

COPY --from=rust --link /root/.rustup/ /root/.rustup/

ENV PATH="/root/.cargo/bin:${PATH}"

FROM dep AS dev

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

COPY --chown=${USER} --from=rust --link /root/.cargo/ .cargo/

COPY --chown=${USER} --from=rust --link /root/.rustup/ .rustup/

ENV PATH="/home/${USER}/.cargo/bin:${PATH}"
