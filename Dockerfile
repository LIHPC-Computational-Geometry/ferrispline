FROM ubuntu:noble-20260410 AS base

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
    libgl1 \
    pip \
    pre-commit \
    python3 \
    python3-dev   \
    python3-matplotlib \
    python3-pytest \
    python3-scipy && \
    rm -rf /var/lib/apt/lists/* && \
    pip install --break-system-packages \
    pyvista \
    maturin

FROM dep AS latest

COPY --from=rust --link /root/.cargo/ /root/.cargo/
COPY --from=rust --link /root/.rustup/ /root/.rustup/

ENV CARGO_HOME=/root/.cargo
ENV RUSTUP_HOME=/root/.rustup
ENV PATH="${CARGO_HOME}/bin:${PATH}"

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
