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

ENV PATH="/root/.cargo/bin:${PATH}"

FROM latest AS build

ADD core_rust/ /opt/core_rust/

ADD nurbslib/ /opt/nurbslib/

WORKDIR /opt/nurbslib/

RUN maturin build --release --out dist && \
    pip install --break-system-packages dist/*.whl --force-reinstall

FROM dep AS test

COPY --from=build --link /opt/nurbslib/ /opt/nurbslib/

ADD sandbox_python/ /opt/sandbox_python/

WORKDIR /opt/sandbox_python/

RUN pip install --break-system-packages -e . && nurbs-convert -f /opt/sandbox_python/vtk/curve_test.vtk

FROM base AS bot

COPY --from=build --link /opt/nurbslib/ /opt/nurbslib/

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
