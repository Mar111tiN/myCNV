# syntax=docker/dockerfile:1
# mycnv: PON-CNV steps of kumo-variant-calling-nf (plan 3b; named stages, runtime last, plan D9/D11).
# Build from the repo root:
#   docker build --target test .                 # runs the test suite, not shipped
#   docker build -t mycnv:2.0.0.dev0 .           # the runtime image (last stage)
# Use: every step is a `mycnv <subcommand>` call, e.g. `mycnv sample --pileup ... --pon-dir ...`.

# ---- build: solve the env once, at its final prefix, and install the package ----
FROM mambaorg/micromamba:2.10.0-debian13-slim AS build
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean -a -y
COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp/mycnv
RUN /opt/conda/bin/python -m pip install --no-cache-dir --no-deps --no-build-isolation /tmp/mycnv \
 && micromamba env export -n base > /opt/conda/env.lock.yml

# ---- test: docker build --target test . ; not part of the shipped image ----
FROM build AS test
WORKDIR /tmp/mycnv
RUN /opt/conda/bin/python -m pytest -q -p no:cacheprovider tests

# ---- runtime: the shipped image (must stay the last stage) ----
FROM debian:13-slim AS runtime
RUN apt-get update \
 && apt-get install -y --no-install-recommends bash procps gzip \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/conda /opt/conda
ENV PATH=/opt/conda/bin:$PATH \
    MPLBACKEND=Agg
### Create image user
RUN useradd -m mycnv_user
USER mycnv_user
WORKDIR /home/mycnv_user
