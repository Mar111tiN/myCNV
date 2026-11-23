# syntax=docker/dockerfile:1
# mycnv: PON-CNV steps of kumo-variant-calling-nf (plan 3b; named stages, runtime last, plan D9/D11).
# Build from the repo root:
#   docker build --target test .                 # runs the public test suite, not shipped
#   docker build -t mycnv:2.0.0.dev0 .           # the runtime image (last stage)
# Use: every step is a `mycnv <subcommand>` call, e.g. `mycnv sample --pileup ... --ponmatrix ...`.
# Only pyproject.toml, README.md, src/ (build) and tests/ (test) enter the image; .dockerignore keeps the git-ignored,
# patient-derived tests/data/restricted/ out of the build context.

# ---- build: solve the env once, at its final prefix, and install the package ----
FROM mambaorg/micromamba:2.10.0-debian13-slim AS build
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean -a -y
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md /tmp/mycnv/
COPY --chown=$MAMBA_USER:$MAMBA_USER src /tmp/mycnv/src
# pip is needed only for this install. It is removed afterwards: its vendored libraries (pip/_vendor/vendor.txt:
# msgpack, urllib3, setuptools/pkg_resources) carried 5 fixed HIGH CVEs in the first CI Trivy gate (2026-11-23), and
# nothing calls pip at runtime. The lock file is written after the removal, so it lists what ships.
RUN /opt/conda/bin/python -m pip install --no-cache-dir --no-deps --no-build-isolation /tmp/mycnv \
 && /opt/conda/bin/mycnv --version \
 && micromamba remove -y -n base --force pip \
 && micromamba env export -n base > /opt/conda/env.lock.yml \
 && micromamba clean -a -y \
 && /opt/conda/bin/mycnv --version

# ---- test: docker build --target test . ; not part of the shipped image ----
FROM build AS test
COPY --chown=$MAMBA_USER:$MAMBA_USER tests /tmp/mycnv/tests
WORKDIR /tmp/mycnv
RUN /opt/conda/bin/python -m pytest -q -p no:cacheprovider tests

# ---- runtime: the shipped image (must stay the last stage) ----
FROM debian:13.7-slim AS runtime
RUN apt-get update \
 && apt-get install -y --no-install-recommends bash procps gzip \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/conda /opt/conda
# MPLCONFIGDIR: Nextflow runs tasks as `-u uid:gid` with a HOME that is not writable; matplotlib then warns and
# rebuilds its cache in a random temp dir on every plot task
ENV PATH=/opt/conda/bin:$PATH \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/mycnv-matplotlib
### Create image user
RUN useradd -m mycnv_user
USER mycnv_user
WORKDIR /home/mycnv_user
