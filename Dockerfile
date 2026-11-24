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

# ---- slim: the tested env minus what nothing runs (plan 3b, 2026-11-24; about half the size) ----
# docker/slim_env.py: force-removes the listed packages that are installed (pytest and its helpers; tk; cairo with its
# X11 libs and the conda fonts: only harfbuzz's optional cairo build uses them, while matplotlib's raqm/harfbuzz/glib
# stay and matplotlib renders with its bundled DejaVu fonts), deletes test suites, headers, static libs, docs,
# translations and __pycache__, and strips the ELF files (--strip-unneeded keeps the dynamic symbols). Checked on a
# copy of the env: full test suite (incl. plots) and tier-A smoke unchanged.
FROM build AS slim
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends binutils \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/slim_env.py /tmp/slim_env.py
RUN /opt/conda/bin/python -I /tmp/slim_env.py /opt/conda --remove pytest pygments pluggy colorama exceptiongroup \
        typing_extensions tomli iniconfig tk pip cairo pixman fontconfig xorg-libx11 xorg-libxext xorg-libxrender \
        xorg-libsm xorg-libice font-ttf-ubuntu font-ttf-source-code-pro font-ttf-dejavu-sans-mono font-ttf-inconsolata \
        fonts-conda-forge fonts-conda-ecosystem \
 && micromamba env export -n base > /opt/conda/env.lock.yml \
 && micromamba clean -a -y \
 && /opt/conda/bin/mycnv --version

# ---- runtime: the shipped image (must stay the last stage) ----
FROM debian:13.7-slim AS runtime
RUN apt-get update \
 && apt-get install -y --no-install-recommends bash procps gzip \
 && rm -rf /var/lib/apt/lists/*
COPY --from=slim /opt/conda /opt/conda
# MPLCONFIGDIR: Nextflow runs tasks as `-u uid:gid` with a HOME that is not writable; matplotlib then warns and
# rebuilds its cache in a random temp dir on every plot task. PYTHONDONTWRITEBYTECODE: the slim stage removed
# __pycache__ and /opt/conda is not writable for the task user, so Python must not try to write .pyc files.
ENV PATH=/opt/conda/bin:$PATH \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/mycnv-matplotlib \
    PYTHONDONTWRITEBYTECODE=1
### Create image user
RUN useradd -m mycnv_user
USER mycnv_user
WORKDIR /home/mycnv_user
