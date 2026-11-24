FROM ubuntu:24.04

LABEL org.opencontainers.image.title="cfDNA Precision Diagnostics"
LABEL org.opencontainers.image.description="Reproducible cfDNA bioinformatics environment"
LABEL org.opencontainers.image.version="1.0.0"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    bwa \
    samtools \
    bcftools \
    bedtools \
    fastqc \
    tabix \
    default-jre-headless \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses-dev \
    ca-certificates \
    curl \
    git \
    procps \
    time \
    gzip \
    bash \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/cfdna

ENV PATH="/opt/cfdna/bin:${PATH}"
ENV PYTHONUNBUFFERED=1
ENV MPLBACKEND=Agg
ENV MPLCONFIGDIR=/tmp/matplotlib
ENV PIP_NO_CACHE_DIR=1

COPY containers/requirements.txt /tmp/requirements.txt

RUN python -m pip install --upgrade pip setuptools wheel \
    && python -m pip install -r /tmp/requirements.txt

WORKDIR /workspace

CMD ["/bin/bash"]
