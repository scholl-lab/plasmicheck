# ---- Build stage ----
FROM python:3.12-slim AS build

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gcc wget \
    zlib1g-dev libbz2-dev liblzma-dev libncurses-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Build minimap2 from source (apt version lags behind)
ARG MINIMAP2_VERSION=2.28
RUN wget -qO- "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}.tar.bz2" \
    | tar -xj && cd minimap2-${MINIMAP2_VERSION} && make && cp minimap2 /usr/local/bin/

# Build samtools from source
ARG SAMTOOLS_VERSION=1.21
RUN wget -qO- "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
    | tar -xj && cd samtools-${SAMTOOLS_VERSION} && ./configure --without-curses && make && make install

# Install PlasmiCheck in a virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
COPY . /tmp/plasmicheck
RUN pip install --no-cache-dir /tmp/plasmicheck && rm -rf /tmp/plasmicheck

# ---- Runtime stage ----
FROM python:3.12-slim AS runtime

LABEL org.opencontainers.image.title="PlasmiCheck" \
      org.opencontainers.image.description="Plasmid DNA contamination detection in sequencing data" \
      org.opencontainers.image.source="https://github.com/scholl-lab/plasmicheck" \
      org.opencontainers.image.licenses="MIT"

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4 libbz2-1.0 liblzma5 zlib1g \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY --from=build /usr/local/bin/minimap2 /usr/local/bin/minimap2
COPY --from=build /usr/local/bin/samtools /usr/local/bin/samtools
COPY --from=build /opt/venv /opt/venv

ENV PATH="/opt/venv/bin:/usr/local/bin:$PATH"

RUN useradd --create-home appuser
USER appuser
WORKDIR /data

ENTRYPOINT ["plasmicheck"]
CMD ["--help"]

# ---- Test stage (validates build) ----
FROM runtime AS test
USER root
RUN plasmicheck --version && minimap2 --version && samtools --version
USER appuser
