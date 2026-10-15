# ---- Build stage ----
FROM python:3.12-slim AS build

COPY --from=ghcr.io/astral-sh/uv:0.11.28@sha256:0f36cb9361a3346885ca3677e3767016687b5a170c1a6b88465ec14aefec90aa /uv /usr/local/bin/uv

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

# Install the locked runtime dependencies and a non-editable PlasmiCheck package.
ENV UV_PROJECT_ENVIRONMENT=/opt/venv UV_PYTHON_DOWNLOADS=never
COPY . /tmp/plasmicheck
WORKDIR /tmp/plasmicheck
RUN uv sync --locked --no-dev --no-editable --no-cache && rm -rf /tmp/plasmicheck

# ---- Runtime stage ----
FROM python:3.12-slim AS runtime-base

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
FROM runtime-base AS test
RUN plasmicheck --version && minimap2 --version && samtools --version \
    && python -c 'from plasmicheck.resources import get_resource_path; assert all(get_resource_path(p).is_file() for p in ("config.json", "py.typed", "templates/report_template.html", "templates/summary_template.html", "static/img/plasmicheck_logo_200px.png"))'

# Keep a plain docker build equivalent to --target runtime.
FROM runtime-base AS runtime
