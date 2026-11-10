FROM debian:trixie-slim AS builder
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    gfortran liblapack-dev libblas-dev build-essential curl ca-certificates && \
    curl -fsSL https://github.com/fortran-lang/fpm/releases/download/v0.10.1/fpm-0.10.1-linux-x86_64 -o /usr/local/bin/fpm && \
    chmod +x /usr/local/bin/fpm && \
    rm -rf /var/lib/apt/lists/* && \
    useradd -m -u 1000 -s /bin/bash appuser
WORKDIR /work
COPY --chown=1000:1000 . /work
USER 1000:1000
RUN /bin/sh scripts/test-offline.sh && /bin/sh scripts/demo-offline.sh --assert

FROM debian:trixie-slim AS runtime
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgfortran5 && \
    rm -rf /var/lib/apt/lists/* && \
    useradd -m -u 1000 -s /bin/bash appuser && \
    test ! -x /usr/bin/gfortran && test ! -f /usr/lib/x86_64-linux-gnu/liblapack.so.3
WORKDIR /work
COPY --from=builder /work/dist/scaled-qr-evidence /usr/local/bin/scaled-qr-evidence
USER 1000:1000
ENTRYPOINT ["/usr/local/bin/scaled-qr-evidence"]
CMD ["demo", "--assert"]
