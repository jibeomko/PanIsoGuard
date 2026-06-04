# PanIsoGuard container: the C++ adjudicator (built from source) + panisoguard-report.
# Build:  docker build -t panisoguard .
# Use:    docker run --rm -v "$PWD":/data -w /data panisoguard \
#             adjudicate --classification cls.txt --isoforms-gtf iso.gtf --ref-gtf ref.gtf --out-prefix run
#         docker run --rm panisoguard panisoguard --version
#         docker run --rm -v "$PWD":/data -w /data panisoguard panisoguard-report --prefix run
# (Pass the command through the default entrypoint — do NOT use --entrypoint, which would
#  bypass the conda-env activation that puts htslib/matplotlib on the path.)
#
# --- stage 1: build the C++ binary against conda htslib ----------------------
FROM mambaorg/micromamba:1.5-jammy AS build
USER root
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        cxx-compiler c-compiler cmake make pkg-config 'htslib>=1.18' zlib \
    && micromamba clean -ay
COPY --chown=root:root . /src
WORKDIR /src
RUN micromamba run -n base bash -lc '\
        cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
              -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
        && cmake --build build -j"$(nproc)" \
        && ./build/panisoguard --version'

# --- stage 2: slim runtime (no compiler) -------------------------------------
FROM mambaorg/micromamba:1.5-jammy AS runtime
LABEL org.opencontainers.image.title="PanIsoGuard" \
      org.opencontainers.image.description="Caller-agnostic adjudication of long-read novel isoforms + PDF report" \
      org.opencontainers.image.source="https://github.com/jibeomko/PanIsoGuard" \
      org.opencontainers.image.licenses="MIT"
USER root
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        'htslib>=1.18' zlib 'python>=3.9' 'matplotlib>=3.4' \
    && micromamba clean -ay
COPY --from=build /src/build/panisoguard /opt/conda/bin/panisoguard
COPY --from=build /src/python /opt/panisoguard-python
RUN micromamba run -n base pip install --no-cache-dir /opt/panisoguard-python \
    && micromamba run -n base panisoguard --version \
    && micromamba run -n base panisoguard-report --help >/dev/null
USER $MAMBA_USER
# micromamba's entrypoint activates the base env (so htslib/matplotlib are on the path).
CMD ["panisoguard", "--help"]
