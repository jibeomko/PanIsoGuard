# PanIsoGuard container: the C++ adjudicator, built from source against conda htslib.
# Build:  docker build -t panisoguard .
# Use:    docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data -w /data panisoguard \
#             panisoguard adjudicate --classification cls.txt --isoforms-gtf iso.gtf --ref-gtf ref.gtf --out-prefix run
#         docker run --rm panisoguard panisoguard --version
# (-u runs as you so outputs can be written to the mounted folder; the image user cannot.)
# (Pass the command through the default entrypoint — do NOT use --entrypoint, which would
#  bypass the conda-env activation that puts htslib on the path.)
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
        && cmake --install build --prefix /opt/pig \
        && ./build/panisoguard --version'

# --- stage 2: slim runtime (no compiler) -------------------------------------
FROM mambaorg/micromamba:1.5-jammy AS runtime
LABEL org.opencontainers.image.title="PanIsoGuard" \
      org.opencontainers.image.description="Caller-agnostic adjudication of long-read novel isoforms" \
      org.opencontainers.image.source="https://github.com/jibeomko/PanIsoGuard" \
      org.opencontainers.image.licenses="MIT"
USER root
RUN micromamba install -y -n base -c conda-forge -c bioconda 'htslib>=1.18' zlib \
    && micromamba clean -ay
COPY --from=build /opt/pig/bin/panisoguard /opt/conda/bin/panisoguard
COPY --from=build /opt/pig/share/panisoguard /opt/conda/share/panisoguard
RUN micromamba run -n base panisoguard --version
USER $MAMBA_USER
# micromamba's entrypoint activates the base env (so htslib is on the path).
CMD ["panisoguard", "--help"]
