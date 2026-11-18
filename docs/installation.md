# Installation

[Back to the README quick start](../README.md#quick-start).

Install the Bioconda package, use a container, or build from source. The shell commands
below are intended for Linux or macOS; check the Bioconda
[files page](https://anaconda.org/bioconda/panisoguard/files) for the platforms available
for each package version.

## Bioconda

With conda available and your intended environment active:

```bash
conda install -c conda-forge -c bioconda panisoguard
panisoguard example tiny
bash tiny/run.sh
```

Run the example commands from a writable working directory that does not already contain
a `tiny` folder. The package is published on the
[bioconda channel](https://anaconda.org/bioconda/panisoguard); its recipe is in
[`recipes/bioconda/`](../recipes/bioconda/).

To install into a separate environment instead, replace the installation command with:

```bash
conda create -n panisoguard -c conda-forge -c bioconda panisoguard
conda activate panisoguard
```

Then run the example commands above. Choose another environment name if `panisoguard`
already exists.

## Bundled examples

The examples ship with the package in version 0.0.6 and later. Run `panisoguard example`
to list them:

- `tiny`: one caller, with a known isoform, a short-read-supported novel isoform, and an artifact.
- `multi_caller`: three callers, from `combine` to a consensus verdict.
- `worked_locus`: one toy gene followed through the rules; start with its `README.md`.

`panisoguard example <name> [DIR]` copies an example into a writable folder. Replace
`<name>` with one of the names above; the optional `DIR` is the destination path. If you
omit `DIR`, the destination is `./<name>`. The destination must not already exist.
For example, from your working directory:

```bash
panisoguard example multi_caller
bash multi_caller/run.sh

panisoguard example tiny tiny_copy
bash tiny_copy/run.sh
```

These two examples use bundled data and run offline after installation. They write results
under their own `output/` folders and compare them with the saved `expected/` files. A
successful tiny run ends with `example output matches expected files`; a successful
multi-caller run ends with `output matches expected`.

For a step-by-step walk-through from inputs to verdicts, see the
[worked locus](../examples/worked_locus/README.md).

## Containers

### Pull a released BioContainers image

Bioconda publishes a ready-made image on BioContainers. This command pulls the 0.1.0
image; available tags are listed on
[quay.io](https://quay.io/repository/biocontainers/panisoguard?tab=tags).

```bash
docker pull quay.io/biocontainers/panisoguard:0.1.0--h8471819_0
```

No source checkout is needed. From a writable working directory without a `tiny` folder,
copy and run the bundled example:

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/work -w /work \
    quay.io/biocontainers/panisoguard:0.1.0--h8471819_0 \
    bash -c "panisoguard example tiny && bash tiny/run.sh"
```

### Build an image from source

To build your own image from the repository's [Dockerfile](../Dockerfile), first obtain
the source checkout, then run `docker build` from its root:

```bash
git clone https://github.com/jibeomko/PanIsoGuard.git
cd PanIsoGuard
docker build -t panisoguard .
```

If you already have a checkout, change to that `PanIsoGuard` directory instead of cloning
again. The local image is named `panisoguard`; it contains the code from your checkout.
From a writable directory without `tiny` or `multi_caller` folders, run:

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/work -w /work panisoguard \
    bash -c "panisoguard example tiny && bash tiny/run.sh"
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/work -w /work panisoguard \
    bash -c "panisoguard example multi_caller && bash multi_caller/run.sh"
```

### Run on your own data

Change to the directory containing your input files, then mount it into the container.
The following command uses the locally built `panisoguard` image; to use the released
image, replace that image name with
`quay.io/biocontainers/panisoguard:0.1.0--h8471819_0`.

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data -w /data panisoguard \
    panisoguard adjudicate --classification cls.txt --isoforms-gtf iso.gtf --ref-gtf ref.gtf --out-prefix run
```

`cls.txt`, `iso.gtf`, and `ref.gtf` are placeholders for your SQANTI3 classification,
caller isoforms, and reference annotation. Paths passed to PanIsoGuard must be accessible
inside the mounted directory. `run` is the output prefix; results are written into the
mounted directory. `-u` runs the container with your user and group IDs so output files
are written as your user. Include `panisoguard` before its subcommand, as shown above.

## Build from source

### Use compilers and dependencies from conda

This route requires `git` and a working conda installation, such as Miniconda or
Miniforge. It installs the C and C++ compilers, CMake, make, zlib, and htslib into a new
environment, so you do not need to install a system compiler or htslib first.

```bash
git clone https://github.com/jibeomko/PanIsoGuard.git
cd PanIsoGuard
conda create -y -n panisoguard --override-channels -c conda-forge -c bioconda \
    cxx-compiler c-compiler cmake make zlib "htslib>=1.18"
conda activate panisoguard

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

bash examples/tiny/run.sh
bash examples/multi_caller/run.sh
```

If you already have a checkout, start from its root instead of cloning again. If an
environment named `panisoguard` already exists, choose a different name in both conda
commands. `--override-channels` uses only conda-forge and bioconda, regardless of your
configured channels. Creating the environment downloads the build dependencies; the
examples then use the local bundled data.

The executable is `build/panisoguard`. To call it as `panisoguard` in your current shell,
run this from the repository root:

```bash
export PATH="$PWD/build:$PATH"
```

### Use your own compiler and htslib

You need a C++17 compiler, CMake 3.20 or newer, and htslib 1.18 or newer. Check your
distribution's htslib version before building. CMake searches the active `$CONDA_PREFIX`
for htslib; for another installation, add `-DCMAKE_PREFIX_PATH=/path/to/prefix` to the
configure command, replacing `/path/to/prefix` with the prefix containing htslib's
`include/` and `lib/` directories.

From the root of a source checkout:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build            # unit, integration and example tests
./build/panisoguard --version
```
