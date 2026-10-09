# R 4.5.1, Node.js 20 and the three agent CLIs in one image.
#
#   docker build -t scrna-benchmark .
#   docker run --rm -it -v "$PWD/.env:/benchmark/.env:ro" \
#     -v "$PWD/api_benchmark/results:/benchmark/api_benchmark/results" \
#     -v "$PWD/coding_agent_benchmark/results:/benchmark/coding_agent_benchmark/results" \
#     scrna-benchmark

FROM node:20-bookworm-slim AS node

FROM rocker/r-ver:4.5.1

# CRAN packages are taken from a dated snapshot so every build gets the same versions
ARG CRAN_DATE=2026-08-02
RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/noble/${CRAN_DATE}'))" \
      >> "${R_HOME}/etc/Rprofile.site"

RUN apt-get update && apt-get install -y --no-install-recommends \
      git cmake libcurl4-openssl-dev libssl-dev libxml2-dev libhdf5-dev libglpk-dev \
      libpng-dev libjpeg-dev libtiff-dev libfontconfig1-dev libfreetype6-dev \
      libharfbuzz-dev libfribidi-dev libcairo2-dev libgsl-dev \
      python3 python3-numpy python3-matplotlib \
    && rm -rf /var/lib/apt/lists/*

# Node.js and npm
COPY --from=node /usr/local/bin/node /usr/local/bin/node
COPY --from=node /usr/local/lib/node_modules /usr/local/lib/node_modules
RUN ln -s ../lib/node_modules/npm/bin/npm-cli.js /usr/local/bin/npm \
 && ln -s ../lib/node_modules/npm/bin/npx-cli.js /usr/local/bin/npx

WORKDIR /benchmark

# R packages (Bioconductor 3.21, the release used for the benchmark)
ENV R_BIOC_VERSION=3.21
COPY setup_env.R .
RUN Rscript setup_env.R

# agent CLIs for the coding-agent benchmark
RUN npm install -g @anthropic-ai/claude-code @openai/codex @google/gemini-cli

# Node dependencies, then the benchmark itself
COPY package.json package-lock.json ./
RUN npm ci
COPY . .

# container paths; these take precedence over R_EXEC and BENCH_DATASET_DIR in a mounted .env
ENV R_EXEC=Rscript
ENV BENCH_DATASET_DIR=/benchmark/data/pbmc3k/filtered_gene_bc_matrices/hg19

CMD ["bash"]