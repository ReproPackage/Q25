FROM ubuntu:22.04

LABEL maintainer="repropackage@qsw25.github.com"

WORKDIR /repro

RUN apt update \
    && DEBIAN_FRONTEND=noninteractive apt install -y build-essential wget cmake gfortran gobjc gobjc++ gnustep gnustep-devel libbz2-dev liblzma-dev libpcre2-dev libcurl4-openssl-dev libcairo2-dev libtiff5-dev libreadline-dev libxml2-dev libharfbuzz-dev libfribidi-dev libglpk-dev libgsl-dev libgmp-dev libmpc-dev libudunits2-dev libgdal-dev libmagick++-dev \
    && apt clean \
    && rm -rf /var/lib/apt/lists/*

ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda
	
ENV PATH=$CONDA_DIR/bin:$PATH

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" >> /etc/apt/sources.list 
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y r-base

RUN R -e "install.packages('tidyverse', dependencies=TRUE)"
RUN R -e "install.packages('patchwork', dependencies=TRUE)"
RUN R -e "install.packages('tikzDevice', dependencies=TRUE)"
RUN R -e "install.packages('scales', dependencies=TRUE)"
RUN R -e "install.packages('ggh4x', dependencies=TRUE)"
RUN R -e "install.packages('ggpmisc', dependencies=TRUE)"

RUN conda create -n quark_install python=3.10

SHELL ["conda", "run", "-n", "quark_install", "/bin/bash", "-c"]

RUN conda config --add channels conda-forge
RUN conda update -n base -c defaults conda
RUN conda install -c dlr-sc quark=1.1
RUN conda install conda-forge::scip=9.0
RUN conda install -c conda-forge pathos
RUN pip install pandas
RUN pip install p_tqdm
RUN pip install sympy
RUN pip install scipy
RUN pip install sortedcontainers
RUN pip install multiset

COPY ./src/ ./

#OpenShell
CMD ["/bin/bash"]
