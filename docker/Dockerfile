FROM ubuntu:16.04

SHELL [ "/bin/bash", "--login", "-c" ]

# Create a non-root user
ARG username=carlos-chagas
ARG uid=1000
ARG gid=100
ENV USER $username
ENV UID $uid
ENV GID $gid
ENV HOME /home/$USER
RUN adduser --disabled-password \
    --gecos "Non-root user" \
    --uid $UID \
    --gid $GID \
    --home $HOME \
    $USER

COPY environment.yml /tmp/
RUN chown $UID:$GID /tmp/environment.yml

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git wget build-essential cmake

USER $USER

# install miniconda
ENV MINICONDA_VERSION 4.8.2
ENV CONDA_DIR $HOME/miniconda3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_$MINICONDA_VERSION-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

# make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH

# make conda activate command available from /bin/bash --login shells
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile

# make conda activate command available from /bin/bash --interative shells
RUN conda init bash

ENV CXX /usr/bin/gcc

# create a project directory inside user home
WORKDIR $HOME
RUN git clone --recursive https://github.com/rodtheo/RLZI
WORKDIR $HOME/RLZI
RUN git submodule update --remote
RUN chmod +x install.sh
RUN bash install.sh

ENV PATH=$HOME/RLZI/bin:$PATH

# build the conda environment
ENV ENV_PREFIX $PWD/env
RUN conda update --name base --channel defaults conda && \
    conda env create --prefix $ENV_PREFIX --file /tmp/environment.yml --force && \
    conda clean --all --yes

RUN conda activate rlzi-env
