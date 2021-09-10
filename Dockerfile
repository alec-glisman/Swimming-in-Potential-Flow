# Base image
FROM ubuntu:20.04

# Label directives
LABEL maintainer=alec.glisman@gmail.com
LABEL project=bodies-in-potential-flow
LABEL version=0.1
LABEL environment=dev

# General environment directives
ENV DEBIAN_FRONTEND=noninteractive 


# Update packages on base image and install sudo
RUN apt-get update --fix-missing && \
    apt-get upgrade -y && \
    apt-get install -y sudo build-essential software-properties-common

# Add additional package repositories
RUN add-apt-repository ppa:ubuntu-toolchain-r/test

# Install APT packages
RUN apt-get update && apt-get install -y \
    git wget curl file \
    zsh fonts-powerline
# command line packages
RUN apt-get update && apt-get install -y \
    cmake \
    gcc-11 g++-11 \
    nvidia-cuda-toolkit \
    perl cpanminus \
    python3 python3-dev python3-pip
# coding languages
RUN apt-get update && apt-get install -y \
    libboost-all-dev libspdlog-dev catch \
    locales tzdata
# library depdencencies

# Install eigen3 (v3.3.9-2)
WORKDIR "/tmp"
RUN wget -O libeigen3-dev_3.3.9-2_all.deb http://launchpadlibrarian.net/519614686/libeigen3-dev_3.3.9-2_all.deb
RUN apt-get install -y ./libeigen3-dev_3.3.9-2_all.deb
RUN rm libeigen3-dev_3.3.9-2_all.deb

# Install Intel MKL via Intel OneAPI
WORKDIR "/tmp"
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
RUN apt-get update && apt-get install -y intel-oneapi-mkl
ENV PATH="/opt/intel/oneapi/compiler/latest/linux/bin/intel64:${PATH}"

# Download oh-my-zsh and make it the default terminal
RUN wget https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh -O - | zsh || true
RUN chsh -s $(which zsh)


# Copy the folders that are not ignored into the Docker image
WORKDIR "/"
COPY . /bodies-in-potential-flow/


# Encoding
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Change terminal timezone to PDT 
ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN sudo dpkg-reconfigure -f noninteractive tzdata

# Reduce image size
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get autoclean && apt-get autoremove


# Python requirements
WORKDIR "/bodies-in-potential-flow/requirements/Python"
RUN pip3 install --upgrade pip && \
    pip3 install --upgrade virtualenv
RUN pip3 install -r requirements.txt

# Perl requirements
WORKDIR "/bodies-in-potential-flow/requirements/Perl"
RUN cpanm --installdeps .


# REVIEW: Other options of commands to run in Docker container
#     CMD ["zsh"] # launches zsh terminal to test run commands after build
WORKDIR "/bodies-in-potential-flow"
ENTRYPOINT [ "perl", "scripts/collinear-swimmer-wall.pl" ]