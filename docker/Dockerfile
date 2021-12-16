# SECTION: Preamble and general set-up

# Base Image: Ubuntu 20.04 LTS
FROM ubuntu:20.04

# Label: Image
LABEL maintainer=alec.glisman@gmail.com
LABEL project=bodies-in-potential-flow
LABEL version=0.1
LABEL environment=dev

# Environment: OS
ENV DEBIAN_FRONTEND=noninteractive 

# Environment: Intel OneAPI
ENV ONEAPI_BASE="/opt/intel/oneapi"

# Environment: Intel MKL
# ENV MKL_DEBUG_CPU_TYPE=5
ENV MKLROOT="${ONEAPI_BASE}/mkl/latest"
ENV PATH="${MKLROOT}/bin/intel64:${PATH}"
ENV LDFLAGS="${LDFLAGS} -L${MKLROOT}/lib/intel64"
ENV CPPFLAGS="${CPPFLAGS} -L${MKLROOT}/include/intel64"

# Update packages on base image and general use packages
RUN apt-get update && apt-get install -y apt-utils

RUN apt-get update --fix-missing && \
    apt-get upgrade -y

RUN apt-get update && apt-get install -y \
    build-essential \
    software-properties-common \
    sudo \
    wget

# !SECTION: (Preamble and general set-up)



# SECTION: Add repositories

# Add newer gcc package repositories
RUN add-apt-repository ppa:ubuntu-toolchain-r/test

# Add intel package repositories
WORKDIR "/tmp"
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list

# !SECTION (Add repositories)



# SECTION: APT package manager installs

# OS Date and time
RUN apt-get update && apt-get install -y \
    locales tzdata

# MAKE
RUN apt-get update && apt-get install -y \
    make \
    cmake

# Ninja generator
RUN apt-get update && apt-get install -y \
    ninja-build

# C++ and libraries
RUN apt-get update && apt-get install -y \
    gcc-11 \
    g++-11 \
    lcov \
    gcovr \
    libboost-all-dev

# Perl and package manager
RUN apt-get update && apt-get install -y \
    perl \
    cpanminus

# Python
RUN apt-get update && apt-get install -y \
    python3 \
    python3-dev \
    python3-pip

# Intel MKL (v2021.3.0) via Intel OneAPI repositories
RUN apt-get update && apt-get install -y \
    intel-oneapi-mkl-devel-2021.3.0

# Reduce image size
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get autoclean && apt-get autoremove

# !SECTION (APT package manager installs)



# SECTION: OS Modifications

# Encoding
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Change terminal timezone to America/Los_Angeles 
ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN dpkg-reconfigure -f noninteractive tzdata

# !SECTION (OS Modifications)



# SECTION: Language dependencies

# Copy requirements directory into docker image
WORKDIR "/"
COPY requirements /bodies-in-potential-flow/requirements

# Perl dependencies
WORKDIR "/bodies-in-potential-flow/requirements/Perl"
RUN cpanm --installdeps .

# Python dependencies
WORKDIR "/bodies-in-potential-flow/requirements/Python"
RUN pip3 install --upgrade pip && \
    pip3 install -r requirements.txt

# vcpkg
# RUN apt-get update && apt-get install -y \
#     curl \
#     zip \
#     unzip \
#     tar

# Reduce image size
WORKDIR "/bodies-in-potential-flow"
RUN rm -rf requirements

# !SECTION (Language dependencies)



# SECTION: Repository include package installs

# Copy include directory into docker image
WORKDIR "/"
COPY include /bodies-in-potential-flow/include
# COPY vcpkg /bodies-in-potential-flow/vcpkg

# Install eigen3 (commit 8ce341caf2947e4b5ac4580c20254ae7d828b009 from git repo)
WORKDIR "/bodies-in-potential-flow/include/Eigen"
RUN mkdir build
WORKDIR "/bodies-in-potential-flow/include/Eigen/build"
RUN cmake .. && make install

# Install spdlog (v1.9.1 from git repo)
WORKDIR "/bodies-in-potential-flow/include/spdlog"
RUN mkdir build
WORKDIR "/bodies-in-potential-flow/include/spdlog/build"
RUN cmake .. && make install

# Install catch2 (v2.13.7 from git repo)
WORKDIR "/bodies-in-potential-flow/include/Catch2"
RUN mkdir build
WORKDIR "/bodies-in-potential-flow/include/Catch2/build"
RUN cmake .. && make install

# vcpkg setup
# WORKDIR "/bodies-in-potential-flow/vcpkg"
# RUN ./bootstrap-vcpkg.sh -disableMetrics

# Reduce image size
WORKDIR "/bodies-in-potential-flow"
RUN rm -rf include
COPY include/cmake_modules /bodies-in-potential-flow/include/cmake_modules

# !SECTION (Repository include package installs)



# SECTION: Intel/NVIDA package initialization

# Activate Intel MKL variables
SHELL ["/bin/bash", "-c"] 
RUN source "${MKLROOT}/env/vars.sh"

# Activate VCPKG
COPY vcpkg /bodies-in-potential-flow/vcpkg
COPY vcpkg.json /bodies-in-potential-flow/vcpkg.json
RUN /bodies-in-potential-flow/vcpkg/./bootstrap-vcpkg.sh

# !SECTION (Intel/NVIDA package initialization)



# SECTION: Prepare simulation

# Copy directories for simulation
COPY input /bodies-in-potential-flow/input
COPY python /bodies-in-potential-flow/python
COPY scripts /bodies-in-potential-flow/scripts
COPY src /bodies-in-potential-flow/src
COPY test /bodies-in-potential-flow/test
COPY CMakeLists.txt /bodies-in-potential-flow/CMakeLists.txt

# Make output directories for data IO from container
WORKDIR "/bodies-in-potential-flow"
RUN mkdir data
RUN mkdir temp

# Build simulation
ENTRYPOINT [ "perl", "scripts/collinear-swimmer-internal-dynamics.pl" ]   # starts simulations from perl script

# !SECTION (Prepare simulation)
