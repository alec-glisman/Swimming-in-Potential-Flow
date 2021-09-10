# Base image
FROM ubuntu:20.04

# Label directives
LABEL maintainer=alec.glisman@gmail.com
LABEL project=bodies-in-potential-flow
LABEL version=0.1
LABEL environment=dev

# General environment directives
ENV DEBIAN_FRONTEND=noninteractive 

# Copy the folders that are not ignored into the Docker image
COPY . /bodies-in-potential-flow/
WORKDIR "/bodies-in-potential-flow"

# Update packages on base image
RUN apt-get -y update
RUN apt-get upgrade -y
RUN apt-get install -y sudo


# Install APT packages
RUN apt-get install -y build-essential \        
    software-properties-common \
    git \
    wget \
    curl \
    zsh \
    fonts-powerline \
    cmake \
    perl \
    cpanminus \
    libboost-all-dev \
    libspdlog-dev \
    catch \
    tzdata

# Install gcc-11
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update && apt-get install gcc-11 g++-11
# Set gcc-11 as default
RUN update-alternatives --install \
    /usr/bin/gcc gcc /usr/bin/gcc-9 90 \
    --slave /usr/bin/g++ g++ /usr/bin/g++-9 \
    --slave /usr/bin/gcov gcov /usr/bin/gcov-9 \
    --slave /usr/bin/gcc-ar gcc-ar /usr/bin/gcc-ar-9 \
    --slave /usr/bin/gcc-ranlib gcc-ranlib /usr/bin/gcc-ranlib-9 \
    --slave /usr/bin/cpp cpp /usr/bin/cpp-9 
RUN update-alternatives --install \
    /usr/bin/gcc gcc /usr/bin/gcc-11 110 \
    --slave /usr/bin/g++ g++ /usr/bin/g++-11 \
    --slave /usr/bin/gcov gcov /usr/bin/gcov-11 \
    --slave /usr/bin/gcc-ar gcc-ar /usr/bin/gcc-ar-11 \
    --slave /usr/bin/gcc-ranlib gcc-ranlib /usr/bin/gcc-ranlib-11  \ 
    --slave /usr/bin/cpp cpp /usr/bin/cpp-11

# Install miniconda to /miniconda
WORKDIR "/"
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

# Install Intel MKL via OneAPI
WORKDIR "/tmp"
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
RUN apt-get update
RUN apt-get install -y intel-basekit

# Install homebrew
WORKDIR "$HOME"
RUN git clone https://github.com/Homebrew/brew $HOME/.linuxbrew/Homebrew \
    && mkdir $HOME/.linuxbrew/bin \
    && ln -s ../Homebrew/bin/brew $HOME/.linuxbrew/bin \
    && eval $($HOME/.linuxbrew/bin/brew shellenv) \
    && brew --version
ENV PATH=$HOME/.linuxbrew/bin:$PATH

# Install eigen3 head
RUN brew install eigen --HEAD


# Conda requirements
WORKDIR "/bodies-in-potential-flow"
RUN conda env create -f requirements/Python/environment.yml 

# Perl requirements
WORKDIR "/bodies-in-potential-flow/requirements/Perl"
RUN cpanm --installdeps .


# Download oh-my-zsh and make it the default terminal
RUN wget https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh -O - | zsh || true
RUN chsh -s $(which zsh)

# Change terminal timezone to PDT 
ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN sudo dpkg-reconfigure -f noninteractive tzdata


# REVIEW: Other options of commands to run in Docker container
#     CMD ["zsh"] # launches zsh terminal to test run commands after build
# Run the program to build from CMAKE
ENTRYPOINT ["perl", "scripts/collinear-swimmer-wall.pl"]  
