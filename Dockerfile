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


# Update packages on base image and install sudo
RUN apt-get update && apt-get upgrade -y && apt-get install -y sudo

# Install APT packages
RUN apt-get update && apt-get install -y \
    build-essential \        
    software-properties-common \
    git \
    wget \
    curl \
    file \
    zsh \
    fonts-powerline \
    cmake \
    perl \
    cpanminus \
    libboost-all-dev \
    libspdlog-dev \
    catch \
    locales \
    tzdata

# Install gcc-11
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update && apt-get install -y gcc-11 g++-11

# Install homebrew
WORKDIR "/"
RUN mkdir /homebrew
RUN curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1 -C homebrew \
    && eval $(/homebrew/bin/brew shellenv) \
    && brew --version
ENV PATH="/homebrew/bin:${PATH}"
# Install eigen3 head
RUN brew install eigen --HEAD

# Install miniconda3 to /miniconda
WORKDIR "/"
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/miniconda/bin:${PATH}"
RUN conda update -y conda

# Install Intel MKL via Intel OneAPI
WORKDIR "/tmp"
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
RUN apt-get update && apt-get install -y intel-oneapi-mkl

# Download oh-my-zsh and make it the default terminal
RUN wget https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh -O - | zsh || true
RUN chsh -s $(which zsh)

# Conda requirements
WORKDIR "/bodies-in-potential-flow"
RUN conda env create -f requirements/Python/environment.yml 

# Perl requirements
WORKDIR "/bodies-in-potential-flow/requirements/Perl"
RUN cpanm --installdeps .


# Encoding
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Change terminal timezone to PDT 
ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN sudo dpkg-reconfigure -f noninteractive tzdata

# Reduce image size
RUN rm -rf /var/lib/apt/lists/*


# REVIEW: Other options of commands to run in Docker container
#     CMD ["zsh"] # launches zsh terminal to test run commands after build
WORKDIR "/bodies-in-potential-flow"
ENTRYPOINT ["perl", "scripts/collinear-swimmer-wall.pl"]  # Run the program to build from CMAKE
