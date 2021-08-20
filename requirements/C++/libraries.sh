#!/usr/bin/env bash

# NOTE: Assumes that homebrew is already installed on computer
brew install boost spdlog catch2
brew install eigen --HEAD  # HEAD needed for Eigen::seqN()


# NOTE: Using Assuming DEBIAN linux environment
# Intel GPU Drivers
# @SOURCE: https://dgpu-docs.intel.com/installation-guides/ubuntu/ubuntu-focal.html
sudo apt-get install -y gpg-agent wget
wget -qO - https://repositories.intel.com/graphics/intel-graphics.key | sudo apt-key add -
sudo apt-add-repository 'deb [arch=amd64] https://repositories.intel.com/graphics/ubuntu focal main'

# Intel (oneAPI dependencies)
# @SOURCE: https://dgpu-docs.intel.com/installation-guides/ubuntu/ubuntu-focal.html
sudo apt install -y intel-opencl-icd intel-level-zero-gpu level-zero intel-media-va-driver-non-free libmfx1
sudo apt install -y libigc-dev intel-igc-cm libigdfcl-dev libigfxcmrt-dev level-zero-dev


# Install OneAPI
# @SOURCE: https://software.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top/installation/install-using-package-managers/apt.html

# use wget to fetch the Intel repository public key
cd /tmp
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
# add to your apt sources keyring so that archives signed with this key will be trusted.
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
# remove the public key
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB

# configure the APT client to use Intel's repository
echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list

# Install new packages
sudo apt update
sudo apt install -y intel-basekit
sudo apt install -y intel-hpckit
