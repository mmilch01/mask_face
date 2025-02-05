FROM ubuntu:22.04

ENV MCR_HOME=/usr/local/MATLAB_Runtime \
	MASKFACE_HOME=/usr/local/maskface \
	MCR_CACHE_ROOT=/tmp \
	FSLDIR=/nrgpackages/packages/fsl \
	PATH=/usr/local/miniconda3/bin:/nrgpackages/packages/fsl/bin:/nrgpackages/tools/nrg-improc:$PATH

COPY  mcr/for_redistribution_files_only/* /usr/local/maskface/
COPY  nrg-improc /nrgpackages/tools/nrg-improc/
COPY  fsl	/nrgpackages/packages/fsl/
COPY  lib64/* /lib64/

RUN apt-get update && apt-get install -y curl gpg
RUN curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg; \
  curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
RUN apt-get update && apt install -y imagemagick wget zip unzip bzip2 git bc default-jre

#libX11 libXt java 

RUN mkdir -p /usr/local/MATLAB_Runtime && \
	mkdir -p /docker_mount && \
	mkdir -p /usr/local/maskface && \
	mkdir -p /input && \
	cd /tmp && wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda3-py38_23.1.0-1-Linux-x86_64.sh; chmod +x Miniconda3-py38_23.1.0-1-Linux-x86_64.sh && \
    ./Miniconda3-py38_23.1.0-1-Linux-x86_64.sh -u -b -p /usr/local/miniconda3 && \
    rm -rf /tmp/* && \
    pip install --upgrade pip
    
RUN cd /var/local; git clone https://github.com/MIC-DKFZ/HD-BET; cd HD-BET; git checkout ae16068; pip install -e . ; cd HD_BET; \
    wget -O 0.model https://zenodo.org/record/2540695/files/0.model?download=1; \    
    sed -i "s#folder_with_parameter_files.*#folder_with_parameter_files = \'/var/local/HD-BET/HD_BET\'#g" paths.py; \
    pip install fslpy

COPY  mcr/MATLAB_MCR /usr/local/MATLAB_MCR/
RUN rm -rf /tmp/* && apt clean

WORKDIR /docker_mount
ENTRYPOINT [ "/bin/bash", "-l", "-c" ]

