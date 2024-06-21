#!/bin/bash

if [ -z "$1" ]; then 
	echo "usage: debug_image.sh <container args>"
	exit -1
fi

echo docker run --network="host" -u $(id -u ${USER}):$(id -g ${USER}) -v `pwd`:/docker_mount -it --rm registry.nrg.wustl.edu/docker/nrg-repo/facemasking:latest run_facemasking2_xnat $@
docker run --network="host" -u $(id -u ${USER}):$(id -g ${USER}) -v `pwd`:/docker_mount -it --rm registry.nrg.wustl.edu/docker/nrg-repo/facemasking:latest run_facemasking2_xnat $@

