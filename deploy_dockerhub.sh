docker tag registry.nrg.wustl.edu/docker/nrg-repo/facemasking:latest \
    xnat/facemasking:1.1
docker login docker.io
docker push xnat/facemasking:1.1

