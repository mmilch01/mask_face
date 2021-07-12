sudo docker tag registry.nrg.wustl.edu/docker/nrg-repo/facemasking:latest \
    xnat/facemasking:1.0
sudo docker login docker.io
sudo docker push xnat/facemasking:1.0

