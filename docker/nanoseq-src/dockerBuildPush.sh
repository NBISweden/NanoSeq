#!/bin/bash

IMAGE_VERSION=0.1 # Update this line

IMAGE_NAME="ghcr.io/nbisweden/nanoseq-src"

# Build images

	docker build -t $IMAGE_NAME:$IMAGE_VERSION -t $IMAGE_NAME:latest . || {
		echo "Error building image, exited script"
		exit 1
	}

# Push images

 	docker push $IMAGE_NAME:$IMAGE_VERSION || {
 		echo "Error pushing image, exited script"
 		exit 1
 	}

 	docker push $IMAGE_NAME:latest || {
 		echo "Error pushing image, exited script"
 		exit 1
 	}
