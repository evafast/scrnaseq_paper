## Build the docker image from within folder
```bash
docker build --rm -t signac_seurat:1.0.0 .
```

## This changes the name of the image
```bash
docker tag 9e18d0cf3d23 evafast1/signac_seurat:1.0.0
```

## this does the login
```bash
docker login -u evafast1 -p <password>
```

## this pushes to the repo
```bash
docker push evafast1/signac_seurat:1.0.0
```

## Launch the docker image
```bash
docker run \
--rm \
-d \
--name signac \
-p 8882:8888 \
-e JUPYTER_ENABLE_LAB=YES \
-v /Users/efast/Documents/:/home/jovyan/work \
signac_seurat:1.0.0
```