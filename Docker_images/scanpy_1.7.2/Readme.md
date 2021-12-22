## Build the docker image from within folder
```bash
docker build --rm -t scanpy:1.7.2 .
```

## This changes the name of the image
```bash
docker tag 3360c6e4dabe evafast1/scanpy:1.7.2
```

## this does the login
```bash
docker login -u evafast1 -p <password>
```

## this pushes to the repo
```bash
docker push evafast1/scanpy:1.7.2
```

## Launch the docker image
```bash
docker run \
--rm \
-d \
--name scanpy2 \
-p 8881:8888 \
-e JUPYTER_ENABLE_LAB=YES \
-v /Users/efast/:/home/jovyan/work \
scanpy:1.7.2
```