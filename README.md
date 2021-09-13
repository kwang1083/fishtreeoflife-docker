# fishtreeoflife-docker

## Updating the Docker image

```
docker build -t jonchang/fishtreeoflife-docker:master .
docker push jonchang/fishtreeoflife-docker:master
```

## Inspecting the image interactively

```
docker run -it jonchang/fishtreeoflife-docker:master
```

## Copying a file

```
docker run --rm --entrypoint cat jonchang/fishtreeoflife-docker:master  /path/to/file > path/to/destination
```


## Building without cache

```
docker build -t jonchang/fishtreeoflife-docker:master .
```
