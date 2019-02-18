# fishtreeoflife-docker

## Updating the Docker image

```
docker build -t jonchang/fishtreeoflife-docker:latest .
docker push jonchang/fishtreeoflife-docker:latest
```

## Inspecting the image interactively

```
docker run -it jonchang/fishtreeoflife-docker:latest
```

## Copying a file

```
docker run --rm --entrypoint cat jonchang/fishtreeoflife-docker:latest  /path/to/file > path/to/destination
```


## Building without cache

```
docker build -t jonchang/fishtreeoflife-docker:latest .
```

## See image layers

```

``
