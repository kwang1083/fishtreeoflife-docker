run:
	docker run -it jonchang/fishtreeoflife-docker:latest

build:
	docker build -t jonchang/fishtreeoflife-docker:latest .

push:
	docker push jonchang/fishtreeoflife-docker:latest
