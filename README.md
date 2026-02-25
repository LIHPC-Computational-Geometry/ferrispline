# vtkConverter

## Image for CI

One can run:

```bash
cd /path/to/vtk_converter/
docker login registry.gitlab.com
docker build --network=host --tag registry.gitlab.com/maxime-stauffert/vtk_converter:latest --target latest .devcontainer/
docker push registry.gitlab.com/maxime-stauffert/vtk_converter
```

## Image for dev

One can run:

```bash
cd /path/to/vtk_converter/
docker build --build-arg USER=${USER} --network=host --tag registry.gitlab.com/maxime-stauffert/vtk_converter:dev --target dev .devcontainer/
```

## Run tests

One can run:

```bash
xhost +local:docker
docker run --interactive --network=host --rm --tty --volume ./:/home/${USER}/vtk_converter/ registry.gitlab.com/maxime-stauffert/vtk_converter:dev
cd vtk_converter/
pre-commit run
pytest
```
