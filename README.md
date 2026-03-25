# vtkConverter

## Image for dev

One can run:

```bash
cd /path/to/vtk_converter/
CI_COMMIT_REF_NAME="$(git rev-parse --abbrev-ref HEAD)"
CI_REGISTRY_IMAGE='registry.gitlab.com/maxime-stauffert/vtk_converter'
docker build --build-arg USER=${USER} \
             --network=host \
             --tag $CI_REGISTRY_IMAGE:$CI_COMMIT_REF_NAME \
             --target dev \
             .devcontainer/
```

## Run tests

One can run:

```bash
xhost +local:docker
docker run --interactive \
           --network=host \
           --rm \
           --tty \
           --volume ./:/home/${USER}/vtk_converter/ \
           $CI_REGISTRY_IMAGE:$CI_COMMIT_REF_NAME
cd vtk_converter/
pre-commit run
pytest
cargo test
```
