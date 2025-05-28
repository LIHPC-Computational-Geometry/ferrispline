# vtkConverter

## Image for CI

One can run:

```bash
cd /path/to/vtk_converter/.devcontainer/
docker login registry.gitlab.com
docker build --network=host --tag registry.gitlab.com/maxime-stauffert/vtk_converter --target python .
docker push registry.gitlab.com/maxime-stauffert/vtk_converter
```
