# MMesh
MMesh creates, reads and writes mesh for Marine Energy Resource Assessment Canada.

## Installation
This package was developed,tested and built using conda. MMesh uses mshapely, shapely, numpy, scipy, fiona, matplotlib.
Only tested with python >=3.6

```bash
conda create -n mmesh python=3.8
conda activate mmesh
conda install -c meracan mmesh
```

For developers and debugging:
```bash
conda create -n mmesh python=3.8
conda activate osmgmsh
conda install -c conda-forge numpy scipy fiona shapely pyproj requests geojson tqdm matplotlib gmsh
pip install -e ./mshapely
pip install -e ./slfpy
pip install -e ./mmesh
```

### Usage, user guide and examples
[Docs](doc/doc_mmesh.ipynb)

### Testing
[Docs](test/README.md)

### License
[License](LICENSE)