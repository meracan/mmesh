# Mesh helper





For developers and debugging:
```bash
mkdir ../data
cd osmgmsh
conda activate osmgmsh
PYTHONPATH="../mshapely/:../slfpy/:../mmesh/" python3 test/test_io.py
PYTHONPATH="../mshapely/:../slfpy/:../mmesh/" python3 test/test_mmesh.py

```

```
jupyter notebook --ip=0.0.0.0 --port=8080 --no-browser
```