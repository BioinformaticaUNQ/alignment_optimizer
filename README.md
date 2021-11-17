### Steps
## 1. Config `venv` 

```bash
$ python3 -m venv .venv
$ source .venv/bin/activate
```

## 2. Install project locally 

```bash
# install base dependencies
$ pip install -e .

# install testing dependencies
$ pip install -e .[testing]
```

