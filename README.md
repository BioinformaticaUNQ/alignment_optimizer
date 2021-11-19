# Steps
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

## 3. Run example command
```bash
# path final_proyect_alignment/ run:
$ python src/project_alignment_optimizer/main.py 10
```

## 4. Run tests
```bash
# run tests with python version 3.8
$ tox -e py38
```
