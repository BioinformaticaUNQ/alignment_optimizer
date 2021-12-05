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
$ python src/project_alignment_optimizer/program/main.py

# COMMANDS:
# Run alignment:
$ python src/project_alignment_optimizer/program/main.py align -f {FILE_PATH.fasta} -qs {QUERY_SEQUENCE}

# Run config -h (get config helper):
$ python src/project_alignment_optimizer/program/main.py config -h
# Run config -r (reset config.env file to default):
$ python src/project_alignment_optimizer/program/main.py config -r True
# Run config -k KEY -v VALUE (change value from config.env):
$ python src/project_alignment_optimizer/program/main.py config -k {KEY_NAME} -v {NEW_VALUE}

# Run view_config (get config.env variables):
$ python src/project_alignment_optimizer/program/main.py view_config
```

## 4. Run tests
```bash
# run tests with python version 3.8
$ pytest
```
