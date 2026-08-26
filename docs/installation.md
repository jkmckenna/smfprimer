# Installation

smfprimer requires Python 3.11 or newer. Primer3-py and Biopython are installed
as package dependencies.

## Development installation

From a clone of the repository:

```console
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
smfprimer --version
```

## Documentation tools

Install the documentation dependency group with a package manager that
supports PEP 735 groups, or install the declared dependency directly:

```console
python -m pip install "mkdocs-material==9.*"
mkdocs serve
```

Open `http://127.0.0.1:8000` for live preview. A production build is created
with:

```console
mkdocs build --strict
```

The generated `site/` directory is disposable and should not be committed.

## Bowtie 1

Bowtie is optional. It is needed only when `--bowtie-index` is supplied.
Confirm that the Bowtie 1 executable—not Bowtie 2—is available:

```console
bowtie --version
bowtie-build --version
```
