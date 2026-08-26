# Contributing to the documentation

Documentation source is ordinary Markdown under `docs/`. Navigation is
explicitly defined in `mkdocs.yml`; add every new user-facing page there.

## Preview locally

```console
python -m pip install "mkdocs-material==9.*"
mkdocs serve
```

## Validate before committing

```console
mkdocs build --strict
pytest -q
ruff check .
ruff format --check .
```

Strict builds reject broken internal links and pages omitted from navigation.

## Tutorial conventions

- Run commands from the repository root.
- Prefer files under `docs/examples/` so examples are reproducible.
- Do not commit human genomes, Bowtie indexes, private plasmid maps, or large
  generated outputs.
- State coordinate conventions and expected output artifacts.
- Mark relaxed demonstration constraints clearly so users do not mistake them
  for assay recommendations.

Pushes to `main` build and deploy the site through the documentation GitHub
Actions workflow. Pull requests build the site without deploying it.

## Enable GitHub Pages once

Before the first deployment, a repository administrator must open
**Settings → Pages** and select **GitHub Actions** under **Build and
deployment → Source**. GitHub's default workflow token can deploy an enabled
Pages site, but it cannot enable Pages for a repository that has no Pages site
yet.
