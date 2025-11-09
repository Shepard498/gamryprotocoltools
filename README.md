# Octave Gamry Pipelines

Turnkey structure to parse Gamry `.DTA` files, run processing pipelines (EIS, activation, polarization and OCP), and export plots and spreadsheets.

## Quick start
```bash
# clone and enter (after you create a GitHub repo)
git clone <your-repo-url>.git
cd <repo>

# optional: set up a venv for helpers (not required for Octave)
# python -m venv .venv && source .venv/bin/activate

# run a smoke test (Octave must be installed)
octave -qf --eval "run('octave/scripts/run_all.m')"
```

## Layout
```
octave/           # all .m code lives here
  lib/            # small utilities (e.g., merge_opts)
  pipelines/      # pipeline_*.m entry points
  scripts/        # top-level runner scripts
data/
  raw/            # input .DTA files (never edited)
  interim/        # temporary intermediate exports
  processed/      # final spreadsheets and cleaned data
plots/            # .svg/.png outputs
docs/             # longer docs, images for README
examples/         # sample references / usage snippets
.github/
  workflows/      # CI configuration
```

## Contributing
- Open an issue first for big changes.
- Use a feature branch and open a PR into `main`.
- Add or update tests in `octave/tests` when changing logic.

## License
MIT (see `LICENSE`).

