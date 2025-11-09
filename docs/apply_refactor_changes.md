# Applying the pipeline export refactor

The refactor that introduces `finalize_pipeline_outputs` and `export_svg_safe` is
committed on the `work` branch. Follow the steps below to bring those changes
into another clone of the repository.

## 1. Ensure a clean working tree
Before updating, commit or stash any local edits:

```bash
git status
git add <files>
git commit -m "Save local work"  # or use `git stash`
```

## 2. Fetch the updated branch
Retrieve the latest commits from the remote that hosts the refactor (replace
`origin` if you use another remote name):

```bash
git fetch origin
```

## 3. Check out the target branch
Switch to the branch that carries the refactor. If you do not have it locally
yet, create it to track the remote branch:

```bash
git checkout work                # if the branch already exists locally
# or
git checkout -b work origin/work  # create a local branch that tracks origin/work
```

## 4. Fast-forward to the latest commit
Update your local branch to match the remote history:

```bash
git pull --ff-only
```

If the pull reports that the branch has diverged because of your local commits,
use `git pull --rebase` instead or reset to the remote tip with:

```bash
git reset --hard origin/work
```

## 5. Verify the refactor is present
List the most recent commits and confirm that "Refactor pipeline export helpers"
appears:

```bash
git log --oneline -5
```

Once the commit is visible, you have the centralized export helper and shared
SVG fallback available locally.

## 6. Run the pipelines (optional)
Add the `octave` folder to your Octave path and invoke any pipeline function.
They now rely on the shared helpers automatically.

```octave
addpath('octave');
pipeline_activacion(...);
```

That is all you need to incorporate the proposed changes into your repository.
