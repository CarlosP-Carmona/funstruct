# Publishing funstruct to GitHub — first-time setup

One-time steps, in order. Total time ~10 minutes. Everything happens in the
Terminal (macOS: Applications > Utilities > Terminal, or the Terminal tab in
RStudio).

## 1. Create the empty repository on github.com

1. Log in at github.com (user: CarlosP-Carmona).
2. Click the "+" (top right) > "New repository".
3. Repository name: `funstruct`. Leave EVERYTHING else unticked
   (no README, no .gitignore, no licence — we already have those).
4. Click "Create repository". Keep the page open.

## 2. Tell git who you are (one-time, per computer)

```bash
git config --global user.name  "Carlos P. Carmona"
git config --global user.email "perezcarmonacarlos@gmail.com"
```

## 3. Initialize and push

In Terminal, go to the funstruct folder (adjust the path if you moved it):

```bash
cd "/Users/carlosperezcarmona/Dropbox/Manuscritos en marcha/R PACKAGE for Functional structure/funstruct"
git init -b main
git add .
git commit -m "Scaffold: fspace class, fs_space, as_fspace, fs_reduce, fs_rotate"
git remote add origin https://github.com/CarlosP-Carmona/funstruct.git
git push -u origin main
```

When git asks for a password: GitHub no longer accepts account passwords on
the command line. Two options, easiest first:

- **GitHub Desktop** (recommended if the above feels alien): install from
  desktop.github.com, log in, File > Add local repository, select the
  funstruct folder, then "Publish repository". It handles authentication
  for you.
- **Personal access token**: github.com > your avatar > Settings >
  Developer settings > Personal access tokens > Tokens (classic) >
  Generate new token, tick "repo", copy it, and paste it as the password
  when git asks. macOS keychain remembers it afterwards.

## 4. Check that CI runs

After the push, open github.com/CarlosP-Carmona/funstruct > "Actions" tab.
A workflow called R-CMD-check starts automatically: it regenerates the
documentation and runs the full check + test suite on macOS, Windows and
two Ubuntu versions. Green tick = everything passes. Red cross = tell
Claude; the logs are readable through the GitHub API.

## 5. A note on Dropbox + git

Working with a git repo inside Dropbox is fine for a single-author project
(and gives you an extra backup), but if Dropbox ever shows sync conflicts
in the `.git` folder, move the repo out of Dropbox and keep GitHub as the
single source of truth. GitHub is the real backup from now on.

## Day-to-day afterwards

```bash
git add -A
git commit -m "what changed"
git push
```

Or the equivalent three clicks in GitHub Desktop. That's all the git you
need for this project.
