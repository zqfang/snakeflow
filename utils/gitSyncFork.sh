#!bin/sh


# ref: https://help.github.com/articles/syncing-a-fork/

# step 1: set upstream
git remote -v
git remote add upstream https://github.com/fastai/courses

# step 2: Fetch the branches and their respective commits from the upstream repository.
git fetch upstream

# checkout local master branch
git checkout master
# Merge the changes from upstream/master into your local master branch. 
git merge upstream/master


