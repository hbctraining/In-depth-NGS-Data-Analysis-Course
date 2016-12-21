To update the master branch of this repo with the most updated branch (rebase the master), perform the following steps on your local computer:
```
git pull origin newest_branch
git checkout newest_branch
git merge -s ours master
git checkout master
git merge newest_branch
```
