# Update the master with the most recent branch (rebase)
git pull origin newest_branch
git checkout newest_branch
git merge -s ours master
git checkout master
git merge newest_branch
