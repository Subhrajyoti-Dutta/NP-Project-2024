git pull
git add .
set name=git config user.name
echo %name%
git commit -m "`git config user.name` %DATE% %TIME%"