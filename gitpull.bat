git pull
git add .
set /p name=
git config user.name | (set name)
echo %name%
git commit -m " %DATE% %TIME%"