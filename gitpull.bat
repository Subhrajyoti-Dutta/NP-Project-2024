git pull
git add .
git config user.name | (set /p name= & set name)
echo %name%
git commit -m " %DATE% %TIME%"