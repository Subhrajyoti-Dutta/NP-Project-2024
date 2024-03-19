git pull
git add .
for /f %%i in ('git config user.name') do set VAR=%%i
echo %name%
git commit -m " %DATE% %TIME%"