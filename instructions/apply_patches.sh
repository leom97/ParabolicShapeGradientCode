#!/bin/bash

echo "Insert path to environment site-packages (e.g. \$HOME/anaconda3/envs/article_env/venv_39/lib/python3.9/site-packages, without final slash):"
read site_packages_path

cp "fixes/fenics/tools.py" $site_packages_path"/ffc/uflacs/tools.py"