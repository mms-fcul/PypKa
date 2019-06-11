python setup.py sdist bdist_wheel
latest_wheel=`ls -t dist/pypka*.whl | head -n 1`
pip install --upgrade $latest_wheel --user
