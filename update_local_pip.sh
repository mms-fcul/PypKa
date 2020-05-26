python3 setup.py sdist bdist_wheel
latest_wheel=`ls -t dist/pypka*.whl | head -n 1`
pip3 install --upgrade $latest_wheel --user

#python3 -m twine upload dist/pypka-0.0.6*
