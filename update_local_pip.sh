# Build new whell
rm -rf build pypKa.egg-info
python3 setup.py sdist bdist_wheel
latest_wheel=`ls -t dist/pypKa*.whl | head -n 1`

# Upload to Pypi
python3 -m twine upload ${latest_wheel::-17}*

# Local install
#pip3 install --upgrade $latest_wheel --user