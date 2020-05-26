pytest -s --cov=../ test.py
coverage combine */*/.coverage .coverage
coverage report -i --omit=*__init__.py,*__main__.py,*mc_setup.py,*dist-packages/numpy*,*delphi4py.py
coverage xml -i --omit=*__init__.py,*__main__.py,*mc_setup.py,*dist-packages/numpy*,*delphi4py.py
coverage html -i --omit=*__init__.py,*__main__.py,*mc_setup.py,*dist-packages/numpy*,*delphi4py.py

export CODACY_PROJECT_TOKEN=d28668518b1944cabb9336acf2efe032
python-codacy-coverage -r coverage.xml

google-chrome htmlcov/index.html 
