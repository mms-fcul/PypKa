#pytest -s --cov=../ test.py
#coverage combine */*/.coverage .coverage
#coverage report -i --omit=*__init__.py,*__main__.py,*mc_setup.py,*dist-packages/numpy*,*delphi4py.py
#coverage xml -i --omit=*__init__.py,*__main__.py,*mc_setup.py,*dist-packages/numpy*,*delphi4py.py
#coverage html -i --omit=*__init__.py,*__main__.py,*mc_setup.py,*dist-packages/numpy*,*delphi4py.py

export CODACY_PROJECT_TOKEN=dc26f09bc8074b28828a2ae96d70bf6f
bash <(curl -Ls https://coverage.codacy.com/get.sh) report -r coverage.xml

#google-chrome htmlcov/index.html
