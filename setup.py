#!/usr/bin/env python3


from os.path import dirname, join
from glob import glob
from setuptools import setup, find_packages


DIR = (dirname(__file__) or '.')

setup_args = {}
setup_args.update(
      name='ACE',
      version='0.0.1',
      description='ACE',
      author='Jin Seok (Andy) Lee, Dhuvarakesh Karthikeyan',
      author_email='ajslee@unc.edu, dhuvik@ad.unc.edu',
      packages=find_packages(),
      scripts=glob(join(DIR, 'scripts/*.py'))
)

setup(**setup_args)