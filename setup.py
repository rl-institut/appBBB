#! /usr/bin/env python

from setuptools import setup

setup(name='appBBB',
      version='0.0.1',
      author='Elisa Gaudchau, Birgit Schachler',
      author_email='elisa.gaudchau@rl-institut.de',
      description='A model of the heat and power systems of Brandenburg and Berlin',
      package_dir={'appBBB': 'appBBB'},
      install_requires=['oemof == 0.0.9',
                        'feedinlib == 0.0.10',
                        'oemof.db'],
      dependency_links = ['http://github.com/oemof/oemof.db/tarball/master'] 
      )
