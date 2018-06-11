#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

requirements = [i.strip() for i in open('lib_requirements.txt').readlines()]

setup(
    name='larval_gonad_ovary',
    version='0.0.1',
    description="Local library for the larval gonad ovary project",
    author="Justin M Fear",
    author_email='justin.m.fear@gmail.com',
    url='https://github.com/jfear/larval_gonad_ovary',
    packages=['lib'],
    install_requires=requirements,
    license="MIT license",
)
