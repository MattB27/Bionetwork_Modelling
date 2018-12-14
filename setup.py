#!/usr/bin/env python

from setuptools import setup
from setuptools.command.test import test as TestCommand
from setuptools import Command
import sys
import os
from subprocess import call


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['tests']
        self.test_suite = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


class PyDoc(Command):
    description = 'run documentation'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        cwd = os.getcwd()
        os.chdir('docs/')
        call(['make', 'html'])
        os.chdir(cwd)


setup(
    name='sbml_reader',
    version='0.0.1',
    packages=['sbml_reader'],
    include_package_data=True,
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'jupyter',
        'gekko'
    ],
    tests_require=[
        'pytest'
    ],
    cmdclass={
        'test': PyTest,
        'document': PyDoc
    },
)
