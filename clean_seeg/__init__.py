"""Automated preprocessing of sEEG data"""
__version__ = '0.1'

from .clean_seeg import cleanSEEG
from .clean_PLI import removePLI, zapline, cleanline