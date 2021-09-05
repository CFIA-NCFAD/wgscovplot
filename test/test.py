#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd


def inc(x):
    return x + 1


def test_answer():
    assert inc(4) == 5
