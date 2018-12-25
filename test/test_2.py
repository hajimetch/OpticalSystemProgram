# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 16:28:22 2018

@author: hajimetch
"""

from enum import Enum


class Color(Enum):
    red = (1, '赤')
    green = (2, '緑')
    blue = (3, '青')

    def __new__(cls, value, japanese_name):
        obj = object.__new__(cls)
        print('object.__new__ called')
        obj._value_ = value
        obj.japanese_name = japanese_name
        return obj


print('before Color(1) called')
color = Color(1)
print('after Color(1) called')
color2 = Color(2)
print('after Color(2) called')
