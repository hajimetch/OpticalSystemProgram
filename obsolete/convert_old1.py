# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 09:59:09 2018

@author: hajimetch
"""

import configparser
from enum import Enum


class Projection(Enum):
    SPOT = 'spot'
    PATTERN = 'pattern'


class Reflector(Enum):
    TRUE = 'true'
    FALSE = 'false'


class Telephoto(Enum):
    TRUE = 'true'
    FALSE = 'false'

    def number(self):
        if self == Telephoto.TRUE:
            return 1
        elif self == Telephoto.FALSE:
            return 0
        else:
            return None


class InputDataType(Enum):
    SURFACE = 'surface data'
    POWER = 'power data'
    POWER_LIGHT = 'power and light data'
    GROUP = 'group data'

    def number(self):
        if self == InputDataType.SURFACE:
            return 0
        elif self == InputDataType.POWER:
            return 1
        elif self == InputDataType.POWER_LIGHT:
            return 2
        elif self == InputDataType.GROUP:
            return 3
        else:
            return None


class TracingNumber(Enum):
    SINGLE = 'single'
    MULTI = 'multi'

    def number(self):
        if self == TracingNumber.SINGLE:
            return 0
        elif self == TracingNumber.MULTI:
            return 5
        else:
            return None


class IntervalAdjustment(Enum):
    NONE = 'none'
    FOCAL_LENGTH = 'specify focal length'
    FOCAL_POSITION = 'specify focal position'
    TELEPHOTO_POWER_0 = 'telephoto system power set to 0'
    FOCAL_POSITION_LENGTH = 'specify focal position and length'

    def number(self):
        if self == IntervalAdjustment.NONE:
            return 0
        elif self == IntervalAdjustment.FOCAL_LENGTH:
            return 1
        elif self == IntervalAdjustment.FOCAL_POSITION:
            return 2
        elif self == IntervalAdjustment.TELEPHOTO_POWER_0:
            return 3
        elif self == IntervalAdjustment.FOCAL_POSITION_LENGTH:
            return 4
        else:
            return None


class ObjectDistance(Enum):
    INFINITY = 'to infinity'
    FIRST_POWER = 'to first power'
    FIRST_SURFACE = 'to first surface'
    PRINCIPAL_POSITION = 'to principal position of projection side'

    def number(self):
        if self == ObjectDistance.INFINITY:
            return 0
        elif self == ObjectDistance.FIRST_POWER:
            return 1
        elif self == ObjectDistance.FIRST_SURFACE:
            return 2
        elif self == ObjectDistance.PRINCIPAL_POSITION:
            return 3
        else:
            return None


class IncidentLightRange(Enum):
    SYSTEM_SOURCE = 'by system and light source'
    APERTURE_IRIS_SOURCE = 'by lens aperture iris and light source'
    APERTURE_IRIS = 'by lens aperture and iris'

    def number(self):
        if self == IncidentLightRange.SYSTEM_SOURCE:
            return 1
        elif self == IncidentLightRange.APERTURE_IRIS_SOURCE:
            return 2
        elif self == IncidentLightRange.APERTURE_IRIS:
            return 3
        else:
            return None


class TracingDirection(Enum):
    SCREEN_TO_APERTURE = 'screen to aperture'
    SCREEN_TO_APERTURE_EXT = 'screen to aperture extended'
    SOURCE_TO_APERTURE = 'light source to aperture'
    SOURCE_TO_APERTURE_EXT = 'light source to aperture extended'

    def number(self):
        if self == TracingDirection.SCREEN_TO_APERTURE:
            return 0
        elif self == TracingDirection.SCREEN_TO_APERTURE_EXT:
            return -1
        elif self == TracingDirection.SOURCE_TO_APERTURE:
            return 1
        elif self == TracingDirection.SOURCE_TO_APERTURE_EXT:
            return 2
        else:
            return None


class SplitPattern(Enum):
    NONE = 'none'
    FIRST_INTERVAL = 'first power and interval'
    LAST_INTERVAL = 'last power and interval'
    RATIO_INTERVAL = 'power ratio and interval'
    FIRST_LAST = 'first and last powers'

    def number(self):
        if self == SplitPattern.NONE:
            return 0
        elif self == SplitPattern.FIRST_INTERVAL:
            return 1
        elif self == SplitPattern.LAST_INTERVAL:
            return 2
        elif self == SplitPattern.RATIO_INTERVAL:
            return 3
        elif self == SplitPattern.FIRST_LAST:
            return 4
        else:
            return None


class PTLCurvatureSpecification(Enum):
    FRONT = 'front'
    BACK = 'back'
    FRONT_OVER_BACK = 'front/back'
    BACK_OVER_FRONT = 'back/front'
    FIRST_LENS = 'first lens'
    SECOND_LENS = 'second lens'
    APERTURE_MIRROR = 'aperture or mirror'
    REFLECTOR = 'reflector'

    def number(self):
        if self == PTLCurvatureSpecification.FRONT:
            return 0
        elif self == PTLCurvatureSpecification.BACK:
            return 1
        elif self == PTLCurvatureSpecification.FRONT_OVER_BACK:
            return 2
        elif self == PTLCurvatureSpecification.BACK_OVER_FRONT:
            return 3
        elif self == PTLCurvatureSpecification.FIRST_LENS:
            return 4
        elif self == PTLCurvatureSpecification.SECOND_LENS:
            return 5
        elif self == PTLCurvatureSpecification.APERTURE_MIRROR:
            return 6
        elif self == PTLCurvatureSpecification.REFLECTOR:
            return 7
        else:
            return None


class ScreenDisplayDimension(Enum):
    XYPLANE_COORDINATES = 'XY-plane and ray passing point coordinates'
    XYPLANE_3D = 'XY-plane and 3D'
    XYPLANE_XZPLANE = 'XY-plane and XZ-plane'
    SPLIT = 'split view'

    def number(self):
        if self == ScreenDisplayDimension.XYPLANE_COORDINATES:
            return 0
        elif self == ScreenDisplayDimension.XYPLANE_3D:
            return 1
        elif self == ScreenDisplayDimension.XYPLANE_XZPLANE:
            return 2
        elif self == ScreenDisplayDimension.SPLIT:
            return 3
        else:
            return None


class ElementBoundsException(Exception):
    pass


def writedataline(fout, config, args):
    data = []
    for e in args:
        data.append(getdata(fout, config, e))
    fout.write(','.join(data) + '\n')


def getdata(fout, config, e):
    if len(e) == 1:
        return e[0]
    elif len(e) == 2:
        return config.get(e[0], e[1], fallback='')
    elif len(e) == 3:
        return config.get(e[0], e[1], fallback=e[2])
    else:
        raise ElementBoundsException(Exception)


def validoneof(config, args):
    for e in args:
        if (config.get(e[0], e[1], fallback='') != ''):
            return [e]
    return []


def validsof(config, args):
    ones = []
    for e in args:
        if (config.get(e[0], e[1], fallback='') != ''):
            ones.append(e)
    return ones


config = configparser.ConfigParser()
config.read('example.ini')

Projection = Projection(
    config.get('Parameters', 'projection', fallback='spot'))
Reflector = Reflector(config.get('Parameters', 'reflector', fallback='false'))
Telephoto = Telephoto(config.get('Parameters', 'telephoto', fallback='false'))
InputDataType = InputDataType(
    config.get('Parameters', 'input data type', fallback='surface data'))
TracingNumber = TracingNumber(
    config.get('Parameters', 'tracing number', fallback='single'))
IntervalAdjustment = IntervalAdjustment(
    config.get(
        'Parameters', 'lenses/powers interval adjustment', fallback='none'))
ObjectDistance = ObjectDistance(
    config.get('Parameters', 'object distance', fallback='to infinity'))
IncidentLightRange = IncidentLightRange(
    config.get(
        'Parameters',
        'incident light range definition',
        fallback='by system and light source'))
TracingDirection = TracingDirection(
    config.get(
        'Parameters', 'tracing direction', fallback='screen to aperture'))
ScreenDisplayDimension = ScreenDisplayDimension(
    config.get(
        'View Parameters',
        'screen display dimension',
        fallback='XY-plane and ray passing point coordinates'))

SplitPatterns = []
for e in config.get('System', 'split patterns', fallback='none').split(','):
    SplitPatterns.append(SplitPattern(e.strip()))

PTLCurvatureSpecifications = []
for e in config.get(
        'System', 'power to lens curvature specifications',
        fallback='front').split(','):
    PTLCurvatureSpecifications.append(PTLCurvatureSpecification(e.strip()))

with open('datafile.txt', 'wt') as fout:
    fout.write('*SDT:data ')
    sysname = 'FS' if Projection == Projection.SPOT else ''
    sysname += 'REF' if Reflector == Reflector.TRUE else ''
    sysname += config.get('Parameters', 'system name').strip()
    writedataline(fout, config, [[sysname]])

    fout.write('*IAF:data ')
    writedataline(fout, config, [[str(Telephoto.number())]])

    fout.write('*MSDT:data ')
    args = []
    args.append(['Telephoto', 'distance between systems'])
    args.append(['Telephoto', 'focal length'])
    args.append(['Telephoto', 'principal position of object side'])
    args.append(['Telephoto', 'principal position of projection side'])
    args.append(['Telephoto', 'aperture for display'])
    args.append(['Telephoto', 'refractive index of projection region'])
    args.append(['Telephoto', 'adjusting lenses/powers interval number'])
    writedataline(fout, config, args)

    fout.write('*SCL:data ')
    writedataline(fout, config, [['Parameters', 'scale', '100']])

    fout.write('*NPSL:data ')
    args = []
    args.append(['System', 'number of powers'])
    args.append(['System', 'number of surfaces'])
    args.append(['System', 'number of lenses'])
    if (InputDataType == InputDataType.GROUP):
        args.append(['System', 'number of groups'])
    ones = []
    ones.append(['System', 'adjusting surface interval number'])
    ones.append(['System', 'adjusting power interval number'])
    args.extend(validoneof(config, ones))
    if (IntervalAdjustment == IntervalAdjustment.FOCAL_POSITION or
            IntervalAdjustment == IntervalAdjustment.FOCAL_POSITION_LENGTH):
        args.append(['System', 'focal position'])
    if (IntervalAdjustment == IntervalAdjustment.TELEPHOTO_POWER_0):
        args.append(['Telephoto', 'system position'])
    if (IntervalAdjustment == IntervalAdjustment.FOCAL_LENGTH or
            IntervalAdjustment == IntervalAdjustment.FOCAL_POSITION_LENGTH):
        args.append(['System', 'focal length'])
    if (IntervalAdjustment == IntervalAdjustment.FOCAL_POSITION_LENGTH):
        args.append(['System', 'iris position'])
    if (IntervalAdjustment == IntervalAdjustment.TELEPHOTO_POWER_0):
        ones = []
        ones.append(['Telephoto', 'number of surfaces'])
        ones.append(['Telephoto', 'number of powers'])
        args.extend(validoneof(config, ones))
    writedataline(fout, config, args)

    fout.write('*NSPR:data ')
    args = []
    args.append(['System', 'surfaces of powers'])
    args.append(['System', 'division ratios of powers'])
    writedataline(fout, config, args)

    fout.write('*INC:data ')
    args = []
    args.append([str(ObjectDistance.number())])
    args.append([str(IncidentLightRange.number())])
    ones = []
    ones.append(['System', 'aperture radius'])
    ones.append(['System', 'light source radius'])
    args.extend(validoneof(config, ones))
    args.append(['System', 'object distance'])
    args.append([str(TracingDirection.number())])
    if (TracingDirection == TracingDirection.SCREEN_TO_APERTURE_EXT):
        args.append(['System', 'aperture to light source'])
        args.extend(validsof(config, [['System', 'mirror to valve']]))
    elif (TracingDirection == TracingDirection.SOURCE_TO_APERTURE or
          TracingDirection == TracingDirection.SOURCE_TO_APERTURE_EXT):
        ones = []
        ones.append(['System', 'aperture position'])
        ones.append(['System', 'aperture radius real'])
        args.extend(validsof(config, ones))
    writedataline(fout, config, args)

    fout.write('*INCA:data ')
    args = []
    if (IncidentLightRange == IncidentLightRange.SYSTEM_SOURCE or
            IncidentLightRange == IncidentLightRange.APERTURE_IRIS_SOURCE):
        args.append(['Light Source', 'light source file'])
        args.append(['Light Source', 'diffusion radius'])
        args.append(['Light Source', 'filament number'])
        args.append(['Light Source', 'filament position'])
        args.append(['Light Source', 'filament length'])
        args.append(['Light Source', 'number of aperture segments'])
        args.append(['Light Source', 'angle range'])
        args.append(['Light Source', 'number of tracing rays'])
        args.append(['Light Source', 'iris surface number'])
        args.append(['Light Source', 'limit surface number'])
    writedataline(fout, config, args)

    fout.write('*SCRN:data ')
    args = []
    if (TracingDirection == TracingDirection.SOURCE_TO_APERTURE):
        args.append(['System', 'light source to screen'])
        args.append(['System', 'screen radius'])
    writedataline(fout, config, args)

    fout.write('*EF:data ')
    writedataline(fout, config, [[
        str(InputDataType.number() +
            (TracingNumber.number() + IntervalAdjustment.number()) * 10)
    ]])

    fout.write('*DSCV:data ')
    args = []
    if (InputDataType == InputDataType.SURFACE):
        args.append(['System', 'surface intervals'])
        args.append(['System', 'surface curvature radiuses'])
    writedataline(fout, config, args)

    fout.write('*DSPW:data ')
    args = []
    if (InputDataType == InputDataType.POWER):
        args.append(['System', 'power intervals'])
        args.append(['System', 'power values'])
    elif (InputDataType == InputDataType.POWER_LIGHT):
        args = []
        args.append(['System', 'incident heights at powers start'])
        args.append(['System', 'power intervals'])
        args.append(['System', 'incident heights at powers'])
    elif (InputDataType == InputDataType.GROUP):
        args = []
        args.append(['System', 'incident heights at groups start'])
        args.append(['System', 'group intervals'])
        args.append(['System', 'incident heights at groups'])
    writedataline(fout, config, args)

    fout.write('*PSPL:data ')
    args = []
    if (InputDataType == InputDataType.GROUP):
        for e in SplitPatterns:
            args.append([str(e.number())])
        args.append(['System', 'split parameters first'])
        args.append(['System', 'split parameters second'])
    writedataline(fout, config, args)

    fout.write('*GLAS:data ')
    writedataline(fout, config, [['System', 'glass type names']])

    fout.write('*DSTK:data ')
    writedataline(fout, config, [['System', 'surface intervals of powers']])

    fout.write('*CURV:data ')
    args = []
    for e in PTLCurvatureSpecifications:
        args.append([str(e.number())])
    args.append(['System', 'curvature specification values'])
    writedataline(fout, config, args)

    fout.write('*LRAD:data ')
    writedataline(fout, config, [['System', 'lens radiuses']])

    fout.write('*LAR:data ')
    writedataline(fout, config, [['System', 'lens effective radiuses']])

    fout.write('*ALR:data ')
    writedataline(fout, config, [['System', 'surfaces axis length ratios']])

    fout.write('*ASPH:data ')
    writedataline(fout, config, [['System', 'surfaces aspheric coefficient']])

    fout.write('*DISP:data ')
    args = []
    args.append(['View Parameters', 'X-axis length'])
    args.append(['View Parameters', 'screen edge to first surface'])
    args.append(['View Parameters', 'Y-coordinate at screen center'])
    args.append(['View Parameters', 'X-axis grid pitch'])
    args.append(['View Parameters', 'Y-axis grid pitch'])
    args.append([str(ScreenDisplayDimension.number())])
    writedataline(fout, config, args)

    fout.write('*DEND:data ')
    args = []
    args.append(['Light Source', 'number of mirror splitting'])
    args.append(['Light Source', 'number of mirror facets'])
    args.append(['Light Source', 'heights of mirror facet edges'])
    args.append(
        ['Light Source', 'angles between mirror facets and optical axis'])
    writedataline(fout, config, args)
